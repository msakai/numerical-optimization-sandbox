{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
-- translated from https://github.com/chokkan/liblbfgs
module LineSearch
  ( Params (..)
  , defaultParams
  , lineSearch
  , lineSearchMoreThuente
  ) where

import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra ((<.>))


clip :: Ord a => a -> a -> a -> a
clip lo hi x
  | hi < x = hi
  | x < lo = lo
  | otherwise = x

midpoint :: Fractional a => a -> a -> a
midpoint x y = x + 0.5 * (y - x)

signdiff :: (Fractional a, Ord a) => a -> a -> Bool
signdiff x y = x * (y / abs y) < 0


-- | The minimizer of the interpolated cubic.
cubicMinimizer
  :: (Ord a, Floating a)
  => a -- ^ The value of one point, @u@.
  -> a -- ^ The value of @f(u)@.
  -> a -- ^ The value of @f'(u)@.
  -> a -- ^ The value of another point, @v@.
  -> a -- ^ The value of @f(v)@.
  -> a -- ^ The value of @f'(v)@.
  -> a
cubicMinimizer u fu du v fv dv = u + r * d
  where
    d = v - u
    theta = (fu - fv) * 3 / d + du + dv
    s = maximum [p, q, r]
      where
        p = abs theta
        q = abs du
        r = abs dv
    a = theta / s
    gamma = (if v < u then negate else id) $
            s * sqrt (a * a - (du / s) * (dv / s))
    p = gamma - du + theta
    q = gamma - du + gamma + dv
    r = p / q


-- | The minimizer of the interpolated cubic.
cubicMinimizer2
  :: (Ord a, Floating a)
  => a -- ^ The value of one point, @u@.
  -> a -- ^ The value of @f(u)@.
  -> a -- ^ The value of @f'(u)@.
  -> a -- ^ The value of another point, @v@.
  -> a -- ^ The value of @f(v)@.
  -> a -- ^ The value of @f'(v)@.
  -> a -- ^ The minimum value.
  -> a -- ^ The maximum value.
  -> a
cubicMinimizer2 u fu du v fv dv xmin xmax
  | r < 0.0 && gamma /= 0.0 = v - r * d
  | a < 0 = xmax
  | otherwise = xmin
  where
    d = v - u
    theta = (fu - fv) * 3 / d + du + dv
    s = maximum [p, q, r]
      where
        p = abs theta
        q = abs du
        r = abs dv
    a = theta / s
    gamma = (if u < v then negate else id) $
            s * sqrt (max 0 (a * a - (du / s) * (dv / s)))
    p = gamma - dv + theta
    q = gamma - dv + gamma + du
    r = p / q


quardMinimizer
  :: (Ord a, Fractional a)
  => a -- ^ The value of one point, @u@.
  -> a -- ^ The value of @f(u)@.
  -> a -- ^ The value of @f'(u)@.
  -> a -- ^ The value of another point, @v@.
  -> a -- ^ The value of @f(v)@.
  -> a
quardMinimizer u fu du v fv = u + du / ((fu - fv) / a + du) / 2 * a
  where
    a = v - u


quardMinimizer2
  :: (Ord a, Fractional a)
  => a -- ^ The value of one point, @u@.
  -> a -- ^ The value of @f'(u)@.
  -> a -- ^ The value of another point, @v@.
  -> a -- ^ The value of @f'(v)@.
  -> a
quardMinimizer2 u du v dv = v + dv / (dv - du) * a
  where
    a = u - v


data Error
  = ERR_OUTOFINTERVAL
  | ERR_INCREASEGRADIENT
  | ERR_INCORRECT_TMINMAX
  | ERR_INVALIDPARAMETERS
  | ERR_MAXIMUMSTEP
  | ERR_MINIMUMSTEP
  | ERR_ROUNDING_ERROR
  | ERR_WIDTHTOOSMALL
  | ERR_MAXIMUMLINESEARCH
  deriving (Eq, Ord, Enum, Bounded, Show)


{- | Update a safeguarded trial value and interval for line search.

The parameter @x@ represents the step with the least function value.
The parameter @t@ represents the current step. This function assumes
that the derivative at the point of @x@ in the direction of the step.
If the bracket is set to true, the minimizer has been bracketed in
an interval of uncertainty with endpoints between @x@ and @y@.
-}
updateTrialInterval
  :: (Ord a, Floating a)
  => (a, a, a) -- ^ The value of one endpoint @x@, the value of @f(x)@ and the value of @f'(x)@
  -> (a, a, a) -- ^ The value of another endpoint @y@, the value of @f(y)@ and the value of @f'(y)@
  -> (a, a, a) -- ^ The value of the trial value @t@, the value of @f(t)@ and the value of @f'(t)@
  -> a         -- ^ The minimum value for the trial value, @t@.
  -> a         -- ^ The maximum value for the trial value, @t@.
  -> Bool      -- ^ The predicate if the trial value is bracketed.
  -> Either Error (a, (a, a, a), (a, a, a), Bool)
     -- ^ 'Error' or new trial value, updated @(x, f(x), f'(x))@, updated @(y, f(y), f'(y))@, and updated bracketed predicate.
updateTrialInterval (x, fx, dx) (y, fy, dy) (t, ft, dt) tmin tmax brackt
  | brackt && (t <= min x y || max x y <= t) =
      -- The trival value t is out of the interval.
      Left ERR_OUTOFINTERVAL
  | brackt && (0 <= dx * (t - x)) =
      -- The function must decrease from x.
      Left ERR_INCREASEGRADIENT
  | brackt && tmax < tmin =
      -- Incorrect tmin and tmax specified.
      Left ERR_INCORRECT_TMINMAX
  | otherwise = Right (newt3, (x', fx', dx'), (y', fy', dy'), brackt')
  where
    dsign = signdiff dt dx

    (newt1, brackt', bound)
      | fx < ft =
          {-
            Case 1: a higher function value.
            The minimum is brackt. If the cubic minimizer is closer
            to x than the quadratic one, the cubic one is taken, else
            the average of the minimizers is taken.
          -}
          let mc = cubicMinimizer x fx dx t ft dt
              mq = quardMinimizer x fx dx t ft
          in ( if abs (mc - x) < abs (mq - x)
               then mc
               else midpoint mc mq
             , True
             , True
             )
      | dsign =
          {-
            Case 2: a lower function value and derivatives of
            opposite sign. The minimum is brackt. If the cubic
            minimizer is closer to x than the quadratic (secant) one,
            the cubic one is taken, else the quadratic one is taken.
           -}
          let mc = cubicMinimizer x fx dx t ft dt
              mq = quardMinimizer2 x dx t dt
          in ( if abs (mc - t) > abs (mq - t)
               then mc
               else mq
             , True
             , False
             )
      | abs dt < abs dx =
          {-
            Case 3: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative decreases.
            The cubic minimizer is only used if the cubic tends to
            infinity in the direction of the minimizer or if the minimum
            of the cubic is beyond t. Otherwise the cubic minimizer is
            defined to be either tmin or tmax. The quadratic (secant)
            minimizer is also computed and if the minimum is brackt
            then the the minimizer closest to x is taken, else the one
            farthest away is taken.
           -}
          let mc = cubicMinimizer2 x fx dx t ft dt tmin tmax
              mq = quardMinimizer2 x dx t dt
          in ( if brackt then
                 if abs (t - mc) < abs (t - mq)
                 then mc
                 else mq
               else
                 if abs (t - mc) > abs (t - mq)
                 then mc
                 else mq
             , brackt
             , True
             )
      | otherwise =
          {-
            Case 4: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative does
            not decrease. If the minimum is not brackt, the step
            is either tmin or tmax, else the cubic minimizer is taken.
           -}
          ( if brackt then
              cubicMinimizer t ft dt y fy dy
            else if x < t then
              tmax
            else
              tmin
          , brackt
          , False
          )

    {-
        Update the interval of uncertainty. This update does not
        depend on the new step or the case analysis above.

        - Case a: if f(x) < f(t),
            x <- x, y <- t.
        - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
            x <- t, y <- y.
        - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0,
            x <- t, y <- x.
     -}
    (x', fx', dx', y', fy', dy')
      | fx < ft   = (x, fx, dx, t, ft, dt) -- Case a
      | dsign     = (t, ft, dt, x, fx, dx) -- Case c
      | otherwise = (t, ft, dt, y, fy, dy) -- Case b

    -- Clip the new trial value in [tmin, tmax].
    newt2 = clip tmin tmax newt1

    -- Redefine the new trial value if it is close to the upper bound of the interval.
    newt3
      | brackt' && bound && (if x' < y' then mq < newt2 else newt2 < mq) = mq
      | otherwise = newt2
      where
        delta = 0.66
        mq = x' + delta * (y' - x')


data Params a
  = Params
  { paramsMinStep :: a
    -- ^ The minimum step of the line search routine.
    --
    -- The default value is @1e-20@. This value need not be modified unless
    -- the exponents are too large for the machine being used, or unless the
    -- problem is extremely badly scaled (in which case the exponents should
    -- be increased).

  , paramsMaxStep :: a
    -- ^ The maximum step of the line search.
    --
    -- The default value is @1e+20@. This value need not be modified unless
    -- the exponents are too large for the machine being used, or unless the
    -- problem is extremely badly scaled (in which case the exponents should
    -- be increased).

  , paramsFTol :: a
    -- ^ A parameter to control the accuracy of the line search routine.
    --
    -- The default value is @1e-4@. This parameter should be greater
    -- than zero and smaller than @0.5@.
    -- "μ" in [MoreThuente1994].

  , paramsGTol :: a
    -- ^ A parameter to control the accuracy of the line search routine.
    --
    -- The default value is @0.9@. If the function and gradient
    -- evaluations are inexpensive with respect to the cost of the
    -- iteration (which is sometimes the case when solving very large
    -- problems) it may be advantageous to set this parameter to a small
    -- value. A typical small value is @0.1@. This parameter should be
    -- greater than the 'paramsFTol' parameter (@1e-4@) and smaller than
    -- @1.0@.
    -- "η" in [MoreThuente1994].

  , paramsXTol :: a
    -- ^ The machine precision for floating-point values.
    --
    -- This parameter must be a positive value set by a client program to
    -- estimate the machine precision. The line search routine will terminate
    -- with the status code ('ERR_ROUNDING_ERROR') if the relative width
    -- of the interval of uncertainty is less than this parameter.

  , paramsMaxLineSearch :: Int
    -- ^ The maximum number of trials for the line search.
    --
    -- This parameter controls the number of function and gradients evaluations
    -- per iteration for the line search routine. The default value is @40@.
  }

defaultParams :: Floating a => Params a
defaultParams
  = Params
  { paramsMinStep = 1e-20
  , paramsMaxStep = 1e+20
  , paramsFTol = 1e-4
  , paramsGTol = 0.9
  , paramsXTol = 1.0e-16
  , paramsMaxLineSearch = 40
  }

lineSearch
  :: forall a. (Ord a, Floating a, LA.Numeric a)
  => Params a
  -> (LA.Vector a -> (a, LA.Vector a))
  -> (LA.Vector a, a, LA.Vector a)
  -> LA.Vector a
  -> a
  -> (Maybe Error, a, (LA.Vector a, a, LA.Vector a))
lineSearch = lineSearchMoreThuente

lineSearchMoreThuente
  :: forall a. (Ord a, Floating a, LA.Numeric a)
  => Params a
  -> (LA.Vector a -> (a, LA.Vector a))
  -> (LA.Vector a, a, LA.Vector a)
  -> LA.Vector a
  -> a
  -> (Maybe Error, a, (LA.Vector a, a, LA.Vector a))
lineSearchMoreThuente params evaluate (x0, f0, g0) s step0
  | step0 < 0 = (Just ERR_INVALIDPARAMETERS, 0, (x0, f0, g0))
  | 0 < dg0   = (Just ERR_INCREASEGRADIENT, 0, (x0, f0, g0))
  | otherwise = seq dgtest $ go 0 Nothing False True width0 prevWidth0 (0, f0, dg0) (0, f0, dg0) step0
  where
    -- φ(step) = f(x0 + step * s) = f(x)
    dg0 = g0 <.> s                   -- φ'(0)
    dgtest = paramsFTol params * dg0 -- μ φ'(0)
    width0 = paramsMaxStep params - paramsMinStep params
    prevWidth0 = 2.0 * width0

    go
      :: Int
      -> Maybe Error
      -> Bool -> Bool
      -> a
      -> a
      -> (a, a, a)
      -> (a, a, a)
      -> a
      -> (Maybe Error, a, (LA.Vector a, a, LA.Vector a))
    go count uinfo brackt stage1 width prev_width (stx, fx, dgx) (sty, fy, dgy) step_
      | brackt && (step <= stmin || stmax <= step || uinfo /= Nothing) =
          -- Rounding errors prevent further progress.
          (Just ERR_ROUNDING_ERROR, step, (x, f, g))
      | step == paramsMaxStep params && sufficientDecrease && dg <= dgtest =
          -- The step is the maximum value.
          (Just ERR_MAXIMUMSTEP, step, (x, f, g))
      | step == paramsMinStep params && (not sufficientDecrease || dgtest <= dg) =
          -- The step is the minimum value.
          (Just ERR_MINIMUMSTEP, step, (x, f, g))
      | brackt && stmax - stmin <= paramsXTol params * stmax =
          -- Relative width of the interval of uncertainty is at most xtol.
          (Just ERR_WIDTHTOOSMALL, step, (x, f, g))
      | paramsMaxLineSearch params <= count + 1 =
          -- Maximum number of iteration
          (Just ERR_MAXIMUMLINESEARCH, step, (x, f, g))
      | sufficientDecrease && abs dg <= paramsGTol params * (-dg0) =
          -- The sufficient decrease condition and the directional derivative condition hold.
          (Nothing, step, (x, f, g))
      | otherwise =
          go (count + 1) uinfo' brackt' stage1' width' prev_width' (stx', fx', dgx') (sty', fy', dgy') step''
      where
        {-
          Set the minimum and maximum steps to correspond to the
          present interval of uncertainty.
         -}
        (stmin, stmax)
          | brackt    = (min stx sty, max stx sty)
          | otherwise = (stx, step_ + 4 * (step_ - stx))

        -- Clip the step in the range of [stepmin, stepmax].
        step = f $ clip (paramsMinStep params) (paramsMaxStep params) $ step_
          where
            {- If an unusual termination is to occur then let
               step be the lowest point obtained so far.
             -}
            f step
              | brackt && (step <= stmin || stmax <= step || paramsMaxLineSearch params <= count + 1 || uinfo /= Nothing) = stx
              | brackt && (stmax - stmin <= paramsXTol params * stmax) = stx
              | otherwise = step

        -- Compute the current value of x
        x = x0 `LA.add` LA.scale step s

        -- Evaluate the function and gradient values.
        (f, g) = evaluate x
        dg = g <.> s

        -- φ(α) <= φ(0) + α μ φ'(0)
        sufficientDecrease :: Bool
        sufficientDecrease = f <= f0 + step * dgtest

        {-
          In the first stage we seek a step for which the modified
          function has a nonpositive value and nonnegative derivative.

          ψ(α) <= 0
          ⇔ φ(α) - φ(0) - μ φ'(0) <= 0
          ⇔ φ(α) <= φ(0) + μ φ'(0)
          ⇔ sufficientDecrease

          ψ'(α) >= 0
          ⇔ φ'(α) >= 0
          で、min(μ,η) φ'(0) <= φ'(α) の左辺は少し余裕を持たせている?

          ftol (μ) の方が gtol (η) より小さい前提のはずなのに min をとっているのは何故?
         -}
        stage1'
          | stage1 && sufficientDecrease && min (paramsFTol params) (paramsGTol params) * dg0 <= dg = False
          | otherwise = stage1

        (uinfo', step', (stx', fx', dgx'), (sty', fy', dgy'), brackt')
            | stage1' && not sufficientDecrease && f <= fx =
                {- A modified function is used to predict the step only if
                   we have not obtained a step for which the modified
                   function has a nonpositive function value and nonnegative
                   derivative, and if a lower function value has been
                   obtained but the decrease is not sufficient.
                 -}
                let -- Define the modified function and derivative values.
                    -- dgtest = μ φ'(0) で
                    -- ψ(α) = φ(α) - φ(0) - μ φ'(0) α だとすると - φ(0) の項は定数差なので無視されている?
                    fm = f - step * dgtest
                    fxm = fx - stx * dgtest
                    fym = fy - sty * dgtest
                    dgm = dg - dgtest
                    dgxm = dgx - dgtest
                    dgym = dgy - dgtest
                in  case updateTrialInterval (stx, fxm, dgxm) (sty, fym, dgym) (step, fm, dgm) stmin stmax brackt of
                      Left err ->
                        ( Just err
                        , step
                        , (stx, fx, dgx)
                        , (sty, fy, dgy)
                        , brackt
                        )
                      Right (step', (stx', fxm', dgxm'), (sty', fym', dgym'), brackt') ->
                        ( Nothing
                        , step'
                        , (stx', fxm' + stx' * dgtest, dgxm' + dgtest)
                        , (sty', fym' + sty' * dgtest, dgym' + dgtest)
                        , brackt'
                        )
            | otherwise =
                case updateTrialInterval (stx, fx, dgx) (sty, fy, dgy) (step, f, dg) stmin stmax brackt of
                  Left err ->
                    ( Just err
                    , step
                    , (stx, fx, dgx)
                    , (sty, fy, dgy)
                    , brackt
                    )
                  Right (step', (stx', fx', dgx'), (sty', fy', dgy'), brackt') ->
                    ( Nothing
                    , step'
                    , (stx', fx', dgx')
                    , (sty', fy', dgy')
                    , brackt'
                    )

        (step'', prev_width', width')
            | brackt' =
                ( if 0.66 * prev_width <= abs (sty' - stx')
                  then midpoint stx' sty'
                  else step'
                , width
                , abs (sty' - stx')
                )
            | otherwise = (step', prev_width, width)
