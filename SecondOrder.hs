{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module SecondOrder
  ( newtonMethod
  , gaussNewton
  , levenbergMarquardt
  , bfgs
  , bfgsV
  , lbfgs
  , lbfgsV

  , rosenbrock
  ) where

import qualified Data.Foldable as F
import Data.Reflection (Reifies)
import Data.Sequence (Seq, ViewL (..), (<|))
import qualified Data.Sequence as Seq
import qualified Data.Traversable as T
import qualified Data.Vector.Generic as VG
import Foreign.Storable
import Numeric.AD
import Numeric.AD.Internal.Reverse (Reverse, Tape)
import Numeric.AD.Rank1.Sparse (Sparse)
import Numeric.LinearAlgebra
import qualified LineSearch as LS


zipWithTV :: (Traversable f, Storable b) => (a -> b -> c) -> f a -> Vector b -> f c
zipWithTV u x v = snd $ T.mapAccumL (\i x_i -> (i+1, u x_i (v VG.! i))) 0 x


-- | Perform a newton method.
--
-- >>> let sq x = x * x
-- >>> let rosenbrock [x,y] = sq (1 - x) + 100 * sq (y - sq x)
-- >>> rosenbrock [0,0]
-- 1
-- >>> rosenbrock (newtonMethod rosenbrock [0, 0] !! 2) < 0.1
-- True
newtonMethod
  :: forall f a. (Traversable f, Field a)
  => (forall s. f (AD s (Sparse a)) -> AD s (Sparse a))
  -> f a -> [f a]
newtonMethod f x0 = go x0
  where
    n = length x0

    go :: f a -> [f a]
    go x = x : go (zipWithTV (-) x (h <\> g))
      where
        _y :: a
        gh :: f (a, f a)
        (_y, gh) = hessian' f x

        g :: Vector a
        g = fromList $ map fst $ F.toList gh

        h :: Matrix a
        h = (n >< n) $ concat $ map (F.toList . snd) $ F.toList gh


gaussNewton
  :: forall f g a. (Traversable f, Traversable g, Field a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> g (Reverse s a))
  -> f a -> [f a]
gaussNewton f x0 = go x0
  where
    m = length x0

    go :: f a -> [f a]
    go x = x : go (zipWithTV (-) x (pinv j #> r))
      where
        rj :: g (a, f a)
        rj = jacobian' f x

        r :: Vector a
        r = fromList $ map fst $ F.toList rj

        j :: Matrix a
        j = (length rj >< m) $ concat $ map (F.toList . snd) (F.toList rj)


-- example from https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_gaussNewton = gaussNewton f x0
  where
    f :: Fractional a => [a] -> [a]
    f [beta1,beta2] = [y - beta1*x / (beta2 + x) | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]


-- | Levenberg–Marquardt algorithm with Tikhonov Dampling
levenbergMarquardt
  :: forall f g a. (Traversable f, Traversable g, Field a, Ord a)
  => a
  -> (forall s. Reifies s Tape => f (Reverse s a) -> g (Reverse s a))
  -> f a -> [f a]
levenbergMarquardt lambda0 f x0 = go lambda0 x0
  where
    m = length x0

    go :: a -> f a -> [f a]
    go lambda x = x : go lambda' x'
      where
        rj :: g (a, f a)
        rj = jacobian' f x

        n = length rj

        r :: Vector a
        r = fromList $ map fst $ F.toList rj

        j :: Matrix a
        j = (n >< m) $ concat $ map (F.toList . snd) (F.toList rj)

        gnMat :: Matrix a
        gnMat = unSym (mTm j) `add` diag (VG.replicate m lambda)

        delta :: Vector a
        delta = gnMat <\> scale (-1) (r <# j)

        x' :: f a
        x' = zipWithTV (+) x delta

        loss, loss' :: a
        loss  = sum [r_j * r_j | r_j <- toList r] / fromIntegral n
        loss' = sum [r_j * r_j | (r_j, _) <- F.toList (jacobian' f x')] / fromIntegral n

        approx :: Vector a -> a
        approx delta
          = loss
          + (scale (2 / fromIntegral n) r <# j) `dot` delta
          + (delta `dot` (scale (2 / fromIntegral n) gnMat #> delta)) / 2

        rho :: a
        rho = (loss' - loss) / (approx delta - approx (VG.replicate m 0))

        lambda' :: a
        lambda'
          | rho > 3/4 = lambda * 2 / 3
          | rho < 1/4 = lambda * 3 / 2
          | otherwise = lambda


-- example from https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_levenbergMarquardt = levenbergMarquardt 1.0 f x0
  where
    f :: Fractional a => [a] -> [a]
    f [beta1,beta2] = [y - beta1*x / (beta2 + x) | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]


-- | Broyden–Fletcher–Goldfarb–Shanno algorithm
--
-- https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm
bfgs
  :: forall f a. (Traversable f, Field a, Ord a, Normed (Vector a), Show a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> Reverse s a)
  -> f a -> [f a]
bfgs f x0 = map fromVector $ bfgsV evaluate (toVector x0)
  where
    fromVector :: Vector a -> f a
    fromVector = zipWithTV (\_ x -> x) x0

    toVector :: f a -> Vector a
    toVector = fromList . F.toList

    evaluate :: Vector a -> (a, Vector a)
    evaluate x =
      case grad' f (fromVector x) of
        (obj, g) -> (obj, toVector g)


-- | Limited-memory BFGS
--
-- https://en.wikipedia.org/wiki/Limited-memory_BFGS
lbfgs
  :: forall f a. (Traversable f, Field a, Ord a, Normed (Vector a), Show a)
  => Int
  -> (forall s. Reifies s Tape => f (Reverse s a) -> Reverse s a)
  -> f a -> [f a]
lbfgs m f x0 = map fromVector $ lbfgsV m evaluate (toVector x0)
  where
    fromVector :: Vector a -> f a
    fromVector = zipWithTV (\_ x -> x) x0

    toVector :: f a -> Vector a
    toVector = fromList . F.toList

    evaluate :: Vector a -> (a, Vector a)
    evaluate x =
      case grad' f (fromVector x) of
        (obj, g) -> (obj, toVector g)


bfgsV
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => (Vector a -> (a, Vector a))
  -> Vector a -> [Vector a]
bfgsV f x0 = go (ident n) alpha0 (x0, o0, g0)
  where
    n = VG.length x0
    (o0, g0) = f x0

    alpha0 :: a
    alpha0 = realToFrac $ 1 / norm_2 g0

    epsilon :: Double
    epsilon = 1e-5

    go :: Matrix a -> a -> (Vector a, a, Vector a) -> [Vector a]
    go bInv alpha_ (x, o, g) = x :
      if converged then
        []
      else
        case err of
          Just e -> error (show e)
          Nothing
            | sy > 0    -> go bInv' 1.0 (x', o', g')
            | otherwise -> error ("curvature condition failed: " ++ show sy)
      where
        converged :: Bool
        converged = norm_2 g / max (norm_2 x) 1 <= epsilon

        p :: Vector a
        p = scale (-1) $ bInv #> g

        (err, alpha, (x', o', g')) = LS.lineSearch LS.defaultParams f (x, o, g) p alpha_

        s, y :: Vector a
        s = scale alpha p
        y = g' `add` scale (-1) g
        sy :: a
        sy = s <.> y

        bInv' :: Matrix a
        bInv' =
            bInv
            `add` scale ((sy + y <.> (bInv #> y)) / sy**2) (s `outer` s)
            `add` scale (-1 / sy) (((bInv #> y) `outer` s) `add` (s `outer` (y <# bInv)))


lbfgsV
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Int
  -> (Vector a -> (a, Vector a))
  -> Vector a -> [Vector a]
lbfgsV m f x0 = go Seq.empty alpha0 (x0, o0, g0)
  where
    (o0, g0) = f x0

    alpha0 :: a
    alpha0 = realToFrac $ 1 / norm_2 g0

    epsilon :: Double
    epsilon = 1e-5

    go :: Seq (Vector a, Vector a, a) -> a -> (Vector a, a, Vector a) -> [Vector a]
    go hist alpha_ (x, o, g) = x :
      if converged then
        []
      else
        case err of
          Just e -> error (show e)
          Nothing
            | sy > 0    -> go (Seq.take m ((s,y,rho) <| hist)) 1.0 (x', o', g')
            | otherwise -> error ("curvature condition failed: " ++ show sy)
      where
        converged :: Bool
        converged = norm_2 g / max (norm_2 x) 1 <= epsilon

        p :: Vector a
        p = scale (-1) (f (F.toList hist) g)
          where
            f :: [(Vector a, Vector a, a)] -> Vector a -> Vector a
            f ((s,y,rho) : xs) q = z `add` scale (alpha - beta) s
              where
                alpha = rho * (s <.> q)
                z = f xs (q `add` scale (- alpha) y)
                beta = rho * (y <.> z)
            f [] q =
              case Seq.viewl hist of
                EmptyL -> q
                (s, y, _rho) :< _ -> scale (s <.> y / y <.> y) q

        (err, alpha, (x', o', g')) = LS.lineSearch LS.defaultParams f (x, o, g) p alpha_

        s, y :: Vector a
        s = scale alpha p
        y = g' `add` scale (-1) g
        sy, rho :: a
        sy = s <.> y
        rho = 1 / sy


-- example from https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_BFGS = bfgs f x0
  where
    f :: (Fractional a, Floating a) => [a] -> a
    f [beta1,beta2] = sum [(y - beta1*x / (beta2 + x))**2 | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]


-- example from https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_LBFGS = lbfgs 10 f x0
  where
    f :: (Fractional a, Floating a) => [a] -> a
    f [beta1,beta2] = sum [(y - beta1*x / (beta2 + x))**2 | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]


rosenbrock [x,y] = sq (1 - x) + 100 * sq (y - sq x)
  where
    -- Note that 'sq x = x * x' did not work with Kahn mode.
    -- https://github.com/ekmett/ad/pull/84
    sq x = x ** 2
