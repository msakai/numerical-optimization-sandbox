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

import Control.Exception (assert)
import qualified Data.Foldable as F
import Data.IntSet (IntSet)
import qualified Data.IntSet as IntSet
import Data.List (sortBy)
import Data.Ord
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
import qualified Numeric.LinearAlgebra as LA
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
            | sy > 0    -> go (updateBFGSHessianInv s y sy bInv) 1.0 (x', o', g')
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


updateBFGSHessian
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Vector a -> Vector a -> a -> Matrix a -> Matrix a
updateBFGSHessian s y sy b =
  b
  `add` scale (1 / sy) (y `outer` y)
  `add` scale (-1 / (s <.> bs)) (bs `outer` bs)
  where
    bs = b #> s


updateBFGSHessianInv
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Vector a -> Vector a -> a -> Matrix a -> Matrix a
updateBFGSHessianInv s y sy bInv =
  bInv
  `add` scale ((sy + y <.> (bInv #> y)) / sy**2) (s `outer` s)
  `add` scale (-1 / sy) (((bInv #> y) `outer` s) `add` (s `outer` (y <# bInv)))


type LBFGSState a = (Int, Int, Seq (Vector a, Vector a, a))


updateLBFGSState
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Vector a -> Vector a -> a -> LBFGSState a -> LBFGSState a
updateLBFGSState s y sy (n, m, hist) = (n, m, Seq.take m ((s,y,sy) <| hist))


lbfgsTheta
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> a
lbfgsTheta (_n, _m, hist) =
  case Seq.viewl hist of
    EmptyL -> 1
    (_s, y, sy) :< _ -> y <.> y / sy


lbfgsMultiplyHessianInv
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Vector a -> Vector a
lbfgsMultiplyHessianInv state@(_n, _m, hist) g = f (F.toList hist) g
  where
    theta :: a
    theta = lbfgsTheta state

    f :: [(Vector a, Vector a, a)] -> Vector a -> Vector a
    f ((s,y,sy) : xs) q = z `add` scale (alpha - beta) s
      where
        rho = 1 / sy
        alpha = rho * (s <.> q)
        z = f xs (q `add` scale (- alpha) y)
        beta = rho * (y <.> z)
    f [] q = scale (1 / theta) q


lbfgsHessianInv
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Matrix a
lbfgsHessianInv state@(n, _m, hist) = F.foldr f h0 hist
  where
    theta :: a
    theta = lbfgsTheta state

    h0 = scale (1 / theta) (ident n)
    f (s,y,sy) h = updateBFGSHessianInv s y sy h


-- http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf
lbfgsHessian'
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Matrix a
lbfgsHessian' state@(n, _m, hist)
  | Seq.null hist = scale theta (ident n)
  | otherwise =
      assert (LA.size matY == (n,m)) $
      assert (LA.size matS == (n,m)) $
      assert (LA.size matW == (n,2*m)) $
      assert (LA.size matL == (m,m)) $
      assert (LA.size matD == (m,m)) $
      assert (LA.size matM == (2*m,2*m)) $
        scale theta (ident n) `sub` (matW LA.<> matM LA.<> tr matW)
  where
    theta :: a
    theta = lbfgsTheta state

    m :: Int
    m = Seq.length hist

    matY, matS, matW, matL, matD, matM :: Matrix a
    matY = fromColumns [y | (_s,y,_sy) <- F.toList hist]
    matS = fromColumns [s | (s,_y,_sy) <- F.toList hist]
    matW = matY ||| scale theta matS
    matL = (m >< m)
         [ if i > j then s `dot` y else 0
         | (i, (s, _y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         , (j, (_s, y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         ]
    matD = diag $ LA.fromList [sy | (_s, _y, sy) <- F.toList hist]
    matM = inv $ fromBlocks
           [ [scale (-1) matD, tr matL]
           , [matL, scale theta (tr matS LA.<> matS)]
           ]


-- http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf
lbfgsMultiplyHessian'
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Vector a -> Vector a
lbfgsMultiplyHessian' state@(n, _m, hist) x
  | Seq.null hist = scale theta x
  | otherwise =
      assert (LA.size matY == (n,m)) $
      assert (LA.size matS == (n,m)) $
      assert (LA.size matW == (n,2*m)) $
      assert (LA.size matL == (m,m)) $
      assert (LA.size matD == (m,m)) $
      assert (LA.size matMInv == (2*m,2*m)) $
        scale theta x `sub` (matW #> (matMInv <\> (tr matW #> x)))
  where
    theta :: a
    theta = lbfgsTheta state

    m :: Int
    m = Seq.length hist

    matY, matS, matW, matL, matD, matMInv :: Matrix a
    matY = fromColumns [y | (_s,y,_sy) <- F.toList hist]
    matS = fromColumns [s | (s,_y,_sy) <- F.toList hist]
    matW = matY ||| scale theta matS
    matL = (m >< m)
         [ if i > j then s `dot` y else 0
         | (i, (s, _y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         , (j, (_s, y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         ]
    matD = diag $ LA.fromList [sy | (_s, _y, sy) <- F.toList hist]
    matMInv =
      fromBlocks
      [ [scale (-1) matD, tr matL]
      , [matL, scale theta (tr matS LA.<> matS)]
      ]


-- http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf
lbfgsHessianInv'
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Matrix a
lbfgsHessianInv' state@(n, _m, hist)
  | Seq.null hist = scale (1 / theta) (ident n)
  | otherwise =
      assert (LA.size matY == (n,m)) $
      assert (LA.size matS == (n,m)) $
      assert (LA.size matW == (n,2*m)) $
      assert (LA.size matR == (m,m)) $
      assert (LA.size matD == (m,m)) $
      assert (LA.size matM == (2*m,2*m)) $
        scale (1 / theta) (ident n) `add` (matW LA.<> matM LA.<> tr matW)
  where
    theta :: a
    theta = lbfgsTheta state

    m :: Int
    m = Seq.length hist

    matY, matS, matW, matR, matRInv, matD, matM :: Matrix a
    matY = fromColumns [y | (_s,y,_sy) <- F.toList hist]
    matS = fromColumns [s | (s,_y,_sy) <- F.toList hist]
    matW = scale (1 / theta) matY ||| matS
    matR = (m >< m)
         [ if i <= j then s `dot` y else 0
         | (i, (s, _y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         , (j, (_s, y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         ]
    matRInv = inv matR
    matD = diag $ LA.fromList [sy | (_s, _y, sy) <- F.toList hist]
    matM = fromBlocks $
           [ [konst 0 (m,m), scale (-1) matRInv]
           , [scale (-1) (tr matRInv), tr matRInv LA.<> (matD `add` scale (1 / theta) (tr matY LA.<> matY)) LA.<> matRInv]
           ]


-- http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf
lbfgsMultiplyHessianInv'
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => LBFGSState a -> Vector a -> Vector a
lbfgsMultiplyHessianInv' state@(n, _m, hist) g
  | Seq.null hist = scale (1 / theta) g
  | otherwise =
      assert (LA.size matY == (n,m)) $
      assert (LA.size matS == (n,m)) $
      assert (LA.size matR == (m,m)) $
      assert (LA.size vecD == m) $
        scale (1 / theta) g `add` vecWMWtg
  where
    theta :: a
    theta = lbfgsTheta state

    m :: Int
    m = Seq.length hist

    vecYtg = scale (1 / theta) (tr matY #> g)
    vecStg = tr matS #> g
    vecRInvStg = matRInv #> vecStg
    vecMWtg1 = scale (-1) vecRInvStg
    vecMWtg2 = scale (-1) (tr matRInv #> vecYtg)
               `add`
               (tr matRInv #> ((vecD `mul` vecRInvStg) `add` scale (1 / theta) (tr matY #> matY #> vecRInvStg)))
    vecWMWtg = scale (1 / theta) (matY #> vecMWtg1) `add` (matS #> vecMWtg2)

    mul = VG.zipWith (*)

    matY, matS, matR, matRInv :: Matrix a
    matY = fromColumns [y | (_s,y,_sy) <- F.toList hist]
    matS = fromColumns [s | (s,_y,_sy) <- F.toList hist]
    matR = (m >< m)
         [ if i <= j then s `dot` y else 0
         | (i, (s, _y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         , (j, (_s, y, _sy)) <- zip [m-1,m-2..] (F.toList hist)
         ]
    matRInv = inv matR
    vecD :: Vector a
    vecD = LA.fromList [sy | (_s, _y, sy) <- F.toList hist]


lbfgsV
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Int
  -> (Vector a -> (a, Vector a))
  -> Vector a -> [Vector a]
lbfgsV m f x0 = go (n, m, Seq.empty) alpha0 (x0, o0, g0)
  where
    n = VG.length x0
    (o0, g0) = f x0

    alpha0 :: a
    alpha0 = realToFrac $ 1 / norm_2 g0

    epsilon :: Double
    epsilon = 1e-5

    go :: LBFGSState a -> a -> (Vector a, a, Vector a) -> [Vector a]
    go state alpha_ (x, o, g) = x :
      if converged then
        []
      else
        case err of
          Just e -> error (show e)
          Nothing
            | sy > 0    -> go (updateLBFGSState s y sy state) 1.0 (x', o', g')
            | otherwise -> error ("curvature condition failed: " ++ show sy)
      where
        converged :: Bool
        converged = norm_2 g / max (norm_2 x) 1 <= epsilon

        p :: Vector a
        p = scale (-1) (lbfgsMultiplyHessianInv state g)

        (err, alpha, (x', o', g')) = LS.lineSearch LS.defaultParams f (x, o, g) p alpha_

        s, y :: Vector a
        s = scale alpha p
        y = g' `add` scale (-1) g
        sy :: a
        sy = s <.> y


-- | Compute generalized Cauchy point of a function f(x0) + g (x - x0) + (1/2) (x - x0)^T B (x - x0)
generalizedCauchyPoint
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Vector a -- ^ x0
  -> a        -- ^ f(x0)
  -> Vector a -- ^ g
  -> (Vector a -> Vector a) -- ^ @(B '#>')@
  -> Vector a -- ^ lower bounds
  -> Vector a -- ^ upper bounds
  -> (Vector a, IntSet) -- ^ generalized cauchy point and its active set
generalizedCauchyPoint x0 f0 g multiplyB lb ub = go 0 x0 d0 IntSet.empty breakpoints
  where
    breakpoints :: [(a, Int, a)]
    breakpoints =
      sortBy (comparing (\(ti, _xi, _i) -> ti)) $
      concatMap
        (\(i, gi) ->
           case gi `compare` 0 of
             LT -> [(((x0 VG.! i) - (ub VG.! i)) / gi, i, ub VG.! i)]
             GT -> [(((x0 VG.! i) - (lb VG.! i)) / gi, i, lb VG.! i)]
             EQ -> []
        ) $
      zip [(0::Int)..] (VG.toList g)

    d0 :: Vector a
    d0 = scale (-1) g

    {-
    Let Δt := t - tj and x := xj + Δt dj, then:
    f(x0) + g (x - x0) + (1/2) (x - x0)^T B (x - x0)
    = f(x0) + g (xj + Δt dj - x0) + (1/2) (xj + Δt dj - x0)^T B (xj + Δt dj - x0)
    = f(x0) + g (Δt dj + (xj - x0)) + (1/2) (Δt dj + (xj - x0))^T B (Δt dj + (xj - x0))
    = f(x0) + g (Δt dj) + g (xj - x0) + (1/2) ((Δt dj) B (Δt dj) + (Δt dj) B (xj - x0) + (xj - x0)^T B (Δt dj) + (xj - x0)^T B (xj - x0))
    = (f(x0) + g (xj - x0) + (1/2) (xj - x0)^T B (xj - x0))
      + (g dj + dj B (xj - x0) Δt
      + (1/2) (dj B dj) Δt²
    = a0 + a1 Δt + (1/2) a2 Δt²
    -}
    go :: a -> Vector a -> Vector a -> IntSet -> [(a, Int, a)] -> (Vector a, IntSet)
    go tj xj dj as bps =
      case bps of
        [] -> (xj `add` scale dt_opt dj, as)
        (tj', i, val) : bps'
          | tj + dt_opt < tj' -> (xj `add` scale dt_opt dj, as)
          | otherwise -> go tj' (xj VG.// [(i, val)]) (dj VG.// [(i, 0)]) (IntSet.insert i as) bps'
      where
        z = xj `sub` x0
        _a0 = f0 + (g <.> z) + (z <.> multiplyB z) / 2
        a1 = g <.> dj + dj <.> multiplyB z
        a2 = dj <.> multiplyB dj
        dt_opt = - a1 / a2


sub :: (Additive (c t), Linear t c, Num t) => c t -> c t -> c t
sub x y = x `add` scale (-1) y


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
