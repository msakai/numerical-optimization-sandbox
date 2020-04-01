{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module SecondOrder
  ( newtonMethod
  , gaussNewton
  , levenbergMarquardt
  ) where

import qualified Data.Foldable as F
import Data.Reflection (Reifies)
import qualified Data.Traversable as T
import qualified Data.Vector.Generic as VG
import Foreign.Storable

import Numeric.AD
import Numeric.AD.Internal.Reverse (Reverse, Tape)
import Numeric.AD.Rank1.Sparse (Sparse)

import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra as HMat


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
  :: forall f a. (Traversable f, Storable a, Field a)
  => (forall s. f (AD s (Sparse a)) -> AD s (Sparse a))
  -> f a -> [f a]
newtonMethod f x0 = go x0
  where
    n = length x0

    go :: f a -> [f a]
    go x = x : go (zipWithTV (-) x (h <\> g))
      where
        (_y, gh) = hessian' f x

        g :: Vector a
        g = fromList $ map fst  $ F.toList gh

        h :: Matrix a
        h = (n >< n) $ concat $ map (F.toList . snd) $ F.toList gh


gaussNewton
  :: forall f g a. (Traversable f, Traversable g, Storable a, Field a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> g (Reverse s a))
  -> f a -> [f a]
gaussNewton f x0 = go x0
  where
    m = length x0

    go :: f a -> [f a]
    go x = x : go (zipWithTV (-) x (pinv j #> r))
      where
        rj = jacobian' f x

        r :: Vector a
        r = fromList $ map fst $ F.toList rj

        j :: Matrix a
        j = (length rj >< m) $ concat $ map (F.toList . snd) (F.toList rj)


-- https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_gaussNewton :: [[Double]]
test_gaussNewton = gaussNewton f (fmap auto x0)
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
  :: forall f g a. (Traversable f, Traversable g, Storable a, Field a, Ord a, Show a, Show (f a))
  => a
  -> (forall s. Reifies s Tape => f (Reverse s a) -> g (Reverse s a))
  -> f a -> [f a]
levenbergMarquardt lambda0 f x0 = go lambda0 x0
  where
    m = length x0

    go :: a -> f a -> [f a]
    go lambda x = x : go lambda' x'
      where
        rj = jacobian' f x
        n = length rj

        r :: Vector a
        r = fromList $ map fst $ F.toList rj

        j :: Matrix a
        j = (n >< m) $ concat $ map (F.toList . snd) (F.toList rj)

        gnMat :: Matrix a
        gnMat = (tr j HMat.<> j) `add` diag (VG.replicate m lambda)

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


test_levenbergMarquardt :: [[Double]]
test_levenbergMarquardt = levenbergMarquardt 1.0 f (fmap auto x0)
  where
    f :: Fractional a => [a] -> [a]
    f [beta1,beta2] = [y - beta1*x / (beta2 + x) | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]
