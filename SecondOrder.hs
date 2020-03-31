{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module SecondOrder
  ( newtonMethod
  , gaussNewtonMethod
  ) where

import qualified Data.Foldable as F
import Data.Reflection (Reifies)
import qualified Data.Traversable as T
import qualified Data.Vector.Generic as VG
import Foreign.Storable

import Numeric.AD
import Numeric.AD.Internal.Reverse (Reverse, Tape)
import Numeric.AD.Rank1.Sparse (Sparse (..))

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
    go x = x : go (zipWithTV (-) x d)
      where
        (_y, gh) = hessian' f x

        g :: Vector a
        g = fromList $ map fst  $ F.toList gh
       
        h :: Matrix a
        h = (n >< n) $ concat $ map (F.toList . snd) $ F.toList gh

        d :: Vector a
        d = h <\> g


gaussNewtonMethod
  :: forall f g a. (Traversable f, Traversable g, Storable a, Field a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> g (Reverse s a))
  -> f a -> [f a]
gaussNewtonMethod f x0 = go x0
  where
    m = length x0

    go :: f a -> [f a]
    go x = x : go (zipWithTV (-) x (pinv j #> r))
      where
        rj = jacobian' f x
        r = fromList $ map fst $ F.toList rj
        j = (length rj >< m) $ concat $ map (F.toList . snd) (F.toList rj)


-- https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
test_gaussNewtonMethod :: [[Double]]
test_gaussNewtonMethod = gaussNewtonMethod f (fmap auto x0)
  where
    f :: Fractional a => [a] -> [a]
    f [beta1,beta2] = [y - beta1*x / (beta2 + x) | (x,y) <- zip xs ys]
      where
        xs = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
        ys = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

    x0 :: [Double]
    x0 = [0.9, 0.2]
