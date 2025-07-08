{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module ProximalGradientMethod where

import Data.Foldable
import Data.Reflection (Reifies)
import qualified Data.Vector.Storable as V
import Numeric.AD
import Numeric.AD.Mode.Reverse
import Numeric.AD.Internal.Reverse (Tape)
import Numeric.LinearAlgebra (norm_2)
import Test.QuickCheck

-- ------------------------------------------------------------------------
-- https://qiita.com/msekino/items/9f217fcd735513627f65

proximalGradientMethod
  :: (Traversable f, Ord a, Fractional a)
  => a
  -> a
  -> (forall s. Reifies s Tape => f (Reverse s a) -> Reverse s a)
  -> (f a -> a, a -> f a -> f a)
  -> f a -> [f a]
proximalGradientMethod eta0 beta f (g, prox) x = map fst $ iterate h (x, eta0)
  where
    f_hat eta x y = fy + sum (zipWith (*) (toList gfy) zs) + (1 / (2*eta)) * sum (map (^(2::Int)) zs)
      where
        (fy, gfy) = grad' f y
        zs = zipWith (-) (toList x) (toList y)

    h (x, eta)
      | fst (grad' f x') <= f_hat eta x' x = (x', eta)
      | otherwise = h (x, eta * beta)
      where
        x' = prox eta (zipWithTF (\xk gk -> xk - eta*gk) x (grad f x))

l2reg :: (Functor f, Foldable f, Fractional a) => a -> (f a -> a, a -> f a -> f a)
l2reg lam = (g, prox)
  where
    g x = (lam / 2) * sum [xk^(2::Int) | xk <- toList x]
    prox eta = fmap (/ (1 + eta * lam))

l1reg :: (Functor f, Foldable f, Ord a, Fractional a) => a -> (f a -> a, a -> f a -> f a)
l1reg lam = (g, prox)
  where
    g x = lam * sum [abs xk | xk <- toList x]
    prox eta = fmap h
      where
        h xk
          | xk >  lam * eta = xk - lam * eta
          | xk < -lam * eta = xk + lam * eta
          | otherwise = 0

-- https://qiita.com/AnchorBlues/items/4e50d3b98a40c8b3086e
glassoReg
  :: (Functor f, Foldable f)
  => (f Double -> [V.Vector Double])
  -> (f Double -> [V.Vector Double] -> f Double)
  -> Double
  -> (f Double -> Double, Double -> f Double -> f Double)
glassoReg group update lam = (g, prox)
  where
    g x = lam * sum [norm_2 xg | xg <- group x]
    prox eta x = update x [if norm_2 xg >= eta * lam then xg - V.map ((eta * lam / norm_2 xg) *) xg else 0 | xg <- group x]

isProximalOperator :: (Foldable f, Show (f a), Ord a, Fractional a, Show a, Arbitrary a) => (f a -> a, a -> f a -> f a)-> Gen (f a) -> Property
isProximalOperator (g, prox) gen = 
  forAll arbitrary $ \(Positive eta) ->
    forAll gen $ \y ->
      let obj x = eta * g x + (1/2) * sum (zipWith (\xk yk -> (xk - yk)^(2::Int)) (toList x) (toList y))
          x_opt = prox eta y
          obj_opt = obj x_opt
       in counterexample (show (x_opt, obj_opt)) $
            forAll gen $ \x -> counterexample (show (obj x)) $ obj_opt <= obj x

prop_l2reg_prox :: Property
prop_l2reg_prox =
  forAll arbitrary $ \(Positive n, Positive (lam :: Rational)) ->
    isProximalOperator (l2reg lam) (vectorOf n arbitrary)

prop_l1reg_prox :: Property
prop_l1reg_prox =
  forAll arbitrary $ \(Positive n, Positive (lam :: Rational)) ->
    isProximalOperator (l1reg lam) (vectorOf n arbitrary)

-- ------------------------------------------------------------------------
-- https://wiki.haskell.org/Foldable_and_Traversable#Generalising_zipWith

data Supply s v = Supply { unSupply :: [s] -> ([s],v) }

instance Functor (Supply s) where 
  fmap f av = Supply (\l -> let (l',v) = unSupply av l in (l',f v))

instance Applicative (Supply s) where
  pure v    = Supply (\l -> (l,v))
  af <*> av = Supply (\l -> let (l',f)  = unSupply af l
                                (l'',v) = unSupply av l'
                            in (l'',f v))

runSupply :: (Supply s v) -> [s] -> v
runSupply av l = snd $ unSupply av l

supply :: Supply s s
supply = Supply (\(x:xs) -> (xs,x))

zipTF :: (Traversable t, Foldable f) => t a -> f b -> t (a,b)
zipTF t f = runSupply (traverse (\a -> (,) a <$> supply) t) (toList f)

zipWithTF :: (Traversable t,Foldable f) => (a -> b -> c) -> t a -> f b -> t c
zipWithTF g t f = runSupply  (traverse (\a -> g a <$> supply) t) (toList f)

zipWithTFM :: (Traversable t,Foldable f,Monad m) => 
              (a -> b -> m c) -> t a -> f b -> m (t c)
zipWithTFM g t f = sequence (zipWithTF g t f)

zipWithTFA :: (Traversable t,Foldable f,Applicative m) => 
              (a -> b -> m c) -> t a -> f b -> m (t c)
zipWithTFA g t f = sequenceA (zipWithTF g t f)

-- ------------------------------------------------------------------------
