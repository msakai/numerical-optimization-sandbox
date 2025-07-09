module Projection where

import qualified Data.Vector.Storable as V
import Numeric.LinearAlgebra ((<.>), dot, norm_2, scale)


class Member a where
  member :: V.Vector Double -> a -> Bool

notMember :: Member a => V.Vector Double -> a -> Bool
notMember a c = not $ a `member` c

class Member a => Projectable a where
  project :: a -> V.Vector Double -> V.Vector Double


data Box = Box !(V.Vector Double) !(V.Vector Double)
  deriving (Eq, Show)

instance Member Box where
  member xs (Box lb ub) = V.and $ V.zipWith3 f lb ub xs
    where
      f l u x = l <= x && x <= u

instance Projectable Box where
  project (Box lb ub) x = V.zipWith3 f lb ub x
    where
      f l u = max l . min u


data HalfSpace = HalfSpace !(V.Vector Double) Double
  deriving (Eq, Show)

instance Member HalfSpace where
  member x (HalfSpace a b) = a `dot` x >= b

instance Projectable HalfSpace where
  project c@(HalfSpace a b) x
    | x `member` c = x
    | otherwise = x - scale ((a <.> x - b) / (norm_2 a) ** 2) a


newtype Ball = Ball Double
  deriving (Eq, Show)

instance Member Ball where
  member x (Ball r) = norm_2 x <= r

instance Projectable Ball where
  project (Ball r) x
    | d <= r = x
    | otherwise = scale (r / d) x
    where
      d = norm_2 x


_testHalfSpace = (x `notMember` c, y `member` c)
  where
    c = HalfSpace (V.fromList [-2, -3]) (-5)
    x = V.fromList [2, 2]
    y = project c x

_testBall = (x `notMember` c, y `member` c)
  where
    c = Ball 2
    x = V.fromList [2, 2]
    y = project c x
