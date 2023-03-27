{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module QuadProg
  ( QuadProg (..)
  , evalQuadProgObj
  , evalQuadProg
  , solveQuadProg

  , test_unconstrained_qp
  , test_unconstrained_qp_2
  , test_no_minimum
  , test_nonbinding_constraints
  , test_some_binding_constraints
  , test_some_binding_constraints_2
  ) where

import Control.Monad
import Control.Exception (assert)
import Data.IntSet (IntSet)
import qualified Data.IntSet as IntSet
import qualified Data.Vector.Generic as VG
import Numeric.LinearAlgebra


-- | Quadratic Programming problem: $min \{\frac{1}{2} x^T Q x + c^T x \mid A x \le b\}$
data QuadProg a =
  QuadProg
    !(Matrix a) -- ^ $Q \in R^{n\times n}$
    !(Vector a) -- ^ $c \in R^n$
    !(Matrix a) -- ^ $A \in R^{m\times n}$
    !(Vector a) -- ^ $b \in R^m$
  deriving (Show)


evalQuadProgObj
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => QuadProg a -> Vector a -> a
evalQuadProgObj (QuadProg qs c _ _) x = (x <.> (qs #> x)) / 2 + (c <.> x)


evalQuadProg
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => QuadProg a -> Vector a -> a -> Maybe a
evalQuadProg qp@(QuadProg _ _ as b) x tol = do
  guard $ VG.all (>= - tol) (b `sub` (as #> x))
  return $ evalQuadProgObj qp x


-- | Solve Quadratic Preogramming (QP) problem using Active Set Method
--
-- http://www.fujilab.dnj.ynu.ac.jp/lecture/system5.pdf
solveQuadProg
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => QuadProg a -- ^ QP problem
  -> Vector a   -- ^ initial solution
  -> [Vector a]
solveQuadProg (QuadProg qs c as b) x0
  | not (and [size qs == (n,n), size c == n, size as == (m,n), size b == m, size x0 == n]) = error "dimention mismatch"
  -- TODO: check positive-definiteness of Q
  | VG.any (< -tol) slack0 = error "infeasible initial solution"
  | otherwise = go ws0 x0
  where
    tol = 1e-8

    n = VG.length x0
    m = rows as
    slack0 = b `sub` (as #> x0)
    wsAll = IntSet.fromList [0 .. m-1]
    ws0 = IntSet.fromList [i | i <- [0..m-1], slack0 VG.! i < tol]

    go :: IntSet -> Vector a -> [Vector a]
    go !ws !x = assert (size as' == (m',n)) $ assert (size x' == n) $ assert (size y' == m') $ (x :) $
      if VG.all (\z -> abs z < tol) (x' `sub` x) then
        if VG.all (> - tol) y' then
          -- converged
          let _y :: Vector a
              _y = VG.replicate m 0 VG.// zip (IntSet.toAscList ws) (VG.toList y')
           in []
        else
          go (IntSet.delete (wsV VG.! VG.minIndex y') ws) x
      else
        let d = x' `sub` x
            alphas =
              [ (i, alpha)
              | i <- IntSet.toList (wsAll IntSet.\\ ws)
              , let as_i = as ! i
              , let v_i = as_i <.> d
              , v_i > 0
              , let alpha = ((b VG.! i) - as_i <.> x) / v_i
              ]
            alpha = minimum (1 : map snd alphas)
            x''  = x `add` scale alpha d
            ws'' = ws `IntSet.union` IntSet.fromList [i | (i, alpha') <- alphas, alpha' <= alpha]
         in go ws'' x''
      where
        m' = IntSet.size ws
        wsV :: Vector Int
        wsV = VG.fromListN m' (IntSet.toAscList ws)
        as' = as ? IntSet.toAscList ws
        b' = VG.fromListN m' [b VG.! i | i <- IntSet.toAscList ws]
        (x',y') = VG.splitAt n (mat <\> (scale (-1) c VG.++ b'))
          where
            mat =
              qs ||| tr' as'
              ===
              as' ||| konst 0 (m',m')


sub :: (Additive (c t), Linear t c, Num t) => c t -> c t -> c t
sub x y = x `add` scale (-1) y


-- https://www.fsb.miamioh.edu/lij14/400_slide_qp.pdf
-- Example 1: unconstrained QP
test_unconstrained_qp = (x, evalQuadProgObj prob x)  -- ((4,2), -32)
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,0,0,8])
             (VG.fromList [-8,-16])
             ((1 >< 2) [0,0])
             (VG.fromList [0])
    x0 = VG.fromList [0,0]
    x = last $ solveQuadProg prob x0


test_unconstrained_qp_2 = (x, evalQuadProgObj prob x)  -- ((4,2), -32)
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,0,0,8])
             (VG.fromList [-8,-16])
             ((0 >< 2) [])
             (VG.fromList [])
    x0 = VG.fromList [0,0]
    x = last $ solveQuadProg prob x0


-- Example 2: a QP problem that has no minimum
test_no_minimum = (x, evalQuadProgObj prob x)  -- should be error, but ...
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,4,4,8])
             (VG.fromList [-8,-16])
             ((1 >< 2) [0,0])
             (VG.fromList [0])
    x0 = VG.fromList [0,0]
    x = last $ solveQuadProg prob x0


-- Example 3: Constrained QP with non-binding constraints
test_nonbinding_constraints = (x, evalQuadProgObj prob x) -- ((4,2), -32)
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,0,0,8])
             (VG.fromList [-8,-16])
             ((2 >< 2) [-1,-1,-1,0])
             (VG.fromList [-5,-3])
    x0 = VG.fromList [10,10]
    x = last $ solveQuadProg prob x0


-- Example 4: Some constraints are binding
test_some_binding_constraints = (x, evalQuadProgObj prob x)  -- ((4.5, 2), -31.75)
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,0,0,8])
             (VG.fromList [-8,-16])
             ((2 >< 2) [-1,-1,-1,0])
             (VG.fromList [-5, -4.5])
    x0 = VG.fromList [10, 10]
    x = last $ solveQuadProg prob x0


-- Example 5: Some constraints are binding
test_some_binding_constraints_2 = (x, evalQuadProgObj prob x)  -- ((4.8, 2.2), -31.2)
  where
    prob :: QuadProg Double
    prob = QuadProg
             ((2 >< 2) [2,0,0,8])
             (VG.fromList [-8,-16])
             ((2 >< 2) [-1,-1,-1,0])
             (VG.fromList [-7, -3])
    x0 = VG.fromList [10, 10]
    x = last $ solveQuadProg prob x0
