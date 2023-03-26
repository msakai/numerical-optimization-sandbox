{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
module QP where

import Control.Exception (assert)
import Data.IntSet (IntSet)
import qualified Data.IntSet as IntSet
import qualified Data.Vector.Generic as VG
import Numeric.LinearAlgebra


-- | Solve $min \{\frac{1}{2} x^T Q x + c^T x \mid A x \le b\}$ using active set method
--
-- http://www.fujilab.dnj.ynu.ac.jp/lecture/system5.pdf
quadprog
  :: forall a. (Field a, Ord a, Normed (Vector a), Show a)
  => Matrix a -- ^ $Q \in R^{n\times n}$
  -> Vector a -- ^ $c \in R^n$
  -> Matrix a -- ^ $A \in R^{m\times n}$
  -> Vector a -- ^ $b \in R^m$
  -> Vector a -- ^ $x_0 \in R^n$
  -> [Vector a]
quadprog qs c as b x0
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

    sub x y = x `add` scale (-1) y

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
        b' = VG.generate m' (b VG.!)
        (x',y') = VG.splitAt n (mat <\> (scale (-1) c VG.++ b'))
          where
            mat =
              qs ||| tr' as'
              ===
              as' ||| konst 0 (m',m')


-- https://www.fsb.miamioh.edu/lij14/400_slide_qp.pdf
-- Example 1: unconstrained QP
test_unconstrained_qp = (x, f x)  -- ((4,2), -32)
  where
    qs = (2 >< 2) [2,0,0,8]
    c = VG.fromList [-8,-16]
    as = (1 >< 2) [0,0]
    b = VG.fromList [0]
    x0 = VG.fromList [0,0]
    f x = (x <.> (qs #> x)) / 2 + (c <.> x)

    x :: Vector Double
    x = last $ quadprog qs c as b x0


-- Example 2: a QP problem that has no minimum
test_no_minimum = (x, f x)  -- should be error, but ...
  where
    qs = (2 >< 2) [2,4,4,8]
    c = VG.fromList [-8,-16]
    as = (1 >< 2) [0,0]
    b = VG.fromList [0]
    x0 = VG.fromList [0,0]
    f x = (x <.> (qs #> x)) / 2 + (c <.> x)

    x :: Vector Double
    x = last $ quadprog qs c as b x0


-- Example 3: Constrained QP with non-binding constraints
test_nonbinding_constraints = (x, f x) -- ((4,2), -32)
  where
    qs = (2 >< 2) [2,0,0,8]
    c = VG.fromList [-8,-16]
    as = (2 >< 2) [-1,-1,-1,0]
    b = VG.fromList [-5,-3]
    x0 = VG.fromList [10,10]
    f x = (x <.> (qs #> x)) / 2 + (c <.> x)

    x :: Vector Double
    x = last $ quadprog qs c as b x0


-- Example 4: Some constraints are binding
test_some_binding_constraints = (x, f x)  -- should be ((4.5, 2), -31.75), but ...
  where
    qs = (2 >< 2) [2,0,0,8]
    c = VG.fromList [-8,-16]
    as = (2 >< 2) [-1,-1,-1,0]
    b = VG.fromList [-5, -4.5]
    x0 = VG.fromList [10, 10]
    f x = (x <.> (qs #> x)) / 2 + (c <.> x)

    x :: Vector Double
    x = last $ quadprog qs c as b x0


-- Example 5: Some constraints are binding
test_some_binding_constraints_2 = (xs, f x)  -- ((4.8, 4.2), -31.2)
  where
    qs = (2 >< 2) [2,0,0,8]
    c = VG.fromList [-8,-16]
    as = (2 >< 2) [-1,-1,-1,0]
    b = VG.fromList [-7, -3]
    x0 = VG.fromList [10, 10]
    f x = (x <.> (qs #> x)) / 2 + (c <.> x)
    xs = quadprog qs c as b x0

    x :: Vector Double
    x = last xs
