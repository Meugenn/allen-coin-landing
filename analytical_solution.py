#!/usr/bin/env python3
"""
Derive the exact analytical solution for the optimal Y-partition.

Numerical result: P=(1/2, 11/24), cuts to (1/2,0), (1,7/8), (0,7/8).
Let's verify these are exact rational values and derive the diameter.
"""
from fractions import Fraction
import numpy as np

print("=" * 65)
print("  EXACT ANALYTICAL SOLUTION")
print("=" * 65)

# ── Check if numerical values are simple fractions ──
py_num = 0.4583333334
a_num = 0.8750000000

py_frac = Fraction(py_num).limit_denominator(100)
a_frac = Fraction(a_num).limit_denominator(100)
print(f"\n  py ≈ {py_num} = {py_frac} = {float(py_frac)}")
print(f"  a  ≈ {a_num} = {a_frac} = {float(a_frac)}")

# py = 11/24, a = 7/8
py = Fraction(11, 24)
a = Fraction(7, 8)
print(f"\n  Exact: py = {py} = {float(py):.10f}")
print(f"  Exact: a  = {a} = {float(a):.10f}")

# ── The partition ──
print(f"\n  Center point P = (1/2, {py})")
print(f"  Cut 1: P → (1/2, 0)    [bottom center]")
print(f"  Cut 2: P → (1, {a})    [right edge]")
print(f"  Cut 3: P → (0, {a})    [left edge]")

# ── Verify areas ──
# Piece 1 (right): (1/2, 11/24), (1/2, 0), (1, 0), (1, 7/8)
# Shoelace:
p1_verts = [(Fraction(1,2), Fraction(11,24)),
            (Fraction(1,2), Fraction(0)),
            (Fraction(1,1), Fraction(0)),
            (Fraction(1,1), Fraction(7,8))]

def exact_area(verts):
    n = len(verts)
    s = Fraction(0)
    for i in range(n):
        j = (i+1) % n
        s += verts[i][0] * verts[j][1] - verts[j][0] * verts[i][1]
    return abs(s) / 2

# Piece 1 (right): P, (1/2,0), (1,0), (1,7/8)
A1 = exact_area(p1_verts)

# Piece 2 (top): P, (1,7/8), (1,1), (0,1), (0,7/8)
p2_verts = [(Fraction(1,2), Fraction(11,24)),
            (Fraction(1,1), Fraction(7,8)),
            (Fraction(1,1), Fraction(1,1)),
            (Fraction(0,1), Fraction(1,1)),
            (Fraction(0,1), Fraction(7,8))]
A2 = exact_area(p2_verts)

# Piece 3 (left): P, (0,7/8), (0,0), (1/2,0)
p3_verts = [(Fraction(1,2), Fraction(11,24)),
            (Fraction(0,1), Fraction(7,8)),
            (Fraction(0,1), Fraction(0,1)),
            (Fraction(1,2), Fraction(0,1))]
A3 = exact_area(p3_verts)

print(f"\n  Areas (exact fractions):")
print(f"    Piece 1 (right): {A1} = {float(A1):.10f}")
print(f"    Piece 2 (top):   {A2} = {float(A2):.10f}")
print(f"    Piece 3 (left):  {A3} = {float(A3):.10f}")
print(f"    Sum: {A1+A2+A3}")
assert A1 == Fraction(1,3) and A2 == Fraction(1,3) and A3 == Fraction(1,3)
print("    ✓ All areas exactly 1/3")

# ── Compute exact diameters ──
# Piece 1: diameter between (1/2, 0) and (1, 7/8)
d1_sq = (Fraction(1,1) - Fraction(1,2))**2 + (Fraction(7,8))**2
d1_sq_simplified = d1_sq
print(f"\n  Diameter² of piece 1: (1-1/2)² + (7/8-0)² = {d1_sq} = {float(d1_sq):.10f}")

# Also check (1/2,11/24) to (1,0)
d1b_sq = (Fraction(1,2))**2 + (Fraction(11,24))**2
print(f"  Also: (1/2,11/24)↔(1,0): d² = {d1b_sq} = {float(d1b_sq):.10f}")

# And (1/2,11/24) to (1,7/8)
d1c_sq = (Fraction(1,2))**2 + (Fraction(7,8) - Fraction(11,24))**2
print(f"  Also: P↔(1,7/8): d² = {d1c_sq} = {float(d1c_sq):.10f}")

# Piece 2: max distance among (1/2,11/24), (1,7/8), (1,1), (0,1), (0,7/8)
# Key pairs: (1,7/8)↔(0,1) and (0,7/8)↔(1,1) — both should be the same by symmetry
d2a_sq = Fraction(1)**2 + (Fraction(1) - Fraction(7,8))**2
d2b_sq = Fraction(1)**2 + (Fraction(1) - Fraction(7,8))**2
# Also (1,7/8)↔(0,7/8) = 1
d2c_sq = Fraction(1)**2
# And (0,1)↔(1,7/8)
d2d_sq = Fraction(1)**2 + (Fraction(1) - Fraction(7,8))**2
print(f"\n  Piece 2 distances²:")
print(f"    (1,7/8)↔(0,1): {d2a_sq} = {float(d2a_sq):.10f}")
print(f"    (0,7/8)↔(1,1): {d2b_sq} = {float(d2b_sq):.10f}")
print(f"    (1,7/8)↔(0,7/8): {d2c_sq} = {float(d2c_sq)}")
print(f"    P↔(1,1): {(Fraction(1,2)**2 + (1-Fraction(11,24))**2)} = {float(Fraction(1,2)**2 + (1-Fraction(11,24))**2):.10f}")
print(f"    P↔(0,1): same by symmetry")

# The diameter of all three pieces
print(f"\n  All piece diameters:")
d_exact_sq = d1_sq  # = 1/4 + 49/64
print(f"    d² = 1/4 + 49/64 = {Fraction(1,4) + Fraction(49,64)} = {d1_sq}")
print(f"    d² = {d1_sq.numerator}/{d1_sq.denominator}")
print(f"    d  = √({d1_sq.numerator}/{d1_sq.denominator}) = √{d1_sq.numerator}/√{d1_sq.denominator}")
print(f"       = √{d1_sq.numerator}/{int(d1_sq.denominator**0.5)}")
print(f"       = {float(d1_sq)**0.5:.10f}")

# Verify piece 2 diameter
p2_pairs = []
for i in range(len(p2_verts)):
    for j in range(i+1, len(p2_verts)):
        dx = p2_verts[i][0] - p2_verts[j][0]
        dy = p2_verts[i][1] - p2_verts[j][1]
        d_sq = dx**2 + dy**2
        p2_pairs.append((i, j, d_sq))
max_p2 = max(p2_pairs, key=lambda x: x[2])
print(f"\n  Piece 2 max distance: v{max_p2[0]}↔v{max_p2[1]}, d²={max_p2[2]} = {float(max_p2[2]):.10f}")

# ── Final summary ──
d_opt = float(d1_sq)**0.5
d_strips = float(Fraction(10,9))**0.5

print(f"\n{'='*65}")
print(f"  FINAL RESULT")
print(f"{'='*65}")
print(f"\n  Optimal Y-partition of unit square into 3 equal-area convex pieces:")
print(f"    Center: P = (1/2, 11/24)")
print(f"    Cuts to: (1/2, 0), (1, 7/8), (0, 7/8)")
print(f"")
print(f"    Piece 1 (right):  (1/2, 11/24) → (1/2, 0) → (1, 0) → (1, 7/8)")
print(f"    Piece 2 (top):    (1/2, 11/24) → (1, 7/8) → (1, 1) → (0, 1) → (0, 7/8)")
print(f"    Piece 3 (left):   (1/2, 11/24) → (0, 7/8) → (0, 0) → (1/2, 0)")
print(f"")
print(f"    All areas = 1/3 ✓")
print(f"    All diameters equal = √({d1_sq.numerator}/{d1_sq.denominator}) ≈ {d_opt:.10f}")
print(f"")
print(f"    Strips baseline: √(10/9) ≈ {d_strips:.10f}")
print(f"    Improvement: {(d_strips - d_opt)/d_strips*100:.4f}%")
print(f"    Ratio: {d_opt/d_strips:.10f}")

# Cross-check: is √(65/64) correct?
print(f"\n  Cross-check: 65/64 = {Fraction(65,64)} = {float(Fraction(65,64)):.10f}")
print(f"  √(65/64) = √65 / 8 = {65**0.5/8:.10f}")
print(f"  Our d   = {d_opt:.10f}")
print(f"  Match: {abs(d_opt - 65**0.5/8) < 1e-10}")
