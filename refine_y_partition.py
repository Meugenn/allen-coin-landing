#!/usr/bin/env python3
"""
Refine the Y-partition solution with higher precision and deeper analysis.
The initial result found: P=(0.5, 0.4583), cuts to t=0.5, 1.875, 3.125
giving max diameter ≈ 1.0078, beating strips (√10/3 ≈ 1.0541) by ~4.4%.
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly

# ─── Geometry (same as before) ───

CORNERS_T  = [0.0, 1.0, 2.0, 3.0]
CORNERS_XY = [np.array([0.,0.]), np.array([1.,0.]),
              np.array([1.,1.]), np.array([0.,1.])]

def bpt(t):
    t = t % 4.0
    if t <= 1: return np.array([t, 0.])
    if t <= 2: return np.array([1., t-1.])
    if t <= 3: return np.array([3.-t, 1.])
    return np.array([0., 4.-t])

def arc_corners(t0, t1):
    arc_len = (t1 - t0) % 4.0
    if arc_len < 1e-10: return []
    out = []
    for ct, cxy in zip(CORNERS_T, CORNERS_XY):
        off = (ct - t0) % 4.0
        if 1e-10 < off < arc_len - 1e-10:
            out.append((off, cxy.copy()))
    out.sort()
    return [xy for _, xy in out]

def poly_area(v):
    n = len(v)
    a = sum(v[i][0]*v[(i+1)%n][1] - v[(i+1)%n][0]*v[i][1] for i in range(n))
    return abs(a)/2.

def poly_diam(v):
    n = len(v)
    return max(np.sqrt((v[i][0]-v[j][0])**2+(v[i][1]-v[j][1])**2)
               for i in range(n) for j in range(i+1,n))

def is_convex(v):
    n = len(v)
    if n < 3: return False
    sign = 0
    for i in range(n):
        e1 = v[(i+1)%n] - v[i]
        e2 = v[(i+2)%n] - v[(i+1)%n]
        c = e1[0]*e2[1] - e1[1]*e2[0]
        if abs(c) < 1e-12: continue
        s = 1 if c > 0 else -1
        if sign == 0: sign = s
        elif s != sign: return False
    return True

def y_pieces(px, py, t1, t2, t3):
    P = np.array([px, py])
    p = [bpt(t) for t in [t1, t2, t3]]
    pc1 = [P, p[0]] + arc_corners(t1, t2) + [p[1]]
    pc2 = [P, p[1]] + arc_corners(t2, t3) + [p[2]]
    pc3 = [P, p[2]] + arc_corners(t3, t1) + [p[0]]
    return pc1, pc2, pc3

# ─── High-precision objective ───

def y_obj_precise(x):
    px, py = x[0], x[1]
    if not (0.001 < px < 0.999 and 0.001 < py < 0.999):
        return 1000.0
    ts = sorted([x[2] % 4, x[3] % 4, x[4] % 4])
    gaps = [ts[1]-ts[0], ts[2]-ts[1], ts[0]+4-ts[2]]
    if any(g < 0.01 for g in gaps):
        return 1000.0
    try:
        pcs = y_pieces(px, py, *ts)
        areas = [poly_area(p) for p in pcs]
        if any(a < 0.01 for a in areas): return 1000.0
        if not all(is_convex(p) for p in pcs): return 1000.0
        pen = sum((a - 1.0/3)**2 for a in areas) * 100000
        return max(poly_diam(p) for p in pcs) + pen
    except:
        return 1000.0

# ─── Constrained formulation using symmetry ───
# By symmetry about x=0.5, optimal Y-partition has:
#   P = (0.5, py)
#   t1 = 0.5  (midpoint of bottom edge)
#   t2 = 2 - a  (on right edge at (1, 1-a))... actually let me figure out the symmetry.
#
# From the initial result: P=(0.5, 0.4583), t=(0.5, 1.875, 3.125)
# t=0.5 → (0.5, 0) bottom center
# t=1.875 → (1.0, 0.875) right edge
# t=3.125 → (0.0, 0.875) left edge
# The two side points are symmetric about x=0.5!
#
# So the symmetric parameterization is: py and a, where:
#   P = (0.5, py)
#   t1 = 0.5  (bottom center)
#   t2 = 1 + a  (right edge at y=a)
#   t3 = 3 + (1-a)  = 4-a  (left edge at y=a)

def y_sym_pieces(py, a):
    """Symmetric Y-partition: P=(0.5, py), cuts to (0.5,0), (1,a), (0,a)."""
    P = np.array([0.5, py])
    p1 = np.array([0.5, 0.0])     # bottom center
    p2 = np.array([1.0, a])       # right edge
    p3 = np.array([0.0, a])       # left edge

    # Determine which corners are in each piece based on value of a
    # Piece 1 (right): P → (0.5,0) → arc to (1,a) → back to P
    # The arc from t1=0.5 to t2=1+a goes through corner (1,0) at t=1 if a>0
    # Piece 2 (top): P → (1,a) → arc to (0,a) → back to P
    # The arc from t2=1+a to t3=4-a goes through corners at t=2 and t=3
    # if 1+a < 2 (i.e. a<1) we include (1,1), and if 4-a > 3 (i.e. a<1) we include (0,1)
    # Piece 3 (left): P → (0,a) → arc to (0.5,0) → back to P
    # The arc from t3=4-a to t1=0.5(+4)=4.5 goes through corner (0,0) at t=0(=4)

    t1 = 0.5
    t2 = 1.0 + a
    t3 = 4.0 - a  # = 0 + (4-a), on left edge at y=a

    pc1 = [P, p1] + arc_corners(t1, t2) + [p2]
    pc2 = [P, p2] + arc_corners(t2, t3) + [p3]
    pc3 = [P, p3] + arc_corners(t3, t1) + [p1]
    return pc1, pc2, pc3

def y_sym_obj(x):
    py, a = x
    if not (0.01 < py < 0.99 and 0.01 < a < 0.99):
        return 1000.0
    try:
        pcs = y_sym_pieces(py, a)
        areas = [poly_area(p) for p in pcs]
        if any(ar < 0.01 for ar in areas): return 1000.0
        if not all(is_convex(p) for p in pcs): return 1000.0
        pen = sum((ar - 1.0/3)**2 for ar in areas) * 100000
        return max(poly_diam(p) for p in pcs) + pen
    except:
        return 1000.0


def main():
    print("=" * 65)
    print("  REFINING Y-PARTITION SOLUTION")
    print("=" * 65)

    # ── Phase 1: High-res global search (general Y) ──
    print("\n─── Phase 1: General Y-partition, high-res DE ───")
    best_gen = None
    for seed in range(20):
        res = differential_evolution(
            y_obj_precise,
            [(0.1, 0.9), (0.1, 0.9), (0, 4), (0, 4), (0, 4)],
            maxiter=3000, seed=seed, tol=1e-14,
            polish=True, popsize=40,
            mutation=(0.5, 1.5), recombination=0.9
        )
        if res.fun < 50 and (best_gen is None or res.fun < best_gen.fun):
            best_gen = res

    if best_gen:
        px, py = best_gen.x[:2]
        ts = sorted([best_gen.x[i] % 4 for i in [2,3,4]])
        pcs = y_pieces(px, py, *ts)
        areas = [poly_area(p) for p in pcs]
        diams = [poly_diam(p) for p in pcs]
        print(f"  P = ({px:.8f}, {py:.8f})")
        print(f"  t = ({ts[0]:.8f}, {ts[1]:.8f}, {ts[2]:.8f})")
        print(f"  Boundary points:")
        for i, t in enumerate(ts):
            bp = bpt(t)
            print(f"    p{i+1} = ({bp[0]:.8f}, {bp[1]:.8f})  [t={t:.8f}]")
        print(f"  Areas:     {[f'{a:.8f}' for a in areas]}")
        print(f"  Diameters: {[f'{d:.8f}' for d in diams]}")
        print(f"  Max diam:  {max(diams):.8f}")

    # ── Phase 2: Symmetric Y-partition (2D search — very precise) ──
    print("\n─── Phase 2: Symmetric Y-partition (2 parameters) ───")
    best_sym = None
    for seed in range(50):
        res = differential_evolution(
            y_sym_obj,
            [(0.1, 0.9), (0.1, 0.99)],
            maxiter=5000, seed=seed, tol=1e-15,
            polish=True, popsize=50,
            mutation=(0.5, 1.5), recombination=0.9
        )
        if res.fun < 50 and (best_sym is None or res.fun < best_sym.fun):
            best_sym = res

    if best_sym:
        py, a = best_sym.x
        pcs = y_sym_pieces(py, a)
        areas = [poly_area(p) for p in pcs]
        diams = [poly_diam(p) for p in pcs]
        print(f"  py = {py:.10f}")
        print(f"  a  = {a:.10f}")
        print(f"  P = (0.5, {py:.10f})")
        print(f"  Boundary points: (0.5, 0), (1, {a:.10f}), (0, {a:.10f})")
        print(f"  Areas:     {[f'{ar:.10f}' for ar in areas]}")
        print(f"  Diameters: {[f'{d:.10f}' for d in diams]}")
        print(f"  Max diam:  {max(diams):.10f}")
        print(f"  √10/3:    {np.sqrt(10)/3:.10f}")
        print(f"  Improvement over strips: {(np.sqrt(10)/3 - max(diams))/(np.sqrt(10)/3)*100:.6f}%")

    # ── Phase 3: Analytical check ──
    print("\n─── Phase 3: Analytical structure ───")
    if best_sym:
        py, a = best_sym.x
        pcs = y_sym_pieces(py, a)
        print("\n  Piece details:")
        for i, pc in enumerate(pcs):
            print(f"\n  Piece {i+1} ({len(pc)} vertices):")
            for j, v in enumerate(pc):
                print(f"    v{j}: ({v[0]:.8f}, {v[1]:.8f})")
            print(f"    Area: {poly_area(pc):.8f}")
            print(f"    Diam: {poly_diam(pc):.8f}")
            print(f"    Convex: {is_convex(pc)}")

            # Find the pair achieving the diameter
            n = len(pc)
            for ii in range(n):
                for jj in range(ii+1, n):
                    d = np.sqrt((pc[ii][0]-pc[jj][0])**2 + (pc[ii][1]-pc[jj][1])**2)
                    if abs(d - poly_diam(pc)) < 1e-8:
                        print(f"    Diam achieved by: v{ii}↔v{jj}")

    # ── Phase 4: Fine grid search for verification ──
    print("\n─── Phase 4: Grid search verification ───")
    best_grid = (None, 1e10)
    for py_i in np.linspace(0.3, 0.6, 300):
        for a_i in np.linspace(0.7, 0.99, 300):
            val = y_sym_obj([py_i, a_i])
            if val < best_grid[1]:
                best_grid = ((py_i, a_i), val)

    py_g, a_g = best_grid[0]
    pcs_g = y_sym_pieces(py_g, a_g)
    areas_g = [poly_area(p) for p in pcs_g]
    diams_g = [poly_diam(p) for p in pcs_g]
    print(f"  Grid best: py={py_g:.6f}, a={a_g:.6f}")
    print(f"  Max diam:  {max(diams_g):.8f}")
    print(f"  Areas:     {[f'{ar:.6f}' for ar in areas_g]}")

    # ── Visualization ──
    if best_sym:
        py, a = best_sym.x
        pcs = y_sym_pieces(py, a)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

        # Left: the partition
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        labels = ['Right (bottom-right)', 'Top', 'Left (bottom-left)']
        for i, (pc, col, lab) in enumerate(zip(pcs, colors, labels)):
            poly = MplPoly(pc, closed=True, facecolor=col,
                           edgecolor='black', lw=2, alpha=0.7)
            ax1.add_patch(poly)
            cx = np.mean([v[0] for v in pc])
            cy = np.mean([v[1] for v in pc])
            ax1.text(cx, cy,
                     f'{lab}\nd={poly_diam(pc):.4f}\nA={poly_area(pc):.4f}',
                     ha='center', va='center', fontsize=9, fontweight='bold')

        # Mark center point and cuts
        P = np.array([0.5, py])
        ax1.plot(*P, 'ko', ms=8, zorder=5)
        ax1.annotate(f'P=({P[0]:.3f},{P[1]:.3f})', P, textcoords='offset points',
                     xytext=(10, -15), fontsize=9)
        for t_val in [0.5, 1+a, 4-a]:
            bp = bpt(t_val)
            ax1.plot([P[0], bp[0]], [P[1], bp[1]], 'k-', lw=1.5, zorder=4)
            ax1.plot(*bp, 'rs', ms=6, zorder=5)

        ax1.set_xlim(-0.08, 1.08)
        ax1.set_ylim(-0.08, 1.08)
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        ax1.set_title(f'Optimal Y-partition\nmax diameter = {max(poly_diam(p) for p in pcs):.6f}\n'
                       f'(strips: {np.sqrt(10)/3:.6f}, improvement: '
                       f'{(np.sqrt(10)/3 - max(poly_diam(p) for p in pcs))/np.sqrt(10)/3*100:.2f}%)',
                       fontsize=12)

        # Right: heatmap of objective in (py, a) space
        py_range = np.linspace(0.2, 0.7, 200)
        a_range = np.linspace(0.5, 0.99, 200)
        Z = np.full((len(a_range), len(py_range)), np.nan)
        for i, ai in enumerate(a_range):
            for j, pyj in enumerate(py_range):
                val = y_sym_obj([pyj, ai])
                if val < 50:
                    Z[i, j] = val

        im = ax2.imshow(Z, extent=[py_range[0], py_range[-1], a_range[0], a_range[-1]],
                        aspect='auto', origin='lower', cmap='viridis_r',
                        vmin=0.95, vmax=1.2)
        ax2.plot(py, a, 'r*', ms=15, zorder=5, label=f'Optimum ({py:.4f}, {a:.4f})')
        ax2.set_xlabel('py (center point height)')
        ax2.set_ylabel('a (side boundary point height)')
        ax2.set_title('Objective landscape (symmetric Y-partition)')
        ax2.legend(fontsize=10)
        plt.colorbar(im, ax=ax2, label='max diameter + area penalty')

        plt.tight_layout()
        plt.savefig('square_partition_refined.png', dpi=150, bbox_inches='tight')
        print("\n  Saved refined visualization → square_partition_refined.png")


if __name__ == '__main__':
    main()
