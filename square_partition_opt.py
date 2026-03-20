#!/usr/bin/env python3
"""
Computationally optimize: partition a unit square into 3 equal-area convex
pieces, minimizing the maximum diameter of any piece.

Explores three partition topologies:
  A) Nested chords (generalizes parallel strips)
  B) Adjacent chords (corner-cutting)
  C) Y-partition (three cuts from an interior point)
"""

import time
import numpy as np
from scipy.optimize import differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly

# ─────────────────────────────────────────────────
# Geometry utilities
# ─────────────────────────────────────────────────

CORNERS_T  = [0.0, 1.0, 2.0, 3.0]
CORNERS_XY = [np.array([0.,0.]), np.array([1.,0.]),
              np.array([1.,1.]), np.array([0.,1.])]


def bpt(t):
    """Map perimeter parameter t ∈ [0,4) → point on unit square boundary.
    CCW: t=0→(0,0), t=1→(1,0), t=2→(1,1), t=3→(0,1)."""
    t = t % 4.0
    if t <= 1: return np.array([t, 0.])
    if t <= 2: return np.array([1., t - 1.])
    if t <= 3: return np.array([3. - t, 1.])
    return np.array([0., 4. - t])


def arc_corners(t0, t1):
    """Square corners strictly inside the CCW arc from t0 to t1, in order."""
    arc_len = (t1 - t0) % 4.0
    if arc_len < 1e-10:
        return []
    out = []
    for ct, cxy in zip(CORNERS_T, CORNERS_XY):
        off = (ct - t0) % 4.0
        if 1e-10 < off < arc_len - 1e-10:
            out.append((off, cxy.copy()))
    out.sort()
    return [xy for _, xy in out]


def poly_area(verts):
    """Shoelace formula for polygon area."""
    n = len(verts)
    if n < 3:
        return 0.0
    a = 0.0
    for i in range(n):
        j = (i + 1) % n
        a += verts[i][0] * verts[j][1] - verts[j][0] * verts[i][1]
    return abs(a) / 2.0


def poly_diam(verts):
    """Diameter = max pairwise distance among vertices of convex polygon."""
    mx = 0.0
    n = len(verts)
    for i in range(n):
        for j in range(i + 1, n):
            d2 = (verts[i][0] - verts[j][0])**2 + (verts[i][1] - verts[j][1])**2
            if d2 > mx:
                mx = d2
    return np.sqrt(mx)


def is_convex(verts):
    """Check if polygon vertices form a convex polygon."""
    n = len(verts)
    if n < 3:
        return False
    sign = 0
    for i in range(n):
        e1 = verts[(i + 1) % n] - verts[i]
        e2 = verts[(i + 2) % n] - verts[(i + 1) % n]
        cross = e1[0] * e2[1] - e1[1] * e2[0]
        if abs(cross) < 1e-10:
            continue
        s = 1 if cross > 0 else -1
        if sign == 0:
            sign = s
        elif s != sign:
            return False
    return True


# ─────────────────────────────────────────────────
# Topology A: Nested chords
# Points s1 < s2 < s3 < s4 on boundary.
# Chord_outer: s1↔s4,  Chord_inner: s2↔s3
# Piece1 = arc(s4→s1),  Piece2 = between chords,  Piece3 = arc(s2→s3)
# (This generalizes parallel strips)
# ─────────────────────────────────────────────────

def nested_pieces(s1, s2, s3, s4):
    p = [bpt(s) for s in [s1, s2, s3, s4]]
    pc1 = [p[3]] + arc_corners(s4, s1) + [p[0]]
    pc2 = [p[0]] + arc_corners(s1, s2) + [p[1], p[2]] + arc_corners(s3, s4) + [p[3]]
    pc3 = [p[1]] + arc_corners(s2, s3) + [p[2]]
    return pc1, pc2, pc3


def nested_obj(x):
    s = sorted([xi % 4 for xi in x])
    gaps = [s[1]-s[0], s[2]-s[1], s[3]-s[2], s[0]+4-s[3]]
    if any(g < 0.05 for g in gaps):
        return 100.0
    try:
        pcs = nested_pieces(*s)
        areas = [poly_area(p) for p in pcs]
        if any(a < 0.05 for a in areas):
            return 100.0
        if not all(is_convex(p) for p in pcs):
            return 100.0
        pen = sum((a - 1.0/3)**2 for a in areas) * 10000
        return max(poly_diam(p) for p in pcs) + pen
    except Exception:
        return 100.0


# ─────────────────────────────────────────────────
# Topology B: Adjacent chords
# s1 < s2 < s3 < s4. Chord1: s1↔s2, Chord2: s3↔s4.
# Piece1 = arc(s1→s2), Piece2 = middle, Piece3 = arc(s3→s4)
# ─────────────────────────────────────────────────

def adjacent_pieces(s1, s2, s3, s4):
    p = [bpt(s) for s in [s1, s2, s3, s4]]
    pc1 = [p[0]] + arc_corners(s1, s2) + [p[1]]
    pc2 = [p[1]] + arc_corners(s2, s3) + [p[2], p[3]] + arc_corners(s4, s1) + [p[0]]
    pc3 = [p[2]] + arc_corners(s3, s4) + [p[3]]
    return pc1, pc2, pc3


def adjacent_obj(x):
    s = sorted([xi % 4 for xi in x])
    gaps = [s[1]-s[0], s[2]-s[1], s[3]-s[2], s[0]+4-s[3]]
    if any(g < 0.05 for g in gaps):
        return 100.0
    try:
        pcs = adjacent_pieces(*s)
        areas = [poly_area(p) for p in pcs]
        if any(a < 0.05 for a in areas):
            return 100.0
        if not all(is_convex(p) for p in pcs):
            return 100.0
        pen = sum((a - 1.0/3)**2 for a in areas) * 10000
        return max(poly_diam(p) for p in pcs) + pen
    except Exception:
        return 100.0


# ─────────────────────────────────────────────────
# Topology C: Y-partition
# Interior point P=(px,py), three boundary points t1 < t2 < t3
# Three cuts: P→p1, P→p2, P→p3
# ─────────────────────────────────────────────────

def y_pieces(px, py, t1, t2, t3):
    P = np.array([px, py])
    p = [bpt(t) for t in [t1, t2, t3]]
    pc1 = [P, p[0]] + arc_corners(t1, t2) + [p[1]]
    pc2 = [P, p[1]] + arc_corners(t2, t3) + [p[2]]
    pc3 = [P, p[2]] + arc_corners(t3, t1) + [p[0]]
    return pc1, pc2, pc3


def y_obj(x):
    px, py = x[0], x[1]
    if not (0.01 < px < 0.99 and 0.01 < py < 0.99):
        return 100.0
    ts = sorted([x[2] % 4, x[3] % 4, x[4] % 4])
    gaps = [ts[1]-ts[0], ts[2]-ts[1], ts[0]+4-ts[2]]
    if any(g < 0.05 for g in gaps):
        return 100.0
    try:
        pcs = y_pieces(px, py, *ts)
        areas = [poly_area(p) for p in pcs]
        if any(a < 0.05 for a in areas):
            return 100.0
        if not all(is_convex(p) for p in pcs):
            return 100.0
        pen = sum((a - 1.0/3)**2 for a in areas) * 10000
        return max(poly_diam(p) for p in pcs) + pen
    except Exception:
        return 100.0


# ─────────────────────────────────────────────────
# Run optimizations
# ─────────────────────────────────────────────────

def run_topology(name, obj_fn, bounds, seeds):
    """Run differential_evolution with multiple seeds, return best results."""
    best = None
    for seed in seeds:
        t0 = time.time()
        res = differential_evolution(
            obj_fn, bounds,
            maxiter=1500, seed=seed, tol=1e-12,
            polish=True, mutation=(0.5, 1.5),
            recombination=0.9, popsize=25
        )
        dt = time.time() - t0
        if res.fun < 50:
            if best is None or res.fun < best.fun:
                best = res
            print(f"  seed={seed:4d}: obj={res.fun:.6f}  ({dt:.1f}s)")
        else:
            print(f"  seed={seed:4d}: infeasible  ({dt:.1f}s)")
    return best


def main():
    print("=" * 65)
    print("  SQUARE → 3 EQUAL-AREA CONVEX PIECES: MINIMIZE MAX DIAMETER")
    print("=" * 65)

    results = []  # (name, max_diam, pieces)
    seeds = [42, 7, 123, 999, 2024, 314, 271, 161]

    # ── Reference: vertical strips ──
    strips = [
        [np.array([0,0]), np.array([1/3,0]), np.array([1/3,1]), np.array([0,1])],
        [np.array([1/3,0]), np.array([2/3,0]), np.array([2/3,1]), np.array([1/3,1])],
        [np.array([2/3,0]), np.array([1,0]), np.array([1,1]), np.array([2/3,1])],
    ]
    d_ref = max(poly_diam(p) for p in strips)
    print(f"\nReference: parallel strips  →  max diam = {d_ref:.6f}")
    print(f"  (√10/3 = {np.sqrt(10)/3:.6f})")
    results.append(("Parallel strips", d_ref, strips))

    # ── Topology A: Nested chords ──
    print("\n─── Topology A: Nested Chords (generalizes strips) ───")
    res_a = run_topology("Nested", nested_obj, [(0, 4)] * 4, seeds)
    if res_a and res_a.fun < 50:
        s = sorted([xi % 4 for xi in res_a.x])
        pcs = nested_pieces(*s)
        areas = [poly_area(p) for p in pcs]
        if all(abs(a - 1/3) < 0.01 for a in areas):
            md = max(poly_diam(p) for p in pcs)
            results.append(("Nested chords", md, pcs))
            print(f"  ★ Best: max_d={md:.6f}  params=({', '.join(f'{si:.4f}' for si in s)})")

    # ── Topology B: Adjacent chords ──
    print("\n─── Topology B: Adjacent Chords (corner cuts) ───")
    res_b = run_topology("Adjacent", adjacent_obj, [(0, 4)] * 4, seeds)
    if res_b and res_b.fun < 50:
        s = sorted([xi % 4 for xi in res_b.x])
        pcs = adjacent_pieces(*s)
        areas = [poly_area(p) for p in pcs]
        if all(abs(a - 1/3) < 0.01 for a in areas):
            md = max(poly_diam(p) for p in pcs)
            results.append(("Adjacent chords", md, pcs))
            print(f"  ★ Best: max_d={md:.6f}  params=({', '.join(f'{si:.4f}' for si in s)})")

    # ── Topology C: Y-partition ──
    print("\n─── Topology C: Y-Partition (3 cuts from interior point) ───")
    res_c = run_topology("Y-part", y_obj,
                         [(0.05, 0.95), (0.05, 0.95), (0, 4), (0, 4), (0, 4)],
                         seeds)
    if res_c and res_c.fun < 50:
        px, py = res_c.x[0], res_c.x[1]
        ts = sorted([res_c.x[2] % 4, res_c.x[3] % 4, res_c.x[4] % 4])
        pcs = y_pieces(px, py, *ts)
        areas = [poly_area(p) for p in pcs]
        if all(abs(a - 1/3) < 0.01 for a in areas):
            md = max(poly_diam(p) for p in pcs)
            results.append(("Y-partition", md, pcs))
            print(f"  ★ Best: max_d={md:.6f}  P=({px:.4f},{py:.4f})  t=({ts[0]:.4f},{ts[1]:.4f},{ts[2]:.4f})")

    # ── Summary ──
    print("\n" + "=" * 65)
    print("  RESULTS SUMMARY")
    print("=" * 65)
    results.sort(key=lambda r: r[1])
    for name, d, pcs in results:
        areas = [poly_area(p) for p in pcs]
        diams = [poly_diam(p) for p in pcs]
        print(f"  {d:.6f}  {name}")
        for i, (p, a, dm) in enumerate(zip(pcs, areas, diams)):
            print(f"      piece {i+1}: area={a:.5f}  diam={dm:.5f}  verts={len(p)}")

    best_name, best_d, best_pcs = results[0]
    print(f"\n  ★ OPTIMAL: {best_name}  →  max diameter = {best_d:.6f}")
    print(f"    √10/3 = {np.sqrt(10)/3:.6f}")
    print(f"    Improvement over strips: {(d_ref - best_d)/d_ref * 100:.4f}%")

    for i, p in enumerate(best_pcs):
        verts_str = " → ".join(f"({v[0]:.4f},{v[1]:.4f})" for v in p)
        print(f"    Piece {i+1}: {verts_str}")

    # ── Visualization ──
    n_plots = min(len(results), 4)
    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 6))
    if n_plots == 1:
        axes = [axes]
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']

    seen = set()
    idx = 0
    for name, d, pcs in results:
        key = f"{d:.3f}"
        if key in seen:
            continue
        if idx >= n_plots:
            break
        seen.add(key)
        ax = axes[idx]
        for i, (pc, col) in enumerate(zip(pcs, colors)):
            poly = MplPoly(pc, closed=True, facecolor=col,
                           edgecolor='black', lw=2, alpha=0.7)
            ax.add_patch(poly)
            cx = np.mean([v[0] for v in pc])
            cy = np.mean([v[1] for v in pc])
            ax.text(cx, cy,
                    f'd={poly_diam(pc):.3f}\nA={poly_area(pc):.3f}',
                    ha='center', va='center', fontsize=10, fontweight='bold')
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{name}\nmax d = {d:.4f}', fontsize=12)
        idx += 1

    plt.tight_layout()
    plt.savefig('square_partition.png', dpi=150, bbox_inches='tight')
    print("\n  Saved visualization → square_partition.png")


if __name__ == '__main__':
    main()
