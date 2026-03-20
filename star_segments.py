"""
Key insight: the problem says "lines" but in a drawn star figure,
lines are SEGMENTS, not infinite lines. Infinite lines create
"phantom" intersection points outside the star that break face adjacency.
With segments clipped to the star boundary, we avoid these phantom points
and may get 10 triangles.

Also trying non-regular pentagrams where segment clipping is different.
"""

import math
from itertools import combinations
import random

EPS = 1e-9

def make_line(p1, p2):
    a = p2[1] - p1[1]
    b = p1[0] - p2[0]
    c = p2[0]*p1[1] - p1[0]*p2[1]
    return (a, b, c)

def ix_lines(l1, l2):
    d = l1[0]*l2[1] - l2[0]*l1[1]
    if abs(d) < EPS:
        return None
    return ((l1[1]*l2[2] - l2[1]*l1[2]) / d, (l2[0]*l1[2] - l1[0]*l2[2]) / d)

def peq(p1, p2):
    return abs(p1[0]-p2[0]) < EPS and abs(p1[1]-p2[1]) < EPS

def wn(pt, poly):
    x, y = pt
    w = 0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i+1)%n]
        if y1 <= y:
            if y2 > y and (x2-x1)*(y-y1) - (x-x1)*(y2-y1) > 0:
                w += 1
        else:
            if y2 <= y and (x2-x1)*(y-y1) - (x-x1)*(y2-y1) < 0:
                w -= 1
    return w

def point_on_seg(p, seg_start, seg_end):
    """Check if point p is on segment [seg_start, seg_end]."""
    dx = seg_end[0] - seg_start[0]
    dy = seg_end[1] - seg_start[1]
    len_sq = dx*dx + dy*dy
    if len_sq < EPS*EPS:
        return peq(p, seg_start)
    t = ((p[0]-seg_start[0])*dx + (p[1]-seg_start[1])*dy) / len_sq
    return -EPS < t < 1+EPS

def seg_seg_intersect(s1_start, s1_end, s2_start, s2_end):
    """Find intersection point of two segments, or None."""
    l1 = make_line(s1_start, s1_end)
    l2 = make_line(s2_start, s2_end)
    p = ix_lines(l1, l2)
    if p is None:
        return None
    if point_on_seg(p, s1_start, s1_end) and point_on_seg(p, s2_start, s2_end):
        return p
    return None


def count_triangles_with_segments(V, extra_segments):
    """
    Count triangular faces in the planar subdivision created by
    the pentagram SEGMENTS plus additional SEGMENTS.

    Uses the line-triple approach but only considering intersection
    points that are actually within the drawn segments.
    """
    # Pentagram segments
    pent_segs = [(V[i], V[(i+2)%5]) for i in range(5)]
    pent_lines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

    # All segments and their corresponding infinite lines
    all_segs = pent_segs + extra_segments
    all_lines = pent_lines + [make_line(s[0], s[1]) for s in extra_segments]
    n = len(all_lines)

    star = [V[0], V[2], V[4], V[1], V[3]]

    # Compute pairwise intersections (only those within BOTH segments)
    ixd = {}
    for i, j in combinations(range(n), 2):
        p = seg_seg_intersect(all_segs[i][0], all_segs[i][1],
                              all_segs[j][0], all_segs[j][1])
        if p is not None:
            ixd[(i,j)] = p
            ixd[(j,i)] = p

    # Sorted intersection points on each segment
    lp = {}
    for i in range(n):
        pts = []
        for j in range(n):
            if i != j and (i,j) in ixd:
                pts.append(ixd[(i,j)])
        # Sort along segment direction
        dx = all_segs[i][1][0] - all_segs[i][0][0]
        dy = all_segs[i][1][1] - all_segs[i][0][1]
        pts.sort(key=lambda p: p[0]*dx + p[1]*dy)
        m = []
        for p in pts:
            if not m or not peq(p, m[-1]):
                m.append(p)
        lp[i] = m

    # Count triangular faces
    count = 0
    triangles = []
    for i, j, k in combinations(range(n), 3):
        if (i,j) not in ixd or (i,k) not in ixd or (j,k) not in ixd:
            continue
        pij, pik, pjk = ixd[(i,j)], ixd[(i,k)], ixd[(j,k)]
        if peq(pij, pik) or peq(pij, pjk) or peq(pik, pjk):
            continue

        ok = True
        for li, pa, pb in [(i, pij, pik), (j, pij, pjk), (k, pik, pjk)]:
            pts = lp[li]
            ia = ib = -1
            for idx, p in enumerate(pts):
                if peq(p, pa): ia = idx
                if peq(p, pb): ib = idx
            if ia == -1 or ib == -1 or abs(ia - ib) != 1:
                ok = False
                break

        if ok:
            c = ((pij[0]+pik[0]+pjk[0])/3, (pij[1]+pik[1]+pjk[1])/3)
            if wn(c, star) != 0:
                count += 1
                triangles.append((i, j, k))
    return count, triangles


def clip_line_to_star_segments(line, V):
    """Clip an infinite line to the star boundary (pentagram segments).
    Returns segment (p1, p2) or None."""
    pts = []
    for i in range(5):
        p1 = V[i]
        p2 = V[(i+2)%5]
        p = seg_seg_intersect(p1, p2,
                               # Create a long segment along the line
                               (line[1]*100, -line[0]*100),  # hack
                               (-line[1]*100, line[0]*100))
        # Better: find where line intersects segment
        seg_line = make_line(p1, p2)
        ix_pt = ix_lines(line, seg_line)
        if ix_pt and point_on_seg(ix_pt, p1, p2):
            is_dup = False
            for existing in pts:
                if peq(ix_pt, existing):
                    is_dup = True
                    break
            if not is_dup:
                pts.append(ix_pt)

    if len(pts) < 2:
        return None
    # Sort along line direction
    a, b, c = line
    dx, dy = -b, a
    pts.sort(key=lambda p: p[0]*dx + p[1]*dy)
    return (pts[0], pts[-1])


# === TEST: Regular pentagram ===
print("="*60)
print("REGULAR PENTAGRAM - INFINITE LINES vs SEGMENTS")
print("="*60)

V_reg = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(V_reg[i], V_reg[(i+2)%5]) for i in range(5)]
star_reg = [V_reg[0], V_reg[2], V_reg[4], V_reg[1], V_reg[3]]

# Inner vertices
inner_reg = []
for i, j in combinations(range(5), 2):
    p = ix_lines(pent_lines[i], pent_lines[j])
    if p and not any(peq(p, v) for v in V_reg) and not any(peq(p, v) for v in inner_reg):
        inner_reg.append(p)

all_special = V_reg + inner_reg
pt_names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner_reg))]

# Base
c, _ = count_triangles_with_segments(V_reg, [])
print(f"Base (segments): {c} triangles")

# Generate candidate lines through special points
def same_line(l1, l2):
    def norm(l):
        a, b, c = l
        n = math.sqrt(a*a + b*b)
        if n < EPS: return l
        a, b, c = a/n, b/n, c/n
        if abs(a) > EPS:
            if a < 0: a, b, c = -a, -b, -c
        elif b < 0:
            a, b, c = -a, -b, -c
        return (a, b, c)
    n1, n2 = norm(l1), norm(l2)
    return abs(n1[0]-n2[0]) < 1e-6 and abs(n1[1]-n2[1]) < 1e-6 and abs(n1[2]-n2[2]) < 1e-6

cand_lines = []
cand_names = []
for i, j in combinations(range(len(all_special)), 2):
    line = make_line(all_special[i], all_special[j])
    if any(same_line(line, pl) for pl in pent_lines):
        continue
    if any(same_line(line, cl) for cl in cand_lines):
        continue
    cand_lines.append(line)
    cand_names.append(f"{pt_names[i]}--{pt_names[j]}")

print(f"Candidate lines: {len(cand_lines)}")

# Test all pairs with SEGMENT clipping
print("\nTesting with segments (clipped to star):")
results_seg = []
for ci, cj in combinations(range(len(cand_lines)), 2):
    # Clip each line to the star boundary
    seg1 = clip_line_to_star_segments(cand_lines[ci], V_reg)
    seg2 = clip_line_to_star_segments(cand_lines[cj], V_reg)
    if seg1 is None or seg2 is None:
        continue

    c, tris = count_triangles_with_segments(V_reg, [seg1, seg2])
    results_seg.append((c, ci, cj))

results_seg.sort(reverse=True)
print("Top 15 (segments):")
for c, ci, cj in results_seg[:15]:
    print(f"  {c}: {cand_names[ci]} + {cand_names[cj]}")

from collections import Counter
dist = Counter(c for c, _, _ in results_seg)
print(f"\nDistribution: {dict(sorted(dist.items()))}")

# === Now try with NON-REGULAR pentagram + segments ===
print("\n" + "="*60)
print("NON-REGULAR PENTAGRAM + SEGMENTS")
print("="*60)

random.seed(42)
best = 0
best_config = None

for trial in range(10000):
    # Generate non-regular pentagram
    V = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        angle = base_angle + random.gauss(0, 0.2)
        radius = 1.0 + random.gauss(0, 0.2)
        V.append((radius * math.cos(angle), radius * math.sin(angle)))

    plines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

    # Compute special points
    sp = list(V)
    for i, j in combinations(range(5), 2):
        p = ix_lines(plines[i], plines[j])
        if p and all(not peq(p, s) for s in sp):
            sp.append(p)
    if len(sp) < 10:
        continue

    # Candidate lines
    clns = []
    for i, j in combinations(range(len(sp)), 2):
        line = make_line(sp[i], sp[j])
        if any(same_line(line, pl) for pl in plines):
            continue
        if any(same_line(line, cl) for cl in clns):
            continue
        clns.append(line)

    # Try all pairs with segment clipping
    for ci, cj in combinations(range(len(clns)), 2):
        seg1 = clip_line_to_star_segments(clns[ci], V)
        seg2 = clip_line_to_star_segments(clns[cj], V)
        if seg1 is None or seg2 is None:
            continue

        c, _ = count_triangles_with_segments(V, [seg1, seg2])
        if c > best:
            best = c
            best_config = (V, seg1, seg2)
            print(f"  Trial {trial}: {c} triangles!")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  V = {V}")
                print(f"  seg1 = {seg1}")
                print(f"  seg2 = {seg2}")

    if trial % 2000 == 0:
        print(f"  ... trial {trial}, best: {best}")

print(f"\nBest non-regular + segments: {best}")

# === Random segments (not through special points) ===
print("\n" + "="*60)
print("NON-REGULAR + RANDOM SEGMENTS")
print("="*60)

random.seed(789)
best_rand = 0

for trial in range(50000):
    V = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        angle = base_angle + random.gauss(0, 0.3)
        radius = 1.0 + random.gauss(0, 0.3)
        V.append((radius * math.cos(angle), radius * math.sin(angle)))

    # Random lines, clipped to star
    def rline():
        ang = random.uniform(0, math.pi)
        off = random.uniform(-0.8, 0.8)
        return (math.cos(ang), math.sin(ang), -off)

    l1, l2 = rline(), rline()
    seg1 = clip_line_to_star_segments(l1, V)
    seg2 = clip_line_to_star_segments(l2, V)
    if seg1 is None or seg2 is None:
        continue

    c, _ = count_triangles_with_segments(V, [seg1, seg2])
    if c > best_rand:
        best_rand = c
        print(f"  Trial {trial}: {c} triangles!")
        if c >= 10:
            print(f"  *** FOUND 10! ***")
            print(f"  V = {V}")
            print(f"  seg1 = {seg1}")
            print(f"  seg2 = {seg2}")

    if trial % 10000 == 0:
        print(f"  ... trial {trial}, best: {best_rand}")

print(f"\nBest random: {best_rand}")

# === SUMMARY ===
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Regular pentagram (segments): max {results_seg[0][0] if results_seg else 'N/A'}")
print(f"Non-regular + special points (segments): {best}")
print(f"Non-regular + random (segments): {best_rand}")
