"""
Different interpretation: count ALL triangles visible in the figure,
not just atomic (minimal) triangular faces.

A "visible triangle" = 3 intersection points A,B,C where each pair
lies on a common drawn line, and the interior of triangle ABC lies
inside the star.

This counts both small (atomic) and large (composite) triangles.
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

def ix(l1, l2):
    d = l1[0]*l2[1] - l2[0]*l1[1]
    if abs(d) < EPS: return None
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
            if y2 > y and (x2-x1)*(y-y1) - (x-x1)*(y2-y1) > 0: w += 1
        else:
            if y2 <= y and (x2-x1)*(y-y1) - (x-x1)*(y2-y1) < 0: w -= 1
    return w

def point_on_seg(p, a, b):
    dx, dy = b[0]-a[0], b[1]-a[1]
    l2 = dx*dx + dy*dy
    if l2 < EPS*EPS: return peq(p, a)
    t = ((p[0]-a[0])*dx + (p[1]-a[1])*dy) / l2
    if t < -EPS or t > 1+EPS: return False
    proj = (a[0]+t*dx, a[1]+t*dy)
    return (p[0]-proj[0])**2 + (p[1]-proj[1])**2 < 1e-12

def point_in_triangle(p, t1, t2, t3):
    """Check if point p is inside triangle t1-t2-t3."""
    def sign(p1, p2, p3):
        return (p1[0]-p3[0])*(p2[1]-p3[1]) - (p2[0]-p3[0])*(p1[1]-p3[1])
    d1 = sign(p, t1, t2)
    d2 = sign(p, t2, t3)
    d3 = sign(p, t3, t1)
    has_neg = (d1 < -EPS) or (d2 < -EPS) or (d3 < -EPS)
    has_pos = (d1 > EPS) or (d2 > EPS) or (d3 > EPS)
    return not (has_neg and has_pos)


def count_all_triangles(V, extra_segs):
    """
    Count ALL visible triangles in the figure.
    A triangle = 3 intersection points where each pair lies on a common line/segment,
    and the triangle's centroid is inside the star.
    """
    pent_segs = [(V[i], V[(i+2)%5]) for i in range(5)]
    all_segs = pent_segs + list(extra_segs)
    all_lines = [make_line(s[0], s[1]) for s in all_segs]
    n = len(all_lines)
    star = [V[0], V[2], V[4], V[1], V[3]]

    # Find ALL intersection points (within segments)
    points = []
    point_lines = {}  # point_idx -> set of line indices it's on

    for i, j in combinations(range(n), 2):
        p = ix(all_lines[i], all_lines[j])
        if p and point_on_seg(p, all_segs[i][0], all_segs[i][1]) and \
               point_on_seg(p, all_segs[j][0], all_segs[j][1]):
            # Check if duplicate
            found = -1
            for k, existing in enumerate(points):
                if peq(p, existing):
                    found = k
                    break
            if found >= 0:
                point_lines[found].add(i)
                point_lines[found].add(j)
            else:
                idx = len(points)
                points.append(p)
                point_lines[idx] = {i, j}

    np = len(points)

    # For each pair of points, check which lines they share
    pair_line = {}  # (i,j) -> line index if they share a common line
    for i in range(np):
        for j in range(i+1, np):
            common = point_lines[i] & point_lines[j]
            if common:
                pair_line[(i,j)] = min(common)  # any common line
                pair_line[(j,i)] = min(common)

    # Count triangles
    triangles = []
    for a, b, c in combinations(range(np), 3):
        if (a,b) not in pair_line or (a,c) not in pair_line or (b,c) not in pair_line:
            continue
        # Check that centroid is inside star
        cx = (points[a][0] + points[b][0] + points[c][0]) / 3
        cy = (points[a][1] + points[b][1] + points[c][1]) / 3
        if wn((cx, cy), star) != 0:
            triangles.append((a, b, c))

    return len(triangles), triangles, points


def clip_to_star(p1, p2, V):
    line = make_line(p1, p2)
    pts = []
    for i in range(5):
        s1, s2 = V[i], V[(i+2)%5]
        sl = make_line(s1, s2)
        p = ix(line, sl)
        if p and point_on_seg(p, s1, s2):
            if not any(peq(p, q) for q in pts):
                pts.append(p)
    if len(pts) < 2: return None
    a, b, c = line
    pts.sort(key=lambda p: -b*p[0] + a*p[1])
    return (pts[0], pts[-1])


# === BASE PENTAGRAM ===
print("=" * 60)
print("ALL-TRIANGLES COUNTING")
print("=" * 60)

V = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
plines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

# Inner vertices
inner = []
for i, j in combinations(range(5), 2):
    p = ix(plines[i], plines[j])
    if p and all(not peq(p, v) for v in V) and all(not peq(p, v) for v in inner):
        inner.append(p)

c, tris, pts = count_all_triangles(V, [])
print(f"Base pentagram: {c} all-triangles")
print(f"  Intersection points: {len(pts)}")

# === Test known 9-atomic-triangle config ===
seg1 = clip_to_star(V[2], inner[4], V)
seg2 = clip_to_star(V[4], inner[0], V)
c, tris, pts = count_all_triangles(V, [seg1, seg2])
print(f"\n9-atomic config (V2-I4 + V4-I0): {c} all-triangles")
print(f"  Intersection points: {len(pts)}")


# === SEARCH for config with exactly 10 all-triangles ===
# Starting from base of 10, we want configs where the total becomes exactly 10
# OR where a different base gives 5 and the modified gives 10
print("\n" + "=" * 60)
print("SEARCHING FOR 10 ALL-TRIANGLES")
print("=" * 60)

# First: find all configs and their all-triangle counts
sp = V + inner
best_10 = None

print("\nLines through special points (regular pentagram):")
cand_segs = []
cand_names = []
for i, j in combinations(range(len(sp)), 2):
    line = make_line(sp[i], sp[j])
    # Skip pentagram lines
    is_pent = False
    for pl in plines:
        a1,b1,c1 = line; a2,b2,c2 = pl
        n1 = math.sqrt(a1*a1+b1*b1); n2 = math.sqrt(a2*a2+b2*b2)
        if n1 > EPS and n2 > EPS:
            a1,b1,c1 = a1/n1,b1/n1,c1/n1
            a2,b2,c2 = a2/n2,b2/n2,c2/n2
            if abs(a1) > EPS:
                if a1 < 0: a1,b1,c1 = -a1,-b1,-c1
            elif b1 < 0: a1,b1,c1 = -a1,-b1,-c1
            if abs(a2) > EPS:
                if a2 < 0: a2,b2,c2 = -a2,-b2,-c2
            elif b2 < 0: a2,b2,c2 = -a2,-b2,-c2
            if abs(a1-a2)<1e-4 and abs(b1-b2)<1e-4 and abs(c1-c2)<1e-4:
                is_pent = True; break
    if is_pent: continue

    seg = clip_to_star(sp[i], sp[j], V)
    if seg:
        dup = False
        for existing, _ in cand_segs:
            if (peq(seg[0], existing[0]) and peq(seg[1], existing[1])) or \
               (peq(seg[0], existing[1]) and peq(seg[1], existing[0])):
                dup = True; break
        if not dup:
            names = [f"V{k}" for k in range(5)] + [f"I{k}" for k in range(5)]
            cand_segs.append((seg, f"{names[i]}-{names[j]}"))

print(f"Candidate segments: {len(cand_segs)}")

from collections import Counter
results = []
for ci, cj in combinations(range(len(cand_segs)), 2):
    c, _, _ = count_all_triangles(V, [cand_segs[ci][0], cand_segs[cj][0]])
    results.append((c, ci, cj))
    if c == 10:
        print(f"  *** 10 all-triangles: {cand_segs[ci][1]} + {cand_segs[cj][1]} ***")

results.sort(reverse=True)
print(f"\nTop results (all-triangle count):")
for c, ci, cj in results[:20]:
    print(f"  {c}: {cand_segs[ci][1]} + {cand_segs[cj][1]}")

dist = Counter(c for c, _, _ in results)
print(f"\nDistribution: {dict(sorted(dist.items()))}")


# === Also try non-regular pentagrams ===
print("\n" + "=" * 60)
print("NON-REGULAR PENTAGRAMS - ALL TRIANGLES")
print("=" * 60)

random.seed(42)
target_counts = {}

for trial in range(3000):
    radii = [1.0 + random.uniform(-0.3, 0.3) for _ in range(5)]
    angles = [math.pi/2 + 2*math.pi*k/5 + random.uniform(-0.2, 0.2) for k in range(5)]
    Vn = [(r*math.cos(a), r*math.sin(a)) for r, a in zip(radii, angles)]

    pln = [make_line(Vn[i], Vn[(i+2)%5]) for i in range(5)]
    inn = []
    for i, j in combinations(range(5), 2):
        p = ix(pln[i], pln[j])
        if p and all(not peq(p, v) for v in Vn) and all(not peq(p, v) for v in inn):
            inn.append(p)
    if len(inn) != 5: continue

    # First check base count
    base_c, _, _ = count_all_triangles(Vn, [])

    spn = Vn + inn

    # Try lines through special points
    csn = []
    for i, j in combinations(range(len(spn)), 2):
        line = make_line(spn[i], spn[j])
        is_pent = False
        for pl in pln:
            a1,b1,c1 = line; a2,b2,c2 = pl
            n1 = math.sqrt(a1*a1+b1*b1); n2 = math.sqrt(a2*a2+b2*b2)
            if n1 > EPS and n2 > EPS:
                a1,b1,c1 = a1/n1,b1/n1,c1/n1; a2,b2,c2 = a2/n2,b2/n2,c2/n2
                if abs(a1) > EPS:
                    if a1 < 0: a1,b1,c1 = -a1,-b1,-c1
                elif b1 < 0: a1,b1,c1 = -a1,-b1,-c1
                if abs(a2) > EPS:
                    if a2 < 0: a2,b2,c2 = -a2,-b2,-c2
                elif b2 < 0: a2,b2,c2 = -a2,-b2,-c2
                if abs(a1-a2)<1e-4 and abs(b1-b2)<1e-4 and abs(c1-c2)<1e-4:
                    is_pent = True; break
        if is_pent: continue
        seg = clip_to_star(spn[i], spn[j], Vn)
        if seg:
            dup = False
            for existing in csn:
                if (peq(seg[0], existing[0]) and peq(seg[1], existing[1])) or \
                   (peq(seg[0], existing[1]) and peq(seg[1], existing[0])):
                    dup = True; break
            if not dup:
                csn.append(seg)

    for ci, cj in combinations(range(min(len(csn), 30)), 2):
        c, _, _ = count_all_triangles(Vn, [csn[ci], csn[cj]])
        key = (base_c, c)
        if key not in target_counts:
            target_counts[key] = 0
        target_counts[key] += 1

        if base_c == 5 and c == 10:
            print(f"  Trial {trial}: BASE=5, WITH_LINES=10 !!!")
            print(f"  Radii: {radii}")
            print(f"  Angles: {angles}")

        if c == 10 and base_c <= 10:
            if trial < 100 or c not in [v for (_,v) in target_counts if target_counts[(_,v)] > 5]:
                print(f"  Trial {trial}: base={base_c}, with_lines={c}")

    if trial % 500 == 0:
        print(f"  ... trial {trial}")

print(f"\n(base, with_lines) distribution:")
for key in sorted(target_counts.keys()):
    if target_counts[key] > 0:
        print(f"  base={key[0]}, lines={key[1]}: {target_counts[key]} configs")
