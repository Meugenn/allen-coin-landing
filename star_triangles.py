"""
Задача: пятиконечная звезда (пентаграмма) содержит 5 треугольников (кончики).
Добавить 2 прямые линии так, чтобы треугольных областей стало 10.

Подход:
1. Моделируем пентаграмму как arrangement of lines
2. Для каждой пары дополнительных прямых считаем треугольные грани (faces)
3. Считаем только грани внутри звезды (non-zero winding number)
"""

import math
from itertools import combinations
import random

EPS = 1e-9

def make_line(p1, p2):
    """Returns (a, b, c) for ax + by + c = 0"""
    a = p2[1] - p1[1]
    b = p1[0] - p2[0]
    c = p2[0]*p1[1] - p1[0]*p2[1]
    return (a, b, c)

def intersect_lines(l1, l2):
    det = l1[0]*l2[1] - l2[0]*l1[1]
    if abs(det) < EPS:
        return None
    x = (l1[1]*l2[2] - l2[1]*l1[2]) / det
    y = (l2[0]*l1[2] - l1[0]*l2[2]) / det
    return (x, y)

def pt_eq(p1, p2):
    return abs(p1[0]-p2[0]) < EPS and abs(p1[1]-p2[1]) < EPS

def winding_number(point, polygon):
    """Compute winding number of point w.r.t. polygon."""
    x, y = point
    wn = 0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i+1) % n]
        if y1 <= y:
            if y2 > y:
                # Upward crossing
                cross = (x2-x1)*(y-y1) - (x-x1)*(y2-y1)
                if cross > 0:
                    wn += 1
        else:
            if y2 <= y:
                # Downward crossing
                cross = (x2-x1)*(y-y1) - (x-x1)*(y2-y1)
                if cross < 0:
                    wn -= 1
    return wn

def triangle_centroid(p1, p2, p3):
    return ((p1[0]+p2[0]+p3[0])/3, (p1[1]+p2[1]+p3[1])/3)

def count_triangles(lines, star_polygon=None):
    """
    Count triangular faces in arrangement of lines.
    If star_polygon given, only count faces whose centroid is inside the star.
    """
    n = len(lines)

    # All pairwise intersections
    ix = {}
    for i, j in combinations(range(n), 2):
        p = intersect_lines(lines[i], lines[j])
        if p:
            ix[(i,j)] = p
            ix[(j,i)] = p

    # Sorted unique points on each line
    lpts = {}
    for i in range(n):
        pts = []
        for j in range(n):
            if i != j and (i,j) in ix:
                pts.append(ix[(i,j)])
        a, b, c = lines[i]
        dx, dy = -b, a
        pts.sort(key=lambda p: p[0]*dx + p[1]*dy)

        # Merge coincident points
        merged = []
        for p in pts:
            if not merged or not pt_eq(p, merged[-1]):
                merged.append(p)
        lpts[i] = merged

    # Count triangular faces
    count = 0
    triangles = []

    for i, j, k in combinations(range(n), 3):
        if (i,j) not in ix or (i,k) not in ix or (j,k) not in ix:
            continue

        pij, pik, pjk = ix[(i,j)], ix[(i,k)], ix[(j,k)]

        if pt_eq(pij, pik) or pt_eq(pij, pjk) or pt_eq(pik, pjk):
            continue

        ok = True
        for line_idx, pa, pb in [(i, pij, pik), (j, pij, pjk), (k, pik, pjk)]:
            pts = lpts[line_idx]
            ia = ib = -1
            for idx, p in enumerate(pts):
                if pt_eq(p, pa): ia = idx
                if pt_eq(p, pb): ib = idx
            if ia == -1 or ib == -1 or abs(ia - ib) != 1:
                ok = False
                break

        if ok:
            # Check if inside star
            if star_polygon is not None:
                centroid = triangle_centroid(pij, pik, pjk)
                if winding_number(centroid, star_polygon) == 0:
                    continue
            count += 1
            triangles.append((i, j, k, pij, pik, pjk))

    return count, triangles


# === PENTAGRAM SETUP ===
outer = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(outer[i], outer[(i+2)%5]) for i in range(5)]

# Compute inner vertices
inner = []
for i in range(5):
    for j in range(i+1, 5):
        p = intersect_lines(pent_lines[i], pent_lines[j])
        if p:
            is_outer = any(pt_eq(p, v) for v in outer)
            if not is_outer:
                is_dup = any(pt_eq(p, v) for v in inner)
                if not is_dup:
                    inner.append(p)

print(f"Outer vertices: {len(outer)}")
print(f"Inner vertices: {len(inner)}")

# Star polygon boundary (for inside check)
# The star {5/2} visits: V0, V2, V4, V1, V3
star_polygon = [outer[0], outer[2], outer[4], outer[1], outer[3]]

# Verify base
count0, _ = count_triangles(pent_lines, star_polygon)
print(f"Base pentagram: {count0} triangles inside star")
count0_all, _ = count_triangles(pent_lines)
print(f"Base pentagram: {count0_all} triangles total (all bounded)")

# === STRATEGY 1: Lines through special points ===
print("\n=== Strategy 1: Lines through pairs of 10 special points ===")
all_points = outer + inner
point_names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner))]

candidate_lines_sp = []
candidate_descs_sp = []

for i, j in combinations(range(len(all_points)), 2):
    line = make_line(all_points[i], all_points[j])

    # Skip if same as existing pentagram line
    is_existing = False
    for pl in pent_lines:
        p = intersect_lines(line, pl)
        if p is None:  # parallel
            # Check if same line (a point on one satisfies the other)
            a, b, c = pl
            tx, ty = all_points[i]
            if abs(a*tx + b*ty + c) < EPS:
                is_existing = True
                break
    if is_existing:
        continue

    # Skip duplicates
    is_dup = False
    for cl, _ in candidate_lines_sp:
        p = intersect_lines(line, cl)
        if p is None:
            a, b, c = cl
            tx, ty = all_points[i]
            if abs(a*tx + b*ty + c) < EPS:
                is_dup = True
                break
    if not is_dup:
        candidate_lines_sp.append((line, f"{point_names[i]}--{point_names[j]}"))

print(f"Candidate lines through special points: {len(candidate_lines_sp)}")

# Test all pairs
results_sp = []
for ci, cj in combinations(range(len(candidate_lines_sp)), 2):
    all_lines = pent_lines + [candidate_lines_sp[ci][0], candidate_lines_sp[cj][0]]
    c_in, _ = count_triangles(all_lines, star_polygon)
    c_all, _ = count_triangles(all_lines)
    results_sp.append((c_in, c_all, ci, cj))

results_sp.sort(key=lambda x: x[0], reverse=True)
print("Top results (inside star / total):")
for c_in, c_all, ci, cj in results_sp[:15]:
    print(f"  {c_in} / {c_all}: {candidate_lines_sp[ci][1]} + {candidate_lines_sp[cj][1]}")

# === STRATEGY 2: Pentagon sides ===
print("\n=== Strategy 2: Sides of outer pentagon ===")
pent_sides = [make_line(outer[i], outer[(i+1)%5]) for i in range(5)]

for ci, cj in combinations(range(5), 2):
    all_lines = pent_lines + [pent_sides[ci], pent_sides[cj]]
    c_in, _ = count_triangles(all_lines, star_polygon)
    c_all, _ = count_triangles(all_lines)
    side_i = f"side V{ci}-V{(ci+1)%5}"
    side_j = f"side V{cj}-V{(cj+1)%5}"
    print(f"  {c_in} / {c_all}: {side_i} + {side_j}")

# === STRATEGY 3: Mix special point lines with pentagon sides ===
print("\n=== Strategy 3: One special-point line + one pentagon side ===")
results_mix = []
for ci in range(len(candidate_lines_sp)):
    for si in range(5):
        all_lines = pent_lines + [candidate_lines_sp[ci][0], pent_sides[si]]
        c_in, _ = count_triangles(all_lines, star_polygon)
        c_all, _ = count_triangles(all_lines)
        results_mix.append((c_in, c_all, ci, si))

results_mix.sort(key=lambda x: x[0], reverse=True)
print("Top results (inside star / total):")
for c_in, c_all, ci, si in results_mix[:10]:
    print(f"  {c_in} / {c_all}: {candidate_lines_sp[ci][1]} + side V{si}-V{(si+1)%5}")

# === STRATEGY 4: Random lines ===
print("\n=== Strategy 4: Random line search ===")
random.seed(42)
best_in = 0
best_config = None

for trial in range(300000):
    def rline():
        angle = random.uniform(0, math.pi)
        offset = random.uniform(-0.8, 0.8)
        return (math.cos(angle), math.sin(angle), -offset)

    l1, l2 = rline(), rline()
    all_lines = pent_lines + [l1, l2]
    c_in, _ = count_triangles(all_lines, star_polygon)

    if c_in > best_in:
        best_in = c_in
        best_config = (l1, l2)
        c_all, _ = count_triangles(all_lines)
        print(f"  Trial {trial}: {c_in} inside / {c_all} total")

print(f"\nBest random: {best_in} triangles inside star")

# === STRATEGY 5: Systematic grid of horizontal/vertical lines ===
print("\n=== Strategy 5: Systematic horizontal + vertical lines ===")
best_sys = 0
best_sys_config = None

# Generate horizontal and vertical lines at various offsets
h_lines = [(0, 1, -y) for y in [i*0.05 for i in range(-18, 19)]]  # y = const
v_lines = [(1, 0, -x) for x in [i*0.05 for i in range(-18, 19)]]  # x = const
# Also add diagonal lines
d_lines = []
for angle_deg in range(0, 180, 5):
    angle = math.radians(angle_deg)
    a, b = math.cos(angle), math.sin(angle)
    for offset in [i*0.05 for i in range(-18, 19)]:
        d_lines.append((a, b, -offset))

all_candidate_lines = h_lines + v_lines + d_lines

print(f"Systematic candidates: {len(all_candidate_lines)}")

# Sample pairs (too many for exhaustive search)
results_sys = {}
random.seed(123)
sample_size = 500000

for _ in range(sample_size):
    i = random.randint(0, len(all_candidate_lines)-1)
    j = random.randint(0, len(all_candidate_lines)-1)
    if i == j:
        continue

    l1, l2 = all_candidate_lines[i], all_candidate_lines[j]
    all_lines = pent_lines + [l1, l2]
    c_in, _ = count_triangles(all_lines, star_polygon)

    if c_in > best_sys:
        best_sys = c_in
        best_sys_config = (l1, l2, i, j)
        c_all, _ = count_triangles(all_lines)
        print(f"  New best: {c_in} inside / {c_all} total, lines: {l1}, {l2}")

print(f"\nBest systematic: {best_sys} triangles inside star")

# === SUMMARY ===
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Base pentagram: {count0} triangles")
print(f"Best from special points: {results_sp[0][0]} inside / {results_sp[0][1]} total")
print(f"Best from random: {best_in}")
print(f"Best from systematic: {best_sys}")
