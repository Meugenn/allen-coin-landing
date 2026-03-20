"""
Расширенный анализ треугольников в пентаграмме.

Два способа подсчёта:
1. Атомарные треугольные ГРАНИ (regions) — подразделение плоскости
2. ВСЕ видимые треугольники — любые 3 точки пересечения, образующие треугольник,
   все стороны которого лежат на линиях фигуры и внутри звезды

Это важно, потому что в головоломках "сколько треугольников"
считают ВСЕ треугольники, а не только атомарные области.
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

def intersect_lines(l1, l2):
    det = l1[0]*l2[1] - l2[0]*l1[1]
    if abs(det) < EPS:
        return None
    x = (l1[1]*l2[2] - l2[1]*l1[2]) / det
    y = (l2[0]*l1[2] - l1[0]*l2[2]) / det
    return (x, y)

def pt_eq(p1, p2):
    return abs(p1[0]-p2[0]) < EPS and abs(p1[1]-p2[1]) < EPS

def dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def winding_number(point, polygon):
    x, y = point
    wn = 0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i+1) % n]
        if y1 <= y:
            if y2 > y:
                cross = (x2-x1)*(y-y1) - (x-x1)*(y2-y1)
                if cross > 0:
                    wn += 1
        else:
            if y2 <= y:
                cross = (x2-x1)*(y-y1) - (x-x1)*(y2-y1)
                if cross < 0:
                    wn -= 1
    return wn

def point_in_convex_polygon(point, polygon):
    """Check if point is inside convex polygon (the outer pentagon)."""
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i+1) % n]
        cross = (x2-x1)*(point[1]-y1) - (point[0]-x1)*(y2-y1)
        if cross < -EPS:
            return False
    return True

def segment_inside_star(p1, p2, star_polygon, n_samples=20):
    """Check if segment p1->p2 stays inside the star (non-zero winding)."""
    for i in range(1, n_samples):
        t = i / n_samples
        x = p1[0] + t * (p2[0] - p1[0])
        y = p1[1] + t * (p2[1] - p1[1])
        if winding_number((x, y), star_polygon) == 0:
            return False
    return True

def count_atomic_triangles(lines, star_polygon):
    """Count triangular FACES of the arrangement that are inside the star."""
    n = len(lines)
    ix = {}
    for i, j in combinations(range(n), 2):
        p = intersect_lines(lines[i], lines[j])
        if p:
            ix[(i,j)] = p
            ix[(j,i)] = p

    lpts = {}
    for i in range(n):
        pts = []
        for j in range(n):
            if i != j and (i,j) in ix:
                pts.append(ix[(i,j)])
        a, b, c = lines[i]
        dx, dy = -b, a
        pts.sort(key=lambda p: p[0]*dx + p[1]*dy)
        merged = []
        for p in pts:
            if not merged or not pt_eq(p, merged[-1]):
                merged.append(p)
        lpts[i] = merged

    count = 0
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
            centroid = ((pij[0]+pik[0]+pjk[0])/3, (pij[1]+pik[1]+pjk[1])/3)
            if winding_number(centroid, star_polygon) != 0:
                count += 1
    return count


def count_all_visible_triangles(lines, star_polygon):
    """
    Count ALL visible triangles: any triple of lines forming a triangle
    where all 3 vertices are inside/on the star and all 3 sides stay inside.
    """
    n = len(lines)
    ix = {}
    for i, j in combinations(range(n), 2):
        p = intersect_lines(lines[i], lines[j])
        if p:
            ix[(i,j)] = p
            ix[(j,i)] = p

    count = 0
    triangles = []
    for i, j, k in combinations(range(n), 3):
        if (i,j) not in ix or (i,k) not in ix or (j,k) not in ix:
            continue
        pij, pik, pjk = ix[(i,j)], ix[(i,k)], ix[(j,k)]
        if pt_eq(pij, pik) or pt_eq(pij, pjk) or pt_eq(pik, pjk):
            continue

        # Check all vertices are inside star (or on boundary - use relaxed check)
        all_inside = True
        for p in [pij, pik, pjk]:
            w = winding_number(p, star_polygon)
            if w == 0:
                # Check if on boundary (within EPS of star boundary)
                on_boundary = False
                for si in range(len(star_polygon)):
                    p1 = star_polygon[si]
                    p2 = star_polygon[(si+1)%len(star_polygon)]
                    # Distance from point to segment
                    dx, dy = p2[0]-p1[0], p2[1]-p1[1]
                    len_sq = dx*dx + dy*dy
                    if len_sq < EPS:
                        continue
                    t = max(0, min(1, ((p[0]-p1[0])*dx + (p[1]-p1[1])*dy) / len_sq))
                    closest = (p1[0] + t*dx, p1[1] + t*dy)
                    d = dist(p, closest)
                    if d < 1e-6:
                        on_boundary = True
                        break
                if not on_boundary:
                    all_inside = False
                    break
        if not all_inside:
            continue

        # Check all 3 sides stay inside star
        sides_ok = True
        for pa, pb in [(pij, pik), (pij, pjk), (pik, pjk)]:
            if not segment_inside_star(pa, pb, star_polygon, 30):
                sides_ok = False
                break

        if sides_ok:
            count += 1
            triangles.append((i, j, k))

    return count, triangles


# === PENTAGRAM ===
outer = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(outer[i], outer[(i+2)%5]) for i in range(5)]

# Star polygon {5/2}: visits V0, V2, V4, V1, V3
star_polygon = [outer[0], outer[2], outer[4], outer[1], outer[3]]

# Inner vertices
inner = []
for i in range(5):
    for j in range(i+1, 5):
        p = intersect_lines(pent_lines[i], pent_lines[j])
        if p and not any(pt_eq(p, v) for v in outer) and not any(pt_eq(p, v) for v in inner):
            inner.append(p)

all_special = outer + inner
pt_names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner))]

print("="*60)
print("БАЗОВАЯ ПЕНТАГРАММА (5 линий)")
print("="*60)
at = count_atomic_triangles(pent_lines, star_polygon)
vt, vt_list = count_all_visible_triangles(pent_lines, star_polygon)
print(f"Атомарные треугольные грани: {at}")
print(f"Все видимые треугольники:    {vt}")
print(f"C(5,3) = {len(list(combinations(range(5), 3)))} всего троек линий")

# Show which triples are visible
for i, j, k in vt_list:
    pij = intersect_lines(pent_lines[i], pent_lines[j])
    pik = intersect_lines(pent_lines[i], pent_lines[k])
    pjk = intersect_lines(pent_lines[j], pent_lines[k])
    vertices = []
    for p in [pij, pik, pjk]:
        name = "?"
        for idx, sp in enumerate(all_special):
            if pt_eq(p, sp):
                name = pt_names[idx]
                break
        vertices.append(name)
    print(f"  L{i},L{j},L{k}: вершины {vertices}")

# === Lines through special points ===
print("\n" + "="*60)
print("ЛИНИИ ЧЕРЕЗ ОСОБЫЕ ТОЧКИ")
print("="*60)

def same_line(l1, l2):
    def normalize(l):
        a, b, c = l
        norm = math.sqrt(a*a + b*b)
        if norm < EPS: return l
        a, b, c = a/norm, b/norm, c/norm
        if abs(a) > EPS:
            if a < 0: a, b, c = -a, -b, -c
        elif b < 0:
            a, b, c = -a, -b, -c
        return (a, b, c)
    n1, n2 = normalize(l1), normalize(l2)
    return abs(n1[0]-n2[0]) < 1e-7 and abs(n1[1]-n2[1]) < 1e-7 and abs(n1[2]-n2[2]) < 1e-7

cand_lines = []
cand_descs = []
for i, j in combinations(range(len(all_special)), 2):
    line = make_line(all_special[i], all_special[j])
    if any(same_line(line, pl) for pl in pent_lines):
        continue
    if any(same_line(line, cl) for cl in cand_lines):
        continue
    cand_lines.append(line)
    cand_descs.append(f"{pt_names[i]}--{pt_names[j]}")

print(f"Кандидатов: {len(cand_lines)}")

# Test all pairs - both counts
results = []
for ci, cj in combinations(range(len(cand_lines)), 2):
    lines7 = pent_lines + [cand_lines[ci], cand_lines[cj]]
    at = count_atomic_triangles(lines7, star_polygon)
    vt, _ = count_all_visible_triangles(lines7, star_polygon)
    results.append((at, vt, ci, cj))

results.sort(key=lambda x: x[1], reverse=True)
print("\nТоп-15 по ВСЕМ видимым треугольникам (атомарные / все видимые):")
for at, vt, ci, cj in results[:15]:
    print(f"  {at} / {vt}: {cand_descs[ci]} + {cand_descs[cj]}")

results.sort(key=lambda x: x[0], reverse=True)
print("\nТоп-10 по атомарным граням:")
for at, vt, ci, cj in results[:10]:
    print(f"  {at} / {vt}: {cand_descs[ci]} + {cand_descs[cj]}")

# Check for exactly 10 visible
ten_vis = [(at, vt, ci, cj) for at, vt, ci, cj in results if vt == 10]
print(f"\nКонфигурации с ровно 10 видимыми треугольниками: {len(ten_vis)}")
for at, vt, ci, cj in ten_vis:
    print(f"  {at} atomic / {vt} visible: {cand_descs[ci]} + {cand_descs[cj]}")

# Distribution
from collections import Counter
at_dist = Counter(at for at, vt, ci, cj in results)
vt_dist = Counter(vt for at, vt, ci, cj in results)
print(f"\nРаспределение атомарных: {dict(sorted(at_dist.items()))}")
print(f"Распределение видимых:   {dict(sorted(vt_dist.items()))}")

# === Pentagon sides ===
print("\n" + "="*60)
print("СТОРОНЫ ПЯТИУГОЛЬНИКА")
print("="*60)
pent_sides = [make_line(outer[i], outer[(i+1)%5]) for i in range(5)]
side_names = [f"side({i}-{(i+1)%5})" for i in range(5)]

for ci, cj in combinations(range(5), 2):
    lines7 = pent_lines + [pent_sides[ci], pent_sides[cj]]
    at = count_atomic_triangles(lines7, star_polygon)
    vt, _ = count_all_visible_triangles(lines7, star_polygon)
    print(f"  {at} / {vt}: {side_names[ci]} + {side_names[cj]}")

# Mix: special + pentagon side
print("\nМикс (особая точка + сторона пятиугольника):")
mix_results = []
for ci in range(len(cand_lines)):
    for si in range(5):
        lines7 = pent_lines + [cand_lines[ci], pent_sides[si]]
        at = count_atomic_triangles(lines7, star_polygon)
        vt, _ = count_all_visible_triangles(lines7, star_polygon)
        mix_results.append((at, vt, ci, si))

mix_results.sort(key=lambda x: x[1], reverse=True)
print("Топ-10 по видимым:")
for at, vt, ci, si in mix_results[:10]:
    print(f"  {at} / {vt}: {cand_descs[ci]} + {side_names[si]}")

# === Random search for ALL visible ===
print("\n" + "="*60)
print("СЛУЧАЙНЫЙ ПОИСК (оптимизация видимых треугольников)")
print("="*60)
random.seed(42)
best_vt = 0
best_config = None

for trial in range(100000):
    angle1 = random.uniform(0, math.pi)
    offset1 = random.uniform(-0.8, 0.8)
    angle2 = random.uniform(0, math.pi)
    offset2 = random.uniform(-0.8, 0.8)
    l1 = (math.cos(angle1), math.sin(angle1), -offset1)
    l2 = (math.cos(angle2), math.sin(angle2), -offset2)

    lines7 = pent_lines + [l1, l2]
    vt, _ = count_all_visible_triangles(lines7, star_polygon)

    if vt > best_vt:
        best_vt = vt
        at = count_atomic_triangles(lines7, star_polygon)
        best_config = (l1, l2)
        print(f"  Trial {trial}: {at} atomic / {vt} visible")

print(f"\nЛучший случайный: {best_vt} видимых треугольников")

# === SUMMARY ===
print("\n" + "="*60)
print("ИТОГИ")
print("="*60)
print(f"Базовая пентаграмма: 5 атомарных, {count_all_visible_triangles(pent_lines, star_polygon)[0]} видимых")
print(f"Максимум атомарных с 2 линиями: {max(at for at, vt, ci, cj in results)}")
print(f"Максимум видимых с 2 линиями:   {max(vt for at, vt, ci, cj in results)} (через особые точки)")
print(f"Максимум видимых случайный:      {best_vt}")
