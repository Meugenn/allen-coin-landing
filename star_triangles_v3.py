"""
Финальная версия: считаем треугольники в пентаграмме двумя способами.

Способ 1: АТОМАРНЫЕ грани (regions) — неделимые треугольные области
Способ 2: ВСЕ ВИДИМЫЕ — любой треугольник из 3 линий фигуры,
           вершины которого внутри звезды (в пределах пятиугольника).
           Это "головоломочный" подсчёт: "сколько треугольников видишь?"
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

def point_in_pentagon(point, pentagon_vertices):
    """Check if point is inside or on the convex hull (regular pentagon)."""
    n = len(pentagon_vertices)
    # Check using cross products - all should be same sign
    pos = neg = False
    for i in range(n):
        x1, y1 = pentagon_vertices[i]
        x2, y2 = pentagon_vertices[(i+1)%n]
        cross = (x2-x1)*(point[1]-y1) - (point[0]-x1)*(y2-y1)
        if cross > EPS:
            pos = True
        if cross < -EPS:
            neg = True
        if pos and neg:
            return False
    return True

# === PENTAGRAM ===
outer = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(outer[i], outer[(i+2)%5]) for i in range(5)]

# Star polygon {5/2} for atomic face check
star_polygon = [outer[0], outer[2], outer[4], outer[1], outer[3]]

# Convex pentagon for visible triangle check (vertices in order)
# outer vertices are at angles 90, 162, 234, 306, 378=18 degrees
# Sort by angle for convex hull
pentagon = sorted(outer, key=lambda p: math.atan2(p[1], p[0]))

# Inner vertices
inner = []
for i in range(5):
    for j in range(i+1, 5):
        p = intersect_lines(pent_lines[i], pent_lines[j])
        if p and not any(pt_eq(p, v) for v in outer) and not any(pt_eq(p, v) for v in inner):
            inner.append(p)

all_special = outer + inner
pt_names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner))]


def count_atomic(lines):
    """Count triangular FACES inside star."""
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


def count_visible(lines):
    """Count ALL visible triangles: triples of lines forming triangles
    with all vertices inside the pentagon (the figure's extent)."""
    n = len(lines)
    ix = {}
    for i, j in combinations(range(n), 2):
        p = intersect_lines(lines[i], lines[j])
        if p:
            ix[(i,j)] = p
            ix[(j,i)] = p

    count = 0
    tris = []
    for i, j, k in combinations(range(n), 3):
        if (i,j) not in ix or (i,k) not in ix or (j,k) not in ix:
            continue
        pij, pik, pjk = ix[(i,j)], ix[(i,k)], ix[(j,k)]
        if pt_eq(pij, pik) or pt_eq(pij, pjk) or pt_eq(pik, pjk):
            continue

        # All 3 vertices must be inside the pentagon
        if (point_in_pentagon(pij, pentagon) and
            point_in_pentagon(pik, pentagon) and
            point_in_pentagon(pjk, pentagon)):
            count += 1
            tris.append((i, j, k))
    return count, tris


# === BASE PENTAGRAM ===
print("="*60)
print("БАЗОВАЯ ПЕНТАГРАММА")
print("="*60)
at = count_atomic(pent_lines)
vt, vt_list = count_visible(pent_lines)
print(f"Атомарные грани:      {at}")
print(f"Все видимые:          {vt}")
print(f"(C(5,3) = 10 троек)")
print()
for i, j, k in vt_list:
    pij = intersect_lines(pent_lines[i], pent_lines[j])
    pik = intersect_lines(pent_lines[i], pent_lines[k])
    pjk = intersect_lines(pent_lines[j], pent_lines[k])
    vnames = []
    for p in [pij, pik, pjk]:
        name = "?"
        for idx, sp in enumerate(all_special):
            if pt_eq(p, sp):
                name = pt_names[idx]
                break
        vnames.append(name)
    print(f"  L{i},L{j},L{k}: {vnames}")


# === LINES THROUGH SPECIAL POINTS ===
print("\n" + "="*60)
print("ПОИСК: ЛИНИИ ЧЕРЕЗ ОСОБЫЕ ТОЧКИ")
print("="*60)

def same_line_check(l1, l2):
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
    return abs(n1[0]-n2[0]) < 1e-7 and abs(n1[1]-n2[1]) < 1e-7 and abs(n1[2]-n2[2]) < 1e-7

cand = []
cand_names = []
for i, j in combinations(range(len(all_special)), 2):
    line = make_line(all_special[i], all_special[j])
    if any(same_line_check(line, pl) for pl in pent_lines):
        continue
    if any(same_line_check(line, cl) for cl in cand):
        continue
    cand.append(line)
    cand_names.append(f"{pt_names[i]}--{pt_names[j]}")

print(f"Кандидатов: {len(cand)}")

results = []
for ci, cj in combinations(range(len(cand)), 2):
    lines7 = pent_lines + [cand[ci], cand[cj]]
    a = count_atomic(lines7)
    v, _ = count_visible(lines7)
    results.append((a, v, ci, cj))

# Sort by visible
results_by_vis = sorted(results, key=lambda x: x[1], reverse=True)
print("\nТоп-20 по видимым (атомарные / видимые):")
for a, v, ci, cj in results_by_vis[:20]:
    print(f"  {a:2d} / {v:2d}: {cand_names[ci]} + {cand_names[cj]}")

# Sort by atomic
results_by_at = sorted(results, key=lambda x: x[0], reverse=True)
print("\nТоп-10 по атомарным:")
for a, v, ci, cj in results_by_at[:10]:
    print(f"  {a:2d} / {v:2d}: {cand_names[ci]} + {cand_names[cj]}")

# Distribution
from collections import Counter
print(f"\nРаспределение видимых: {dict(sorted(Counter(v for _,v,_,_ in results).items()))}")
print(f"Распределение атомарных: {dict(sorted(Counter(a for a,_,_,_ in results).items()))}")

# 10 visible
ten = [r for r in results if r[1] == 10]
print(f"\nС ровно 10 видимыми: {len(ten)}")
for a, v, ci, cj in ten:
    print(f"  {a} / {v}: {cand_names[ci]} + {cand_names[cj]}")


# === RANDOM SEARCH ===
print("\n" + "="*60)
print("СЛУЧАЙНЫЙ ПОИСК")
print("="*60)

random.seed(42)
best_vis = 0
best_at = 0
best_cfg = None

for trial in range(200000):
    def rline():
        ang = random.uniform(0, math.pi)
        off = random.uniform(-0.7, 0.7)
        return (math.cos(ang), math.sin(ang), -off)

    l1, l2 = rline(), rline()
    lines7 = pent_lines + [l1, l2]
    v, _ = count_visible(lines7)

    if v > best_vis:
        best_vis = v
        a = count_atomic(lines7)
        best_at = a
        best_cfg = (l1, l2)
        print(f"  Trial {trial}: {a} atomic / {v} visible")

print(f"\nЛучший: {best_at} atomic / {best_vis} visible")

# === Also try one line of each type ===
print("\n" + "="*60)
print("ОДНА ЛИНИЯ (для понимания)")
print("="*60)
for ci in range(len(cand)):
    lines6 = pent_lines + [cand[ci]]
    a = count_atomic(lines6)
    v, _ = count_visible(lines6)
    print(f"  {a:2d} / {v:2d}: {cand_names[ci]}")

# === FINAL SUMMARY ===
print("\n" + "="*60)
print("ИТОГ")
print("="*60)
base_vis = count_visible(pent_lines)[0]
base_at = count_atomic(pent_lines)
print(f"Базовая пентаграмма: {base_at} атомарных граней, {base_vis} видимых треугольников")
print(f"Максимум атомарных с 2 линиями через вершины: {max(a for a,v,ci,cj in results)}")
print(f"Максимум видимых с 2 линиями через вершины:   {max(v for a,v,ci,cj in results)}")
print(f"Максимум видимых случайный поиск:              {best_vis}")
