"""
Search for non-regular pentagrams where 2 additional lines create 10 triangular faces.

Key insight: in a regular pentagram, lines through vertices create concurrencies
(3+ lines through same point), limiting faces to 9 max. In a non-regular pentagram,
these concurrencies can be broken, potentially allowing 10.

Approach:
1. Parametrize pentagram by varying vertex positions
2. For each shape, try lines through pairs of special points + random lines
3. Find configurations with 10 atomic triangular faces inside star
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

def count_atomic_triangles(lines, star):
    """Count triangular faces of line arrangement inside star."""
    n = len(lines)
    ixd = {}
    for i, j in combinations(range(n), 2):
        p = ix_lines(lines[i], lines[j])
        if p:
            ixd[(i,j)] = p
            ixd[(j,i)] = p

    lp = {}
    for i in range(n):
        pts = [ixd[(i,j)] for j in range(n) if i != j and (i,j) in ixd]
        a, b, _ = lines[i]
        dx, dy = -b, a
        pts.sort(key=lambda p: p[0]*dx + p[1]*dy)
        m = []
        for p in pts:
            if not m or not peq(p, m[-1]):
                m.append(p)
        lp[i] = m

    count = 0
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
    return count


def get_special_points(V, pent):
    """Get all intersection points of pentagram lines."""
    sp = list(V)
    for i, j in combinations(range(5), 2):
        p = ix_lines(pent[i], pent[j])
        if p and all(not peq(p, s) for s in sp):
            sp.append(p)
    return sp

def get_candidate_lines(sp, pent):
    """Get candidate lines through pairs of special points, excluding pentagram lines."""
    cands = []
    for i, j in combinations(range(len(sp)), 2):
        line = make_line(sp[i], sp[j])
        if any(same_line(line, pl) for pl in pent):
            continue
        if any(same_line(line, cl) for cl in cands):
            continue
        cands.append(line)
    return cands


# === APPROACH 1: Perturb regular pentagram, use special-point lines ===
print("="*60)
print("APPROACH 1: Perturbed pentagram + special point lines")
print("="*60)

random.seed(42)
best_overall = 0
best_config = None

for trial in range(5000):
    # Generate pentagram with perturbed vertices
    V = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        angle = base_angle + random.gauss(0, 0.15)
        radius = 1.0 + random.gauss(0, 0.15)
        V.append((radius * math.cos(angle), radius * math.sin(angle)))

    pent = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
    star = [V[0], V[2], V[4], V[1], V[3]]

    # Check that the star is valid (all pentagram lines intersect properly)
    try:
        sp = get_special_points(V, pent)
        if len(sp) < 10:
            continue
    except:
        continue

    cands = get_candidate_lines(sp, pent)

    for ci, cj in combinations(range(len(cands)), 2):
        lines7 = pent + [cands[ci], cands[cj]]
        c = count_atomic_triangles(lines7, star)
        if c > best_overall:
            best_overall = c
            best_config = (V, cands[ci], cands[cj], sp)
            print(f"  Trial {trial}: {c} triangles!")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  Vertices: {V}")

    if trial % 1000 == 0:
        print(f"  ... trial {trial}, best so far: {best_overall}")

print(f"\nApproach 1 best: {best_overall}")


# === APPROACH 2: Perturb pentagram + random lines ===
print("\n" + "="*60)
print("APPROACH 2: Perturbed pentagram + random lines")
print("="*60)

random.seed(123)
best2 = 0

for trial in range(50000):
    V = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        angle = base_angle + random.gauss(0, 0.2)
        radius = 1.0 + random.gauss(0, 0.2)
        V.append((radius * math.cos(angle), radius * math.sin(angle)))

    pent = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
    star = [V[0], V[2], V[4], V[1], V[3]]

    # 2 random lines
    def rline():
        ang = random.uniform(0, math.pi)
        off = random.uniform(-0.8, 0.8)
        return (math.cos(ang), math.sin(ang), -off)

    l1, l2 = rline(), rline()
    lines7 = pent + [l1, l2]
    c = count_atomic_triangles(lines7, star)
    if c > best2:
        best2 = c
        print(f"  Trial {trial}: {c} triangles!")
        if c >= 10:
            print(f"  *** FOUND 10! ***")
            print(f"  Vertices: {V}")
            print(f"  L1: {l1}")
            print(f"  L2: {l2}")

    if trial % 10000 == 0:
        print(f"  ... trial {trial}, best so far: {best2}")

print(f"\nApproach 2 best: {best2}")


# === APPROACH 3: Very deformed stars (elongated, skewed) ===
print("\n" + "="*60)
print("APPROACH 3: Highly deformed stars")
print("="*60)

random.seed(456)
best3 = 0

for trial in range(20000):
    # Try very different shapes: elongated, compressed, skewed
    V = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        # Large perturbations
        angle = base_angle + random.uniform(-0.5, 0.5)
        radius = random.uniform(0.3, 2.0)
        V.append((radius * math.cos(angle), radius * math.sin(angle)))

    pent = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
    star = [V[0], V[2], V[4], V[1], V[3]]

    try:
        sp = get_special_points(V, pent)
        if len(sp) < 10:
            continue
    except:
        continue

    # Try a subset of special-point pairs (to save time)
    cands = get_candidate_lines(sp, pent)
    if len(cands) < 2:
        continue

    # Try up to 50 random pairs of candidates
    pairs_to_try = list(combinations(range(len(cands)), 2))
    if len(pairs_to_try) > 50:
        random.shuffle(pairs_to_try)
        pairs_to_try = pairs_to_try[:50]

    for ci, cj in pairs_to_try:
        lines7 = pent + [cands[ci], cands[cj]]
        c = count_atomic_triangles(lines7, star)
        if c > best3:
            best3 = c
            print(f"  Trial {trial}: {c} triangles!")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  Vertices: {V}")

    if trial % 5000 == 0:
        print(f"  ... trial {trial}, best so far: {best3}")

print(f"\nApproach 3 best: {best3}")


# === APPROACH 4: Fix good 2-line config, deform star to break concurrency ===
print("\n" + "="*60)
print("APPROACH 4: Fix lines from regular, deform star to break concurrencies")
print("="*60)

# Start with regular pentagram's best config
V_reg = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_reg = [make_line(V_reg[i], V_reg[(i+2)%5]) for i in range(5)]
sp_reg = get_special_points(V_reg, pent_reg)
cands_reg = get_candidate_lines(sp_reg, pent_reg)

# Find the 9-triangle configs
best_pairs = []
for ci, cj in combinations(range(len(cands_reg)), 2):
    lines7 = pent_reg + [cands_reg[ci], cands_reg[cj]]
    star_reg = [V_reg[0], V_reg[2], V_reg[4], V_reg[1], V_reg[3]]
    c = count_atomic_triangles(lines7, star_reg)
    if c == 9:
        best_pairs.append((ci, cj, cands_reg[ci], cands_reg[cj]))

print(f"Found {len(best_pairs)} 9-triangle configs in regular pentagram")

# For each 9-triangle config, try perturbing the star
best4 = 0
random.seed(789)

for pi, (ci, cj, l1_fixed, l2_fixed) in enumerate(best_pairs):
    print(f"\nConfig {pi}: lines {ci} and {cj}")

    for trial in range(10000):
        # Small perturbation of vertices
        eps = 0.05 + trial * 0.0001  # gradually increase perturbation
        V = []
        for k in range(5):
            base_angle = math.pi/2 + 2*math.pi*k/5
            angle = base_angle + random.gauss(0, eps)
            radius = 1.0 + random.gauss(0, eps)
            V.append((radius * math.cos(angle), radius * math.sin(angle)))

        pent = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
        star = [V[0], V[2], V[4], V[1], V[3]]

        # Use the FIXED lines from regular pentagram
        lines7 = pent + [l1_fixed, l2_fixed]
        c = count_atomic_triangles(lines7, star)
        if c > best4:
            best4 = c
            print(f"  Trial {trial}: {c} triangles! (eps={eps:.4f})")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  Vertices: {V}")
                break

        # Also try NEW optimal lines for this deformed pentagram
        sp = get_special_points(V, pent)
        cands = get_candidate_lines(sp, pent)
        if len(cands) >= 2:
            # Try a few random pairs
            for _ in range(5):
                ci2 = random.randint(0, len(cands)-1)
                cj2 = random.randint(0, len(cands)-1)
                if ci2 == cj2:
                    continue
                lines7b = pent + [cands[ci2], cands[cj2]]
                c2 = count_atomic_triangles(lines7b, star)
                if c2 > best4:
                    best4 = c2
                    print(f"  Trial {trial}: {c2} triangles (new lines)!")
                    if c2 >= 10:
                        print(f"  *** FOUND 10! ***")
                        print(f"  Vertices: {V}")

    if best4 >= 10:
        break

print(f"\nApproach 4 best: {best4}")


# === FINAL SUMMARY ===
print("\n" + "="*60)
print("FINAL SUMMARY")
print("="*60)
print(f"Approach 1 (perturbed + special pts): {best_overall}")
print(f"Approach 2 (perturbed + random lines): {best2}")
print(f"Approach 3 (highly deformed + special pts): {best3}")
print(f"Approach 4 (fixed lines + deformed star): {best4}")
