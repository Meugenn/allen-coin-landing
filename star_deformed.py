"""
Extreme deformations: try stars where the central region changes shape.

Key idea: in a very deformed star, the "central pentagon" might become
a triangle or quadrilateral, fundamentally changing the face structure.

Also try: very elongated stars, asymmetric stars, stars with coincident
inner vertices.
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


class PG:
    def __init__(self):
        self.v = []
        self.adj = {}

    def av(self, p):
        for i, v in enumerate(self.v):
            if abs(v[0]-p[0]) < EPS and abs(v[1]-p[1]) < EPS: return i
        idx = len(self.v)
        self.v.append(p)
        self.adj[idx] = []
        return idx

    def ae(self, i, j):
        if i == j: return
        if j not in self.adj[i]: self.adj[i].append(j)
        if i not in self.adj[j]: self.adj[j].append(i)

    def faces(self):
        for v in self.adj:
            vx, vy = self.v[v]
            self.adj[v].sort(key=lambda u: math.atan2(self.v[u][1]-vy, self.v[u][0]-vx))
        visited = set()
        fs = []
        for u in self.adj:
            for v in self.adj[u]:
                if (u, v) in visited: continue
                face = []
                cu, cv = u, v
                steps = 0
                while True:
                    visited.add((cu, cv))
                    face.append(cu)
                    nb = self.adj[cv]
                    if cu not in nb: break
                    idx = nb.index(cu)
                    nv = nb[(idx - 1) % len(nb)]
                    cu, cv = cv, nv
                    steps += 1
                    if cu == u and cv == v: break
                    if steps > 200: face = []; break
                if len(face) >= 3: fs.append(face)
        return fs

    def area(self, face):
        a = 0
        n = len(face)
        for i in range(n):
            x1, y1 = self.v[face[i]]
            x2, y2 = self.v[face[(i+1)%n]]
            a += x1*y2 - x2*y1
        return a / 2


def count_tris(V, extra_segs):
    segs = [(V[i], V[(i+2)%5]) for i in range(5)] + list(extra_segs)
    star = [V[0], V[2], V[4], V[1], V[3]]

    G = PG()
    for s1, s2 in segs:
        G.av(s1); G.av(s2)

    for i, j in combinations(range(len(segs)), 2):
        l1 = make_line(segs[i][0], segs[i][1])
        l2 = make_line(segs[j][0], segs[j][1])
        p = ix(l1, l2)
        if p and point_on_seg(p, segs[i][0], segs[i][1]) and point_on_seg(p, segs[j][0], segs[j][1]):
            G.av(p)

    for s1, s2 in segs:
        dx, dy = s2[0]-s1[0], s2[1]-s1[1]
        l2 = dx*dx + dy*dy
        pts = []
        for vi, v in enumerate(G.v):
            if point_on_seg(v, s1, s2):
                t = ((v[0]-s1[0])*dx + (v[1]-s1[1])*dy) / l2 if l2 > EPS else 0
                pts.append((t, vi))
        pts.sort()
        filt = []
        for t, vi in pts:
            if not filt or filt[-1][1] != vi: filt.append((t, vi))
        for k in range(len(filt)-1):
            G.ae(filt[k][1], filt[k+1][1])

    faces = G.faces()
    tri = 0
    face_info = []
    for face in faces:
        n = len(face)
        cx = sum(G.v[v][0] for v in face) / n
        cy = sum(G.v[v][1] for v in face) / n
        a = G.area(face)
        inside = wn((cx, cy), star) != 0
        if inside and abs(a) > 1e-8:
            face_info.append(n)
            if n == 3: tri += 1
    return tri, sorted(face_info)


def clip(p1, p2, V):
    line = make_line(p1, p2)
    pts = []
    for i in range(5):
        s1, s2 = V[i], V[(i+2)%5]
        sl = make_line(s1, s2)
        p = ix(line, sl)
        if p and point_on_seg(p, s1, s2):
            if not any(peq(p, q) for q in pts): pts.append(p)
    if len(pts) < 2: return None
    a, b, c = line
    pts.sort(key=lambda p: -b*p[0] + a*p[1])
    return (pts[0], pts[-1])


# ============================================================
# APPROACH 1: Try stars where central region is a triangle/quad
# ============================================================
print("=" * 60)
print("EXTREME DEFORMATIONS")
print("=" * 60)

random.seed(2024)
best = 0
best_info = None

# Strategy: generate stars with extreme aspect ratios
for trial in range(50000):
    # Various extreme deformation strategies
    strategy = trial % 5

    if strategy == 0:
        # Very elongated vertically
        V = []
        stretch = random.uniform(2.0, 5.0)
        for k in range(5):
            a = math.pi/2 + 2*math.pi*k/5 + random.gauss(0, 0.1)
            r = 1.0 + random.gauss(0, 0.1)
            V.append((r*math.cos(a), r*math.sin(a)*stretch))

    elif strategy == 1:
        # Move one vertex far away
        idx = random.randint(0, 4)
        V = []
        for k in range(5):
            a = math.pi/2 + 2*math.pi*k/5
            r = 1.0
            if k == idx:
                r = random.uniform(2.0, 5.0)
                a += random.gauss(0, 0.3)
            V.append((r*math.cos(a), r*math.sin(a)))

    elif strategy == 2:
        # Two vertices very close together
        V = []
        close_pair = random.randint(0, 4)
        for k in range(5):
            a = math.pi/2 + 2*math.pi*k/5
            r = 1.0
            if k == (close_pair + 1) % 5:
                a = math.pi/2 + 2*math.pi*close_pair/5 + random.uniform(0.05, 0.3)
            V.append((r*math.cos(a), r*math.sin(a)))

    elif strategy == 3:
        # Completely random positions
        V = [(random.uniform(-2, 2), random.uniform(-2, 2)) for _ in range(5)]

    else:
        # Parametric with large perturbations
        V = []
        for k in range(5):
            a = math.pi/2 + 2*math.pi*k/5 + random.uniform(-0.8, 0.8)
            r = random.uniform(0.3, 3.0)
            V.append((r*math.cos(a), r*math.sin(a)))

    # Compute inner vertices
    plines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
    inner = []
    for i, j in combinations(range(5), 2):
        p = ix(plines[i], plines[j])
        if p and all(not peq(p, v) for v in V) and all(not peq(p, v) for v in inner):
            inner.append(p)
    if len(inner) < 3: continue

    sp = V + inner

    # Check base structure
    base_c, base_faces = count_tris(V, [])

    # Generate candidate lines
    cand = []
    for i, j in combinations(range(len(sp)), 2):
        seg = clip(sp[i], sp[j], V)
        if seg:
            dup = False
            for existing in cand:
                if (peq(seg[0], existing[0]) and peq(seg[1], existing[1])) or \
                   (peq(seg[0], existing[1]) and peq(seg[1], existing[0])):
                    dup = True; break
            if not dup:
                # Also check it's not a pentagram line
                line = make_line(seg[0], seg[1])
                is_pent = False
                for pl in plines:
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
                if not is_pent:
                    cand.append(seg)

    # Limit to avoid combinatorial explosion
    if len(cand) > 40:
        random.shuffle(cand)
        cand = cand[:40]

    for ci, cj in combinations(range(len(cand)), 2):
        c, faces = count_tris(V, [cand[ci], cand[cj]])
        if c > best:
            best = c
            best_info = (V, cand[ci], cand[cj], faces, base_faces)
            print(f"  Trial {trial} (strat {strategy}): {c} triangles! base={base_faces} faces={faces}")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  V = {V}")
                print(f"  seg1 = {cand[ci]}")
                print(f"  seg2 = {cand[cj]}")

    if trial % 5000 == 0:
        print(f"  ... trial {trial}, best: {best}")

print(f"\nBest: {best}")
if best_info:
    print(f"Base faces: {best_info[4]}")
    print(f"With lines: {best_info[3]}")


# ============================================================
# APPROACH 2: Lines through arbitrary points on star edges
# ============================================================
print("\n" + "=" * 60)
print("LINES THROUGH EDGE POINTS (not just vertices)")
print("=" * 60)

V_reg = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]

random.seed(777)
best2 = 0

for trial in range(100000):
    # Pick 2 random points on star edges, draw line through them
    # Repeat for second line
    def random_edge_point(V):
        # Pick random pentagram edge and random point on it
        i = random.randint(0, 4)
        t = random.uniform(0.05, 0.95)
        p1, p2 = V[i], V[(i+2)%5]
        return (p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1]))

    # Two lines, each defined by 2 points on star edges
    ep1 = random_edge_point(V_reg)
    ep2 = random_edge_point(V_reg)
    seg1 = clip(ep1, ep2, V_reg)

    ep3 = random_edge_point(V_reg)
    ep4 = random_edge_point(V_reg)
    seg2 = clip(ep3, ep4, V_reg)

    if seg1 is None or seg2 is None: continue

    c, faces = count_tris(V_reg, [seg1, seg2])
    if c > best2:
        best2 = c
        print(f"  Trial {trial}: {c} triangles! faces={faces}")
        if c >= 10:
            print(f"  *** FOUND 10! ***")

    if trial % 20000 == 0:
        print(f"  ... trial {trial}, best: {best2}")

print(f"\nBest edge-points: {best2}")


# ============================================================
# APPROACH 3: Non-regular star + random lines (massive search)
# ============================================================
print("\n" + "=" * 60)
print("MASSIVE RANDOM SEARCH")
print("=" * 60)

random.seed(314)
best3 = 0

for trial in range(200000):
    # Random star
    V = []
    for k in range(5):
        a = math.pi/2 + 2*math.pi*k/5 + random.gauss(0, 0.4)
        r = 1.0 + random.gauss(0, 0.4)
        if r < 0.3: r = 0.3
        V.append((r*math.cos(a), r*math.sin(a)))

    # Random line through star
    def rline(V):
        a = random.uniform(0, math.pi)
        o = random.uniform(-1.5, 1.5)
        line = (math.cos(a), math.sin(a), -o)
        pts = []
        for i in range(5):
            s1, s2 = V[i], V[(i+2)%5]
            sl = make_line(s1, s2)
            p = ix(line, sl)
            if p and point_on_seg(p, s1, s2):
                if not any(peq(p, q) for q in pts): pts.append(p)
        if len(pts) < 2: return None
        a, b, c = line
        pts.sort(key=lambda p: -b*p[0] + a*p[1])
        return (pts[0], pts[-1])

    seg1 = rline(V)
    seg2 = rline(V)
    if seg1 is None or seg2 is None: continue

    c, faces = count_tris(V, [seg1, seg2])
    if c > best3:
        best3 = c
        print(f"  Trial {trial}: {c} triangles! faces={faces}")
        if c >= 10:
            print(f"  *** FOUND 10! ***")
            print(f"  V = {V}")
            print(f"  seg1 = {seg1}")
            print(f"  seg2 = {seg2}")

    if trial % 50000 == 0:
        print(f"  ... trial {trial}, best: {best3}")

print(f"\nBest random: {best3}")


print("\n" + "=" * 60)
print("FINAL RESULTS")
print("=" * 60)
print(f"Extreme deformations: {best}")
print(f"Edge-point lines: {best2}")
print(f"Massive random: {best3}")
