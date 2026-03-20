"""
Final comprehensive analysis.

Key question: IS 10 even theoretically possible?

Using Euler's formula: V - E + F = 2
For a planar graph bounded by the star region.

Let's analyze the constraints mathematically.
"""

import math
from itertools import combinations

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

def point_on_seg(p, a, b):
    dx, dy = b[0]-a[0], b[1]-a[1]
    l2 = dx*dx + dy*dy
    if l2 < EPS*EPS: return peq(p, a)
    t = ((p[0]-a[0])*dx + (p[1]-a[1])*dy) / l2
    if t < -EPS or t > 1+EPS: return False
    proj = (a[0]+t*dx, a[1]+t*dy)
    return (p[0]-proj[0])**2 + (p[1]-proj[1])**2 < 1e-12


V = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
plines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

inner = []
for i, j in combinations(range(5), 2):
    p = ix(plines[i], plines[j])
    if p and all(not peq(p, v) for v in V) and all(not peq(p, v) for v in inner):
        inner.append(p)

print("=" * 60)
print("THEORETICAL ANALYSIS: IS 10 POSSIBLE?")
print("=" * 60)

print("""
BASE PENTAGRAM:
- 5 pentagram lines (segments from V[i] to V[(i+2)%5])
- 10 intersection points: 5 outer (V0-V4) + 5 inner (I0-I4)
- Each line has 4 points on it -> 3 segments per line -> 15 edges total
- But edges on the boundary are shared -> internal edges + boundary edges
- Faces inside star: 5 triangles + 1 pentagon = 6 faces

Euler check for the star region:
V=10, E=15 (edges of faces inside), F_inner=6, F_outer=1
V-E+F = 10-15+7 = 2 ✓

ADDING 2 LINES:
Each new line (segment) enters the star, intersects existing lines,
creating new vertices and edges.

A line crossing the star boundary at 2 points adds:
- 2 boundary vertices (on existing edges, splitting them)
- Interior intersection points with each existing line it crosses
- The new segment itself gets split into sub-segments

For MAXIMUM triangulation:
Best case: each new line crosses as many existing lines as possible
inside the star, creating maximum new intersection points.

A line can cross at most 5 pentagram lines + 1 other new line = 6 lines.
But it can only cross lines that it actually intersects INSIDE the star.

Typically, a line crossing through the star hits 4 pentagram segments
(entering through one edge, crossing 3-4 interior segments, exiting).

Let's compute exactly for the 9-triangle configs.
""")

# Detailed analysis for V2->I4 + V4->I0 config
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

seg1 = clip(V[2], inner[4], V)
seg2 = clip(V[4], inner[0], V)

# Count intersections of each new segment with pentagram lines
print("Config V2->I4 + V4->I0:")
for name, seg in [("Seg1 (V2-I4)", seg1), ("Seg2 (V4-I0)", seg2)]:
    line = make_line(seg[0], seg[1])
    hits = 0
    for i in range(5):
        s1, s2 = V[i], V[(i+2)%5]
        sl = make_line(s1, s2)
        p = ix(line, sl)
        if p and point_on_seg(p, s1, s2) and point_on_seg(p, seg[0], seg[1]):
            if not peq(p, seg[0]) and not peq(p, seg[1]):
                hits += 1
                print(f"  {name} crosses pentagram line {i} at ({p[0]:.4f}, {p[1]:.4f})")
    print(f"  {name}: {hits} interior crossings with pentagram lines")

# Check if seg1 and seg2 intersect
l1 = make_line(seg1[0], seg1[1])
l2 = make_line(seg2[0], seg2[1])
p = ix(l1, l2)
if p and point_on_seg(p, seg1[0], seg1[1]) and point_on_seg(p, seg2[0], seg2[1]):
    print(f"  Seg1-Seg2 intersection: ({p[0]:.4f}, {p[1]:.4f})")
else:
    print(f"  Seg1 and Seg2 do NOT intersect inside star!")

# Euler formula analysis
print("""
EULER FORMULA ANALYSIS for 9-triangle config:
Faces: 9 tri + 1 quad + 1 pent = 11 internal + 1 outer = 12
""")

# Let's enumerate V, E precisely
all_segs = [(V[i], V[(i+2)%5]) for i in range(5)] + [seg1, seg2]
all_lines = [make_line(s[0], s[1]) for s in all_segs]

# All vertices (intersections)
all_pts = []
for s in all_segs:
    for p in [s[0], s[1]]:
        if not any(peq(p, q) for q in all_pts):
            all_pts.append(p)

for i, j in combinations(range(7), 2):
    l1 = make_line(all_segs[i][0], all_segs[i][1])
    l2 = make_line(all_segs[j][0], all_segs[j][1])
    p = ix(l1, l2)
    if p and point_on_seg(p, all_segs[i][0], all_segs[i][1]) and \
            point_on_seg(p, all_segs[j][0], all_segs[j][1]):
        if not any(peq(p, q) for q in all_pts):
            all_pts.append(p)

# Count edges: for each segment, count sub-segments between consecutive points
total_edges = 0
for seg in all_segs:
    dx, dy = seg[1][0]-seg[0][0], seg[1][1]-seg[0][1]
    l2 = dx*dx + dy*dy
    pts_on = []
    for p in all_pts:
        if point_on_seg(p, seg[0], seg[1]):
            t = ((p[0]-seg[0][0])*dx + (p[1]-seg[0][1])*dy) / l2 if l2 > EPS else 0
            pts_on.append((t, p))
    pts_on.sort()
    # Deduplicate
    filt = [pts_on[0]]
    for t, p in pts_on[1:]:
        if not peq(p, filt[-1][1]):
            filt.append((t, p))
    total_edges += len(filt) - 1

V_count = len(all_pts)
E_count = total_edges
F_count = 2 - V_count + E_count  # Euler: V - E + F = 2

print(f"Vertices: {V_count}")
print(f"Edges: {E_count}")
print(f"Faces (by Euler): {F_count} (= {F_count-1} bounded + 1 outer)")
print(f"Expected bounded faces: 11 (9 tri + 1 quad + 1 pent)")

print(f"""
FOR 10 TRIANGLES TO EXIST:
The total number of faces is determined by Euler's formula.
With V={V_count}, E={E_count}: F={F_count} ({F_count-1} bounded).

To get 10 triangles, we need:
- Either MORE total faces (need different V,E)
- Or convert the quad to a triangle (merge 2 vertices)

If quad → triangle: V decreases by 1, E decreases by 1, F stays same.
V={V_count-1}, E={E_count-1}, F={F_count}: V-E+F = {V_count-1}-{E_count-1}+{F_count} = 2 ✓

So the quad COULD theoretically become a triangle if we can find
a configuration where the right intersection points coincide.

The quad in the 9-tri config has vertices: I4, I1, I0, X(center).
For it to become a triangle, we need one pair to merge:
- I4=I1: impossible (these are fixed pentagram intersections)
- I4=X: X is intersection of the 2 new lines; I4 is pentagram intersection
  → Need the 2 new lines to cross at pentagram intersection point I4
  → This means BOTH new lines pass through I4
  → But then they both pass through the same inner vertex!
- I1=X: similarly, both new lines pass through I1
- I0=X: both new lines pass through I0
- I4=I0: impossible (fixed pentagram)
- I1=I0: impossible (fixed pentagram)

So for quad → triangle: BOTH new lines must pass through the SAME
inner vertex of the pentagram!

Let's test this!
""")

# TEST: both lines through same inner vertex
print("=" * 60)
print("TESTING: BOTH LINES THROUGH SAME INNER VERTEX")
print("=" * 60)

from itertools import combinations
import random

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


def count_tris_full(V, extra_segs):
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
    return tri, sorted(face_info), G, faces


# For each inner vertex, try all pairs of lines through it
for iv_idx, iv in enumerate(inner):
    print(f"\nInner vertex I{iv_idx} = ({iv[0]:.6f}, {iv[1]:.6f})")

    # Generate lines through this inner vertex
    # Line through iv and various other points
    targets = []
    for i, v in enumerate(V):
        targets.append((v, f"V{i}"))
    for i, v in enumerate(inner):
        if i != iv_idx:
            targets.append((v, f"I{i}"))

    # Also add points on pentagram edges
    for i in range(5):
        for t in [0.25, 0.5, 0.75]:
            p1, p2 = V[i], V[(i+2)%5]
            mp = (p1[0]+t*(p2[0]-p1[0]), p1[1]+t*(p2[1]-p1[1]))
            targets.append((mp, f"L{i}@{t:.2f}"))

    # Also add random directions
    random.seed(42)
    for ri in range(20):
        angle = random.uniform(0, 2*math.pi)
        far = (iv[0] + 5*math.cos(angle), iv[1] + 5*math.sin(angle))
        targets.append((far, f"dir{angle:.2f}"))

    segs_through_iv = []
    seg_names = []
    for target, tname in targets:
        seg = clip(iv, target, V)
        if seg:
            # Check it's not a pentagram line
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
                # Deduplicate
                dup = False
                for existing, _ in segs_through_iv:
                    if (peq(seg[0], existing[0]) and peq(seg[1], existing[1])) or \
                       (peq(seg[0], existing[1]) and peq(seg[1], existing[0])):
                        dup = True; break
                if not dup:
                    segs_through_iv.append((seg, tname))

    print(f"  Unique lines through I{iv_idx}: {len(segs_through_iv)}")

    best_iv = 0
    for ci, cj in combinations(range(len(segs_through_iv)), 2):
        c, faces, _, _ = count_tris_full(V, [segs_through_iv[ci][0], segs_through_iv[cj][0]])
        if c > best_iv:
            best_iv = c
            print(f"  {c} tri: {segs_through_iv[ci][1]} + {segs_through_iv[cj][1]} -> faces={faces}")
            if c >= 10:
                print(f"  *** FOUND 10! ***")

    print(f"  Best through I{iv_idx}: {best_iv}")


# ============================================================
# Also test for non-regular pentagrams with lines through inner vertex
# ============================================================
print("\n" + "=" * 60)
print("NON-REGULAR: LINES THROUGH SAME INNER VERTEX")
print("=" * 60)

random.seed(2025)
best_nr = 0

for trial in range(20000):
    radii = [1.0 + random.gauss(0, 0.3) for _ in range(5)]
    angles = [math.pi/2 + 2*math.pi*k/5 + random.gauss(0, 0.2) for k in range(5)]
    Vn = [(r*math.cos(a), r*math.sin(a)) for r, a in zip(radii, angles)]

    pln = [make_line(Vn[i], Vn[(i+2)%5]) for i in range(5)]
    inn = []
    for i, j in combinations(range(5), 2):
        p = ix(pln[i], pln[j])
        if p and all(not peq(p, v) for v in Vn) and all(not peq(p, v) for v in inn):
            inn.append(p)
    if len(inn) != 5: continue

    # For each inner vertex, try pairs of lines through it
    for iv in inn:
        # Generate random directions through iv
        line_segs = []
        for _ in range(15):
            angle = random.uniform(0, math.pi)
            far = (iv[0] + 5*math.cos(angle), iv[1] + 5*math.sin(angle))
            seg = clip(iv, far, Vn)
            if seg:
                line = make_line(seg[0], seg[1])
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
                if not is_pent:
                    line_segs.append(seg)

        for ci, cj in combinations(range(len(line_segs)), 2):
            c, faces, _, _ = count_tris_full(Vn, [line_segs[ci], line_segs[cj]])
            if c > best_nr:
                best_nr = c
                print(f"  Trial {trial}: {c} tri! faces={faces}")
                if c >= 10:
                    print(f"  *** FOUND 10! ***")
                    print(f"  V={Vn}")

    if trial % 5000 == 0:
        print(f"  ... trial {trial}, best: {best_nr}")

print(f"\nBest non-regular (same vertex): {best_nr}")
