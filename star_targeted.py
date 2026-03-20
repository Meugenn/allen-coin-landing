"""
Targeted approach: In the 9-triangle config, there's 1 quad + 1 pentagon.
For 10 triangles, we need the quad to collapse to a triangle.
This happens when 2 of the quad's vertices merge (a concurrency).

Strategy: identify the quad's vertices, then find pentagram deformations
that cause 2 vertices to coincide.
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

def winding(pt, poly):
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

def point_on_segment(p, s1, s2):
    dx, dy = s2[0]-s1[0], s2[1]-s1[1]
    l2 = dx*dx + dy*dy
    if l2 < EPS*EPS:
        return abs(p[0]-s1[0]) < EPS and abs(p[1]-s1[1]) < EPS
    t = ((p[0]-s1[0])*dx + (p[1]-s1[1])*dy) / l2
    if t < -EPS or t > 1+EPS:
        return False
    proj = (s1[0]+t*dx, s1[1]+t*dy)
    return (p[0]-proj[0])**2 + (p[1]-proj[1])**2 < 1e-14

def seg_intersect(a1, a2, b1, b2):
    l1 = make_line(a1, a2)
    l2 = make_line(b1, b2)
    p = ix_lines(l1, l2)
    if p and point_on_segment(p, a1, a2) and point_on_segment(p, b1, b2):
        return p
    return None

class PlanarGraph:
    def __init__(self):
        self.verts = []
        self.adj = {}

    def add_vertex(self, p):
        for i, v in enumerate(self.verts):
            if abs(v[0]-p[0]) < EPS and abs(v[1]-p[1]) < EPS:
                return i
        idx = len(self.verts)
        self.verts.append(p)
        self.adj[idx] = []
        return idx

    def add_edge(self, i, j):
        if i == j: return
        if j not in self.adj[i]: self.adj[i].append(j)
        if i not in self.adj[j]: self.adj[j].append(i)

    def sort_adjacency(self):
        for v in self.adj:
            vx, vy = self.verts[v]
            self.adj[v].sort(key=lambda u: math.atan2(
                self.verts[u][1]-vy, self.verts[u][0]-vx))

    def find_faces(self):
        self.sort_adjacency()
        visited = set()
        faces = []
        for u in self.adj:
            for v in self.adj[u]:
                if (u, v) in visited: continue
                face = []
                cu, cv = u, v
                steps = 0
                while True:
                    visited.add((cu, cv))
                    face.append(cu)
                    neighbors = self.adj[cv]
                    if cu not in neighbors: break
                    idx = neighbors.index(cu)
                    next_v = neighbors[(idx - 1) % len(neighbors)]
                    cu, cv = cv, next_v
                    steps += 1
                    if cu == u and cv == v: break
                    if steps > 200: face = []; break
                if len(face) >= 3:
                    faces.append(face)
        return faces

    def face_area(self, face):
        area = 0
        n = len(face)
        for i in range(n):
            x1, y1 = self.verts[face[i]]
            x2, y2 = self.verts[face[(i+1)%n]]
            area += x1*y2 - x2*y1
        return area / 2


def build_graph(segments):
    G = PlanarGraph()
    for s1, s2 in segments:
        G.add_vertex(s1)
        G.add_vertex(s2)
    for i, j in combinations(range(len(segments)), 2):
        p = seg_intersect(segments[i][0], segments[i][1],
                          segments[j][0], segments[j][1])
        if p: G.add_vertex(p)
    for s1, s2 in segments:
        dx, dy = s2[0]-s1[0], s2[1]-s1[1]
        l2 = dx*dx + dy*dy
        pts = []
        for vi, v in enumerate(G.verts):
            if point_on_segment(v, s1, s2):
                t = ((v[0]-s1[0])*dx + (v[1]-s1[1])*dy) / l2 if l2 > EPS else 0
                pts.append((t, vi))
        pts.sort()
        filtered = []
        for t, vi in pts:
            if not filtered or filtered[-1][1] != vi:
                filtered.append((t, vi))
        for k in range(len(filtered)-1):
            G.add_edge(filtered[k][1], filtered[k+1][1])
    return G


def count_tris(V, extra_segs):
    segs = [(V[i], V[(i+2)%5]) for i in range(5)] + extra_segs
    star = [V[0], V[2], V[4], V[1], V[3]]
    G = build_graph(segs)
    faces = G.find_faces()
    count = 0
    face_info = []
    for face in faces:
        n = len(face)
        cx = sum(G.verts[v][0] for v in face) / n
        cy = sum(G.verts[v][1] for v in face) / n
        area = G.face_area(face)
        inside = winding((cx, cy), star) != 0
        if inside and n == 3:
            count += 1
        face_info.append((n, area, inside))
    return count, face_info, G, faces


# === ANALYZE THE 9-TRIANGLE CONFIG ===
print("="*60)
print("DETAILED ANALYSIS OF 9-TRIANGLE CONFIG")
print("="*60)

V = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

# Inner vertices
inner = []
for i, j in combinations(range(5), 2):
    p = ix_lines(pent_lines[i], pent_lines[j])
    if p and all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in V):
        if all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in inner):
            inner.append(p)

sp = V + inner
names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner))]

# Use V2--I4 + V4--I0 config
# First, find I4 and I0
print("Special points:")
for i, p in enumerate(sp):
    print(f"  {names[i]}: ({p[0]:.6f}, {p[1]:.6f})")

# Create segment through V2 and I4
def clip_to_star(p1, p2, V):
    """Clip infinite line through p1, p2 to star boundary."""
    line = make_line(p1, p2)
    pts = []
    for i in range(5):
        s1, s2 = V[i], V[(i+2)%5]
        sl = make_line(s1, s2)
        p = ix_lines(line, sl)
        if p and point_on_segment(p, s1, s2):
            is_dup = any(abs(p[0]-q[0])<EPS and abs(p[1]-q[1])<EPS for q in pts)
            if not is_dup:
                pts.append(p)
    if len(pts) < 2:
        return None
    a, b, c = line
    pts.sort(key=lambda p: -b*p[0] + a*p[1])
    return (pts[0], pts[-1])

# V2--I4 and V4--I0
seg1 = clip_to_star(sp[2], sp[9], V)  # V2--I4
seg2 = clip_to_star(sp[4], sp[5], V)  # V4--I0

print(f"\nConfig: V2--I4 + V4--I0")
print(f"Seg1: {seg1}")
print(f"Seg2: {seg2}")

c, face_info, G, faces = count_tris(V, [seg1, seg2])
print(f"\nTriangles: {c}")
print(f"Total faces: {len(faces)}")
print(f"Vertices: {len(G.verts)}")

# Show ALL vertices
print("\nAll vertices:")
for i, v in enumerate(G.verts):
    # Find name
    name = "NEW"
    for j, s in enumerate(sp):
        if abs(v[0]-s[0])<EPS and abs(v[1]-s[1])<EPS:
            name = names[j]
            break
    print(f"  v{i}: ({v[0]:.6f}, {v[1]:.6f}) [{name}]")

# Show ALL faces with vertex details
print("\nAll faces:")
for i, face in enumerate(faces):
    n = len(face)
    area = G.face_area(face)
    cx = sum(G.verts[v][0] for v in face) / n
    cy = sum(G.verts[v][1] for v in face) / n
    star = [V[0], V[2], V[4], V[1], V[3]]
    inside = winding((cx, cy), star) != 0
    kind = "TRI" if n == 3 else f"{n}-gon"
    loc = "INSIDE" if inside else "outside"
    print(f"  Face {i}: {kind}, area={area:.6f}, {loc}")
    for vi in face:
        v = G.verts[vi]
        name = "new"
        for j, s in enumerate(sp):
            if abs(v[0]-s[0])<EPS and abs(v[1]-s[1])<EPS:
                name = names[j]
                break
        print(f"    v{vi}: ({v[0]:.6f}, {v[1]:.6f}) [{name}]")


# === NOW: try to make the quad collapse to triangle ===
# The quad has 4 vertices. If we adjust the pentagram so 2 merge...
# First identify the quad face
print("\n" + "="*60)
print("SEARCHING FOR PENTAGRAM THAT COLLAPSES QUAD TO TRIANGLE")
print("="*60)

# For a parametric pentagram, we vary the outer vertex positions
# and track when the quad face becomes a triangle

# Let's parametrize: V0 at (0, r0), others at normal angles but radius r_k
# and search for r_k values that give 10 triangles

random.seed(42)
best = 0
best_cfg = None

# Strategy: move vertices individually and check
for trial in range(100000):
    # Parametric pentagram: move each vertex radially
    radii = [1.0 + random.gauss(0, 0.3) for _ in range(5)]
    # Also try angle perturbation
    angles = [math.pi/2 + 2*math.pi*k/5 + random.gauss(0, 0.15) for k in range(5)]
    Vn = [(r*math.cos(a), r*math.sin(a)) for r, a in zip(radii, angles)]

    pl = [make_line(Vn[i], Vn[(i+2)%5]) for i in range(5)]

    # Compute inner vertices
    inn = []
    for i, j in combinations(range(5), 2):
        p = ix_lines(pl[i], pl[j])
        if p and all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in Vn):
            if all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in inn):
                inn.append(p)
    if len(inn) != 5:
        continue

    spn = Vn + inn

    # Try ALL pairs of "axis" type lines (outer vertex through non-adjacent inner vertex)
    # These are the ones that gave 9 in regular case
    axis_pairs = []
    for vi in range(5):
        for ii in range(5):
            # Check that inner vertex ii is NOT on a line through outer vertex vi
            # (i.e., it's a "non-adjacent" inner vertex)
            on_line = False
            for pl_idx in range(5):
                if (abs(pl[pl_idx][0]*spn[vi][0] + pl[pl_idx][1]*spn[vi][1] + pl[pl_idx][2]) < 0.01 and
                    abs(pl[pl_idx][0]*inn[ii][0] + pl[pl_idx][1]*inn[ii][1] + pl[pl_idx][2]) < 0.01):
                    on_line = True
                    break
            if not on_line:
                seg = clip_to_star(spn[vi], inn[ii], Vn)
                if seg:
                    axis_pairs.append((vi, ii, seg))

    # Also try inner-inner pairs
    inner_pairs = []
    for i, j in combinations(range(5), 2):
        line = make_line(inn[i], inn[j])
        is_pent = False
        for pli in pl:
            a1, b1, c1 = line
            a2, b2, c2 = pli
            n1 = math.sqrt(a1*a1+b1*b1)
            n2 = math.sqrt(a2*a2+b2*b2)
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
                    is_pent = True
                    break
        if not is_pent:
            seg = clip_to_star(inn[i], inn[j], Vn)
            if seg:
                inner_pairs.append((i, j, seg))

    all_candidates = [(s, f"V{vi}->I{ii}") for vi, ii, s in axis_pairs] + \
                     [(s, f"I{i}->I{j}") for i, j, s in inner_pairs]

    for ci in range(len(all_candidates)):
        for cj in range(ci+1, len(all_candidates)):
            seg_a = all_candidates[ci][0]
            seg_b = all_candidates[cj][0]
            c, _, _, _ = count_tris(Vn, [seg_a, seg_b])
            if c > best:
                best = c
                best_cfg = (Vn, seg_a, seg_b, all_candidates[ci][1], all_candidates[cj][1])
                print(f"  Trial {trial}: {c} triangles! ({all_candidates[ci][1]} + {all_candidates[cj][1]})")
                if c >= 10:
                    print(f"  *** FOUND 10! ***")
                    print(f"  V = {Vn}")
                    print(f"  seg1 = {seg_a}")
                    print(f"  seg2 = {seg_b}")

    if trial % 5000 == 0 and trial > 0:
        print(f"  ... trial {trial}, best: {best}")

print(f"\nBest: {best}")
if best_cfg:
    print(f"Config: {best_cfg[3]} + {best_cfg[4]}")
