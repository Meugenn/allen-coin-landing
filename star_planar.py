"""
Proper planar graph face detection for pentagram + 2 line segments.

The triple-checking approach misses faces that are bounded by segment
endpoints (not line-line intersections). This approach builds the full
planar graph and traverses faces properly.
"""

import math
from itertools import combinations
import random

EPS = 1e-8

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

def winding(pt, poly):
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

def point_on_segment(p, s1, s2):
    dx, dy = s2[0]-s1[0], s2[1]-s1[1]
    l2 = dx*dx + dy*dy
    if l2 < EPS*EPS:
        return abs(p[0]-s1[0]) < EPS and abs(p[1]-s1[1]) < EPS
    t = ((p[0]-s1[0])*dx + (p[1]-s1[1])*dy) / l2
    if t < -EPS or t > 1+EPS:
        return False
    proj = (s1[0]+t*dx, s1[1]+t*dy)
    return (p[0]-proj[0])**2 + (p[1]-proj[1])**2 < EPS*EPS*100


def seg_intersect(a1, a2, b1, b2):
    l1 = make_line(a1, a2)
    l2 = make_line(b1, b2)
    p = ix_lines(l1, l2)
    if p and point_on_segment(p, a1, a2) and point_on_segment(p, b1, b2):
        return p
    return None


class PlanarGraph:
    """Build planar graph from segments and find faces."""

    def __init__(self):
        self.verts = []  # (x, y)
        self.adj = {}    # vertex_idx -> sorted list of neighbor indices

    def add_vertex(self, p):
        for i, v in enumerate(self.verts):
            if abs(v[0]-p[0]) < EPS and abs(v[1]-p[1]) < EPS:
                return i
        idx = len(self.verts)
        self.verts.append(p)
        self.adj[idx] = []
        return idx

    def add_edge(self, i, j):
        if i == j:
            return
        if j not in self.adj[i]:
            self.adj[i].append(j)
        if i not in self.adj[j]:
            self.adj[j].append(i)

    def sort_adjacency(self):
        """Sort neighbors of each vertex by angle."""
        for v in self.adj:
            vx, vy = self.verts[v]
            self.adj[v].sort(key=lambda u: math.atan2(
                self.verts[u][1]-vy, self.verts[u][0]-vx))

    def find_faces(self):
        """Find all faces using clockwise face traversal."""
        self.sort_adjacency()
        visited = set()
        faces = []

        # Collect all directed edges
        directed_edges = []
        for u in self.adj:
            for v in self.adj[u]:
                directed_edges.append((u, v))

        for start_u, start_v in directed_edges:
            if (start_u, start_v) in visited:
                continue

            face = []
            u, v = start_u, start_v
            steps = 0
            while True:
                visited.add((u, v))
                face.append(u)

                # At vertex v, find u in neighbor list, take previous (clockwise)
                neighbors = self.adj[v]
                if u not in neighbors:
                    break
                idx = neighbors.index(u)
                next_v = neighbors[(idx - 1) % len(neighbors)]

                u, v = v, next_v
                steps += 1
                if u == start_u and v == start_v:
                    break
                if steps > 200:
                    face = []
                    break

            if len(face) >= 3:
                faces.append(face)

        return faces

    def face_area(self, face):
        """Signed area of a face (positive = CCW, negative = CW)."""
        area = 0
        n = len(face)
        for i in range(n):
            x1, y1 = self.verts[face[i]]
            x2, y2 = self.verts[face[(i+1)%n]]
            area += x1*y2 - x2*y1
        return area / 2

    def count_triangles_inside_star(self, star_polygon):
        """Count triangular faces inside the star."""
        faces = self.find_faces()
        count = 0
        tri_faces = []
        for face in faces:
            if len(face) != 3:
                continue
            # Centroid
            cx = sum(self.verts[v][0] for v in face) / 3
            cy = sum(self.verts[v][1] for v in face) / 3
            if winding((cx, cy), star_polygon) != 0:
                count += 1
                tri_faces.append(face)
        return count, tri_faces, faces


def build_graph_from_segments(segments):
    """Build planar graph from a list of segments."""
    G = PlanarGraph()

    # Add all endpoints
    for s1, s2 in segments:
        G.add_vertex(s1)
        G.add_vertex(s2)

    # Find all intersection points between segments
    for i, j in combinations(range(len(segments)), 2):
        p = seg_intersect(segments[i][0], segments[i][1],
                          segments[j][0], segments[j][1])
        if p:
            G.add_vertex(p)

    # For each segment, find all vertices on it and create edges
    for s1, s2 in segments:
        dx, dy = s2[0]-s1[0], s2[1]-s1[1]
        l2 = dx*dx + dy*dy

        pts_on_seg = []
        for vi, v in enumerate(G.verts):
            if point_on_segment(v, s1, s2):
                t = ((v[0]-s1[0])*dx + (v[1]-s1[1])*dy) / l2 if l2 > EPS else 0
                pts_on_seg.append((t, vi))

        pts_on_seg.sort()
        # Remove duplicates
        filtered = []
        for t, vi in pts_on_seg:
            if not filtered or filtered[-1][1] != vi:
                filtered.append((t, vi))

        for k in range(len(filtered)-1):
            G.add_edge(filtered[k][1], filtered[k+1][1])

    return G


def count_triangles_for_config(V, extra_segs):
    """Count triangular faces for a pentagram with extra segments."""
    # Pentagram segments
    pent_segs = [(V[i], V[(i+2)%5]) for i in range(5)]
    all_segs = pent_segs + extra_segs
    star = [V[0], V[2], V[4], V[1], V[3]]

    G = build_graph_from_segments(all_segs)
    count, tri_faces, all_faces = G.count_triangles_inside_star(star)
    return count, G, tri_faces, all_faces


# === TEST WITH REGULAR PENTAGRAM ===
print("="*60)
print("REGULAR PENTAGRAM - PLANAR GRAPH APPROACH")
print("="*60)

V = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
star = [V[0], V[2], V[4], V[1], V[3]]

# Base pentagram
c, G, tri, faces = count_triangles_for_config(V, [])
print(f"Base: {c} triangles, {len(faces)} total faces")
print(f"  Vertices: {len(G.verts)}, Edges: {sum(len(v) for v in G.adj.values())//2}")
face_sizes = [len(f) for f in faces]
from collections import Counter
print(f"  Face sizes: {Counter(face_sizes)}")

# Inner vertices
inner = []
for i, j in combinations(range(5), 2):
    p = ix_lines(pent_lines[i], pent_lines[j])
    if p and all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in V):
        if all(abs(p[0]-v[0])>EPS or abs(p[1]-v[1])>EPS for v in inner):
            inner.append(p)

all_sp = V + inner
pt_names = [f"V{i}" for i in range(5)] + [f"I{i}" for i in range(len(inner))]
print(f"\nSpecial points: {len(all_sp)}")

# Candidate lines
cand_lines_inf = []
cand_names = []
for i, j in combinations(range(len(all_sp)), 2):
    line = make_line(all_sp[i], all_sp[j])
    is_pent = False
    for pl in pent_lines:
        # Check if same line
        def norm_line(l):
            a, b, c = l
            n = math.sqrt(a*a+b*b)
            if n < EPS: return l
            a, b, c = a/n, b/n, c/n
            if abs(a) > EPS:
                if a < 0: a,b,c = -a,-b,-c
            elif b < 0: a,b,c = -a,-b,-c
            return (a,b,c)
        if all(abs(a-b) < 1e-6 for a, b in zip(norm_line(line), norm_line(pl))):
            is_pent = True
            break
    if is_pent:
        continue
    is_dup = False
    for cl in cand_lines_inf:
        if all(abs(a-b) < 1e-6 for a, b in zip(norm_line(line), norm_line(cl))):
            is_dup = True
            break
    if not is_dup:
        cand_lines_inf.append(line)
        cand_names.append(f"{pt_names[i]}--{pt_names[j]}")

print(f"Candidate lines: {len(cand_lines_inf)}")

# Clip lines to star boundary and create segments
def clip_to_star(line, V):
    """Clip infinite line to star boundary."""
    pts = []
    for i in range(5):
        p1, p2 = V[i], V[(i+2)%5]
        seg_line = make_line(p1, p2)
        p = ix_lines(line, seg_line)
        if p and point_on_segment(p, p1, p2):
            is_dup = any(abs(p[0]-q[0])<EPS and abs(p[1]-q[1])<EPS for q in pts)
            if not is_dup:
                pts.append(p)
    if len(pts) < 2:
        return None
    a, b, c = line
    pts.sort(key=lambda p: -b*p[0] + a*p[1])
    return (pts[0], pts[-1])

cand_segs = []
cand_seg_names = []
for i, line in enumerate(cand_lines_inf):
    seg = clip_to_star(line, V)
    if seg:
        cand_segs.append(seg)
        cand_seg_names.append(cand_names[i])

print(f"Candidate segments (clipped): {len(cand_segs)}")

# Test all pairs
print("\nTesting all pairs:")
results = []
for ci, cj in combinations(range(len(cand_segs)), 2):
    c, _, _, faces = count_triangles_for_config(V, [cand_segs[ci], cand_segs[cj]])
    results.append((c, ci, cj))

results.sort(reverse=True)
print("Top 15:")
for c, ci, cj in results[:15]:
    print(f"  {c}: {cand_seg_names[ci]} + {cand_seg_names[cj]}")

dist = Counter(c for c, _, _ in results)
print(f"Distribution: {dict(sorted(dist.items()))}")

# === ALSO TRY SEGMENTS TO MIDPOINTS OF EDGES ===
print("\n" + "="*60)
print("EXTENDED CANDIDATE SET: midpoints of edges, etc.")
print("="*60)

# Add midpoints of pentagram segments as additional special points
extended_sp = list(all_sp)
ext_names = list(pt_names)
for i in range(5):
    p1, p2 = V[i], V[(i+2)%5]
    mid = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
    extended_sp.append(mid)
    ext_names.append(f"M{i}")

# Also add midpoints of inner pentagon edges
for i in range(len(inner)):
    for j in range(i+1, len(inner)):
        # Check if adjacent in inner pentagon
        p1, p2 = inner[i], inner[j]
        # Check if this segment is a side of inner pentagon (on a pentagram line)
        is_side = False
        for pl in pent_lines:
            a, b, c = pl
            if abs(a*p1[0]+b*p1[1]+c) < EPS and abs(a*p2[0]+b*p2[1]+c) < EPS:
                is_side = True
                break
        if is_side:
            mid = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
            extended_sp.append(mid)
            ext_names.append(f"MI{i}{j}")

print(f"Extended special points: {len(extended_sp)}")

# Generate more candidates
ext_cand_segs = []
ext_cand_names = []
for i, j in combinations(range(len(extended_sp)), 2):
    line = make_line(extended_sp[i], extended_sp[j])
    is_pent = any(all(abs(a-b) < 1e-6 for a, b in
        zip(norm_line(line), norm_line(pl))) for pl in pent_lines)
    if is_pent:
        continue
    is_dup = any(all(abs(a-b) < 1e-6 for a, b in
        zip(norm_line(line), norm_line(make_line(s[0], s[1]))))
        for s in ext_cand_segs) if ext_cand_segs else False
    if is_dup:
        continue
    seg = clip_to_star(line, V)
    if seg:
        ext_cand_segs.append(seg)
        ext_cand_names.append(f"{ext_names[i]}--{ext_names[j]}")

print(f"Extended candidates: {len(ext_cand_segs)}")

# Test all pairs
ext_results = []
for ci, cj in combinations(range(len(ext_cand_segs)), 2):
    c, _, _, _ = count_triangles_for_config(V, [ext_cand_segs[ci], ext_cand_segs[cj]])
    ext_results.append((c, ci, cj))

ext_results.sort(reverse=True)
print("Top 15 (extended):")
for c, ci, cj in ext_results[:15]:
    print(f"  {c}: {ext_cand_names[ci]} + {ext_cand_names[cj]}")

ext_dist = Counter(c for c, _, _ in ext_results)
print(f"Distribution: {dict(sorted(ext_dist.items()))}")


# === NON-REGULAR PENTAGRAMS ===
print("\n" + "="*60)
print("NON-REGULAR PENTAGRAM SEARCH")
print("="*60)

random.seed(42)
best_nr = 0

for trial in range(3000):
    Vn = []
    for k in range(5):
        base_angle = math.pi/2 + 2*math.pi*k/5
        angle = base_angle + random.gauss(0, 0.25)
        radius = 0.5 + random.random() * 1.5
        Vn.append((radius * math.cos(angle), radius * math.sin(angle)))

    plines = [make_line(Vn[i], Vn[(i+2)%5]) for i in range(5)]

    # Special points
    sp = list(Vn)
    for i, j in combinations(range(5), 2):
        p = ix_lines(plines[i], plines[j])
        if p and all(abs(p[0]-s[0])>EPS or abs(p[1]-s[1])>EPS for s in sp):
            sp.append(p)
    if len(sp) < 10:
        continue

    # Candidate segments
    csegs = []
    for i, j in combinations(range(len(sp)), 2):
        line = make_line(sp[i], sp[j])
        is_pent = any(all(abs(a-b) < 1e-6 for a, b in
            zip(norm_line(line), norm_line(pl))) for pl in plines)
        if is_pent:
            continue
        is_dup = any(all(abs(a-b) < 1e-6 for a, b in
            zip(norm_line(line), norm_line(make_line(s[0], s[1]))))
            for s in csegs) if csegs else False
        if is_dup:
            continue
        seg = clip_to_star(line, Vn)
        if seg:
            csegs.append(seg)

    if len(csegs) < 2:
        continue

    for ci, cj in combinations(range(min(len(csegs), 20)), 2):
        c, _, _, _ = count_triangles_for_config(Vn, [csegs[ci], csegs[cj]])
        if c > best_nr:
            best_nr = c
            print(f"  Trial {trial}: {c} triangles!")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  V = {Vn}")
                print(f"  seg1 = {csegs[ci]}")
                print(f"  seg2 = {csegs[cj]}")

    if trial % 500 == 0:
        print(f"  ... trial {trial}, best: {best_nr}")

print(f"\nBest non-regular: {best_nr}")


# === ANALYSIS: Why is 9 the max? ===
print("\n" + "="*60)
print("ANALYSIS OF 9-TRIANGLE CONFIGURATION")
print("="*60)

# Get the best config
c_best, ci_best, cj_best = results[0]
print(f"Best config: {cand_seg_names[ci_best]} + {cand_seg_names[cj_best]}")
c, G, tri_faces, all_faces = count_triangles_for_config(
    V, [cand_segs[ci_best], cand_segs[cj_best]])

print(f"Triangles: {c}")
print(f"Total faces: {len(all_faces)}")
print(f"Vertices: {len(G.verts)}")
print(f"Edges: {sum(len(v) for v in G.adj.values())//2}")

# Show all faces with their sizes and whether inside star
for i, face in enumerate(all_faces):
    n = len(face)
    cx = sum(G.verts[v][0] for v in face) / n
    cy = sum(G.verts[v][1] for v in face) / n
    area = G.face_area(face)
    inside = winding((cx, cy), star) != 0
    is_tri = "TRI" if n == 3 else "   "
    in_star = "IN" if inside else "  "
    print(f"  Face {i}: {n} edges, area={area:.4f}, {is_tri} {in_star}")
