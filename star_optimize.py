"""
Fresh approach: Gradient-based optimization on non-regular pentagrams.

Key insight we missed: instead of random search, use hill-climbing.
Start from 9-triangle configs and systematically deform the star
to break the concurrency that creates the quadrilateral.

Also: try lines NOT through special points - parameterize lines
by angle and offset, and optimize those too.
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


class PlanarGraph:
    def __init__(self):
        self.verts = []
        self.adj = {}

    def add_v(self, p):
        for i, v in enumerate(self.verts):
            if abs(v[0]-p[0]) < EPS and abs(v[1]-p[1]) < EPS:
                return i
        idx = len(self.verts)
        self.verts.append(p)
        self.adj[idx] = []
        return idx

    def add_e(self, i, j):
        if i == j: return
        if j not in self.adj[i]: self.adj[i].append(j)
        if i not in self.adj[j]: self.adj[j].append(i)

    def find_faces(self):
        for v in self.adj:
            vx, vy = self.verts[v]
            self.adj[v].sort(key=lambda u: math.atan2(
                self.verts[u][1]-vy, self.verts[u][0]-vx))

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
                    nb = self.adj[cv]
                    if cu not in nb: break
                    idx = nb.index(cu)
                    nv = nb[(idx - 1) % len(nb)]
                    cu, cv = cv, nv
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


def build_and_count(V, extra_segs):
    """Build planar graph from pentagram + extra segments, count triangular faces inside star."""
    segs = [(V[i], V[(i+2)%5]) for i in range(5)] + list(extra_segs)
    star = [V[0], V[2], V[4], V[1], V[3]]

    G = PlanarGraph()

    # Add segment endpoints
    for s1, s2 in segs:
        G.add_v(s1)
        G.add_v(s2)

    # Find all pairwise segment intersections
    for i, j in combinations(range(len(segs)), 2):
        l1 = make_line(segs[i][0], segs[i][1])
        l2 = make_line(segs[j][0], segs[j][1])
        p = ix(l1, l2)
        if p and point_on_seg(p, segs[i][0], segs[i][1]) and point_on_seg(p, segs[j][0], segs[j][1]):
            G.add_v(p)

    # Build edges along each segment
    for s1, s2 in segs:
        dx, dy = s2[0]-s1[0], s2[1]-s1[1]
        l2 = dx*dx + dy*dy
        pts = []
        for vi, v in enumerate(G.verts):
            if point_on_seg(v, s1, s2):
                t = ((v[0]-s1[0])*dx + (v[1]-s1[1])*dy) / l2 if l2 > EPS else 0
                pts.append((t, vi))
        pts.sort()
        filtered = []
        for t, vi in pts:
            if not filtered or filtered[-1][1] != vi:
                filtered.append((t, vi))
        for k in range(len(filtered)-1):
            G.add_e(filtered[k][1], filtered[k+1][1])

    faces = G.find_faces()
    tri_count = 0
    face_sizes = []
    for face in faces:
        n = len(face)
        cx = sum(G.verts[v][0] for v in face) / n
        cy = sum(G.verts[v][1] for v in face) / n
        area = G.face_area(face)
        inside = wn((cx, cy), star) != 0
        if inside and area > 1e-8:
            face_sizes.append(n)
            if n == 3:
                tri_count += 1
    return tri_count, face_sizes


def clip_to_star(p1, p2, V):
    """Clip a line (defined by 2 points) to the star boundary segments."""
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


def make_star(radii, angles):
    """Create pentagram vertices from radii and angles."""
    return [(r*math.cos(a), r*math.sin(a)) for r, a in zip(radii, angles)]


def line_from_angle_offset(angle, offset):
    """Create line (a,b,c) from angle and perpendicular offset."""
    return (math.cos(angle), math.sin(angle), -offset)


# ============================================================
# APPROACH 1: Hill-climbing from known 9-triangle configs
# ============================================================
print("=" * 60)
print("HILL-CLIMBING OPTIMIZATION")
print("=" * 60)

random.seed(42)
best_global = 0
best_global_cfg = None

# Start from regular pentagram
base_angles = [math.pi/2 + 2*math.pi*k/5 for k in range(5)]
base_radii = [1.0] * 5

# Get inner vertices of regular pentagram
V0 = make_star(base_radii, base_angles)
plines0 = [make_line(V0[i], V0[(i+2)%5]) for i in range(5)]
inner0 = []
for i, j in combinations(range(5), 2):
    p = ix(plines0[i], plines0[j])
    if p and all(not peq(p, v) for v in V0) and all(not peq(p, v) for v in inner0):
        inner0.append(p)

# Known good line configs: "axis" type (outer vertex to opposite inner vertex)
# V2->I4 and V4->I0 gave 9
# Let's parameterize lines by angle+offset and optimize

def evaluate(radii, angles, line1_params, line2_params):
    """Evaluate a configuration. Returns triangle count."""
    V = make_star(radii, angles)

    l1 = line_from_angle_offset(line1_params[0], line1_params[1])
    l2 = line_from_angle_offset(line2_params[0], line2_params[1])

    seg1 = clip_to_star_line(l1, V)
    seg2 = clip_to_star_line(l2, V)

    if seg1 is None or seg2 is None:
        return 0, []

    return build_and_count(V, [seg1, seg2])


def clip_to_star_line(line, V):
    """Clip line (a,b,c) to star boundary."""
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


# First: find which angle+offset parameterization corresponds to V2->I4 and V4->I0
# V2->I4 line
p1, p2 = V0[2], inner0[4]
# angle of the line
dx, dy = p2[0]-p1[0], p2[1]-p1[1]
line_angle_1 = math.atan2(-dx, dy)  # normal angle
line_offset_1 = math.cos(line_angle_1)*p1[0] + math.sin(line_angle_1)*p1[1]

# V4->I0 line
p1, p2 = V0[4], inner0[0]
dx, dy = p2[0]-p1[0], p2[1]-p1[1]
line_angle_2 = math.atan2(-dx, dy)
line_offset_2 = math.cos(line_angle_2)*p1[0] + math.sin(line_angle_2)*p1[1]

# Verify starting config
c0, fs0 = evaluate(base_radii, base_angles,
                    (line_angle_1, line_offset_1),
                    (line_angle_2, line_offset_2))
print(f"Starting config: {c0} triangles, faces: {sorted(fs0)}")

# Hill-climbing: perturb all 14 parameters (5 radii, 5 angles, 2x2 line params)
def perturb(params, temperature):
    new = list(params)
    idx = random.randint(0, len(new)-1)
    new[idx] += random.gauss(0, temperature)
    # Keep radii positive
    for i in range(5):
        if new[i] < 0.2: new[i] = 0.2
    return new

def pack_params(radii, angles, l1, l2):
    return list(radii) + list(angles) + [l1[0], l1[1], l2[0], l2[1]]

def unpack_params(params):
    radii = params[:5]
    angles = params[5:10]
    l1 = (params[10], params[11])
    l2 = (params[12], params[13])
    return radii, angles, l1, l2

params = pack_params(base_radii, base_angles,
                     (line_angle_1, line_offset_1),
                     (line_angle_2, line_offset_2))

best = c0
best_params = list(params)

# Simulated annealing
for iteration in range(500000):
    temp = 0.3 * (1 - iteration / 500000) + 0.01

    new_params = perturb(params, temp)
    radii, angles, l1, l2 = unpack_params(new_params)

    try:
        c, fs = evaluate(radii, angles, l1, l2)
    except:
        continue

    # Accept if better, or with probability if same
    if c > best or (c == best and random.random() < 0.3):
        params = new_params
        if c > best:
            best = c
            best_params = list(new_params)
            radii, angles, l1, l2 = unpack_params(best_params)
            print(f"  Iter {iteration}: {c} triangles! faces={sorted(fs)} temp={temp:.4f}")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                V = make_star(radii, angles)
                print(f"  V = {V}")
                print(f"  Line1: angle={l1[0]:.6f}, offset={l1[1]:.6f}")
                print(f"  Line2: angle={l2[0]:.6f}, offset={l2[1]:.6f}")
    elif c >= best - 1 and random.random() < math.exp(-1/max(temp, 0.01)):
        params = new_params  # accept slightly worse for exploration

    if iteration % 100000 == 0:
        print(f"  ... iter {iteration}, best: {best}, temp: {temp:.4f}")

print(f"\nHill-climbing best: {best}")


# ============================================================
# APPROACH 2: Exhaustive grid search on line parameters
# For regular pentagram, sweep all line angles/offsets
# ============================================================
print("\n" + "=" * 60)
print("GRID SEARCH ON LINE PARAMETERS")
print("=" * 60)

V_reg = make_star([1.0]*5, base_angles)
best_grid = 0

# Grid of angles and offsets
n_angle = 60
n_offset = 30

angles_grid = [math.pi * i / n_angle for i in range(n_angle)]
offsets_grid = [0.8 * (2*i/(n_offset-1) - 1) for i in range(n_offset)]

# Precompute all valid line segments
valid_lines = []
for ai, ang in enumerate(angles_grid):
    for oi, off in enumerate(offsets_grid):
        line = line_from_angle_offset(ang, off)
        seg = clip_to_star_line(line, V_reg)
        if seg:
            valid_lines.append((ang, off, seg))

print(f"Valid line segments: {len(valid_lines)}")

# Test all pairs
count_checked = 0
for i in range(len(valid_lines)):
    for j in range(i+1, len(valid_lines)):
        seg1 = valid_lines[i][2]
        seg2 = valid_lines[j][2]
        c, fs = build_and_count(V_reg, [seg1, seg2])
        count_checked += 1
        if c > best_grid:
            best_grid = c
            print(f"  {c} triangles: angle1={valid_lines[i][0]:.3f} off1={valid_lines[i][1]:.3f} "
                  f"angle2={valid_lines[j][0]:.3f} off2={valid_lines[j][1]:.3f} faces={sorted(fs)}")
            if c >= 10:
                print(f"  *** FOUND 10! ***")

    if i % 100 == 0 and i > 0:
        print(f"  ... line {i}/{len(valid_lines)}, checked {count_checked}, best: {best_grid}")

print(f"\nGrid search best: {best_grid}")
print(f"Total pairs checked: {count_checked}")


# ============================================================
# APPROACH 3: Non-regular pentagrams with grid search on lines
# ============================================================
print("\n" + "=" * 60)
print("NON-REGULAR + GRID SEARCH")
print("=" * 60)

random.seed(999)
best_nr = 0

for trial in range(2000):
    # Generate non-regular pentagram
    radii = [1.0 + random.uniform(-0.4, 0.4) for _ in range(5)]
    angles = [base_angles[k] + random.uniform(-0.3, 0.3) for k in range(5)]
    V = make_star(radii, angles)

    # Compute inner vertices
    plines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]
    inner = []
    for i, j in combinations(range(5), 2):
        p = ix(plines[i], plines[j])
        if p and all(not peq(p, v) for v in V) and all(not peq(p, v) for v in inner):
            inner.append(p)
    if len(inner) != 5:
        continue

    sp = V + inner

    # Try lines through ALL pairs of special points
    cand_segs = []
    for i, j in combinations(range(len(sp)), 2):
        line = make_line(sp[i], sp[j])
        # Skip if same as pentagram line
        is_pent = False
        for pl in plines:
            a1, b1, c1 = line
            a2, b2, c2 = pl
            n1 = math.sqrt(a1*a1+b1*b1)
            n2 = math.sqrt(a2*a2+b2*b2)
            if n1 > EPS and n2 > EPS:
                a1,b1,c1 = a1/n1,b1/n1,c1/n1
                a2,b2,c2 = a2/n2,b2/n2,c2/n2
                if a1*n1 < 0 or (abs(a1) < EPS and b1 < 0): a1,b1,c1 = -a1,-b1,-c1
                if a2*n2 < 0 or (abs(a2) < EPS and b2 < 0): a2,b2,c2 = -a2,-b2,-c2
                # Normalize sign
                if abs(a1) > EPS:
                    if a1 < 0: a1,b1,c1 = -a1,-b1,-c1
                elif b1 < 0: a1,b1,c1 = -a1,-b1,-c1
                if abs(a2) > EPS:
                    if a2 < 0: a2,b2,c2 = -a2,-b2,-c2
                elif b2 < 0: a2,b2,c2 = -a2,-b2,-c2
                if abs(a1-a2)<1e-4 and abs(b1-b2)<1e-4 and abs(c1-c2)<1e-4:
                    is_pent = True
                    break
        if is_pent:
            continue

        seg = clip_to_star(sp[i], sp[j], V)
        if seg:
            # Deduplicate
            dup = False
            for existing in cand_segs:
                if (peq(seg[0], existing[0]) and peq(seg[1], existing[1])) or \
                   (peq(seg[0], existing[1]) and peq(seg[1], existing[0])):
                    dup = True
                    break
            if not dup:
                cand_segs.append(seg)

    # Test all pairs
    for ci, cj in combinations(range(len(cand_segs)), 2):
        c, fs = build_and_count(V, [cand_segs[ci], cand_segs[cj]])
        if c > best_nr:
            best_nr = c
            print(f"  Trial {trial}: {c} triangles! faces={sorted(fs)}")
            if c >= 10:
                print(f"  *** FOUND 10! ***")
                print(f"  Radii: {radii}")
                print(f"  Angles: {angles}")
                print(f"  Seg1: {cand_segs[ci]}")
                print(f"  Seg2: {cand_segs[cj]}")

    if trial % 500 == 0:
        print(f"  ... trial {trial}, best: {best_nr}")

print(f"\nNon-regular + grid best: {best_nr}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("OPTIMIZATION SUMMARY")
print("=" * 60)
print(f"Hill-climbing: {best}")
print(f"Grid search (regular): {best_grid}")
print(f"Non-regular + special points: {best_nr}")
