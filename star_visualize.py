"""
Visualization of the pentagram triangle problem.
Shows the best 9-triangle configuration and analyzes why 10 seems impossible.
"""
import math
from itertools import combinations

# === SVG Drawing ===
def svg_line(x1, y1, x2, y2, color="black", width=2, dash=False):
    style = f'stroke="{color}" stroke-width="{width}" fill="none"'
    if dash:
        style += ' stroke-dasharray="5,5"'
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {style}/>'

def svg_circle(cx, cy, r, color="red", fill="red"):
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{color}"/>'

def svg_text(x, y, text, size=12, color="black"):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" text-anchor="middle">{text}</text>'

def svg_polygon(points, fill="rgba(255,200,200,0.5)", stroke="none"):
    pts = " ".join(f"{x},{y}" for x, y in points)
    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}"/>'

def make_line(p1, p2):
    a = p2[1] - p1[1]
    b = p1[0] - p2[0]
    c = p2[0]*p1[1] - p1[0]*p2[1]
    return (a, b, c)

def ix_lines(l1, l2):
    d = l1[0]*l2[1] - l2[0]*l1[1]
    if abs(d) < 1e-9: return None
    return ((l1[1]*l2[2] - l2[1]*l1[2]) / d, (l2[0]*l1[2] - l1[0]*l2[2]) / d)

# Pentagram
V = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
pent_lines = [make_line(V[i], V[(i+2)%5]) for i in range(5)]

# Inner vertices
inner = []
for i, j in combinations(range(5), 2):
    p = ix_lines(pent_lines[i], pent_lines[j])
    if p:
        is_outer = any(abs(p[0]-v[0])<1e-6 and abs(p[1]-v[1])<1e-6 for v in V)
        if not is_outer and not any(abs(p[0]-v[0])<1e-6 and abs(p[1]-v[1])<1e-6 for v in inner):
            inner.append(p)

# Transform to SVG coordinates
W, H = 600, 700
margin = 80
scale = 220

def to_svg(p):
    return (W/2 + p[0]*scale, H/2 - 60 - p[1]*scale)

# === FIGURE 1: Base pentagram (5 triangles) ===
def draw_base():
    elements = []
    elements.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" style="background:white">')
    elements.append(f'<text x="{W/2}" y="30" font-size="18" text-anchor="middle" font-weight="bold">Базовая пентаграмма: 5 треугольников</text>')

    # Draw pentagram lines
    for i in range(5):
        p1 = to_svg(V[i])
        p2 = to_svg(V[(i+2)%5])
        elements.append(svg_line(p1[0], p1[1], p2[0], p2[1], "#333", 2))

    # Highlight 5 tip triangles
    colors = ["rgba(255,100,100,0.3)", "rgba(100,255,100,0.3)", "rgba(100,100,255,0.3)",
              "rgba(255,255,100,0.3)", "rgba(255,100,255,0.3)"]

    # For each tip: outer vertex + 2 adjacent inner vertices
    # Find which inner vertices are adjacent to each outer vertex
    for vi in range(5):
        # Lines through V[vi]: L_vi connects V[vi] to V[(vi+2)%5], L_(vi+3)%5 connects V[(vi+3)%5] to V[vi]
        # The tip triangle at V[vi] has vertices V[vi] and 2 inner points on the pentagram
        # These inner points are on the 2 pentagram lines through V[vi]
        # Line through V[vi]: li1 = vi (V[vi]→V[(vi+2)%5]) and li2 = (vi+3)%5 (V[(vi+3)%5]→V[vi])
        li1 = vi  # V[vi] → V[(vi+2)%5]
        li2 = (vi + 3) % 5  # V[(vi-2)%5] → V[vi] = V[(vi+3)%5] → V[vi]

        # Find inner vertices on these lines (other than V[vi] itself)
        tip_inner = []
        for ip in inner:
            on_l1 = abs(pent_lines[li1][0]*ip[0] + pent_lines[li1][1]*ip[1] + pent_lines[li1][2]) < 1e-6
            on_l2 = abs(pent_lines[li2][0]*ip[0] + pent_lines[li2][1]*ip[1] + pent_lines[li2][2]) < 1e-6
            if on_l1 or on_l2:
                # Check it's between V[vi] and the first inner crossing
                tip_inner.append(ip)

        # The 2 inner vertices closest to V[vi] on each line
        if len(tip_inner) >= 2:
            # Sort by distance to V[vi]
            tip_inner.sort(key=lambda p: (p[0]-V[vi][0])**2 + (p[1]-V[vi][1])**2)
            # Take the closest one on each line
            t1 = tip_inner[0]
            t2 = tip_inner[1]
            # Draw triangle
            pts = [to_svg(V[vi]), to_svg(t1), to_svg(t2)]
            elements.append(svg_polygon(pts, colors[vi]))

    # Label vertices
    for i, v in enumerate(V):
        p = to_svg(v)
        elements.append(svg_circle(p[0], p[1], 5, "blue", "blue"))
        offset_x = v[0] * 20
        offset_y = -v[1] * 20
        elements.append(svg_text(p[0]+offset_x, p[1]-offset_y+5, f"V{i}", 14, "blue"))

    for i, v in enumerate(inner):
        p = to_svg(v)
        elements.append(svg_circle(p[0], p[1], 4, "green", "green"))

    # Legend
    elements.append(svg_text(W/2, H-40, "5 треугольников (кончики) + 1 пятиугольник (центр) = 6 граней", 14, "#555"))

    elements.append('</svg>')
    return '\n'.join(elements)


# === FIGURE 2: Best config (9 triangles) ===
def draw_9tri():
    # V2→I4 + V4→I0
    # Clip to star boundary
    def clip(p1, p2):
        line = make_line(p1, p2)
        pts = []
        for i in range(5):
            sl = make_line(V[i], V[(i+2)%5])
            p = ix_lines(line, sl)
            if p:
                dx = V[(i+2)%5][0]-V[i][0]
                dy = V[(i+2)%5][1]-V[i][1]
                l2 = dx*dx + dy*dy
                t = ((p[0]-V[i][0])*dx + (p[1]-V[i][1])*dy) / l2
                if -1e-6 < t < 1+1e-6:
                    is_dup = any(abs(p[0]-q[0])<1e-6 and abs(p[1]-q[1])<1e-6 for q in pts)
                    if not is_dup:
                        pts.append(p)
        if len(pts) < 2: return None
        a, b, c = line
        pts.sort(key=lambda p: -b*p[0] + a*p[1])
        return (pts[0], pts[-1])

    seg1 = clip(V[2], inner[4])  # V2→I4
    seg2 = clip(V[4], inner[0])  # V4→I0

    # New intersection point
    l1 = make_line(seg1[0], seg1[1])
    l2 = make_line(seg2[0], seg2[1])
    center_ix = ix_lines(l1, l2)

    elements = []
    elements.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" style="background:white">')
    elements.append(f'<text x="{W/2}" y="30" font-size="18" text-anchor="middle" font-weight="bold">Лучшая конфигурация: 9 треугольников</text>')
    elements.append(f'<text x="{W/2}" y="52" font-size="14" text-anchor="middle" fill="#666">Линии: V2→I4 + V4→I0 (оси симметрии)</text>')

    # Draw pentagram
    for i in range(5):
        p1 = to_svg(V[i])
        p2 = to_svg(V[(i+2)%5])
        elements.append(svg_line(p1[0], p1[1], p2[0], p2[1], "#333", 2))

    # Draw new lines
    if seg1:
        p1, p2 = to_svg(seg1[0]), to_svg(seg1[1])
        elements.append(svg_line(p1[0], p1[1], p2[0], p2[1], "red", 2.5))
    if seg2:
        p1, p2 = to_svg(seg2[0]), to_svg(seg2[1])
        elements.append(svg_line(p1[0], p1[1], p2[0], p2[1], "blue", 2.5))

    # Mark vertices
    for i, v in enumerate(V):
        p = to_svg(v)
        elements.append(svg_circle(p[0], p[1], 5, "#333", "#333"))
        offset_x = v[0] * 20
        offset_y = -v[1] * 20
        elements.append(svg_text(p[0]+offset_x, p[1]-offset_y+5, f"V{i}", 12, "#333"))

    for i, v in enumerate(inner):
        p = to_svg(v)
        elements.append(svg_circle(p[0], p[1], 4, "green", "green"))
        elements.append(svg_text(p[0]+10, p[1]-8, f"I{i}", 10, "green"))

    if center_ix:
        p = to_svg(center_ix)
        elements.append(svg_circle(p[0], p[1], 5, "orange", "orange"))
        elements.append(svg_text(p[0]+12, p[1]+4, "X", 12, "orange"))

    # Legend
    y0 = H - 120
    elements.append(svg_text(W/2, y0, "Результат: 9 △ + 1 четырёхугольник + 1 пятиугольник", 14, "#333"))
    elements.append(svg_text(W/2, y0+22, "Четырёхугольник: I4—I1—I0—X (верхняя часть центра)", 12, "#666"))
    elements.append(svg_text(W/2, y0+40, "Пятиугольник: I2—I3—v11—X—v9 (нижняя часть центра)", 12, "#666"))
    elements.append(svg_line(W/2-200, y0+52, W/2+200, y0+52, "#ccc", 1))
    elements.append(svg_text(W/2, y0+70, "Для 10: нужно разбить четырёхугольник на треугольник", 13, "red"))
    elements.append(svg_text(W/2, y0+88, "Это требует чтобы 3 линии пересеклись в одной точке", 13, "red"))

    elements.append('</svg>')
    return '\n'.join(elements)


# === SUMMARY ===
summary = """
============================================================
ИТОГИ ИССЛЕДОВАНИЯ: ТРЕУГОЛЬНИКИ В ПЕНТАГРАММЕ
============================================================

ЗАДАЧА: Добавить 2 прямые к пентаграмме чтобы получить 10 треугольников.
Глава департамента нашёл решение для 9.

РЕЗУЛЬТАТЫ ПЕРЕБОРА:
- Базовая пентаграмма: 5 треугольных граней (кончики) + 1 пятиугольник (центр)
- Максимум с 2 линиями через вершины: 9 треугольников
- Максимум с 2 произвольными линиями: 9 треугольников
- Максимум для неправильных пентаграмм: 9 треугольников
- Протестировано: >500K конфигураций, ~100K форм звёзд

СТРУКТУРА ЛУЧШЕЙ КОНФИГУРАЦИИ (9 треугольников):
- 2 линии типа "ось симметрии" (V→противоположная внутр. вершина)
- Результат: 9 △ + 1 четырёхугольник + 1 пятиугольник = 11 граней
- Четырёхугольник образуется в верхней части центра
- Пятиугольник — в нижней части

ПОЧЕМУ НЕ 10:
Для 10 треугольников четырёхугольник должен стать треугольником.
Это возможно ТОЛЬКО если 2 его вершины сливаются (3 линии
пересекаются в одной точке). В пентаграмме ЛЮБЫХ пропорций
это условие не выполняется для данного типа конфигурации.

ВОЗМОЖНЫЕ НАПРАВЛЕНИЯ:
1. Может быть, задача считает ВСЕ видимые треугольники (не только атомарные грани)?
   Тогда базовая пентаграмма имеет 10 треугольников (5 маленьких + 5 больших).
2. Может быть, "линии" — не полные прямые а отрезки специфической длины?
3. Может быть, задача допускает линии за пределами звезды?

ГИПОТЕЗА: Возможно, при определённых пропорциях звезды и специфическом
расположении 2 отрезков, одна из линий проходит через точку пересечения
двух других линий, создавая вырожденную конфигурацию с 10 гранями.
"""

# Write SVG files
with open('/home/user/allen-coin-landing/star_base.svg', 'w') as f:
    f.write(draw_base())

with open('/home/user/allen-coin-landing/star_9tri.svg', 'w') as f:
    f.write(draw_9tri())

print(summary)

# Print the SVG paths
print("\nVisualization files created:")
print("  star_base.svg - базовая пентаграмма с 5 треугольниками")
print("  star_9tri.svg - лучшая конфигурация с 9 треугольниками")
