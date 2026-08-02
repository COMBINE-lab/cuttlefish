#!/usr/bin/env python3
"""Generate the cuttlefish de Bruijn graph logo.

The cuttlefish's arms are paths in a colored compacted de Bruijn graph:
circular vertices joined by edges, each tinted by which input sources contain
it. Vertices shared by two sources are drawn as two-tone circles, which is the
whole point of a *colored* graph. Two arms branch, because a branch is exactly
where a unitig has to end.
"""

import math
import re
import sys

W = H = 256
MARGIN = 6  # keep all geometry this far inside the canvas

# Colour classes = input sources. Distinguishable on both light and dark
# backgrounds, with differing lightnesses so the mark survives greyscale.
SRC = {
    "a": "#E4572E",  # coral
    "b": "#F2A93B",  # amber
    "c": "#17A398",  # teal
    "d": "#4C6EF5",  # indigo
    "e": "#9B5DE5",  # violet
    "f": "#5FAD56",  # green
}

BODY_LIGHT = "#22303F"
BODY_DARK = "#EAF0F7"


def bezier(p0, p1, p2, p3, t):
    u = 1 - t
    return (
        u**3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t**3 * p3[0],
        u**3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t**3 * p3[1],
    )


def arm_points(start, angle_deg, length, n, hook=0.55):
    """Sample n+1 points along one arm.

    `angle_deg` is the heading where the arm leaves the head, measured from
    straight down. The arm then hooks back *toward* vertical by `hook` of that
    angle, so tips point downward the way a cuttlefish's arms hang, instead of
    splaying outward off the canvas.
    """
    a0 = math.radians(angle_deg)
    a1 = math.radians(angle_deg * (1 - hook))
    p0 = start
    mid = (a0 + a1) / 2
    p3 = (start[0] + length * math.sin(mid), start[1] + length * math.cos(mid))
    p1 = (start[0] + length * 0.42 * math.sin(a0), start[1] + length * 0.42 * math.cos(a0))
    p2 = (p3[0] - length * 0.30 * math.sin(a1), p3[1] - length * 0.30 * math.cos(a1))
    return [bezier(p0, p1, p2, p3, i / n) for i in range(n + 1)]


# attachment x, heading, length, vertex count, colour class per edge.
ARMS = [
    (108, -40, 74, 4, ["c", "c", "d", "d"]),
    (115, -25, 62, 3, ["a", "b", "b"]),
    (122, -11, 56, 3, ["e", "d", "d"]),
    (130, 4, 58, 3, ["b", "b", "a"]),
    (138, 18, 64, 3, ["f", "f", "c"]),
    (145, 32, 76, 4, ["d", "d", "e", "e"]),
]
LONG = [  # the two long feeding tentacles
    (102, -55, 100, 5, ["a", "a", "f", "c", "c"]),
    (152, 48, 100, 5, ["e", "d", "d", "b", "b"]),
]
ALL = ARMS + LONG

# A branch: a second edge leaving an interior vertex, keyed by (arm, vertex).
BRANCH = {
    (0, 2): (-58, 36, 2, ["f", "f"]),
    (5, 2): (54, 34, 2, ["c", "c"]),
}

# Vertices belonging to two sources, drawn two-tone: (arm, vertex) -> class.
SHARED = {(0, 2): "f", (2, 1): "b", (5, 2): "c", (6, 3): "d", (7, 2): "a"}

# The favicon is a different drawing, not the logo scaled down. At 16 px the
# eight-arm version is mush, so this one keeps four arms with two vertices
# each, drawn much heavier, and drops the branches and two-tone vertices --
# detail that cannot survive the size is detail that only muddies it.
FAVICON_ARMS = [
    (110, -46, 62, 2, ["a", "a"]),
    (122, -16, 56, 2, ["c", "c"]),
    (134, 16, 56, 2, ["d", "d"]),
    (146, 46, 62, 2, ["b", "b"]),
]

HEAD_Y = 156.0


def build_arms(favicon=False):
    """Return (edge elements, node elements, bounding box)."""
    edges, nodes = [], []
    xs, ys = [], []

    def emit(pts, classes, n, base_w, base_r, key=None):
        for i in range(n):
            (x0, y0), (x1, y1) = pts[i], pts[i + 1]
            w = base_w - (base_w - 2.8) * (i / max(n - 1, 1))
            edges.append(
                f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
                f'stroke="{SRC[classes[i]]}" stroke-width="{w:.2f}" stroke-linecap="round"/>'
            )
            xs.extend([x0, x1])
            ys.extend([y0, y1])
        for i in range(1, n + 1):
            x, y = pts[i]
            r = base_r - (base_r - 3.6) * ((i - 1) / max(n - 1, 1))
            primary = SRC[classes[i - 1]]
            second = SHARED.get((key, i)) if key is not None else None
            if second:
                nodes.append(
                    f'<g><circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{primary}"/>'
                    f'<path d="M {x:.2f} {y - r:.2f} A {r:.2f} {r:.2f} 0 0 1 '
                    f'{x:.2f} {y + r:.2f} Z" fill="{SRC[second]}"/></g>'
                )
            else:
                nodes.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{primary}"/>')
            xs.extend([x - r, x + r])
            ys.extend([y - r, y + r])

    for ai, (ax, angle, length, n, classes) in enumerate(FAVICON_ARMS if favicon else ALL):
        # Outer arms attach slightly higher, following the curve of the head.
        start = (ax, HEAD_Y - abs(ax - 128) * 0.16)
        pts = arm_points(start, angle, length, n)
        if favicon:
            emit(pts, classes, n, 15.0, 13.0)
            continue
        emit(pts, classes, n, 6.8, 6.6, key=ai)
        if (ai, 2) in BRANCH:
            b_angle, b_len, b_n, b_cls = BRANCH[(ai, 2)]
            emit(arm_points(pts[2], b_angle, b_len, b_n), b_cls, b_n, 4.6, 5.2)

    return edges, nodes, (min(xs), max(xs), min(ys), max(ys))


def render(body, eye, favicon=False):
    """`eye` is the sclera colour; it must contrast with `body`, which is what
    the W-shaped pupil is drawn in. On the dark-page variant the body is light,
    so the sclera has to go dark or the eyes disappear entirely."""
    edges, nodes, (x0, x1, y0, y1) = build_arms(favicon)
    # The favicon crops to the drawing rather than to the logo's canvas: at
    # tab size, empty margin is the difference between a recognisable mark and
    # a smudge. The crop is squared about the body's centreline so the animal
    # still sits centred.
    if favicon:
        top, bottom = 14.0, y1 + 6
        side = bottom - top
        view = f'{128 - side / 2:.0f} {top:.0f} {side:.0f} {side:.0f}'
    else:
        view = f'0 0 {W} {H}'
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{view}" '
        f'width="{W}" height="{H}" role="img" '
        f'aria-label="A cuttlefish whose arms are the colored paths of a de Bruijn graph">',
        "<g>" + "".join(edges) + "</g>",
        "<g>" + "".join(nodes) + "</g>",
        f'<g fill="{body}">',
        # Fin margin: a cuttlefish's fin is a continuous skirt running the whole
        # length of the mantle, so head-on it reads as a soft halo around the
        # body rather than as two appendages. Drawn behind the mantle and
        # slightly wider, with the scalloped waviness the fin actually has.
        '<path d="M 128 20 C 152 20 172 44 176 76 C 179 104 174 124 166 136 '
        'C 162 142 156 140 158 132 C 164 112 166 84 162 62 '
        'C 158 38 146 24 128 24 C 110 24 98 38 94 62 '
        'C 90 84 92 112 98 132 C 100 140 94 142 90 136 '
        'C 82 124 77 104 80 76 C 84 44 104 20 128 20 Z" opacity="0.38"/>',
        # Mantle: broad and shield-like -- cuttlefish are stouter than squid --
        # rounded at the posterior and widest about two thirds down.
        '<path d="M 128 18 C 149 18 162 40 165 70 C 167 98 163 118 156 130 '
        'L 100 130 C 93 118 89 98 91 70 C 94 40 107 18 128 18 Z"/>',
        # Head: a touch narrower than the mantle, carrying the eyes and the
        # crown the arms hang from.
        '<path d="M 101 129 L 155 129 C 157 144 150 156 128 158 '
        'C 106 156 99 144 101 129 Z"/>',
        "</g>",
    ]
    # Eyes, with the W-shaped pupil that is a cuttlefish's signature.
    for cx in (114, 142):
        parts.append(f'<ellipse cx="{cx}" cy="140" rx="10.6" ry="8.2" fill="{eye}"/>')
        parts.append(
            f'<path d="M {cx - 7.2} 135.2 L {cx - 3} 140.2 L {cx} 136.6 '
            f'L {cx + 3} 140.2 L {cx + 7.2} 135.2 L {cx + 7.2} 139.8 '
            f'L {cx} 143.2 L {cx - 7.2} 139.8 Z" fill="{body}"/>'
        )
    parts.append("</svg>")

    fits = x0 >= MARGIN and x1 <= W - MARGIN and y1 <= H - MARGIN
    print(
        f"arm bbox: x {x0:.1f}..{x1:.1f}  y {y0:.1f}..{y1:.1f}  "
        f"{'OK' if fits else 'OVERFLOWS CANVAS'}",
        file=sys.stderr,
    )
    return "\n".join(parts)


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "light"
    if target == "dark":
        print(render(BODY_DARK, "#22303F"))
    elif target == "favicon":
        print(render(BODY_LIGHT, "#FFFFFF", favicon=True))
    else:
        print(render(BODY_LIGHT, "#FFFFFF"))
