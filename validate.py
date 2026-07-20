#!/usr/bin/env python3
"""Light-weight validation of a gerber2kicad conversion, without KiCad.

It re-parses the generated .kicad_pcb with an independent s-expression
reader, rebuilds copper/edge geometry from it, and compares against the
copper geometry computed straight from the original Gerber files:

  * copper area per side (union of all objects), % difference
  * bounding boxes
  * object counts (pads vs flashes, tracks, holes)
  * side-by-side SVG renders in <out>/validation/

Usage: python3 validate.py GERBER_DIR KICAD_PCB
"""

import math
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.simplefilter('ignore')

from gerbonara import LayerStack
from gerbonara import apertures as gn_ap
from gerbonara.graphic_objects import Line as GnLine, Arc as GnArc, \
    Flash as GnFlash, Region as GnRegion
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, LineString, box
from shapely.ops import unary_union


def sample_arc(s, m, e, step=0.2):
    """Sample the circular arc through 3 points into a polyline."""
    n = len([s, m, e])
    sx = sy = sz = sxx = syy = sxy = sxz = syz = 0.0
    for x, y in (s, m, e):
        z = x * x + y * y
        sx += x; sy += y; sz += z
        sxx += x * x; syy += y * y; sxy += x * y
        sxz += x * z; syz += y * z

    def det3(mm):
        return (mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
                - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
                + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]))

    d = det3([[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, 3]])
    if abs(d) < 1e-9:
        return [s, m, e]
    a = det3([[sxz, sxy, sx], [syz, syy, sy], [sz, sy, 3]]) / d
    b = det3([[sxx, sxz, sx], [sxy, syz, sy], [sx, sz, 3]]) / d
    c = det3([[sxx, sxy, sxz], [sxy, syy, syz], [sx, sy, sz]]) / d
    cx, cy = a / 2, b / 2
    rr = c + cx * cx + cy * cy
    if rr <= 0:
        return [s, m, e]
    r = math.sqrt(rr)
    a1 = math.atan2(s[1] - cy, s[0] - cx)
    a2 = math.atan2(m[1] - cy, m[0] - cx)
    a3 = math.atan2(e[1] - cy, e[0] - cx)
    d12 = (a2 - a1) % (2 * math.pi)
    d13 = (a3 - a1) % (2 * math.pi)
    if d12 > d13:
        d13 -= 2 * math.pi
    k = max(3, int(abs(d13) * r / step))
    return [(cx + r * math.cos(a1 + d13 * i / k),
             cy + r * math.sin(a1 + d13 * i / k)) for i in range(k + 1)]


# --------------------------------------------------------------------------
# s-expression parsing of the .kicad_pcb
# --------------------------------------------------------------------------

def sexpr_parse(text):
    tokens = re.findall(r'"(?:[^"\\]|\\.)*"|[()]|[^\s()"]+', text)
    pos = 0

    def parse():
        nonlocal pos
        assert tokens[pos] == '('
        pos += 1
        out = []
        while tokens[pos] != ')':
            if tokens[pos] == '(':
                out.append(parse())
            else:
                t = tokens[pos]
                if t.startswith('"'):
                    t = t[1:-1].replace('\\"', '"').replace('\\\\', '\\')
                out.append(t)
                pos += 1
        pos += 1
        return out

    return parse()


def find_all(node, name):
    for item in node:
        if isinstance(item, list) and item and item[0] == name:
            yield item


def find_one(node, name, default=None):
    return next(find_all(node, name), default)


def get_pts(node):
    pts = find_one(node, 'pts', [])
    return [(float(x[1]), float(x[2])) for x in find_all(pts, 'xy')]


def pcb_geometry(pcb_path):
    """-> dict: copper geoms per side, edge polys, counts."""
    tree = sexpr_parse(Path(pcb_path).read_text())
    copper = {'F': [], 'B': []}
    edge = []
    counts = defaultdict(int)

    for seg in find_all(tree, 'segment'):
        a = find_one(seg, 'start'); b = find_one(seg, 'end')
        w = float(find_one(seg, 'width')[1])
        lay = find_one(seg, 'layer')[1]
        g = LineString([(float(a[1]), float(a[2])),
                        (float(b[1]), float(b[2]))]).buffer(w / 2)
        copper[lay[0]].append(g)
        counts['segments'] += 1

    for arc in find_all(tree, 'arc'):
        pts = [(float(p[1]), float(p[2]))
               for p in (find_one(arc, k) for k in ('start', 'mid', 'end'))]
        w = float(find_one(arc, 'width')[1])
        lay = find_one(arc, 'layer')[1]
        g = LineString(sample_arc(*pts)).buffer(w / 2)
        copper[lay[0]].append(g)
        counts['arcs'] += 1

    for via in find_all(tree, 'via'):
        at = find_one(via, 'at')
        size = float(find_one(via, 'size')[1])
        g = Point(float(at[1]), float(at[2])).buffer(size / 2)
        copper['F'].append(g)
        copper['B'].append(g)
        counts['vias'] += 1

    for zone in find_all(tree, 'zone'):
        lay = find_one(zone, 'layer')[1]
        for fp in find_all(zone, 'filled_polygon'):
            pts = get_pts(fp)
            if len(pts) >= 3:
                copper[lay[0]].append(Polygon(pts).buffer(0))
                counts['filled_polygons'] += 1

    for fp in find_all(tree, 'footprint'):
        at = find_one(fp, 'at')
        fx, fy = float(at[1]), float(at[2])
        for pad in find_all(fp, 'pad'):
            counts['pads'] += 1
            ptype, shape = pad[2], pad[3]
            pat = find_one(pad, 'at')
            px, py = fx + float(pat[1]), fy + float(pat[2])
            size = find_one(pad, 'size')
            w, h = float(size[1]), float(size[2])
            layers = find_one(pad, 'layers')[1:]
            if ptype == 'np_thru_hole':
                counts['npth'] += 1
                continue
            if shape == 'circle':
                g = Point(px, py).buffer(w / 2)
            elif shape == 'roundrect':   # octagon approximation
                g = Point(px, py).buffer(w / 2 / math.cos(math.pi / 8),
                                         quad_segs=2)
                counts['octagon_pads'] += 1
            elif shape == 'oval':
                ang = math.radians(float(pat[3]) if len(pat) > 3 else 0)
                half = (w - h) / 2 if w > h else 0
                dx2, dy2 = half * math.cos(ang), -half * math.sin(ang)
                g = LineString([(px - dx2, py - dy2),
                                (px + dx2, py + dy2)]).buffer(min(w, h) / 2)
                counts['oval_pads'] += 1
            else:
                g = box(px - w / 2, py - h / 2, px + w / 2, py + h / 2)
            sides = []
            for l in layers:
                if l.startswith('*'):
                    sides = ['F', 'B']
                    break
                if l.endswith('.Cu'):
                    sides.append(l[0])
            for s in set(sides):
                copper[s].append(g)

    for gp in find_all(tree, 'gr_poly'):
        lay = find_one(gp, 'layer')[1]
        if lay == 'Edge.Cuts':
            pts = get_pts(gp)
            if len(pts) >= 3:
                edge.append(Polygon(pts))
                counts['edge_polys'] += 1
    for name in ('gr_line', 'gr_arc', 'gr_circle'):
        for gl in find_all(tree, name):
            lay = find_one(gl, 'layer')[1]
            counts[f'{name}_{lay}'] += 1
            if lay == 'Edge.Cuts':
                counts['edge_prims'] += 1

    return {'copper': {s: unary_union(gs) if gs else None
                       for s, gs in copper.items()},
            'edge': edge, 'counts': counts}


# --------------------------------------------------------------------------
# reference geometry straight from the Gerbers
# --------------------------------------------------------------------------

def flash_geom(o):
    ap = o.aperture
    if isinstance(ap, gn_ap.RectangleAperture):
        return box(o.x - ap.w / 2, o.y - ap.h / 2,
                   o.x + ap.w / 2, o.y + ap.h / 2)
    if isinstance(ap, gn_ap.PolygonAperture):
        r = ap.diameter / 2
        pts = [(o.x + r * math.cos(ap.rotation + 2 * math.pi * i
                                   / ap.n_vertices),
                o.y + r * math.sin(ap.rotation + 2 * math.pi * i
                                   / ap.n_vertices))
               for i in range(ap.n_vertices)]
        return Polygon(pts)
    d = getattr(ap, 'diameter', 0) or 0.01
    return Point(o.x, o.y).buffer(d / 2)


def gerber_copper(gerber_dir):
    stack = LayerStack.open(gerber_dir)
    gl = dict(stack.graphic_layers)
    out = {}
    counts = defaultdict(int)
    for side, sname in (('top', 'F'), ('bottom', 'B')):
        geoms = []
        for o in gl[(side, 'copper')].objects:
            dark = getattr(o, 'polarity_dark', True)
            if isinstance(o, GnRegion):
                pts = list(o.outline)
                if len(pts) < 3:
                    continue
                g = Polygon(pts).buffer(0)
                counts[f'{sname}_regions'] += 1
            elif isinstance(o, GnLine):
                w = getattr(o.aperture, 'diameter', 0) or 0.01
                g = LineString([(o.x1, o.y1), (o.x2, o.y2)]).buffer(w / 2)
                counts[f'{sname}_lines'] += 1
            elif isinstance(o, GnArc):
                continue
            elif isinstance(o, GnFlash):
                g = flash_geom(o)
                counts[f'{sname}_flashes'] += 1
            else:
                continue
            if dark:
                geoms.append(g)
            else:
                geoms = [p.difference(g) for p in geoms]
        out[sname] = unary_union(geoms) if geoms else None
    drills = []
    for dl in stack.drill_layers:
        for o in dl.objects:
            if isinstance(o, GnFlash):
                drills.append((o.x, o.y, o.aperture.diameter))
    counts['drills'] = len(drills)
    return out, counts, drills


# --------------------------------------------------------------------------
# SVG rendering of shapely geometry
# --------------------------------------------------------------------------

def geom_to_svg_paths(geom, color, flip_y=None):
    if geom is None or geom.is_empty:
        return ''
    geoms = geom.geoms if hasattr(geom, 'geoms') else [geom]
    d = []
    for g in geoms:
        if g.geom_type != 'Polygon':
            continue
        for ring in [g.exterior] + list(g.interiors):
            coords = list(ring.coords)
            part = []
            for i, (x, y) in enumerate(coords):
                if flip_y is not None:
                    y = flip_y - y
                part.append(f'{"M" if i == 0 else "L"}{x:.3f},{y:.3f}')
            d.append(' '.join(part) + ' Z')
    return (f'<path d="{" ".join(d)}" fill="{color}" fill-rule="evenodd" '
            f'stroke="none"/>')


def write_svg(path, layers, bounds, title):
    minx, miny, maxx, maxy = bounds
    pad = 2
    w, h = maxx - minx + 2 * pad, maxy - miny + 2 * pad
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" '
             f'viewBox="{minx - pad} {miny - pad} {w} {h}" '
             f'width="800" height="{800 * h / w:.0f}">',
             f'<rect x="{minx - pad}" y="{miny - pad}" width="{w}" '
             f'height="{h}" fill="#001023"/>',
             f'<title>{title}</title>']
    parts += layers
    parts.append('</svg>')
    Path(path).write_text('\n'.join(parts))


def main():
    gerber_dir, pcb_path = sys.argv[1], sys.argv[2]
    out_dir = Path(pcb_path).parent / 'validation'
    out_dir.mkdir(exist_ok=True)

    print('Parsing generated .kicad_pcb ...')
    kic = pcb_geometry(pcb_path)
    print('Rebuilding reference geometry from Gerbers ...')
    ref, ref_counts, drills = gerber_copper(gerber_dir)

    ok = True
    print('\n=== numeric comparison ===')
    for s, label in (('F', 'front'), ('B', 'back')):
        a_ref = ref[s].area if ref[s] else 0
        a_kic = kic['copper'][s].area if kic['copper'][s] else 0
        diff = 100 * (a_kic - a_ref) / a_ref if a_ref else 0
        # bbox sizes (translation between the two coordinate frames is
        # expected: the converter shifts/mirrors, sizes must match)
        bb_r = ref[s].bounds if ref[s] else (0, 0, 0, 0)
        bb_k = kic['copper'][s].bounds if kic['copper'][s] else (0, 0, 0, 0)
        sz_r = (bb_r[2] - bb_r[0], bb_r[3] - bb_r[1])
        sz_k = (bb_k[2] - bb_k[0], bb_k[3] - bb_k[1])
        dsz = math.hypot(sz_r[0] - sz_k[0], sz_r[1] - sz_k[1])
        # exact-shape check: mirror+translate the reference onto the
        # converted frame and measure the symmetric difference
        xor_pct = None
        if ref[s] and kic['copper'][s]:
            rm = scale(ref[s], xfact=1, yfact=-1, origin=(0, 0))
            rmb = rm.bounds
            rm = translate(rm, xoff=bb_k[0] - rmb[0],
                           yoff=bb_k[1] - rmb[1])
            xor_pct = 100 * rm.symmetric_difference(kic['copper'][s]).area \
                / a_ref
        status = 'OK' if abs(diff) < 5 and dsz < 0.5 and \
            (xor_pct is None or xor_pct < 2) else 'MISMATCH'
        if status != 'OK':
            ok = False
        print(f'{label} copper: area {a_ref:8.1f} -> {a_kic:8.1f} mm2 '
              f'({diff:+.2f}%), bbox {sz_r[0]:.2f}x{sz_r[1]:.2f} -> '
              f'{sz_k[0]:.2f}x{sz_k[1]:.2f} mm, '
              f'shape XOR {xor_pct:.2f}%   [{status}]')

    c = kic['counts']
    n_flash = ref_counts['F_flashes'] + ref_counts['B_flashes']
    print(f"pads: {c['pads']} (copper flashes in gerbers: {n_flash}, "
          f"shared by through-hole pads)")
    print(f"segments: {c['segments']} "
          f"(gerber draws: {ref_counts['F_lines'] + ref_counts['B_lines']}, "
          f"incl. zero-length -> pads)")
    print(f"vias: {c['vias']}  npth: {c['npth']}  "
          f"drills in gerber: {ref_counts['drills']}")
    print(f"holes accounted: "
          f"{c['vias'] + c['npth'] + c.get('pads_th', 0)} "
          f"(vias + npth; plated pads carry the rest)")
    print(f"edge polygons: {c['edge_polys']}")

    print('\n=== SVG renders (validation/) ===')
    ymax = max(g.bounds[3] for g in ref.values() if g)
    for s, label in (('F', 'front'), ('B', 'back')):
        if ref[s]:
            write_svg(out_dir / f'original_{label}.svg',
                      [geom_to_svg_paths(ref[s], '#c87137', flip_y=ymax)],
                      (ref[s].bounds[0], ymax - ref[s].bounds[3],
                       ref[s].bounds[2], ymax - ref[s].bounds[1]),
                      f'original {label} copper')
        if kic['copper'][s]:
            g = kic['copper'][s]
            layers = [geom_to_svg_paths(g, '#2a9d8f')]
            for e in kic['edge']:
                layers.append(f'<path d="{" ".join("%s%.3f,%.3f" % ("M" if i == 0 else "L", x, y) for i, (x, y) in enumerate(e.exterior.coords))} Z" fill="none" stroke="#ffd166" stroke-width="0.2"/>')
            write_svg(out_dir / f'converted_{label}.svg', layers,
                      g.bounds, f'converted {label} copper')
        print(f'  original_{label}.svg / converted_{label}.svg')

    print('\nRESULT:', 'PASS' if ok else 'CHECK NEEDED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
