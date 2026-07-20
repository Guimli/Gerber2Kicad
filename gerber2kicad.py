#!/usr/bin/env python3
"""gerber2kicad — convert a set of Gerber/Excellon fabrication files into a
KiCad 9-openable project (.kicad_pro + .kicad_pcb + .kicad_sch).

Gerber files carry no netlist, footprint or component information, so this
tool reconstructs what it can:
  * copper draws        -> tracks / arcs
  * copper flashes      -> pads (grouped into footprints, see below)
  * copper regions      -> filled zones
  * drill hits          -> plated pads, vias or NPTH holes (heuristic)
  * silk/mask/paste     -> graphic primitives on the matching layers
  * board outline       -> chained + simplified Edge.Cuts polygons
  * connectivity        -> nets derived from copper geometry (shapely),
                           the net with the largest pour area is named GND
  * optional BOM/CPL    -> SMD pads near CPL positions are grouped into a
    (JLCPCB-style CSV)     footprint per component; those components are
                           also placed in the reconstructed schematic with
                           global net labels on their pins.

Usage:
    python3 gerber2kicad.py INPUT -o OUTDIR [--name NAME]
                            [--bom BOM.csv] [--cpl CPL.csv]

INPUT is a directory or a .zip containing the Gerber files.
"""

import argparse
import csv
import dataclasses
import io
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
import warnings
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

warnings.simplefilter('ignore')  # gerbonara is chatty about 0-width draws

from gerbonara import LayerStack
from gerbonara import apertures as gn_ap
from gerbonara.graphic_objects import Line as GnLine, Arc as GnArc, \
    Flash as GnFlash, Region as GnRegion
from shapely.geometry import Point, Polygon, LineString, box
from shapely.ops import unary_union
from shapely.strtree import STRtree

MM = 'mm'
PCB_VERSION = 20240108      # KiCad 8.0 board format, opens fine in KiCad 9
SCH_VERSION = 20231120      # KiCad 8.0 schematic format
GENERATOR = 'gerber2kicad'
GEN_VERSION = '1.0'

COORD_EPS = 0.005           # mm, endpoint snapping for segment chaining
ARC_FIT_TOL = 0.015         # mm, max residual for arc/circle recovery
ARC_FIT_RMAX = 150          # mm, max plausible arc radius
ARC_MIN_SPAN = 0.35         # rad, min swept angle to accept an arc
HOLE_MATCH_TOL = 0.15       # mm, drill <-> flash center matching
MASK_MATCH_TOL = 0.15       # mm, mask flash <-> pad matching
CPL_GROUP_RADIUS = 8.0      # mm, SMD pad -> CPL component association
OUTLINE_SIMPLIFY = 0.01     # mm, Douglas-Peucker tolerance
REGION_SIMPLIFY = 0.005     # mm

# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def nid():
    return str(uuid.uuid4())

def ocr_image(pil_image):
    """OCR a PIL image with tesseract, return (text, mean confidence)."""
    import tempfile as _tf
    with _tf.NamedTemporaryFile(suffix='.png') as f:
        pil_image.save(f.name)
        try:
            out = subprocess.run(
                ['tesseract', f.name, 'stdout', '--psm', '7', 'tsv'],
                capture_output=True, text=True, timeout=30).stdout
        except Exception:
            return None, 0
    words, confs = [], []
    for row in csv.DictReader(out.splitlines(), delimiter='\t'):
        t = row.get('text', '').strip()
        if t and float(row.get('conf', -1)) > 0:
            words.append(t)
            confs.append(float(row['conf']))
    if not words:
        return None, 0
    return ' '.join(words), sum(confs) / len(confs)


def fit_circle(pts):
    """Least-squares (Kåsa) circle fit. -> (cx, cy, r, max_residual).
    Points are centered on their centroid first: without this the normal
    equations lose precision (coordinates ~1e2 mm vs sagittas ~1e-3 mm)."""
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sx = sy = sz = sxx = syy = sxy = sxz = syz = 0.0
    for px, py in pts:
        x, y = px - mx, py - my
        z = x * x + y * y
        sx += x; sy += y; sz += z
        sxx += x * x; syy += y * y; sxy += x * y
        sxz += x * z; syz += y * z

    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    m = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]]
    d = det3(m)
    if abs(d) < 1e-9:
        return None
    a = det3([[sxz, sxy, sx], [syz, syy, sy], [sz, sy, n]]) / d
    b = det3([[sxx, sxz, sx], [sxy, syz, sy], [sx, sz, n]]) / d
    cc = det3([[sxx, sxy, sxz], [sxy, syy, syz], [sx, sy, sz]]) / d
    cx, cy = a / 2, b / 2
    rr = cc + cx * cx + cy * cy
    if rr <= 0:
        return None
    r = math.sqrt(rr)
    cx += mx
    cy += my
    res = max(abs(math.hypot(x - cx, y - cy) - r) for x, y in pts)
    return cx, cy, r, res


def polyline_to_prims(pts, closed, tol=ARC_FIT_TOL, rmax=ARC_FIT_RMAX,
                      min_pts=5, min_span=ARC_MIN_SPAN):
    """Recover circles / arcs / merged straight lines from a tessellated
    polyline.  Returns a list of primitives:
        ('circle', (cx, cy), r)
        ('arc', p_start, p_mid, p_end)
        ('line', p1, p2)
    Endpoints of the input chain are preserved exactly, so contours stay
    connected."""
    def steps_ok(window, cx, cy, r):
        """The polyline must really hug the circle: every chord must stay
        within tol of the arc (sagitta check — kills rectangles, whose
        vertices sit exactly on their circumcircle) and turn in a single
        direction."""
        sign = 0
        for k in range(len(window) - 1):
            v1 = (window[k][0] - cx, window[k][1] - cy)
            v2 = (window[k + 1][0] - cx, window[k + 1][1] - cy)
            da = math.atan2(v1[0] * v2[1] - v1[1] * v2[0],
                            v1[0] * v2[0] + v1[1] * v2[1])
            if abs(da) > 1e-9:
                s = 1 if da > 0 else -1
                if sign and s != sign:
                    return False
                sign = s
            if r * (1 - math.cos(da / 2)) > tol * 1.5:
                return False
        return True

    prims = []
    n = len(pts)
    i = 0
    while i < n - 1:
        best = None
        j = i + min_pts - 1
        while j < n:
            f = fit_circle(pts[i:j + 1])
            if f is None:
                break
            cx, cy, r, res = f
            if res > tol or r > rmax or not steps_ok(pts[i:j + 1], cx, cy, r):
                break
            best = (j, cx, cy, r)
            j += 1
        if best:
            j, cx, cy, r = best
            span = 0.0
            for k in range(i, j):
                v1 = (pts[k][0] - cx, pts[k][1] - cy)
                v2 = (pts[k + 1][0] - cx, pts[k + 1][1] - cy)
                span += math.atan2(v1[0] * v2[1] - v1[1] * v2[0],
                                   v1[0] * v2[0] + v1[1] * v2[1])
            if abs(span) >= min_span:
                if abs(span) > 2 * math.pi - 0.15 and closed \
                        and i == 0 and j == n - 1:
                    prims.append(('circle', (cx, cy), r))
                    i = j
                    continue
                a1 = math.atan2(pts[i][1] - cy, pts[i][0] - cx)
                am = a1 + span / 2
                mid = (cx + r * math.cos(am), cy + r * math.sin(am))
                prims.append(('arc', pts[i], mid, pts[j]))
                i = j
                continue
        # straight run: merge collinear points within tol
        j = i + 1
        ax, ay = pts[i]
        while j + 1 < n:
            bx, by = pts[j + 1]
            ll = math.hypot(bx - ax, by - ay)
            if ll < 1e-9:
                break
            ok = all(
                abs((bx - ax) * (ay - py) - (ax - px) * (by - ay)) / ll <= tol
                for px, py in pts[i + 1:j + 1])
            if not ok:
                break
            j += 1
        prims.append(('line', pts[i], pts[j]))
        i = j
    return prims


def sample_arc(s, m, e, step=0.2):
    """Sample the circular arc through 3 points into a polyline."""
    f = fit_circle([s, m, e])
    if f is None:
        return [s, m, e]
    cx, cy, r, _ = f
    a1 = math.atan2(s[1] - cy, s[0] - cx)
    a2 = math.atan2(m[1] - cy, m[0] - cx)
    a3 = math.atan2(e[1] - cy, e[0] - cx)
    d12 = (a2 - a1) % (2 * math.pi)
    d13 = (a3 - a1) % (2 * math.pi)
    if d12 > d13:                    # mid not between start/end going ccw
        d13 -= 2 * math.pi
    n = max(3, int(abs(d13) * r / step))
    return [(cx + r * math.cos(a1 + d13 * k / n),
             cy + r * math.sin(a1 + d13 * k / n)) for k in range(n + 1)]


def point_seg_dist(p, a, b):
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    ll = dx * dx + dy * dy
    if ll == 0:
        return dist(p, a)
    t = max(0, min(1, ((p[0] - ax) * dx + (p[1] - ay) * dy) / ll))
    return dist(p, (ax + t * dx, ay + t * dy))


def fnum(v):
    """Format a number the way KiCad likes: fixed decimals, no trailing junk."""
    s = f'{v:.6f}'.rstrip('0').rstrip('.')
    return s if s not in ('-0', '') else '0'

def q(s):
    return '"' + str(s).replace('\\', '\\\\').replace('"', '\\"') + '"'

def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


class SExpr:
    """Tiny indenting s-expression writer."""
    def __init__(self):
        self.lines = []
        self.depth = 0

    def raw(self, text):
        self.lines.append('\t' * self.depth + text)

    def open(self, text):
        self.raw('(' + text)
        self.depth += 1

    def close(self, suffix=''):
        self.depth -= 1
        self.raw(')' + suffix)

    def leaf(self, text):
        self.raw('(' + text + ')')

    def text(self):
        return '\n'.join(self.lines) + '\n'


# ---------------------------------------------------------------------------
# intermediate model
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class Pad:
    x: float                 # transformed (KiCad) coords
    y: float
    shape: str               # 'circle' | 'rect' | 'octagon' | 'oval'
    size: tuple              # (w, h) or (d, d); octagon: across-flats
    side: str                # 'F' | 'B' | 'TH'
    drill: float = None
    plated: bool = True
    mask_margin: float = None
    paste: bool = False
    net: int = 0
    number: str = '1'
    component: str = None    # CPL reference if grouped
    rot: float = 0           # degrees, for oval pads
    seg: tuple = None        # ((x1,y1),(x2,y2),w) capsule of an oval pad


@dataclasses.dataclass
class Track:
    x1: float; y1: float; x2: float; y2: float
    width: float
    layer: str
    net: int = 0


@dataclasses.dataclass
class TrackArc:
    start: tuple; mid: tuple; end: tuple
    width: float
    layer: str
    net: int = 0


@dataclasses.dataclass
class Via:
    x: float; y: float
    size: float
    drill: float
    net: int = 0


@dataclasses.dataclass
class ZonePoly:
    poly: object             # shapely Polygon (may have interior holes)
    layer: str
    net: int = 0
    net_name: str = ''


@dataclasses.dataclass
class Component:
    ref: str
    value: str
    pads: list                # list of Pad
    is_smd: bool
    kind: str = 'U'           # symbol class: R C CP D L J SW U TP H
    in_schematic: bool = True


@dataclasses.dataclass
class SilkWord:
    x: float
    y: float
    text: str
    conf: float
    side: str
    claimed: bool = False


@dataclasses.dataclass
class Graphic:
    kind: str                # 'line' | 'arc' | 'poly' | 'circle' | 'rect'
    layer: str
    width: float = 0
    fill: bool = False
    pts: list = None         # line: [p1,p2]; poly: ring; arc: [start,mid,end]
    center: tuple = None     # circle
    radius: float = None     # circle
    corners: tuple = None    # rect: (x1,y1,x2,y2)


LAYER_MAP = {
    ('top', 'copper'): 'F.Cu',
    ('bottom', 'copper'): 'B.Cu',
    ('top', 'silk'): 'F.SilkS',
    ('bottom', 'silk'): 'B.SilkS',
    ('top', 'mask'): 'F.Mask',
    ('bottom', 'mask'): 'B.Mask',
    ('top', 'paste'): 'F.Paste',
    ('bottom', 'paste'): 'B.Paste',
}


class Converter:
    def __init__(self, input_path, out_dir, name=None, bom=None, cpl=None,
                 offset=(25.0, 25.0), gnd_heuristic=True, ocr=True):
        self.input_path = Path(input_path)
        self.out_dir = Path(out_dir)
        self.name = name
        self.bom_path = bom
        self.cpl_path = cpl
        self.offset = offset
        self.gnd_heuristic = gnd_heuristic
        self.ocr = ocr

        self.tracks = []
        self.vias = []
        self.zones = []
        self.pads = []          # flat list of Pad
        self.graphics = []
        self.edge_prims = []    # Edge.Cuts circle/arc/line primitives
        self.nets = {0: ''}     # net number -> name
        self.components = []    # list of Component
        self.silk_words = []
        self.silk_clusters = {'F': [], 'B': []}
        self.board_thickness = 1.6
        self.stats = Counter()
        self.warnings = []

    # -- loading -----------------------------------------------------------

    def load(self):
        path = self.input_path
        self._tmp = None
        if path.is_file() and path.suffix.lower() == '.zip':
            self._tmp = tempfile.TemporaryDirectory()
            with zipfile.ZipFile(path) as z:
                z.extractall(self._tmp.name)
            path = Path(self._tmp.name)
        if not path.is_dir():
            raise SystemExit(f'input {path} is not a directory or zip')
        self.stack = LayerStack.open(path)
        if self.name is None:
            self.name = self.stack.netlist and 'board' or None
        gl = dict(self.stack.graphic_layers)
        def to_mm(o):
            if str(o.unit) != MM:
                try:
                    return o.converted(MM)
                except TypeError:
                    o.convert_to(MM)
            return o

        self.layers = {k: [to_mm(o) for o in v.objects]
                       for k, v in gl.items()}
        self.drills = []
        for dl in self.stack.drill_layers:
            for o in dl.objects:
                o = to_mm(o)
                if isinstance(o, GnFlash):
                    self.drills.append(o)
        job = getattr(self.stack, 'board_thickness', None)
        # try gbrjob for thickness
        for f in (list(path.rglob('*.gbrjob')) or []):
            try:
                jd = json.loads(f.read_text())
                t = jd.get('Overall', {}).get('BoardThickness')
                if t:
                    self.board_thickness = float(t)
                if self.name is None:
                    self.name = jd.get('Overall', {}).get('Name', {}) \
                                  .get('ProjectId')
            except Exception:
                pass
        if self.name is None:
            self.name = self.input_path.stem
        # clear polarity is handled on copper (pour clearances); warn if it
        # shows up anywhere else
        n_clear = sum(1 for (side, use), objs in self.layers.items()
                      for o in objs
                      if use != 'copper'
                      and not getattr(o, 'polarity_dark', True))
        if n_clear:
            self.warnings.append(
                f'{n_clear} clear-polarity (LPC) objects on non-copper '
                f'layers ignored; those layers may show extra material')

    # -- coordinate transform ---------------------------------------------

    def compute_transform(self):
        xs, ys = [], []
        outline = self.layers.get(('mechanical', 'outline'), [])
        src = [o for o in outline
               if getattr(getattr(o, 'aperture', None), 'diameter', 1) == 0] \
              or outline
        if not src:
            src = [o for objs in self.layers.values() for o in objs]
        for o in src:
            (x0, y0), (x1, y1) = o.bounding_box()
            xs += [x0, x1]; ys += [y0, y1]
        self.minx, self.maxx = min(xs), max(xs)
        self.miny, self.maxy = min(ys), max(ys)

    def tx(self, x, y):
        return (round(x - self.minx + self.offset[0], 6),
                round(self.maxy - y + self.offset[1], 6))

    # -- conversion of graphic layers -------------------------------------

    def convert_layers(self):
        for key, objs in self.layers.items():
            if key == ('mechanical', 'outline'):
                self.convert_outline(objs)
                continue
            layer = LAYER_MAP.get(key)
            if layer is None:
                self.warnings.append(f'layer {key} not mapped, skipped')
                continue
            if layer.endswith('.Cu'):
                continue  # copper handled separately (needs nets)
            by_width = defaultdict(list)
            for o in objs:
                if isinstance(o, GnLine):
                    p1, p2 = self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)
                    w = getattr(o.aperture, 'diameter', 0) or 0
                    if p1 == p2:
                        if w > 0:
                            self.graphics.append(Graphic(
                                'circle', layer, fill=True, center=p1,
                                radius=w / 2))
                            self.stats[f'graphic_{layer}'] += 1
                    else:
                        by_width[round(w, 4)].append((p1, p2))
                    continue
                g = self.obj_to_graphic(o, layer)
                if g:
                    self.graphics.append(g)
                    self.stats[f'graphic_{layer}'] += 1
            # draws: chain, then recover circles/arcs/merged lines
            for w, segs in by_width.items():
                width = max(w, 0.01)
                self.stats['draw_segments_in'] += len(segs)
                for pts, closed in chain_segments(segs, COORD_EPS):
                    for prim in polyline_to_prims(pts, closed):
                        self.graphics.append(
                            self.prim_to_graphic(prim, layer, width))
                        self.stats[f'graphic_{layer}'] += 1
                        self.stats['draw_prims_out'] += 1

    def prim_to_graphic(self, prim, layer, width):
        if prim[0] == 'circle':
            self.stats['circles_recovered'] += 1
            return Graphic('circle', layer, width=width, fill=False,
                           center=prim[1], radius=prim[2])
        if prim[0] == 'arc':
            self.stats['arcs_recovered'] += 1
            return Graphic('arc', layer, width=width,
                           pts=[prim[1], prim[2], prim[3]])
        return Graphic('line', layer, width=width, pts=[prim[1], prim[2]])

    def obj_to_graphic(self, o, layer):
        if isinstance(o, GnLine):
            w = getattr(o.aperture, 'diameter', 0) or 0
            p1, p2 = self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)
            if p1 == p2:
                if w <= 0:
                    return None
                return Graphic('circle', layer, fill=True, center=p1,
                               radius=w / 2)
            return Graphic('line', layer, width=max(w, 0.01), pts=[p1, p2])
        if isinstance(o, GnArc):
            w = getattr(o.aperture, 'diameter', 0) or 0
            s, m, e = self.arc_points(o)
            return Graphic('arc', layer, width=max(w, 0.01), pts=[s, m, e])
        if isinstance(o, GnFlash):
            return self.flash_to_graphic(o, layer)
        if isinstance(o, GnRegion):
            pts = self.region_points(o)
            if len(pts) >= 3:
                return Graphic('poly', layer, fill=True, pts=pts)
        return None

    def flash_to_graphic(self, o, layer):
        ap = o.aperture
        c = self.tx(o.x, o.y)
        if isinstance(ap, gn_ap.CircleAperture):
            return Graphic('circle', layer, fill=True, center=c,
                           radius=ap.diameter / 2)
        if isinstance(ap, gn_ap.RectangleAperture):
            return Graphic('rect', layer, fill=True,
                           corners=(c[0] - ap.w / 2, c[1] - ap.h / 2,
                                    c[0] + ap.w / 2, c[1] + ap.h / 2))
        if isinstance(ap, gn_ap.PolygonAperture):
            pts = self.polygon_aperture_points(ap, c)
            return Graphic('poly', layer, fill=True, pts=pts)
        # fallback: bounding box
        (x0, y0), (x1, y1) = o.bounding_box()
        a, b = self.tx(x0, y0), self.tx(x1, y1)
        return Graphic('rect', layer, fill=True,
                       corners=(min(a[0], b[0]), min(a[1], b[1]),
                                max(a[0], b[0]), max(a[1], b[1])))

    def polygon_aperture_points(self, ap, center):
        n = ap.n_vertices
        rot = ap.rotation  # radians
        r = ap.diameter / 2
        pts = []
        for i in range(n):
            a = rot + 2 * math.pi * i / n
            # note: y axis is mirrored by the transform, mirror the angle too
            pts.append((round(center[0] + r * math.cos(a), 6),
                        round(center[1] - r * math.sin(a), 6)))
        return pts

    def arc_points(self, o):
        # gerbonara: cx, cy are offsets of the arc center relative to start
        cx, cy = o.x1 + o.cx, o.y1 + o.cy
        r = math.hypot(o.x1 - cx, o.y1 - cy)
        a1 = math.atan2(o.y1 - cy, o.x1 - cx)
        a2 = math.atan2(o.y2 - cy, o.x2 - cx)
        if o.clockwise:
            if a2 >= a1:
                a2 -= 2 * math.pi
        else:
            if a2 <= a1:
                a2 += 2 * math.pi
        am = (a1 + a2) / 2
        mid = (cx + r * math.cos(am), cy + r * math.sin(am))
        return self.tx(o.x1, o.y1), self.tx(*mid), self.tx(o.x2, o.y2)

    def region_points(self, o, simplify=REGION_SIMPLIFY):
        pts = []
        outline = list(o.outline)
        centers = list(o.arc_centers) if o.arc_centers else [None] * len(outline)
        for i, p in enumerate(outline):
            c = centers[i - 1] if i > 0 and i - 1 < len(centers) else None
            if c:
                # flatten arc from outline[i-1] to p around center c
                clockwise, (ccx, ccy) = c if isinstance(c, tuple) and \
                    len(c) == 2 and isinstance(c[0], bool) else (None, None)
                if ccx is None:
                    pts.append(p)
                    continue
                p0 = outline[i - 1]
                r = math.hypot(p0[0] - ccx, p0[1] - ccy)
                a1 = math.atan2(p0[1] - ccy, p0[0] - ccx)
                a2 = math.atan2(p[1] - ccy, p[0] - ccx)
                if clockwise and a2 >= a1:
                    a2 -= 2 * math.pi
                if not clockwise and a2 <= a1:
                    a2 += 2 * math.pi
                steps = max(2, int(abs(a2 - a1) * r / 0.05))
                for s in range(1, steps + 1):
                    a = a1 + (a2 - a1) * s / steps
                    pts.append((ccx + r * math.cos(a), ccy + r * math.sin(a)))
            else:
                pts.append(p)
        if len(pts) >= 3 and simplify:
            try:
                ring = Polygon(pts).buffer(0)
                if ring.geom_type == 'MultiPolygon':
                    ring = max(ring.geoms, key=lambda g: g.area)
                ring = ring.simplify(simplify)
                if not ring.is_empty and ring.exterior:
                    pts = list(ring.exterior.coords)[:-1]
            except Exception:
                pass
        return [self.tx(x, y) for x, y in pts]

    # -- outline -----------------------------------------------------------

    def convert_outline(self, objs):
        """The outline layer may mix the board edge, hole/slot cutouts and
        keepout/clearance decorations at several line widths.  Strategy:
          * chain segments per line width, keep closed loops
          * the loop with the largest area is the board edge (centerline)
          * zero-width loops are authoritative cutouts (overlapping circles
            forming slots get unioned)
          * non-zero-width loops duplicating the board edge or a zero-width
            cutout are clearance rings -> dropped
          * cutouts that only re-draw a drill hit are dropped (the hole is
            already represented by an NPTH pad)
        """
        by_width = defaultdict(list)
        for o in objs:
            w = getattr(getattr(o, 'aperture', None), 'diameter', 0) or 0
            if isinstance(o, GnLine):
                by_width[round(w, 4)].append(
                    (self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)))
            elif isinstance(o, GnArc):
                s, m, e = self.arc_points(o)
                by_width[round(w, 4)] += [(s, m), (m, e)]

        loops = []   # (shapely Polygon, width, ring pts)
        for w, segs in by_width.items():
            for pts, closed in chain_segments(segs, COORD_EPS):
                spts = [(round(x, 4), round(y, 4)) for x, y in pts]
                if closed and len(spts) >= 4:
                    ring = spts[:-1] if spts[0] == spts[-1] else spts
                    poly = Polygon(ring).buffer(0)
                    if poly.is_empty:
                        continue
                    if poly.geom_type == 'MultiPolygon':
                        poly = max(poly.geoms, key=lambda g: g.area)
                    loops.append((poly, w, ring))
                else:
                    for prim in polyline_to_prims(spts, False):
                        self.edge_prims.append(prim)
                    self.stats['edge_open_chains'] += 1

        if not loops:
            self.stats['edge_prims'] = len(self.edge_prims)
            return

        board = max(loops, key=lambda t: t[0].area)
        self.add_edge_ring(board[2])
        rest = [l for l in loops if l is not board]
        zero = [l for l in rest if l[1] == 0]
        nonzero = [l for l in rest if l[1] != 0]
        refs = [p for p, _, _ in zero] + [board[0]]
        kept_nz = []
        for p, w, ring in nonzero:
            if any(p.intersects(r) for r in refs):
                self.stats['edge_clearance_rings_dropped'] += 1
            else:
                kept_nz.append((p, w, ring))

        cut_polys = [p for p, _, _ in zero] + [p for p, _, _ in kept_nz]
        merged = unary_union(cut_polys) if cut_polys else None
        drill_pts = [(self.tx(d.x, d.y), d.aperture.diameter)
                     for d in self.drills]
        if merged is not None and not merged.is_empty:
            geoms = merged.geoms if hasattr(merged, 'geoms') else [merged]
            for g in geoms:
                cx, cy = g.centroid.x, g.centroid.y
                minx, miny, maxx, maxy = g.bounds
                dia = max(maxx - minx, maxy - miny)
                redundant = any(dist((cx, cy), c) < 0.3 and
                                abs(dia - d) < max(0.3, d * 0.2)
                                for c, d in drill_pts)
                if redundant:
                    self.stats['edge_drill_circles_dropped'] += 1
                    continue
                ring = [(round(x, 4), round(y, 4))
                        for x, y in g.exterior.coords[:-1]]
                if len(ring) >= 3:
                    self.add_edge_ring(ring)
        self.stats['edge_prims'] = len(self.edge_prims)

    def add_edge_ring(self, ring):
        """Store a closed board-edge ring as recovered circle/arc/line
        primitives."""
        pts = list(ring) + [ring[0]]
        for prim in polyline_to_prims(pts, True):
            self.edge_prims.append(prim)
            if prim[0] == 'circle':
                self.stats['circles_recovered'] += 1
            elif prim[0] == 'arc':
                self.stats['arcs_recovered'] += 1
        self.stats['edge_rings'] += 1

    # -- copper ------------------------------------------------------------

    def convert_copper(self):
        self.copper_flashes = {'F': [], 'B': []}
        for side, layer in (('top', 'F.Cu'), ('bottom', 'B.Cu')):
            sname = layer[0]
            pours = []          # shapely polygons, in draw order
            clear_after_draw = 0
            for o in self.layers.get((side, 'copper'), []):
                dark = getattr(o, 'polarity_dark', True)
                if not dark:
                    # subtractive object: punch it out of the pours drawn
                    # so far (the usual pour-with-clearances construction)
                    g = self.obj_to_shapely(o)
                    if g is not None and not g.is_empty:
                        pours = [p.difference(g) for p in pours]
                    continue
                if isinstance(o, GnLine):
                    w = getattr(o.aperture, 'diameter', 0) or 0.01
                    p1, p2 = self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)
                    if p1 == p2:
                        self.copper_flashes[sname].append(
                            (p1, 'circle', (w, w), None))
                    else:
                        self.tracks.append(Track(*p1, *p2, width=w,
                                                 layer=layer))
                elif isinstance(o, GnArc):
                    w = getattr(o.aperture, 'diameter', 0) or 0.01
                    s, m, e = self.arc_points(o)
                    self.tracks.append(TrackArc(s, m, e, width=w, layer=layer))
                elif isinstance(o, GnFlash):
                    self.copper_flashes[sname].append(
                        self.flash_to_padinfo(o))
                elif isinstance(o, GnRegion):
                    pts = self.region_points(o)
                    if len(pts) >= 3:
                        g = Polygon(pts).buffer(0)
                        if not g.is_empty:
                            pours.append(g)
            # merge overlapping pour fragments, keep connected pieces apart
            merged = unary_union(pours) if pours else None
            if merged is not None and not merged.is_empty:
                geoms = merged.geoms if hasattr(merged, 'geoms') else [merged]
                for g in geoms:
                    if g.area > 1e-6:
                        self.zones.append(ZonePoly(g.simplify(REGION_SIMPLIFY),
                                                   layer))
        self.stats['tracks'] = len(self.tracks)
        self.stats['zones'] = len(self.zones)

    def obj_to_shapely(self, o):
        """Transformed shapely geometry of a gerber object (for subtraction
        and connectivity tests)."""
        if isinstance(o, GnRegion):
            pts = self.region_points(o)
            return Polygon(pts).buffer(0) if len(pts) >= 3 else None
        if isinstance(o, GnLine):
            w = getattr(o.aperture, 'diameter', 0) or 0
            p1, p2 = self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)
            return LineString([p1, p2]).buffer(max(w, 0.001) / 2)
        if isinstance(o, GnArc):
            w = getattr(o.aperture, 'diameter', 0) or 0
            s, m, e = self.arc_points(o)
            return LineString([s, m, e]).buffer(max(w, 0.001) / 2)
        if isinstance(o, GnFlash):
            c, shape, size, ap = self.flash_to_padinfo(o)
            if shape == 'rect':
                return box(c[0] - size[0] / 2, c[1] - size[1] / 2,
                           c[0] + size[0] / 2, c[1] + size[1] / 2)
            return Point(c).buffer(max(size) / 2)
        return None

    def flash_to_padinfo(self, o):
        """-> (center, shape, (w,h), aperture) in KiCad coords."""
        ap = o.aperture
        c = self.tx(o.x, o.y)
        if isinstance(ap, gn_ap.CircleAperture):
            return (c, 'circle', (ap.diameter, ap.diameter), ap)
        if isinstance(ap, gn_ap.RectangleAperture):
            return (c, 'rect', (ap.w, ap.h), ap)
        if isinstance(ap, gn_ap.PolygonAperture):
            rot_deg = math.degrees(ap.rotation) % 45
            if ap.n_vertices == 8 and abs(rot_deg - 22.5) < 0.5:
                flat = ap.diameter * math.cos(math.pi / 8)
                return (c, 'octagon', (flat, flat), ap)
            d = ap.diameter
            return (c, 'circle', (d, d), ap)
        (x0, y0), (x1, y1) = o.bounding_box()
        return (c, 'rect', (x1 - x0, y1 - y0), ap)

    # -- drills, pads, vias ------------------------------------------------

    def build_pads_and_vias(self):
        mask_flashes = {'F': [], 'B': []}
        paste_flashes = {'F': [], 'B': []}
        for side, sname in (('top', 'F'), ('bottom', 'B')):
            for o in self.layers.get((side, 'mask'), []):
                if isinstance(o, GnFlash):
                    mask_flashes[sname].append(self.flash_to_padinfo(o))
            for o in self.layers.get((side, 'paste'), []):
                if isinstance(o, GnFlash):
                    paste_flashes[sname].append(self.flash_to_padinfo(o))

        def nearest(cands, c, tol):
            best, bd = None, tol
            for i, item in enumerate(cands):
                d = dist(item[0], c)
                if d < bd:
                    best, bd = i, d
            return best

        used_flash = {'F': set(), 'B': set()}
        used_mask = {'F': set(), 'B': set()}
        used_paste = {'F': set(), 'B': set()}

        # copper polygons for stitching-via detection
        zone_geoms = {'F': [], 'B': []}
        for z in self.zones:
            zone_geoms[z.layer[0]].append(z.poly)
        zone_union = {s: unary_union(gs) if gs else None
                      for s, gs in zone_geoms.items()}

        def capsule_pad(c, hole):
            """A hole with no flash may sit on a short thick draw: an
            obround (oval) pad.  Returns the capsule track or None."""
            best = None
            for t in self.tracks:
                if not isinstance(t, Track):
                    continue
                length = math.hypot(t.x2 - t.x1, t.y2 - t.y1)
                if length > 6 or t.width < hole * 0.8:
                    continue
                d = point_seg_dist(c, (t.x1, t.y1), (t.x2, t.y2))
                if d > max(t.width / 2 - hole / 4, 0.2):
                    continue
                if best is None or length < best[0]:
                    best = (length, t)
            return best[1] if best else None

        for d in self.drills:
            c = self.tx(d.x, d.y)
            hole = d.aperture.diameter
            fi = nearest(self.copper_flashes['F'], c, HOLE_MATCH_TOL)
            bi = nearest(self.copper_flashes['B'], c, HOLE_MATCH_TOL)
            plated = getattr(d.aperture, 'plated', None)
            if fi is None and bi is None:
                cap = capsule_pad(c, hole)
                if cap is not None:
                    # oval pad: remove all capsule draws of this pad from
                    # the track list (both copper sides), eat mask lines
                    mates = [t for t in self.tracks if isinstance(t, Track)
                             and abs(t.width - cap.width) < 0.4
                             and dist(((t.x1 + t.x2) / 2, (t.y1 + t.y2) / 2),
                                      ((cap.x1 + cap.x2) / 2,
                                       (cap.y1 + cap.y2) / 2)) < 0.3]
                    for t in mates:
                        self.tracks.remove(t)
                    mid = ((cap.x1 + cap.x2) / 2, (cap.y1 + cap.y2) / 2)
                    length = math.hypot(cap.x2 - cap.x1, cap.y2 - cap.y1)
                    rot = math.degrees(math.atan2(-(cap.y2 - cap.y1),
                                                  cap.x2 - cap.x1)) % 180
                    margin = None
                    for g in list(self.graphics):
                        if g.kind == 'line' and g.layer.endswith('.Mask'):
                            gm = ((g.pts[0][0] + g.pts[1][0]) / 2,
                                  (g.pts[0][1] + g.pts[1][1]) / 2)
                            if dist(gm, mid) < 1.0:
                                margin = round((g.width - cap.width) / 2, 3)
                                self.graphics.remove(g)
                    self.pads.append(Pad(
                        mid[0], mid[1], 'oval',
                        (length + cap.width, cap.width), 'TH', drill=hole,
                        mask_margin=margin if margin and margin > 0 else None,
                        rot=round(rot, 2),
                        seg=((cap.x1, cap.y1), (cap.x2, cap.y2), cap.width)))
                    self.stats['oval_pads'] += 1
                    continue
                inpour = all(
                    zone_union[s] is not None and
                    zone_union[s].contains(Point(c))
                    for s in ('F', 'B'))
                if plated or (plated is None and inpour):
                    self.vias.append(Via(c[0], c[1], size=hole + 0.3,
                                         drill=hole))
                    self.stats['vias_in_pour'] += 1
                else:
                    p = Pad(c[0], c[1], 'circle', (hole, hole), 'TH',
                            drill=hole, plated=False)
                    self.pads.append(p)
                    self.stats['npth'] += 1
                continue
            ref = self.copper_flashes['F'][fi] if fi is not None \
                else self.copper_flashes['B'][bi]
            _, shape, size, _ = ref
            if fi is not None:
                used_flash['F'].add(fi)
            if bi is not None:
                used_flash['B'].add(bi)
            mi_f = nearest(mask_flashes['F'], c, max(MASK_MATCH_TOL, 0.15))
            mi_b = nearest(mask_flashes['B'], c, max(MASK_MATCH_TOL, 0.15))
            has_mask = mi_f is not None or mi_b is not None
            if mi_f is not None:
                used_mask['F'].add(mi_f)
            if mi_b is not None:
                used_mask['B'].add(mi_b)
            if not has_mask and shape in ('circle', 'octagon') \
               and plated is not False:
                self.vias.append(Via(c[0], c[1], size=size[0], drill=hole))
                self.stats['vias'] += 1
            else:
                margin = None
                if mi_f is not None:
                    msize = mask_flashes['F'][mi_f][2]
                    margin = round((min(msize) - min(size)) / 2, 3)
                    if margin < 0:
                        margin = None
                self.pads.append(Pad(c[0], c[1], shape, size, 'TH',
                                     drill=hole, mask_margin=margin))
                self.stats['tht_pads'] += 1

        # remaining copper flashes = SMD pads (or spurious via pads)
        for sname in ('F', 'B'):
            for i, (c, shape, size, ap) in \
                    enumerate(self.copper_flashes[sname]):
                if i in used_flash[sname]:
                    continue
                mi = nearest(mask_flashes[sname], c, MASK_MATCH_TOL)
                margin = None
                if mi is not None:
                    used_mask[sname].add(mi)
                    msize = mask_flashes[sname][mi][2]
                    margin = round((min(msize) - min(size)) / 2, 3)
                    if margin < 0:
                        margin = None
                pi = nearest(paste_flashes[sname], c, MASK_MATCH_TOL)
                if pi is not None:
                    used_paste[sname].add(pi)
                self.pads.append(Pad(c[0], c[1], shape, size, sname,
                                     mask_margin=margin,
                                     paste=pi is not None))
                self.stats['smd_pads'] += 1

        # unmatched mask flashes -> plain graphics on the mask layer
        for sname, layer in (('F', 'F.Mask'), ('B', 'B.Mask')):
            for i, (c, shape, size, ap) in enumerate(mask_flashes[sname]):
                if i in used_mask[sname]:
                    continue
                self.graphics.append(self.padinfo_graphic(c, shape, size,
                                                          ap, layer))
                self.stats[f'graphic_{layer}'] += 1
            for i, (c, shape, size, ap) in enumerate(paste_flashes[sname]):
                if i in used_paste[sname]:
                    continue
                lay = sname + '.Paste'
                self.graphics.append(self.padinfo_graphic(c, shape, size,
                                                          ap, lay))
                self.stats[f'graphic_{lay}'] += 1

    def padinfo_graphic(self, c, shape, size, ap, layer):
        if shape == 'circle':
            return Graphic('circle', layer, fill=True, center=c,
                           radius=size[0] / 2)
        if shape == 'rect':
            return Graphic('rect', layer, fill=True,
                           corners=(c[0] - size[0] / 2, c[1] - size[1] / 2,
                                    c[0] + size[0] / 2, c[1] + size[1] / 2))
        if isinstance(ap, gn_ap.PolygonAperture):
            return Graphic('poly', layer, fill=True,
                           pts=self.polygon_aperture_points(ap, c))
        return Graphic('circle', layer, fill=True, center=c,
                       radius=size[0] / 2)

    def merge_copper_arcs(self):
        """Recover arcs / circles / merged lines from tessellated copper
        draws.  Runs after pad extraction so oval-pad capsules are gone."""
        by_key = defaultdict(list)
        others = []
        for t in self.tracks:
            if isinstance(t, Track):
                by_key[(t.layer, round(t.width, 4))].append(
                    ((t.x1, t.y1), (t.x2, t.y2)))
            else:
                others.append(t)
        new = list(others)
        n_in = sum(len(v) for v in by_key.values())
        for (layer, w), segs in by_key.items():
            for pts, closed in chain_segments(segs, COORD_EPS):
                for prim in polyline_to_prims(pts, closed):
                    if prim[0] == 'line':
                        new.append(Track(*prim[1], *prim[2], width=w,
                                         layer=layer))
                    elif prim[0] == 'arc':
                        self.stats['arcs_recovered'] += 1
                        new.append(TrackArc(prim[1], prim[2], prim[3],
                                            width=w, layer=layer))
                    else:                      # full circle: two half arcs
                        self.stats['circles_recovered'] += 1
                        (cx, cy), r = prim[1], prim[2]
                        p0, p1 = (cx + r, cy), (cx - r, cy)
                        new.append(TrackArc(p0, (cx, cy + r), p1,
                                            width=w, layer=layer))
                        new.append(TrackArc(p1, (cx, cy - r), p0,
                                            width=w, layer=layer))
        self.tracks = new
        self.stats['tracks'] = len(self.tracks)
        self.stats['copper_segments_merged'] = \
            n_in - sum(1 for t in new if isinstance(t, Track))

    # -- component grouping ------------------------------------------------

    def load_bom_cpl(self):
        self.cpl = {}
        self.bom = {}
        if self.cpl_path:
            with open(self.cpl_path, newline='', encoding='utf-8-sig') as f:
                for row in csv.DictReader(f):
                    row = {k.strip().lower(): v for k, v in row.items() if k}
                    ref = row.get('designator') or row.get('ref')
                    if not ref:
                        continue
                    try:
                        x = float(row.get('mid x', '').replace('mm', ''))
                        y = float(row.get('mid y', '').replace('mm', ''))
                    except ValueError:
                        continue
                    side = (row.get('layer') or 'Top').strip().lower()
                    rot = float(row.get('rotation') or 0)
                    self.cpl[ref.strip()] = {
                        'pos': self.tx(x, y),
                        'side': 'F' if side.startswith('t') else 'B',
                        'rot': rot}
        if self.bom_path:
            with open(self.bom_path, newline='', encoding='utf-8-sig') as f:
                for row in csv.DictReader(f):
                    row = {k.strip().lower(): v for k, v in row.items() if k}
                    refs = (row.get('designator') or '')
                    val = row.get('comment') or row.get('value') or ''
                    fp = row.get('footprint') or ''
                    for ref in re.split(r'[,;]\s*', refs):
                        if ref.strip():
                            self.bom[ref.strip()] = {'value': val,
                                                     'footprint': fp}

    # -- silkscreen analysis (component outlines + OCR of legend text) ----

    GLYPH_MAX_DIAG = 4.2     # mm, larger connected silk chains are outlines
    WORD_MERGE = 0.7         # mm, glyph blobs closer than this form a word
    PAD_TO_OUTLINE = 3.2     # mm, pad-to-silk-cluster association radius
    WORD_TO_COMP = 4.5       # mm, value word to component bbox margin
    LABEL_TO_PAD = 3.0       # mm, net label word to pad radius
    WORD_STOP = {'REMOVE', 'FLASH', 'TEST', 'GENERATOR', 'GENERA', 'TOR',
                 'TTERN', 'PATTER', 'CRT', 'CRI', 'THE', 'WWW', 'REV'}

    VALUE_RE = re.compile(
        r'^\d+(?:[KMR]\d*)?$|^\d+(?:[.,]\d+)?[UNP][FH]$', re.I)
    PART_RES = [
        (re.compile(r'^1N\d{3,4}', re.I), 'D'),
        (re.compile(r'^(SB|BA[TVS]|SS|SR)\d{2,4}', re.I), 'D'),
        (re.compile(r'^(78|79)\d\d$'), 'U'),
        (re.compile(r'^(MT|MC|LM|NE|TL|MAX|PT)\d{3,5}', re.I), 'U'),
    ]
    LABEL_RE = re.compile(r'^(\+\d+V?|[A-Z][A-Z0-9\-]{1,3})$')
    LABEL_STOP = {'ON', 'OFF', 'TOP', 'BOT', 'GBR', 'REV', 'THE', 'AND'}

    def analyze_silkscreen(self):
        """Build per-side component-outline clusters and OCR the legend."""
        self.silk_clusters = {'F': [], 'B': []}
        self.silk_words = []
        has_tess = self.ocr and shutil.which('tesseract')
        if self.ocr and not has_tess:
            self.warnings.append(
                'tesseract not found: silkscreen values/labels not read '
                '(component values will be empty)')
        for side, key in (('F', ('top', 'silk')), ('B', ('bottom', 'silk'))):
            segs = []
            for o in self.layers.get(key, []):
                if isinstance(o, GnLine):
                    segs.append((self.tx(o.x1, o.y1), self.tx(o.x2, o.y2)))
                elif isinstance(o, GnArc):
                    s, m, e = self.arc_points(o)
                    segs += [(s, m), (m, e)]
            if not segs:
                continue
            chains = chain_segments(segs, COORD_EPS)
            glyphs, outlines = [], []
            for pts, closed in chains:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                diag = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
                # closed rings above letter size are component outlines too
                # (e.g. small capacitor circles), open small chains = glyphs
                is_outline = diag > self.GLYPH_MAX_DIAG or \
                    (closed and diag > 2.8)
                (outlines if is_outline else glyphs).append(pts)
            if outlines:
                u = unary_union([LineString(p).buffer(0.25)
                                 for p in outlines if len(p) > 1])
                clusters = list(u.geoms) if hasattr(u, 'geoms') else [u]
                # merge nested outlines (a body box inside a courtyard box)
                hulls = [cl.convex_hull for cl in clusters]
                merged_into = {}
                for i, hi in enumerate(hulls):
                    for j, hj in enumerate(hulls):
                        if i == j or hj.area <= hi.area:
                            continue
                        inter = hi.intersection(hj).area
                        if hi.area and inter / hi.area > 0.8:
                            merged_into[i] = j
                            break
                def root(i):
                    while i in merged_into:
                        i = merged_into[i]
                    return i

                groups = defaultdict(list)
                for i, cl in enumerate(clusters):
                    groups[root(i)].append(cl)
                self.silk_clusters[side] = [unary_union(g)
                                            for g in groups.values()]
            self._glyphs = getattr(self, '_glyphs', {})
            self._glyphs[side] = glyphs
            if glyphs and has_tess:
                self.silk_words += self.ocr_glyphs(glyphs, side)
        self.stats['silk_words'] = len(self.silk_words)

    def ocr_glyphs(self, glyphs, side):
        """Group glyph stroke chains into words and OCR each crop."""
        try:
            import cairosvg
            from PIL import Image
        except ImportError:
            self.warnings.append('cairosvg/pillow missing: OCR skipped')
            return []
        blobs = unary_union([LineString(p).buffer(self.WORD_MERGE)
                             for p in glyphs if len(p) > 1])
        blobs = list(blobs.geoms) if hasattr(blobs, 'geoms') else [blobs]
        by_blob = defaultdict(list)
        for pts in glyphs:
            c = LineString(pts).centroid if len(pts) > 1 else None
            if c is None:
                continue
            for i, b in enumerate(blobs):
                if b.contains(c):
                    by_blob[i].append(pts)
                    break
        words = []
        scale = 30
        for i, strokes in by_blob.items():
            minx, miny, maxx, maxy = blobs[i].bounds
            w, h = maxx - minx + 1, maxy - miny + 1
            if w * h > 1500 or len(strokes) < 2:
                continue
            paths = []
            for pts in strokes:
                d = ' '.join(f'{"M" if j == 0 else "L"}{x:.2f},{y:.2f}'
                             for j, (x, y) in enumerate(pts))
                paths.append(f'<path d="{d}" stroke="black" '
                             f'stroke-width="0.18" fill="none" '
                             f'stroke-linecap="round"/>')
            # bottom silk is mirrored on the physical board: flip for OCR
            flip = (f'transform="translate({minx + maxx},0) scale(-1,1)"'
                    if side == 'B' else '')
            svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
                   f'viewBox="{minx - .5} {miny - .5} {w} {h}" '
                   f'width="{w * scale:.0f}" height="{h * scale:.0f}">'
                   f'<rect x="{minx - .5}" y="{miny - .5}" width="{w}" '
                   f'height="{h}" fill="white"/>'
                   f'<g {flip}>{"".join(paths)}</g></svg>')
            try:
                png = cairosvg.svg2png(bytestring=svg.encode())
            except Exception:
                continue
            im = Image.open(io.BytesIO(png))
            cands = [im]
            if h > w * 1.4:
                cands += [im.rotate(90, expand=True),
                          im.rotate(270, expand=True)]
            best_txt, best_conf = None, 25
            for cim in cands:
                txt, conf = ocr_image(cim)
                if txt and conf > best_conf:
                    best_txt, best_conf = txt, conf
            if best_txt:
                words.append(SilkWord((minx + maxx) / 2, (miny + maxy) / 2,
                                      best_txt, best_conf, side))
        return words

    # -- component grouping ------------------------------------------------

    def group_footprints(self):
        # 1) SMD pads near a CPL position -> that component's footprint
        comp_pads = defaultdict(list)
        for p in self.pads:
            if p.side in ('F', 'B'):
                best, bd = None, CPL_GROUP_RADIUS
                for ref, info in self.cpl.items():
                    if info['side'] != p.side:
                        continue
                    d = dist((p.x, p.y), info['pos'])
                    if d < bd:
                        best, bd = ref, d
                if best:
                    p.component = best
                    comp_pads[best].append(p)
        for ref, plist in comp_pads.items():
            plist.sort(key=lambda p: (round(p.y, 1), p.x))
            for i, p in enumerate(plist, 1):
                p.number = str(i)
            val = self.bom.get(ref, {}).get('value', '')
            self.components.append(Component(ref, val, plist, True, kind='U'))

        # 2) remaining plated pads -> nearest silkscreen outline cluster;
        #    pads deep inside a big closed outline (connector body) match by
        #    containment in the cluster's convex hull
        hulls = {s: [cl.convex_hull for cl in self.silk_clusters[s]]
                 for s in ('F', 'B')}
        cluster_pads = defaultdict(list)   # (side, idx) -> pads
        loose = []
        for p in self.pads:
            if p.component or not p.plated:
                continue
            pt = Point(p.x, p.y)
            best, bd = None, self.PAD_TO_OUTLINE
            for side in (('F', 'B') if p.side == 'TH' else (p.side,)):
                for i, cl in enumerate(self.silk_clusters[side]):
                    d = cl.distance(pt)
                    if d < bd:
                        best, bd = (side, i), d
            if best is None:
                cand = []
                for side in (('F', 'B') if p.side == 'TH' else (p.side,)):
                    for i, h in enumerate(hulls[side]):
                        if h.contains(pt):
                            cand.append((h.area, (side, i)))
                if cand:
                    best = min(cand)[1]
            if best:
                cluster_pads[best].append(p)
            else:
                loose.append(p)

        anon = []
        for (side, i), plist in cluster_pads.items():
            plist.sort(key=lambda p: (round(p.y, 1), p.x))
            for k, p in enumerate(plist, 1):
                p.number = str(k)
            kind = self.classify(plist)
            # 2 pads inside a closed round-ish outline = radial capacitor
            if kind == 'R' and all(
                    hulls[side][i].contains(Point(p.x, p.y))
                    for p in plist):
                kind = 'CP'
            anon.append(Component('?', '', plist, False, kind=kind))
        for p in loose:
            anon.append(Component('?', '', [p], p.side != 'TH', kind='TP'))

        # 3) values / part names from OCR, then reference assignment
        self.assign_values(anon)
        counters = Counter()
        for comp in sorted(anon, key=lambda c: (
                round(min(p.y for p in c.pads), 0),
                min(p.x for p in c.pads))):
            counters[comp.kind] += 1
            comp.ref = f'{comp.kind}{counters[comp.kind]}'
            self.components.append(comp)

        # 4) NPTH holes: plain mechanical footprints, not in the schematic
        h = 0
        for p in self.pads:
            if p.component or p.plated:
                continue
            h += 1
            p.component = f'H{h}'
            self.components.append(Component(
                f'H{h}', f'NPTH {fnum(p.drill)}mm', [p], False,
                kind='H', in_schematic=False))
        for comp in self.components:
            for p in comp.pads:
                p.component = comp.ref
        self.stats['footprints'] = len(self.components)
        self.stats['components_multi_pad'] = sum(
            1 for c in self.components if len(c.pads) > 1)

    def classify(self, plist):
        """Heuristic component class from pad pattern (no OCR needed)."""
        n = len(plist)
        if n == 1:
            return 'TP'
        if n == 2:
            return 'R'      # refined later if a value word says otherwise
        xs = sorted(p.x for p in plist)
        ys = sorted(p.y for p in plist)
        w, hgt = xs[-1] - xs[0], ys[-1] - ys[0]
        rows = len({round(p.y * 2) / 2 for p in plist})
        cols = len({round(p.x * 2) / 2 for p in plist})
        if n >= 4 and (rows <= 2 or cols <= 2) and max(w, hgt) > 7:
            return 'J'      # single/double row header or connector
        return 'U'

    def assign_values(self, comps):
        """Attach the nearest matching OCR word to each component and refine
        its class; leftover label-looking words rename nets."""
        def norm(t):
            t = t.upper().strip(' .:*()[]{}|').replace('$', 'S')
            return (t.replace('®', '0').replace('O', '0')
                    if re.match(r'^\d', t) else t)

        # candidate (component, word, distance) triples for value words
        cands = []
        for comp in comps:
            if comp.kind in ('TP', 'H'):
                continue
            minx = min(p.x for p in comp.pads) - self.WORD_TO_COMP
            maxx = max(p.x for p in comp.pads) + self.WORD_TO_COMP
            miny = min(p.y for p in comp.pads) - self.WORD_TO_COMP
            maxy = max(p.y for p in comp.pads) + self.WORD_TO_COMP
            cx = (minx + maxx) / 2
            cy = (miny + maxy) / 2
            for wo in self.silk_words:
                if not (minx <= wo.x <= maxx and miny <= wo.y <= maxy):
                    continue
                t = norm(wo.text)
                if t in self.WORD_STOP:
                    continue
                kind = None
                for rx, k in self.PART_RES:      # part numbers first: a
                    if rx.match(t):              # '7805' is not a resistor
                        kind = k
                        break
                if kind is None and self.VALUE_RE.match(t):
                    if t.endswith('F'):
                        kind = 'CP' if 'U' in t else 'C'
                    elif t.endswith('H'):
                        kind = 'L'
                    else:
                        kind = 'R'
                if kind is None:
                    continue
                cands.append((dist((wo.x, wo.y), (cx, cy)), comp, wo, t,
                              kind))
        cands.sort(key=lambda c: c[0])
        for d, comp, wo, t, kind in cands:
            if wo.claimed or comp.value:
                continue
            if len(comp.pads) > 4 and kind in ('R', 'C', 'CP', 'L', 'D'):
                continue    # a 21-pin connector is not a resistor
            wo.claimed = True
            comp.value = t
            if len(comp.pads) <= 4 or kind == 'U':
                comp.kind = kind
        # second chance for small parts without a value: re-OCR the area
        # inside/around the component with a restricted character set
        for comp in comps:
            if comp.value or comp.kind not in ('R', 'C', 'CP', 'D', 'L'):
                continue
            t = self.ocr_component_area(comp)
            if t:
                comp.value = t
                if t.endswith('F'):
                    comp.kind = 'CP' if comp.kind in ('CP', 'R') else 'C'
                elif t.endswith('H'):
                    comp.kind = 'L'

        # ON/OFF legend next to a small part = a switch
        for comp in comps:
            if comp.kind not in ('U', 'J') or len(comp.pads) > 6:
                continue
            cx = sum(p.x for p in comp.pads) / len(comp.pads)
            cy = sum(p.y for p in comp.pads) / len(comp.pads)
            near = {w.text.upper().strip(' .:*') for w in self.silk_words
                    if dist((w.x, w.y), (cx, cy)) < 8}
            if {'ON', 'OFF'} & near:
                comp.kind = 'SW'

        # connectors / unknown parts may carry a plain-word legend (SCART..)
        for comp in comps:
            if comp.value or comp.kind not in ('J', 'U', 'SW'):
                continue
            best = None
            bd = self.WORD_TO_COMP + 3
            minx = min(p.x for p in comp.pads)
            maxx = max(p.x for p in comp.pads)
            miny = min(p.y for p in comp.pads)
            maxy = max(p.y for p in comp.pads)
            for wo in self.silk_words:
                t = wo.text.upper().strip(' .:*')
                if wo.claimed or wo.conf < 55 or t in self.WORD_STOP \
                        or not re.match(r'^[A-Z]{4,10}$', t):
                    continue
                d = math.hypot(max(minx - wo.x, 0, wo.x - maxx),
                               max(miny - wo.y, 0, wo.y - maxy))
                if d < bd:
                    best, bd = wo, d
            if best:
                best.claimed = True
                comp.value = best.text.upper().strip(' .:*')
                if comp.value in ('PATTERN',):
                    comp.kind = 'SW'

    def ocr_component_area(self, comp):
        """Targeted OCR retry: render the glyph strokes lying within the
        component's bbox and read them with a value-only character set."""
        if not (self.ocr and shutil.which('tesseract')):
            return None
        try:
            import cairosvg
            from PIL import Image
        except ImportError:
            return None
        side = 'F' if comp.pads[0].side in ('F', 'TH') else 'B'
        minx = min(p.x for p in comp.pads)
        maxx = max(p.x for p in comp.pads)
        miny = min(p.y for p in comp.pads)
        maxy = max(p.y for p in comp.pads)
        # the value is printed between the two pads: inset the long axis to
        # stay clear of the neighbours, widen the short axis a little
        if maxx - minx >= maxy - miny:
            minx += 1.0; maxx -= 1.0
            miny -= 1.3; maxy += 1.3
        else:
            miny += 1.0; maxy -= 1.0
            minx -= 1.3; maxx += 1.3
        strokes = []
        for pts in self._glyphs.get(side, []):
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
            if minx <= cx <= maxx and miny <= cy <= maxy:
                strokes.append(pts)
        if not strokes:
            return None
        bxs = [p[0] for pts in strokes for p in pts]
        bys = [p[1] for pts in strokes for p in pts]
        minx, maxx = min(bxs) - .5, max(bxs) + .5
        miny, maxy = min(bys) - .5, max(bys) + .5
        w, h = maxx - minx, maxy - miny
        if w * h > 400:
            return None
        paths = []
        for pts in strokes:
            d = ' '.join(f'{"M" if j == 0 else "L"}{x:.2f},{y:.2f}'
                         for j, (x, y) in enumerate(pts))
            paths.append(f'<path d="{d}" stroke="black" stroke-width="0.18" '
                         f'fill="none" stroke-linecap="round"/>')
        flip = (f'transform="translate({minx + maxx},0) scale(-1,1)"'
                if side == 'B' else '')
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
               f'viewBox="{minx} {miny} {w} {h}" '
               f'width="{w * 30:.0f}" height="{h * 30:.0f}">'
               f'<rect x="{minx}" y="{miny}" width="{w}" height="{h}" '
               f'fill="white"/><g {flip}>{"".join(paths)}</g></svg>')
        try:
            png = cairosvg.svg2png(bytestring=svg.encode())
        except Exception:
            return None
        im = Image.open(io.BytesIO(png))
        cands = [im]
        if h > w * 1.4:
            cands += [im.rotate(90, expand=True), im.rotate(270, expand=True)]
        best, best_conf = None, 20
        for cim in cands:
            txt, conf = ocr_image(cim)
            if not txt:
                continue
            for tok in txt.upper().split():
                tok = tok.strip(' .:*()[]{}|').replace('$', 'S') \
                         .replace('®', '0')
                if self.VALUE_RE.match(tok) and conf > best_conf:
                    best, best_conf = tok, conf
                    break
        return best

    def name_nets_from_labels(self):
        """Silk words like GND / +12 / BLU next to a pad name its net."""
        renamed = {}
        used_names = set(self.nets.values())
        for wo in sorted(self.silk_words, key=lambda w: -w.conf):
            if wo.claimed or wo.conf < 60:
                continue
            t = wo.text.upper().strip(' .:*_')
            if not self.LABEL_RE.match(t) or t in self.LABEL_STOP:
                continue
            best, bd = None, self.LABEL_TO_PAD
            for p in self.pads:
                if not p.plated:
                    continue
                d = dist((p.x, p.y), (wo.x, wo.y))
                if d < bd:
                    best, bd = p, d
            if best is None or best.net == 0:
                continue
            cur = self.nets[best.net]
            if not cur.startswith('N$') or t in used_names:
                continue
            self.nets[best.net] = t
            used_names.add(t)
            renamed[t] = cur
        for z in self.zones:
            z.net_name = self.nets.get(z.net, z.net_name)
        self.stats['nets_named_from_silk'] = len(renamed)

    # -- connectivity ------------------------------------------------------

    def derive_nets(self):
        items = []   # (geom, layers(set), object, kind)
        for t in self.tracks:
            if isinstance(t, Track):
                g = LineString([(t.x1, t.y1), (t.x2, t.y2)]) \
                    .buffer(t.width / 2)
            else:
                g = LineString(sample_arc(t.start, t.mid, t.end)) \
                    .buffer(t.width / 2)
            items.append((g, {t.layer[0]}, t, 'track'))
        for z in self.zones:
            items.append((z.poly, {z.layer[0]}, z, 'zone'))
        for v in self.vias:
            g = Point(v.x, v.y).buffer(v.size / 2)
            items.append((g, {'F', 'B'}, v, 'via'))
        for p in self.pads:
            if not p.plated:
                continue
            if p.shape == 'oval' and p.seg:
                g = LineString([p.seg[0], p.seg[1]]).buffer(p.seg[2] / 2)
            elif p.shape == 'rect':
                g = box(p.x - p.size[0] / 2, p.y - p.size[1] / 2,
                        p.x + p.size[0] / 2, p.y + p.size[1] / 2)
            else:
                g = Point(p.x, p.y).buffer(max(p.size) / 2)
            lays = {'F', 'B'} if p.side == 'TH' else {p.side}
            items.append((g, lays, p, 'pad'))

        parent = list(range(len(items)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for side in ('F', 'B'):
            idx = [i for i, it in enumerate(items) if side in it[1]]
            geoms = [items[i][0] for i in idx]
            if not geoms:
                continue
            tree = STRtree(geoms)
            for k, gi in enumerate(idx):
                for m in tree.query(geoms[k], predicate='intersects'):
                    union(gi, idx[m])

        groups = defaultdict(list)
        for i in range(len(items)):
            groups[find(i)].append(i)

        # rank: GND = group with the largest zone area (if pours exist)
        def group_zone_area(members):
            return sum(items[i][0].area for i in members
                       if items[i][3] == 'zone')

        ordered = sorted(groups.values(),
                         key=lambda m: (-group_zone_area(m), -len(m)))
        netno = 0
        for gi, members in enumerate(ordered):
            netno += 1
            if self.gnd_heuristic and gi == 0 and group_zone_area(ordered[0]) \
               > 50 and any(items[i][3] == 'zone' for i in members):
                name = 'GND'
            else:
                name = f'N${netno}'
            self.nets[netno] = name
            for i in members:
                items[i][2].net = netno
                if items[i][3] == 'zone':
                    items[i][2].net_name = name
        self.stats['nets'] = netno

    # -- emit: board -------------------------------------------------------

    def emit_pcb(self):
        s = SExpr()
        s.open(f'kicad_pcb (version {PCB_VERSION}) (generator {q(GENERATOR)}) '
               f'(generator_version {q(GEN_VERSION)})')
        s.open('general')
        s.leaf(f'thickness {fnum(self.board_thickness)}')
        s.leaf('legacy_teardrops no')
        s.close()
        s.leaf('paper "A4"')
        s.open('layers')
        for line in (
                '0 "F.Cu" signal', '31 "B.Cu" signal',
                '32 "B.Adhes" user "B.Adhesive"',
                '33 "F.Adhes" user "F.Adhesive"',
                '34 "B.Paste" user', '35 "F.Paste" user',
                '36 "B.SilkS" user "B.Silkscreen"',
                '37 "F.SilkS" user "F.Silkscreen"',
                '38 "B.Mask" user', '39 "F.Mask" user',
                '40 "Dwgs.User" user "User.Drawings"',
                '41 "Cmts.User" user "User.Comments"',
                '42 "Eco1.User" user "User.Eco1"',
                '43 "Eco2.User" user "User.Eco2"',
                '44 "Edge.Cuts" user', '45 "Margin" user',
                '46 "B.CrtYd" user "B.Courtyard"',
                '47 "F.CrtYd" user "F.Courtyard"',
                '48 "B.Fab" user', '49 "F.Fab" user'):
            s.leaf(line)
        s.close()
        s.open('setup')
        s.leaf('pad_to_mask_clearance 0')
        s.leaf('allow_soldermask_bridges_in_footprints no')
        s.open('pcbplotparams')
        s.leaf('layerselection 0x00010fc_ffffffff')
        s.leaf('plot_on_all_layers_selection 0x0000000_00000000')
        s.leaf('disableapertmacros no')
        s.leaf('usegerberextensions no')
        s.leaf('usegerberattributes yes')
        s.leaf('usegerberadvancedattributes yes')
        s.leaf('creategerberjobfile yes')
        s.leaf('dashed_line_dash_ratio 12.000000')
        s.leaf('dashed_line_gap_ratio 3.000000')
        s.leaf('svgprecision 4')
        s.leaf('plotframeref no')
        s.leaf('viasonmask no')
        s.leaf('mode 1')
        s.leaf('useauxorigin no')
        s.leaf('hpglpennumber 1')
        s.leaf('hpglpenspeed 20')
        s.leaf('hpglpendiameter 15.000000')
        s.leaf('pdf_front_fp_property_popups yes')
        s.leaf('pdf_back_fp_property_popups yes')
        s.leaf('dxfpolygonmode yes')
        s.leaf('dxfimperialunits yes')
        s.leaf('dxfusepcbnewfont yes')
        s.leaf('psnegative no')
        s.leaf('psa4output no')
        s.leaf('plotreference yes')
        s.leaf('plotvalue yes')
        s.leaf('plotfptext yes')
        s.leaf('plotinvisibletext no')
        s.leaf('sketchpadsonfab no')
        s.leaf('subtractmaskfromsilk no')
        s.leaf('outputformat 1')
        s.leaf('mirror no')
        s.leaf('drillshape 1')
        s.leaf('scaleselection 1')
        s.leaf('outputdirectory ""')
        s.close()
        s.close()

        for no, name in sorted(self.nets.items()):
            s.leaf(f'net {no} {q(name)}')

        for comp in self.components:
            self.emit_footprint(s, comp.ref, comp.value, comp.pads,
                                comp.is_smd)

        for g in self.graphics:
            self.emit_graphic(s, g)

        for prim in self.edge_prims:
            if prim[0] == 'circle':
                (cx, cy), r = prim[1], prim[2]
                s.open('gr_circle')
                s.leaf(f'center {fnum(cx)} {fnum(cy)}')
                s.leaf(f'end {fnum(cx + r)} {fnum(cy)}')
                s.leaf('stroke (width 0.05) (type solid)')
                s.leaf('fill none')
            elif prim[0] == 'arc':
                s.open('gr_arc')
                s.leaf(f'start {fnum(prim[1][0])} {fnum(prim[1][1])}')
                s.leaf(f'mid {fnum(prim[2][0])} {fnum(prim[2][1])}')
                s.leaf(f'end {fnum(prim[3][0])} {fnum(prim[3][1])}')
                s.leaf('stroke (width 0.05) (type solid)')
            else:
                s.open('gr_line')
                s.leaf(f'start {fnum(prim[1][0])} {fnum(prim[1][1])}')
                s.leaf(f'end {fnum(prim[2][0])} {fnum(prim[2][1])}')
                s.leaf('stroke (width 0.05) (type solid)')
            s.leaf('layer "Edge.Cuts"')
            s.leaf(f'uuid {q(nid())}')
            s.close()

        for t in self.tracks:
            if isinstance(t, Track):
                s.open('segment')
                s.leaf(f'start {fnum(t.x1)} {fnum(t.y1)}')
                s.leaf(f'end {fnum(t.x2)} {fnum(t.y2)}')
                s.leaf(f'width {fnum(t.width)}')
                s.leaf(f'layer {q(t.layer)}')
                s.leaf(f'net {t.net}')
                s.leaf(f'uuid {q(nid())}')
                s.close()
            else:
                s.open('arc')
                s.leaf(f'start {fnum(t.start[0])} {fnum(t.start[1])}')
                s.leaf(f'mid {fnum(t.mid[0])} {fnum(t.mid[1])}')
                s.leaf(f'end {fnum(t.end[0])} {fnum(t.end[1])}')
                s.leaf(f'width {fnum(t.width)}')
                s.leaf(f'layer {q(t.layer)}')
                s.leaf(f'net {t.net}')
                s.leaf(f'uuid {q(nid())}')
                s.close()

        for v in self.vias:
            s.open('via')
            s.leaf(f'at {fnum(v.x)} {fnum(v.y)}')
            s.leaf(f'size {fnum(v.size)}')
            s.leaf(f'drill {fnum(v.drill)}')
            s.leaf('layers "F.Cu" "B.Cu"')
            s.leaf(f'net {v.net}')
            s.leaf(f'uuid {q(nid())}')
            s.close()

        for z in self.zones:
            outline = [(round(x, 4), round(y, 4))
                       for x, y in z.poly.exterior.coords[:-1]]
            if len(outline) < 3:
                continue
            s.open('zone')
            s.leaf(f'net {z.net}')
            s.leaf(f'net_name {q(z.net_name or self.nets.get(z.net, ""))}')
            s.leaf(f'layer {q(z.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.leaf('hatch edge 0.508')
            s.leaf('connect_pads (clearance 0.2)')
            s.leaf('min_thickness 0.1')
            s.leaf('filled_areas_thickness no')
            s.open('fill yes')
            s.leaf('thermal_gap 0.3')
            s.leaf('thermal_bridge_width 0.3')
            s.close()
            s.open('polygon')
            self.emit_pts(s, outline)
            s.close()
            for piece in fracture(z.poly):
                pts = [(round(x, 4), round(y, 4))
                       for x, y in piece.exterior.coords[:-1]]
                if len(pts) < 3:
                    continue
                s.open('filled_polygon')
                s.leaf(f'layer {q(z.layer)}')
                self.emit_pts(s, pts)
                s.close()
            s.close()

        s.close()
        return s.text()

    def emit_pts(self, s, pts):
        s.open('pts')
        for i in range(0, len(pts), 4):
            chunk = pts[i:i + 4]
            s.raw(' '.join(f'(xy {fnum(x)} {fnum(y)})' for x, y in chunk))
        s.close()

    def emit_graphic(self, s, g):
        if g.kind == 'line':
            s.open('gr_line')
            s.leaf(f'start {fnum(g.pts[0][0])} {fnum(g.pts[0][1])}')
            s.leaf(f'end {fnum(g.pts[1][0])} {fnum(g.pts[1][1])}')
            s.leaf(f'stroke (width {fnum(g.width)}) (type solid)')
            s.leaf(f'layer {q(g.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()
        elif g.kind == 'arc':
            s.open('gr_arc')
            s.leaf(f'start {fnum(g.pts[0][0])} {fnum(g.pts[0][1])}')
            s.leaf(f'mid {fnum(g.pts[1][0])} {fnum(g.pts[1][1])}')
            s.leaf(f'end {fnum(g.pts[2][0])} {fnum(g.pts[2][1])}')
            s.leaf(f'stroke (width {fnum(g.width)}) (type solid)')
            s.leaf(f'layer {q(g.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()
        elif g.kind == 'circle':
            s.open('gr_circle')
            s.leaf(f'center {fnum(g.center[0])} {fnum(g.center[1])}')
            s.leaf(f'end {fnum(g.center[0] + g.radius)} {fnum(g.center[1])}')
            s.leaf(f'stroke (width {fnum(g.width or 0)}) (type solid)')
            s.leaf('fill solid' if g.fill else 'fill none')
            s.leaf(f'layer {q(g.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()
        elif g.kind == 'rect':
            x1, y1, x2, y2 = g.corners
            s.open('gr_rect')
            s.leaf(f'start {fnum(x1)} {fnum(y1)}')
            s.leaf(f'end {fnum(x2)} {fnum(y2)}')
            s.leaf(f'stroke (width 0) (type solid)')
            s.leaf('fill solid' if g.fill else 'fill none')
            s.leaf(f'layer {q(g.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()
        elif g.kind == 'poly':
            if len(g.pts) < 3:
                return
            s.open('gr_poly')
            self.emit_pts(s, g.pts)
            s.leaf(f'stroke (width 0) (type solid)')
            s.leaf('fill solid' if g.fill else 'fill none')
            s.leaf(f'layer {q(g.layer)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()

    def emit_footprint(self, s, ref, val, plist, is_smd):
        cx = sum(p.x for p in plist) / len(plist)
        cy = sum(p.y for p in plist) / len(plist)
        side = plist[0].side
        fp_layer = 'B.Cu' if side == 'B' else 'F.Cu'
        name = f'g2k:{ref.replace("$", "_")}'
        s.open(f'footprint {q(name)}')
        s.leaf(f'layer {q(fp_layer)}')
        s.leaf(f'uuid {q(nid())}')
        s.leaf(f'at {fnum(cx)} {fnum(cy)}')
        s.leaf('attr ' + ('smd' if is_smd else 'through_hole'))
        # references go on the Fab layer: the original gerbers have no
        # reference silkscreen, the converted board must not add any
        s.open(f'property "Reference" {q(ref)}')
        s.leaf(f'at 0 -2.5 0')
        s.leaf(f'layer {q("B.Fab" if side == "B" else "F.Fab")}')
        s.leaf(f'uuid {q(nid())}')
        s.leaf('effects (font (size 1 1) (thickness 0.15))')
        s.close()
        s.open(f'property "Value" {q(val)}')
        s.leaf(f'at 0 2.5 0')
        s.leaf(f'layer {q("B.Fab" if side == "B" else "F.Fab")}')
        s.leaf(f'uuid {q(nid())}')
        s.leaf('effects (font (size 1 1) (thickness 0.15))')
        s.close()
        for p in plist:
            dx, dy = round(p.x - cx, 6), round(p.y - cy, 6)
            if p.side == 'TH':
                layers = '"*.Cu" "*.Mask"'
                ptype = 'thru_hole' if p.plated else 'np_thru_hole'
            else:
                lays = [f'"{p.side}.Cu"', f'"{p.side}.Mask"']
                if p.paste:
                    lays.insert(1, f'"{p.side}.Paste"')
                layers = ' '.join(lays)
                ptype = 'smd'
            if p.shape == 'octagon':
                shape = 'roundrect'
            elif p.shape in ('rect', 'oval'):
                shape = p.shape
            else:
                shape = 'circle'
            s.open(f'pad {q(p.number)} {ptype} {shape}')
            at = f'at {fnum(dx)} {fnum(dy)}'
            if p.rot:
                at += f' {fnum(p.rot)}'
            s.leaf(at)
            s.leaf(f'size {fnum(p.size[0])} {fnum(p.size[1])}')
            if p.drill:
                s.leaf(f'drill {fnum(p.drill)}')
            s.leaf(f'layers {layers}')
            if p.shape == 'octagon':
                s.leaf('roundrect_rratio 0')
                s.leaf('chamfer_ratio 0.2929')
                s.leaf('chamfer top_left top_right bottom_left bottom_right')
            s.leaf('remove_unused_layers no')
            if p.plated:
                s.leaf(f'net {p.net} {q(self.nets.get(p.net, ""))}')
            if p.mask_margin:
                s.leaf(f'solder_mask_margin {fnum(p.mask_margin)}')
            s.leaf(f'uuid {q(nid())}')
            s.close()
        s.close()

    # -- emit: schematic ---------------------------------------------------
    #
    # Each reconstructed component gets a real symbol (R / C / polarized C /
    # D / L / test point / generic box) laid out on a grid, with a global
    # net label on every pin.  The net names come from the copper
    # connectivity analysis (and silkscreen labels when OCR found some), so
    # the schematic *is* the reconstructed netlist.

    TWO_PIN_KINDS = ('R', 'C', 'CP', 'D', 'L')

    def emit_sch(self, root_uuid):
        s = SExpr()
        s.open(f'kicad_sch (version {SCH_VERSION}) (generator {q(GENERATOR)}) '
               f'(generator_version {q(GEN_VERSION)})')
        s.leaf(f'uuid {q(root_uuid)}')
        s.leaf('paper "A3"')

        comps = [c for c in self.components if c.in_schematic]

        def sort_key(c):
            order = {'U': 0, 'J': 1, 'SW': 2, 'D': 3, 'L': 4, 'CP': 5,
                     'C': 6, 'R': 7, 'TP': 8}
            m = re.match(r'([A-Z$]+)(\d+)', c.ref.replace('U$', 'U'))
            num = int(m.group(2)) if m else 0
            return (order.get(c.kind, 9), num)

        comps.sort(key=sort_key)

        used_kinds = sorted({self.symbol_name(c) for c in comps})
        s.open('lib_symbols')
        for name in used_kinds:
            self.emit_symbol_def(s, name)
        s.close()

        x0, y0 = 30.48, 40.64
        x, y = x0, y0
        row_h = 0
        for comp in comps:
            w, h = self.symbol_cell(comp)
            if x + w > 380:
                x = x0
                y += row_h
                row_h = 0
            self.emit_symbol_instance(s, comp, x + w / 2, y, root_uuid)
            x += w
            row_h = max(row_h, h)

        s.open('text')
        s.raw(q('Schematic reconstructed from Gerber files by gerber2kicad.\\n'
                'Nets = copper connectivity; component grouping = silkscreen '
                'outlines;\\nvalues + net names = silkscreen OCR (verify '
                'before use!).\\nPin numbers are geometric (top-left to '
                'bottom-right), not the\\nmanufacturer pinout.'))
        s.leaf('exclude_from_sim no')
        s.leaf('at 25.4 25.4 0')
        s.leaf('effects (font (size 2 2)) (justify left bottom)')
        s.leaf(f'uuid {q(nid())}')
        s.close()

        s.open('sheet_instances')
        s.leaf('path "/" (page "1")')
        s.close()
        s.close()
        return s.text()

    def symbol_name(self, comp):
        if comp.kind in self.TWO_PIN_KINDS and len(comp.pads) == 2:
            return f'G2K-{comp.kind}'
        if len(comp.pads) == 1:
            return 'G2K-TP'
        return f'G2K-BOX{len(comp.pads)}'

    def symbol_cell(self, comp):
        """(width, height) of the grid cell for a component instance."""
        name = self.symbol_name(comp)
        if name.startswith('G2K-BOX'):
            n = len(comp.pads)
            half = (n + 1) // 2
            return (55.88, half * 2.54 + 17.78)
        if name == 'G2K-TP':
            return (33.02, 12.7)
        return (40.64, 12.7)

    def emit_symbol_def(self, s, name):
        pins = []      # (num, x, y, angle, length)
        gfx = []       # raw s-expr strings
        if name == 'G2K-R':
            gfx.append('(rectangle (start -2.54 -1.016) (end 2.54 1.016) '
                       '(stroke (width 0.254) (type default)) '
                       '(fill (type none)))')
            pins = [('1', -3.81, 0, 0, 1.27), ('2', 3.81, 0, 180, 1.27)]
        elif name in ('G2K-C', 'G2K-CP'):
            gfx.append('(polyline (pts (xy -0.762 -2.032) (xy -0.762 2.032))'
                       ' (stroke (width 0.508) (type default)) '
                       '(fill (type none)))')
            gfx.append('(polyline (pts (xy 0.762 -2.032) (xy 0.762 2.032))'
                       ' (stroke (width 0.508) (type default)) '
                       '(fill (type none)))')
            if name == 'G2K-CP':
                gfx.append('(polyline (pts (xy -2.286 1.016) (xy -1.27 1.016)'
                           ') (stroke (width 0.254) (type default)) '
                           '(fill (type none)))')
                gfx.append('(polyline (pts (xy -1.778 0.508) (xy -1.778 '
                           '1.524)) (stroke (width 0.254) (type default)) '
                           '(fill (type none)))')
            pins = [('1', -3.81, 0, 0, 3.048), ('2', 3.81, 0, 180, 3.048)]
        elif name == 'G2K-D':
            gfx.append('(polyline (pts (xy -1.27 1.016) (xy -1.27 -1.016) '
                       '(xy 1.27 0) (xy -1.27 1.016)) '
                       '(stroke (width 0.254) (type default)) '
                       '(fill (type outline)))')
            gfx.append('(polyline (pts (xy 1.27 -1.016) (xy 1.27 1.016)) '
                       '(stroke (width 0.254) (type default)) '
                       '(fill (type none)))')
            pins = [('1', -3.81, 0, 0, 2.54), ('2', 3.81, 0, 180, 2.54)]
        elif name == 'G2K-L':
            for k in range(4):
                x1 = -2.54 + 1.27 * k
                gfx.append(f'(arc (start {fnum(x1)} 0) '
                           f'(mid {fnum(x1 + 0.635)} -0.635) '
                           f'(end {fnum(x1 + 1.27)} 0) '
                           f'(stroke (width 0.254) (type default)) '
                           f'(fill (type none)))')
            pins = [('1', -3.81, 0, 0, 1.27), ('2', 3.81, 0, 180, 1.27)]
        elif name == 'G2K-TP':
            gfx.append('(circle (center 0 0) (radius 0.762) '
                       '(stroke (width 0.254) (type default)) '
                       '(fill (type none)))')
            pins = [('1', -3.81, 0, 0, 3.048)]
        else:                       # G2K-BOXn
            n = int(name[7:])
            half = (n + 1) // 2
            h = half * 1.27 + 2.54
            gfx.append(f'(rectangle (start -7.62 {fnum(h)}) '
                       f'(end 7.62 {fnum(-h)}) '
                       f'(stroke (width 0.254) (type default)) '
                       f'(fill (type background)))')
            for i in range(n):
                left = i < half
                row = i if left else i - half
                py = (half - 1) * 1.27 - row * 2.54
                pins.append((str(i + 1), -10.16 if left else 10.16, py,
                             0 if left else 180, 2.54))
        s.open(f'symbol {q("g2k:" + name)}')
        s.leaf('exclude_from_sim no')
        s.leaf('in_bom yes')
        s.leaf('on_board yes')
        s.open('property "Reference" "U"')
        s.leaf('at 0 3.81 0')
        s.leaf('effects (font (size 1.27 1.27))')
        s.close()
        s.open('property "Value" ""')
        s.leaf('at 0 -3.81 0')
        s.leaf('effects (font (size 1.27 1.27))')
        s.close()
        s.open(f'symbol {q(name + "_0_1")}')
        for g in gfx:
            s.raw(g)
        s.close()
        s.open(f'symbol {q(name + "_1_1")}')
        for num, px, py, ang, plen in pins:
            s.open('pin passive line')
            s.leaf(f'at {fnum(px)} {fnum(py)} {ang}')
            s.leaf(f'length {fnum(plen)}')
            s.leaf(f'name {q("~")} (effects (font (size 1.27 1.27)))')
            s.leaf(f'number {q(num)} (effects (font (size 1.27 1.27)))')
            s.close()
        s.close()
        s.close()

    def symbol_pin_positions(self, comp):
        """[(pad_number, dx, dy, side)] in symbol coords (y up)."""
        name = self.symbol_name(comp)
        if name.startswith('G2K-BOX'):
            n = len(comp.pads)
            half = (n + 1) // 2
            out = []
            for i in range(n):
                left = i < half
                row = i if left else i - half
                py = (half - 1) * 1.27 - row * 2.54
                out.append((str(i + 1), -10.16 if left else 10.16, py,
                            'L' if left else 'R'))
            return out
        if name == 'G2K-TP':
            return [('1', -3.81, 0, 'L')]
        return [('1', -3.81, 0, 'L'), ('2', 3.81, 0, 'R')]

    def emit_symbol_instance(self, s, comp, x, y, root_uuid):
        name = self.symbol_name(comp)
        pinpos = self.symbol_pin_positions(comp)
        half_h = max(abs(py) for _, _, py, _ in pinpos) + 2.54 \
            if pinpos else 2.54
        s.open(f'symbol (lib_id {q("g2k:" + name)})')
        s.leaf(f'at {fnum(x)} {fnum(y)} 0')
        s.leaf('unit 1')
        s.leaf('exclude_from_sim no')
        s.leaf('in_bom yes')
        s.leaf('on_board yes')
        s.leaf('dnp no')
        s.leaf(f'uuid {q(nid())}')
        s.open(f'property "Reference" {q(comp.ref)}')
        s.leaf(f'at {fnum(x)} {fnum(y - half_h - 2.54)} 0')
        s.leaf('effects (font (size 1.27 1.27))')
        s.close()
        s.open(f'property "Value" {q(comp.value)}')
        s.leaf(f'at {fnum(x)} {fnum(y + half_h + 2.54)} 0')
        s.leaf('effects (font (size 1.27 1.27))')
        s.close()
        s.open(f'property "Footprint" '
               f'{q("g2k:" + comp.ref.replace("$", "_"))}')
        s.leaf(f'at {fnum(x)} {fnum(y)} 0')
        s.leaf('effects (font (size 1.27 1.27)) (hide yes)')
        s.close()
        for num, _, _, _ in pinpos:
            s.leaf(f'pin {q(num)} (uuid {q(nid())})')
        s.open('instances')
        s.open(f'project {q(self.name)}')
        s.leaf(f'path {q("/" + root_uuid)} (reference {q(comp.ref)}) '
               f'(unit 1)')
        s.close()
        s.close()
        s.close()
        # one local net label per pin (sheet Y grows downward -> flip dy);
        # the schematic is a single sheet, so local labels are enough to
        # carry the reconstructed netlist
        pads_by_num = {p.number: p for p in comp.pads}
        for num, dx, dy, side in pinpos:
            pad = pads_by_num.get(num)
            if pad is None:
                continue
            net_name = self.nets.get(pad.net) or f'N${pad.net}'
            px, py = x + dx, y - dy
            s.open(f'label {q(net_name)}')
            s.leaf(f'at {fnum(px)} {fnum(py)} {180 if side == "L" else 0}')
            s.leaf('fields_autoplaced yes')
            s.leaf('effects (font (size 1.27 1.27)) '
                   f'(justify {"right" if side == "L" else "left"} bottom)')
            s.leaf(f'uuid {q(nid())}')
            s.close()


    def emit_pro(self, root_uuid):
        return json.dumps({
            'board': {'3dviewports': [], 'design_settings': {},
                      'layer_presets': [], 'viewports': []},
            'boards': [],
            'cvpcb': {'equivalence_files': []},
            'libraries': {'pinned_footprint_libs': [],
                          'pinned_symbol_libs': []},
            'meta': {'filename': f'{self.name}.kicad_pro', 'version': 1},
            'net_settings': {'classes': [{
                'bus_width': 12, 'clearance': 0.2, 'diff_pair_gap': 0.25,
                'diff_pair_via_gap': 0.25, 'diff_pair_width': 0.2,
                'line_style': 0, 'microvia_diameter': 0.3,
                'microvia_drill': 0.1, 'name': 'Default',
                'pcb_color': 'rgba(0, 0, 0, 0.000)',
                'schematic_color': 'rgba(0, 0, 0, 0.000)',
                'track_width': 0.25, 'via_diameter': 0.8,
                'via_drill': 0.4, 'wire_width': 6}],
                'meta': {'version': 3},
                'net_colors': None, 'netclass_assignments': None,
                'netclass_patterns': []},
            'pcbnew': {'last_paths': {}, 'page_layout_descr_file': ''},
            'schematic': {'legacy_lib_dir': '', 'legacy_lib_list': []},
            'sheets': [[root_uuid, 'Root']],
            'text_variables': {},
        }, indent=2)

    # -- driver ------------------------------------------------------------

    def run(self):
        print(f'Loading {self.input_path} ...')
        self.load()
        self.compute_transform()
        print(f'Board bbox: {self.maxx - self.minx:.2f} x '
              f'{self.maxy - self.miny:.2f} mm')
        self.convert_layers()
        self.convert_copper()
        self.build_pads_and_vias()
        self.merge_copper_arcs()
        self.load_bom_cpl()
        self.analyze_silkscreen()
        self.group_footprints()
        self.derive_nets()
        self.name_nets_from_labels()

        self.out_dir.mkdir(parents=True, exist_ok=True)
        root_uuid = nid()

        def write_atomic(path, text):
            tmp = path.with_suffix(path.suffix + '.tmp')
            tmp.write_text(text)
            tmp.replace(path)

        write_atomic(self.out_dir / f'{self.name}.kicad_pcb',
                     self.emit_pcb())
        write_atomic(self.out_dir / f'{self.name}.kicad_sch',
                     self.emit_sch(root_uuid))
        write_atomic(self.out_dir / f'{self.name}.kicad_pro',
                     self.emit_pro(root_uuid))

        print(f'\nWrote project "{self.name}" to {self.out_dir}/')
        for k in sorted(self.stats):
            print(f'  {k}: {self.stats[k]}')
        print(f'  nets: {self.stats["nets"]} '
              f'(names: {", ".join(list(self.nets.values())[1:6])} ...)')
        for w in self.warnings:
            print(f'  WARNING: {w}')
        report = {
            'name': self.name,
            'board_size_mm': [round(self.maxx - self.minx, 3),
                              round(self.maxy - self.miny, 3)],
            'stats': dict(self.stats),
            'nets': self.nets,
            'components': [
                {'ref': c.ref, 'value': c.value, 'kind': c.kind,
                 'pads': len(c.pads),
                 'nets': sorted({self.nets.get(p.net, '') for p in c.pads})}
                for c in self.components if c.in_schematic],
            'warnings': self.warnings,
        }
        (self.out_dir / 'conversion_report.json') \
            .write_text(json.dumps(report, indent=2))
        return report


def fracture(poly, slit=0.0005):
    """Split a shapely Polygon with holes into hole-free polygons by cutting
    thin vertical slits through each hole."""
    todo = [poly]
    result = []
    guard = 0
    while todo and guard < 10000:
        guard += 1
        p = todo.pop()
        if p.is_empty or p.area < 1e-9:
            continue
        if not p.interiors:
            result.append(p)
            continue
        hx = Polygon(p.interiors[0]).representative_point().x
        minx, miny, maxx, maxy = p.bounds
        cut = box(hx - slit, miny - 1, hx + slit, maxy + 1)
        pieces = p.difference(cut)
        pieces = pieces.geoms if hasattr(pieces, 'geoms') else [pieces]
        for g in pieces:
            if g.geom_type == 'Polygon':
                todo.append(g)
    return result


def chain_segments(segs, eps):
    """Chain 2-point segments into polylines. -> [(points, closed)]"""
    def key(p):
        return (round(p[0] / eps), round(p[1] / eps))

    adj = defaultdict(list)
    for i, (a, b) in enumerate(segs):
        adj[key(a)].append((i, a, b))
        adj[key(b)].append((i, b, a))
    used = set()
    chains = []
    for i, (a, b) in enumerate(segs):
        if i in used:
            continue
        used.add(i)
        pts = [a, b]
        # extend forward then backward
        for reverse in (False, True):
            if reverse:
                pts.reverse()
            while True:
                tail = pts[-1]
                nxt = None
                for j, p, o in adj[key(tail)]:
                    if j not in used:
                        nxt = (j, o)
                        break
                if nxt is None:
                    break
                used.add(nxt[0])
                pts.append(nxt[1])
                if key(pts[-1]) == key(pts[0]):
                    break
        closed = key(pts[0]) == key(pts[-1]) and len(pts) > 3
        if closed:
            pts[-1] = pts[0]
        chains.append((pts, closed))
    return chains


def main():
    ap = argparse.ArgumentParser(
        description='Convert Gerber/Excellon files to a KiCad project')
    ap.add_argument('input', help='directory or zip with Gerber files')
    ap.add_argument('-o', '--output', required=True, help='output directory')
    ap.add_argument('--name', help='project name (default: from gbrjob/zip)')
    ap.add_argument('--bom', help='BOM csv (JLCPCB style)')
    ap.add_argument('--cpl', help='CPL/placement csv (JLCPCB style)')
    ap.add_argument('--no-gnd-heuristic', action='store_true',
                    help='do not rename the biggest pour net to GND')
    ap.add_argument('--no-ocr', action='store_true',
                    help='skip silkscreen OCR (component values / net names)')
    args = ap.parse_args()
    conv = Converter(args.input, args.output, name=args.name,
                     bom=args.bom, cpl=args.cpl,
                     gnd_heuristic=not args.no_gnd_heuristic,
                     ocr=not args.no_ocr)
    conv.run()


if __name__ == '__main__':
    main()
