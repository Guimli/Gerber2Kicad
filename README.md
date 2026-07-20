# gerber2kicad

Converts a set of Gerber (RS-274X) + Excellon fabrication files into a
**KiCad** project you can open directly (KiCad 8 and 9): `.kicad_pro`,
`.kicad_pcb` and `.kicad_sch`.

Gerber files carry no netlist, no footprints and no components — the
converter reconstructs what it can, heuristically, from the copper
geometry, the drill file, the silkscreen artwork (including OCR of the
printed text), and an optional BOM/CPL pair.

## What gets reconstructed

| Gerber source | KiCad result |
|---|---|
| Copper draws (D01, round aperture) | Tracks (`segment` / `arc`) |
| Copper flashes (D03) | Pads (circle, rect, octagon via chamfered roundrect) |
| Copper regions (G36/G37, composed dark/clear polarity) | Filled zones (copper pour with clearances) |
| Excellon drill hits | Through-hole pads, vias, or NPTH holes (heuristic, see below) |
| Silkscreen / soldermask / paste | Graphic primitives on the matching layers; **circles and arcs tessellated into short segments are recovered** as real `gr_circle`/`gr_arc` objects (least-squares circle fit with a chord-sagitta check), and collinear segments are merged — about 95% fewer primitives than a literal segment-for-segment import |
| Board outline layer (GKO/GML…) | `Edge.Cuts` rebuilt as lines/arcs/circles (rounded corners become arcs), overlapping-circle slots merged, duplicate clearance rings dropped |
| Copper draws that are really arcs | KiCad `arc` tracks instead of hundreds of tiny straight segments |
| Copper connectivity | **Nets** derived from geometric analysis (shapely); the net with the largest copper pour is named `GND` (disable with `--no-gnd-heuristic`) |
| BOM/CPL CSV pair (JLCPCB-style, optional) | SMD pads near a CPL placement are grouped into one footprint per component, tagged with the BOM reference and value |
| Silkscreen outlines | **Pad-to-component grouping**: every silkscreened component outline (box, circle…) claims the pads that sit inside it, producing multi-pad footprints on the PCB and matching symbols in the schematic |
| Silkscreen text (OCR) | Vectorized text is isolated (stroke chaining) and read with **Tesseract**: component values (`390`, `2K2`, `10uF`, `22uH`…), part numbers (`1N4004`, `SB120`, `MT3608`, `7805` → classified as D/U), connector legends (`SCART`, `AUDIO`), and **net names** (`GND`, `BH`, `RM`… applied to the nearest pad's net) |

## Reconstructed schematic

Every component gets a real symbol matching its class — resistor (R),
capacitor (C) / polarized (CP), diode (D), inductor (L), generic N-pin
connector or IC (J/U/SW), test point (TP) — laid out on a grid, with a
**local net label on every pin**: the schematic *is* the netlist
reconstructed from the PCB copper. Classification comes from the pad
pattern (two pads inside a circle = radial capacitor, pins on a 2.54 mm
row = connector…) and from OCR (a value ending in `uF` → capacitor, `1N…`
/`SB…` → diode, `78xx`/`MTxxxx` → IC, an `ON`/`OFF` legend nearby →
switch). Orphan pads become test points; non-plated holes (`H*`) are kept
out of the schematic.

⚠️ Verify before reuse: pin numbering is purely geometric (not the
manufacturer's pinout), OCR'd values can contain misreads, and two nets
that are only tied together through a component (rather than by copper)
stay separate.

### Drill-hole classification heuristic

- copper flash + matching soldermask opening → plated through-hole pad;
- copper flash with no soldermask opening → via (best effort);
- no flash, but a short, thick copper draw covers the hole → **oblong
  (oval) pad**, the usual Eagle way of drawing connector pins;
- no flash, but the hole sits inside a copper pour on both sides → via
  stitching a ground plane;
- otherwise → non-plated hole (NPTH).

## Installation

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Usage

```bash
.venv/bin/python gerber2kicad.py GERBERS -o OUTDIR \
    [--name NAME] [--bom bom.csv] [--cpl cpl.csv] \
    [--no-gnd-heuristic] [--no-ocr]
```

For silkscreen OCR (component values, net names) also install:
`sudo apt install tesseract-ocr` and `pip install cairosvg pillow`.
Without Tesseract the conversion still runs, just with empty values.

`GERBERS` can be a directory **or a .zip**. Layers are identified
automatically (Protel extensions `.GTL/.GBL/...`, common layer names,
`.gbrjob`). Board thickness and project name are read from the `.gbrjob`
file when present.

Example, using the test fixture from
[Portable_CRT_Test_Pattern_Generator](https://github.com/baritonomarchetto/Portable_CRT_Test_Pattern_Generator):

```bash
.venv/bin/python gerber2kicad.py testdata/gerbers.zip -o output \
    --bom testdata/CRT_pattern_gen_top_bom.csv \
    --cpl testdata/CRT_pattern_gen_top_cpl.csv
```

A `conversion_report.json` file (stats, nets, component list, warnings) is
written to the output directory alongside the KiCad project.

## Validation (no KiCad required)

```bash
.venv/bin/python validate.py testdata/gerbers output/CRT_pattern_gen.kicad_pcb
```

`validate.py` re-parses the generated `.kicad_pcb` with an independent
s-expression reader, rebuilds the copper as shapely geometry, and compares
it against copper computed straight from the original Gerbers: area per
side, bounding box, **symmetric-difference (XOR) of the shapes after
alignment**, and object counts. It also renders side-by-side SVGs into
`OUTDIR/validation/`.

Result on the test project: copper area within −0.05%, shape XOR
0.14–0.18% (tolerances come from polygon simplification and from
approximating octagonal pads as chamfered squares).

## Known limitations

- Nets are purely geometric: two nets that are only tied together off-board
  (or through a component) stay separate; net names not read from
  silkscreen are arbitrary (`N$k`).
- Silkscreen text on the PCB stays vectorized (thousands of tiny
  primitives), not editable as text.
- Component references (`R1`, `C1`…) are generated (top-left to
  bottom-right reading order): Eagle-exported Gerbers only print values on
  silkscreen, not the original references.
- Pad/pin numbering is arbitrary (geometric order) — check before reuse.
- OCR can misread some values (e.g. `100` occasionally read as `108`) —
  check before reuse.
- If you re-fill zones in KiCad (`B` key), the fill is recomputed with
  KiCad's own rules and can differ slightly from the original Gerber (the
  originally-computed fill is what ships in the generated file).
- Clear polarity (LPC) is only handled on copper (pour clearances), not on
  other layers.

## Dependencies

- [gerbonara](https://pypi.org/project/gerbonara/) — RS-274X/Excellon parsing
- [shapely](https://pypi.org/project/shapely/) — geometry (nets, zones, XOR)
- [Tesseract OCR](https://github.com/tesseract-ocr/tesseract) (optional,
  system binary) — reading component values and net names off silkscreen
- `cairosvg`, `pillow` (optional) — rasterizing silkscreen crops for OCR
  and validation SVGs
