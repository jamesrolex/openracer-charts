#!/usr/bin/env python3
"""Extract bundled depth contours for the Abersoch racing area.

**What this makes:** `assets/depth/abersoch-bathymetry.json` — a single
GeoJSON FeatureCollection with depth contour lines, a shallow-water
polygon tint, the drying-ground tint (B-1857) and depth-label anchor
points, bundled into the app as a static asset. At runtime there is NO provider and NO network — the
chart just renders this file. Offline-first by construction.

**Data source:** EMODnet Bathymetry DTM 2024 (open data, CC BY 4.0),
fetched from the EMODnet ERDDAP as a plain CSV grid slice:

    https://erddap.emodnet.eu/erddap/griddap/bathymetry_dtm_2024

Grid resolution is 1/16 arc-minute (~115 m N-S, ~70 m E-W at 52.8 N).
Elevation is metres relative to Lowest Astronomical Tide — the SAME
datum as the Admiralty paper charts, so the contours here read like
the paper chart's. Negative = below LAT (sea), positive = above LAT
(drying / land), blank = no data (land interior).

**Attribution required (CC BY 4.0):** (c) EMODnet Bathymetry Consortium.
DOI: https://doi.org/10.12770/cf51df64-56f9-4a99-b1aa-36b8d7b743a1
The app surfaces this in the chart footer whenever the layer renders.

**Not for navigation.** EMODnet's own terms: the grid "should not be
used for navigation or any purpose relating to safety at sea". This
layer exists to read tidal streams over shoals, not to keep the keel
off the bottom.

Dependencies (documented, build-time only — nothing ships in the app):

    python3 -m venv venv && venv/bin/pip install numpy contourpy
    venv/bin/python extract-depth-contours.py

Regeneration: edit the constants below, re-run, commit the new asset.
The ERDDAP download is cached in `cache/` (gitignored) so re-runs with
new contour settings cost no network.

**Any box, for a chart area the helm cut himself (B-1845):**

    venv/bin/python extract-depth-contours.py --bbox=W,S,E,N --out=path.json

Same contour settings as the bundled coasts, same JSON shape, so the
app reads the file with the parser it already has. EMODnet covers
European waters only. **A box it does not cover writes nothing and
exits 0**, so the cut job never fails for want of depth.

    venv/bin/python extract-depth-contours.py --self-test

builds contours from a tiny synthetic grid. No network.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from contourpy import FillType, LineType, contour_generator

# ── Regions ──────────────────────────────────────────────────────
#
# **One region, one asset, one box.** Boxes match the chart regions in
# `src/domain/chartBasemapAvailability.ts` so a helm who has the
# basemap for a coast also has its depth. Adding a region here means
# adding it to `src/domain/depthAssets.ts` too, or the app will never
# load the file.
#
# `abersoch` is DELIBERATELY narrower than the Cardigan Bay basemap
# box (52.0-53.0 N, -4.85..-4.0). It predates the region catalogue and
# is the racing area, not the cruising ground. Widening it re-cuts a
# shipped asset for no racing gain, so it stays as extracted.

REGIONS: dict[str, dict] = {
    # Abersoch, St Tudwal's Roads + Islands, the Gimblet shoal ground,
    # out to the Tremadog Bay approaches. The original box.
    "abersoch": {
        "bounds": (52.72, 52.95, -4.75, -4.20),
        "out": "abersoch-bathymetry.json",
    },
    # Exe, Teign, Torbay, Dart, Start Bay, Salcombe. Matches
    # SOUTH_DEVON. **The estuary coast**, and the reason this script
    # grew a region table at all: the first non-family tester keeps
    # his boat here and the app drew him no depth whatsoever.
    "south-devon": {
        "bounds": (50.15, 50.75, -3.9, -3.2),
        "out": "south-devon-bathymetry.json",
    },
    # Needles and Hurst round to Selsey, south past St Catherine's,
    # north to Southampton Water and the Hamble. Matches SOLENT.
    #
    # **`band_tol_m` is the B-1103 knob.** The Solent seabed is the
    # busiest of the three coasts (490 contour lines against
    # Abersoch's 146), and its four-band cut at the shared 20 m
    # tolerance is 1,590,175 bytes — over the 1.5 MB budget. Bands
    # are broad tints under a contour line of at least 1 px drawn on
    # the same locus, so band edge slop is invisible where contour
    # slop would not be; the LINES keep the shared 20 m untouched.
    # Raised stepwise until the asset fit: 25 m → 1,570,283 bytes,
    # 30 m → 1,546,410, 35 m → 1,524,713, **40 m → 1,495,529 — under
    # budget**. Chosen 2026-08-27; the four costed options are in
    # B-1103's register row.
    #
    # **B-1857 raised `band_tol_m` again, 40 -> 60 m, to pay for the
    # drying ground.** The Solent dries more than any other coast here
    # (Chichester, Langstone, Portsmouth, the Beaulieu and Lymington
    # flats), and its drying tint costs bytes the 4,471 of headroom at
    # 40 m could never hold. 60 m is the Saronic's value and frees
    # about 105 KB of band geometry. `drying_post_tol_m` then spends
    # fewer points on the drying edge itself; both are measured in
    # B-1857's register row.
    "solent": {
        "bounds": (50.5, 51.0, -1.75, -0.7),
        "out": "solent-bathymetry.json",
        "band_tol_m": 60.0,
        "drying_post_tol_m": 20.0,
    },
    # Corfu and the Corfu channel south past Paxos, Preveza and
    # Lefkas to Meganisi, Ithaca, Kefalonia and Zakynthos. Matches
    # IONIAN. **George's charter water** (B-1664), and the first box
    # in this table that is not British.
    #
    # **`max_depth_m` exists because of this box.** The sanity assert
    # in `fetch_grid` refused anything deeper than 150 m, which is a
    # fact about the English Channel and the Irish Sea and nothing
    # else. The Ionian trench runs past 4,000 m inside these bounds,
    # so the assert had to become a per-region number rather than a
    # hard-coded one. It is still a real guard — a parse bug that
    # doubled a sign or dropped a decimal would blow past 5,000.
    #
    # **The budget trade, measured rather than guessed** (B-1664).
    # 4.14 square degrees against the Solent's 0.53, and at the
    # shared settings the cut is 4,677,300 bytes — three times the
    # 1,500,000 budget. Every figure below is a real run of this
    # script over the 2026-09-09 EMODnet grid:
    #
    #   lines / bands / post   bytes
    #      20 /    20 /    4   4,677,300
    #      60 /   120 /    4   2,993,169
    #      90 /   180 /    4   2,237,020
    #      20 /    20 /   20   1,708,486
    #      20 /    20 /   30   1,339,032  **shipped**
    #      20 /    20 /   40   1,131,716
    #
    # **The line tolerance is the wrong knob, and the sweep proved
    # it.** Trebling it (20 -> 60) bought only a third, because
    # Chaikin re-inflates whatever Douglas-Peucker removes: the
    # tolerance sets the SHAPE and the smoothing then spends four
    # vertices drawing every corner of it. So this cut keeps the
    # shared 20 m shape tolerance every British coast uses — the
    # contours run in the same places, at the same fidelity — and
    # raises `post_tol_m`, which is only how many points are spent
    # saying so. At 30 m that is two screen pixels at z13 (15.1 m/px
    # at 38 N), the zoom a charter skipper pilots in, against the
    # sub-pixel density the 4 m default buys a boat rounding a mark.
    #
    # **Nothing is clipped.** The alternative to a density knob was
    # cutting Corfu or Zakynthos out of a chain George actually
    # sails, and a chart that stops at Ithaca is worse than a chart
    # drawn with fewer points. The bands need no knob here at all —
    # the Ionian trench keeps the 0-30 m ramp pinned to the shore.
    "ionian": {
        "bounds": (37.6, 39.9, 19.5, 21.3),
        "out": "ionian-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 5_000.0,
        "post_tol_m": 30.0,
    },
    # Piraeus and Alimos down the Saronic past Aegina, Poros and
    # Hydra to Spetses, west into the Argolic gulf for Nafplio, east
    # to Sounion. Matches SARONIC.
    #
    # 1.81 square degrees — three times the Solent, half the Ionian —
    # and **the dearer of the two Greek cuts despite being the
    # smaller box.** The Aegean shelf keeps the 0-30 m ramp in play
    # across whole gulfs where the Ionian trench pins it to the
    # shore, so there is more contour per square degree here.
    # Measured, same grid, same day:
    #
    #   lines / bands / post   bytes
    #      20 /    20 /    4   5,670,315
    #      20 /    60 /    4   5,074,455
    #      35 /    70 /    4   4,711,284
    #      45 /    90 /    4   4,217,107
    #      20 /    20 /   12   2,909,950
    #      20 /    20 /   20   2,153,001
    #      20 /    20 /   30   1,701,279   still over budget
    #      20 /    60 /   35   **shipped**
    #
    # The band knob comes back for this coast (B-1103's argument
    # unchanged: a band edge lies under a contour line drawn on the
    # same locus, so band slop is invisible where line slop is not),
    # which is what lets the density knob stop at 35 m rather than
    # the 40 m the bands-untouched column would have needed.
    "saronic": {
        "bounds": (36.85, 38.1, 22.7, 24.15),
        "out": "saronic-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 1_500.0,
        "band_tol_m": 60.0,
        "post_tol_m": 35.0,
    },
    # ── Croatia + the Trieste-Venice leg (B-1743) ───────────────────
    #
    # Owner, 2026-09-14: *"need you to load up croatia and italy maps -
    # for my other brother jonny.. i need him to tets for his motor
    # yacht bluff - its a commerical yacht sunseaker - 29 meters ...
    # depth and draft will be very important"*, *"he will be taking
    # the yacht to venice from croatia this month so i want us
    # ready"*. Bluff: Sunseeker 2014, 28.15 m LOA, 6.5 m beam,
    # **2.05 m draft**, planing, cruise 23 kn. Five non-overlapping
    # boxes, Venice round to Dubrovnik.
    #
    # **`shoal_min_ring_area_m2` is the new knob this row adds.** Every
    # region above drops a filled-band ring under `MIN_RING_AREA_M2`
    # (30,000 m2) as a speckle. For a 2.05 m draft motor yacht that
    # rule is backwards on the SHOALEST rung: a 2,000 m2 patch of 0-5 m
    # water dropped as noise is a real drying bank or a reef a deep-V
    # planing hull can find at 23 kn, and a bigger file is the honest
    # trade. The knob only lowers the floor for the `shallow` (0-5 m)
    # band; the deeper tint rungs (5-10 / 10-20 / 20-30 m) keep the
    # shared 30,000 m2 floor because they are advisory colour, not the
    # line between "sails" and "grounds". See `filled_band()`.
    #
    # **`max_depth_m: 500` is a generous sanity ceiling, not a
    # measurement.** The Adriatic's own deep pits (Jabuka, the South
    # Adriatic Pit off Otranto) sit outside all five boxes; the shelf
    # water these boxes actually cover is shallow throughout. 500 m
    # only exists to catch a parse bug the way the British 150 m and
    # the Greek 1,500/5,000 m ceilings do — it is not a claim about
    # the seabed.
    "venice-trieste": {
        "bounds": (45.10, 45.82, 12.10, 13.40),
        "out": "venice-trieste-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 500.0,
        "shoal_min_ring_area_m2": 2_000.0,
        "post_tol_m": 6.0,
    },
    "istria-kvarner": {
        "bounds": (44.45, 45.82, 13.40, 15.25),
        "out": "istria-kvarner-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 500.0,
        "shoal_min_ring_area_m2": 2_000.0,
        "band_tol_m": 60.0,
        "post_tol_m": 60.0,
    },
    "north-dalmatia": {
        "bounds": (43.65, 44.45, 14.60, 16.05),
        "out": "north-dalmatia-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 500.0,
        "shoal_min_ring_area_m2": 2_000.0,
        "band_tol_m": 90.0,
        "post_tol_m": 90.0,
    },
    # **`latMin` raised 42.85 -> 43.00 (B-1743).** The box as first
    # drafted overlapped `south-dalmatia` by 0.30 deg of latitude and
    # 0.75 deg of longitude — a real rectangle intersection, not a
    # touching edge, and `chartRegions.test.ts` forbids it. 43.00 N
    # is the shared boundary: Split (43.51), Trogir (43.51), Hvar town
    # (43.17) and Vis town (43.06) all sit north of it; Korčula town
    # (42.96) and the Pelješac peninsula's NW tip (~42.97) both sit
    # south of it. Neither named place moves boxes.
    "central-dalmatia": {
        "bounds": (43.00, 43.65, 15.75, 17.35),
        "out": "central-dalmatia-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        "max_depth_m": 500.0,
        "shoal_min_ring_area_m2": 2_000.0,
        "band_tol_m": 90.0,
        "post_tol_m": 90.0,
    },
    # **`latMax` lowered 43.15 -> 43.00 (B-1743)** — the other half of
    # the `central-dalmatia` overlap fix above. Korčula, Lastovo,
    # Mljet, Pelješac, Dubrovnik and Cavtat are all well south of 43.00.
    "south-dalmatia": {
        "bounds": (42.30, 43.00, 16.60, 18.60),
        "out": "south-dalmatia-bathymetry.json",
        # B-1857: microtidal, so no drying tint. See `DRYING_TINT`.
        "drying_tint": False,
        # **Measured, not guessed — the box reaches the South Adriatic
        # Pit.** The first run of this region tripped the sanity assert
        # at a 500 m ceiling: the real grid max is 1,151.2 m, off the
        # Pit's northern flank towards Mljet/Cavtat. 500 m was a British
        # -style guess that the Otranto pit sat outside every box here;
        # it does not. 1,500 m matches the Saronic's own ceiling and
        # clears the measured max with headroom.
        "max_depth_m": 1_500.0,
        "shoal_min_ring_area_m2": 2_000.0,
        "band_tol_m": 90.0,
        "post_tol_m": 90.0,
    },
}

DEFAULT_REGION = "abersoch"

# ── Configuration ────────────────────────────────────────────────

# Bound at import from DEFAULT_REGION; `configure()` rebinds them and
# everything derived from them when a region is named on the command
# line. Kept as module globals rather than threaded through every
# helper — this is a build script, not a library.
LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = REGIONS[DEFAULT_REGION]["bounds"]

# Contour set — matches what a sailor reads on the paper chart.
# Depths in metres below LAT. 0 is the drying line.
DEPTH_LEVELS_M = [0, 2, 5, 10, 20, 30, 50]

# Shallow tint band: 0-5 m. The water that bends the tide.
SHALLOW_BAND_M = (0.0, 5.0)

# Graduated depth tints, metres below LAT (B-1070).
#
# **The marine chart's oldest piece of grammar**: deep water is left
# bare and every step shallower takes one step more blue, so a helm
# reads the seabed as a colour before he reads a number. Until now the
# asset carried exactly ONE filled band (0-5 m above), which is a
# shallow-water WARNING, not cartography - the whole of Cardigan Bay
# outside it painted flat, and the chart read as a street map with
# dots on it.
#
# Each band is a SEPARATE `contourpy` fill between two of the contour
# levels already traced above, so:
#
#  - the bands are DISJOINT (`cg.filled(a, b)` answers a <= z <= b,
#    with holes), which means the app never stacks two tints and the
#    composited colour of every band is exactly one alpha over the
#    water. That is what makes the palette's contrast maths in
#    `src/domain/depthContours.ts` a measurement rather than an
#    estimate.
#  - a band edge IS a contour line the asset already draws, so the
#    tint and the line can never disagree.
#
# **The bands are disjoint in DEPTH; whether they are flush in the
# PLANE was measured, not assumed.** Each band's rings go through the
# smoothing pipeline independently, so a shared edge can diverge. Over
# the 2026-08-27 asset the worst mean seam is 0.53 m and the total
# overlap is 0.048% of the 546 km2 the ramp covers - a fifth of a
# pixel at z17, under a contour line of at least 1 px drawn on the
# same locus. If SIMPLIFY_TOLERANCE_M or the Chaikin settings are ever
# raised, re-measure: the seams scale with them. A region's
# `band_tol_m` (B-1103) raises exactly that seam for its bands —
# the Solent's 40 m knob means its band-band seams and its
# band-under-line offsets scale to the tens-of-metres order. That is
# the trade the knob buys the budget with, accepted because every
# band edge lies under a contour line drawn on the same locus at the
# lines' own untouched 20 m fidelity.
#
# **0-5 m is deliberately NOT re-cut here.** The legacy `shallow`
# polygon above stays byte-identical and remains the shoalest tint;
# splitting it into 0-2 / 2-5 is a fifth close-toned step on a phone
# in sunlight and needs its own contrast pass. See B-1071.
#
# 30 m and deeper is left BARE on purpose. Untinted water is the
# chart's "deep, and nothing here concerns you", and the eye needs a
# rest colour for the ramp to mean anything.
DEPTH_TINT_BANDS_M = [(5.0, 10.0), (10.0, 20.0), (20.0, 30.0)]

# Geometry budget. Tolerance is Douglas-Peucker in metres; raise it
# if the asset creeps over the bundle budget.
SIMPLIFY_TOLERANCE_M = 20.0

# Douglas-Peucker tolerance for `depth-band` polygons ONLY (B-1103).
# Defaults to SIMPLIFY_TOLERANCE_M; a region's `band_tol_m` entry
# raises it for that coast alone. **The contour lines never read
# this** — they keep their fidelity, and a band edge sits under a
# contour line of at least 1 px drawn on the same locus, which is
# what makes band slop invisible where line slop would not be.
# Rebound by `configure()` below.
BAND_SIMPLIFY_TOLERANCE_M = SIMPLIFY_TOLERANCE_M

# Douglas-Peucker tolerance for contour LINES (B-1664). Defaults to
# SIMPLIFY_TOLERANCE_M — every British coast ships at the shared 20 m
# and always has. A region's `line_tol_m` entry raises it for that sea
# alone, and raising it is a REAL loss of fidelity in a way the band
# knob is not: nothing is drawn under a contour line to hide its slop.
# It exists because the Greek boxes are four and three times the
# Solent's area and the alternative was clipping water George sails.
# Rebound by `configure()` below.
LINE_SIMPLIFY_TOLERANCE_M = SIMPLIFY_TOLERANCE_M

# Deepest plausible sounding in the box, metres. **A parse guard, not
# a cartographic limit** — the contour set stops at 50 m either way.
# 150 m is the British default (the Channel and the Irish Sea);
# a region's `max_depth_m` raises it where the sea floor really does
# go deeper, which is every Greek box. Rebound by `configure()`.
MAX_PLAUSIBLE_DEPTH_M = 150.0
COORD_DECIMALS = 5           # ~1 m — matches the grid's honesty
MIN_LINE_LENGTH_M = 250.0    # drop contour crumbs
MIN_RING_AREA_M2 = 30_000.0  # drop shallow-tint speckles (< ~170 m square)

# **The shoal floor, and it is a SEPARATE knob from `MIN_RING_AREA_M2`
# (B-1743).** Every ring below the floor is dropped as a speckle —
# fine for the three deeper tint rungs, which are advisory colour. It
# is the wrong rule for the SHOALEST rung (`kind="shallow"`, 0-5 m):
# for a draught that measures in whole metres, a 2,000-25,000 m2 patch
# of 0-5 m water is a real drying bank or reef, and dropping it as
# noise is worse than a bigger file. Defaults to `MIN_RING_AREA_M2` —
# every British and Greek cut is unaffected — and a region's
# `shoal_min_ring_area_m2` lowers it for that coast's shallow band
# alone. Rebound by `configure()` below.
SHOAL_MIN_RING_AREA_M2 = MIN_RING_AREA_M2
MAX_ASSET_BYTES = 1_500_000

# Smoothing (owner 2026-08-10: "why do the maps look so jagged").
# The raw contourpy output walks the ~115 m grid, so every vertex is
# a cell-edge crossing and the line reads as a polygon mesh, not a
# chart. The pipeline is now: Douglas-Peucker at SIMPLIFY_TOLERANCE_M
# (kills the grid staircase — same tolerance the shipped asset used,
# so the overall SHAPE is unchanged), then CHAIKIN_ITERATIONS of
# corner-cutting (each pass replaces every corner with two points at
# 1/4 and 3/4 along its edges — the polyline converges on a smooth
# quadratic B-spline), then a final light Douglas-Peucker pass at
# POST_SMOOTH_TOLERANCE_M to strip the near-collinear points Chaikin
# leaves on straight runs. Chaikin points are convex combinations of
# existing points, so a smoothed line can NEVER leave the box or
# cross a neighbouring contour it didn't already cross.
CHAIKIN_ITERATIONS = 2
POST_SMOOTH_TOLERANCE_M = 4.0

# **The vertex-density knob, and the one that actually pays** (B-1664).
#
# `SIMPLIFY_TOLERANCE_M` fixes the line's SHAPE; this one fixes how
# many points are spent drawing that shape. Chaikin quadruples the
# vertex count and the 4 m pass above puts back only what is
# collinear, so at 4 m a contour carries a point every few metres —
# right for a boat rounding a mark at 1:10k, and four times more than
# a phone can show at the zoom a charter skipper pilots at.
#
# Measured on the Saronic cut, and this is why the region knobs are
# shaped the way they are: raising the LINE tolerance from 20 m to
# 45 m moved the asset only 5,670,315 -> 4,217,107 bytes, because
# Chaikin re-inflates whatever Douglas-Peucker removed. Raising THIS
# is what brings a Greek cut under budget, and it costs density
# rather than shape. Rebound by `configure()`.
DEFAULT_POST_SMOOTH_TOLERANCE_M = POST_SMOOTH_TOLERANCE_M

# Label anchors along contour lines. The app renders the text.
LABEL_SPACING_M = 1_800.0
MIN_LABEL_LINE_M = 700.0
MAX_LABELS_PER_LINE = 4
MAX_LABELS_TOTAL = 260

# Land / no-data cells become this depth so the 0 m contour closes
# against the coast. Negative depth = above water.
NODATA_DEPTH_M = -5.0

# ── The drying-ground tint (B-1857) ───────────────────────────────
#
# Owner, 2026-10-09, over Abersoch: *"why is teh shore white?"*, *"when
# the middle is deeper and is light ?"*, *"be concistant"*.
#
# **The rule: shallower is never lighter than deeper, right up to the
# high-water line.** The tint ramp starts at the 0 m LAT line. The coast
# the app draws is the HIGH-water line. Every cell between them is
# ground that dries, and it painted as bare water - the same colour as
# 30 m and deeper. At Abersoch the beach dries a long way, so the
# shallowest ground on the chart read as the deepest.
#
# **`kind: "drying"`, a feature of its own - never folded into
# `shallow`.** The app paints it with the 0-5 m tint, in the same
# source. But `shallow` MEANS 0-5 m below LAT: `whatsHere` turns a hit
# into a charted depth of 0 m and the draught check reads that as the
# floor. A beach that dries 3 m must never answer "0-5 m". Kept apart,
# every existing feature stays byte-identical and the planner's
# answers do not move.
#
# **No `depthM`, no `fromM`/`toM`.** An app build older than B-1857
# skips a kind it does not know only when it carries neither, so a
# user cut made by this script stays safe on a phone not yet updated.
#
# **Which cells:** every EMODnet cell above LAT (`depth < 0`), whatever
# its height. Measured over the 2026-10-09 grid: 708 of Abersoch's
# above-LAT cells inside the OSM water are a coastal fringe at 5-8.4 m
# (median 5.18 m) - just over the -5 m no-data convention above - so a
# "-5 < depth < 0" cut left them white: 38% of the near-shore white
# filled against 54% (both stopping at the last data cell, as in the
# table below). Only 17 of the 2,403 above-LAT cells (0.7%) sit under
# the coastline asset's LAND polygons, none higher than 5.2 m.
#
# **Where it stops.** EMODnet's last measured cell often sits short
# of the OSM high-water line, and the grid alone cannot say by how
# much (at Abersoch, 849 no-data cells have their centres in the OSM
# water). Where the coast has a coastline cut, the cut says: a no-data
# cell within `DRYING_SHORE_REACH_CELLS` of real data whose centre is
# in the WATER is unmeasured foreshore, and is tinted. At the land
# edge the tint then runs half a cell further, to the midpoint, and
# stops. A coast with no cut (the Solent, South Devon, every `--bbox`)
# gets the half cell from its last data cell.
#
# Measured at Abersoch against the OSM coastline. "Left" is water
# within 1,500 m of the land that no tint covers; 32.30 km2 before.
# Some of it is deep (30 m and more) and bare on purpose.
#
#   stop at                        filled  left      on land  reach p95/max
#   the last data cell              54%   14.70 km2   0.08%    9 / 24 m
#   half a cell past it             63%   11.90 km2   0.71%   26 / 57 m
#   a whole cell past it            78%    7.24 km2   8.75%   84 / 135 m
#   coastline-aware (shipped)       83%    5.41 km2   2.45%   24 / 57 m
#
# The whole cell would be clipped by the land fill on the vector
# chart, but nothing clips it on the raster basemap, with SAT on, on
# the replays, or on a coast with no coastline cut. The shipped rule
# needs no clip: its reach onto land is the half cell, the same bound
# everywhere. (Drawing the tint below `coastline-land-fill` was ruled
# out: that cut draws the Venice and Grado lagoons as LAND, under 260
# km2 of today's 0-5 m tint.) No-data cells that touch WET water keep
# the far value unless the coastline puts them in the water, so the
# tint never doubles over the 0-5 m band on a steep shore; the ones it
# does put in the water keep the contour grid's own -5 m, so their 0 m
# crossing is the band's own edge.
#
# **The contour LINES are never touched.** They, the labels and every
# `shallow` / `depth-band` polygon come from the untouched grid; the
# drying tint has its own grid. The self-test proves it: strip the
# `drying` features and the file is byte-identical to a run with the
# tint switched off.
#
# **Microtidal seas are off (`drying_tint: False`).** In the Greek and
# Adriatic boxes the above-LAT cells are a coastal fringe at 0.3-0.7 m
# (median), on a sea with tens of centimetres of tide - not drying
# ground. At each box's own line settings the tint would add 98-407 KB
# per box, over the 1.5 MB budget in five of seven, for a change a
# render of the Dubrovnik coast barely shows. Venice is the closest
# call (about 1 m of tide, real lagoon flats); it is off with the rest
# and left to the owner (B-1857). A user's `--bbox` cut
# keeps the default (on): no budget, and a Brittany box has real tide.
DRYING_KIND = "drying"
DRYING_TINT = True

# Post-smoothing density for the drying tint only. Defaults to the
# region's own `POST_SMOOTH_TOLERANCE_M`; a region's
# `drying_post_tol_m` raises it (the Solent, for its budget). The
# shape tolerance is the LINE one, as for `shallow`: the drying edge
# meets the 0-5 m band along the 0 m line and must stay as tight.
DRYING_POST_SMOOTH_TOLERANCE_M = 4.0

# Far below every band: no-data cells the drying tint must not reach.
DRYING_FAR_BELOW_M = 1e6

# How far, in grid cells, the tint may run through no-data cells that
# the coastline cut puts in the WATER (B-1857). One or two cells is
# the rim between EMODnet's last measured cell and the OSM high-water
# line; unbounded, it would flood estuaries and harbours EMODnet never
# surveyed (Pwllheli, the Glaslyn) and call them drying ground.
DRYING_SHORE_REACH_CELLS = 2

# The coastline cut for the named region, or None (`--bbox`, and any
# coast without one). Bound by `configure()` / `configure_bbox()`.
COASTLINE: Path | None = None

HERE = Path(__file__).resolve().parent
M_PER_DEG_LAT = 111_132.0

# Bound by `configure()` below, which runs at import for the default
# region and again when a region is named on the command line.
ERDDAP_URL = ""
CACHE = HERE / "cache" / "unconfigured.csv"
OUT = HERE.parent.parent / "assets" / "depth" / "unconfigured.json"
REGION = DEFAULT_REGION
M_PER_DEG_LON = 111_320.0


def configure(region: str) -> None:
    """Point every box-derived global at `region`.

    **One place computes the box, the URL, the cache path and the
    output path.** They were four independent module constants before
    three regions existed, which is how a re-cut of one coast quietly
    overwrites another coast's asset.
    """
    global LAT_MIN, LAT_MAX, LON_MIN, LON_MAX
    global ERDDAP_URL, CACHE, OUT, REGION, M_PER_DEG_LON
    global BAND_SIMPLIFY_TOLERANCE_M, LINE_SIMPLIFY_TOLERANCE_M
    global MAX_PLAUSIBLE_DEPTH_M, POST_SMOOTH_TOLERANCE_M
    global SHOAL_MIN_RING_AREA_M2
    global DRYING_TINT, DRYING_POST_SMOOTH_TOLERANCE_M, COASTLINE
    if region not in REGIONS:
        known = ", ".join(sorted(REGIONS))
        raise SystemExit(f"unknown region {region!r}. Known regions: {known}")
    REGION = region
    BAND_SIMPLIFY_TOLERANCE_M = REGIONS[region].get(
        "band_tol_m", SIMPLIFY_TOLERANCE_M
    )
    LINE_SIMPLIFY_TOLERANCE_M = REGIONS[region].get(
        "line_tol_m", SIMPLIFY_TOLERANCE_M
    )
    MAX_PLAUSIBLE_DEPTH_M = REGIONS[region].get("max_depth_m", 150.0)
    POST_SMOOTH_TOLERANCE_M = REGIONS[region].get(
        "post_tol_m", DEFAULT_POST_SMOOTH_TOLERANCE_M
    )
    SHOAL_MIN_RING_AREA_M2 = REGIONS[region].get(
        "shoal_min_ring_area_m2", MIN_RING_AREA_M2
    )
    # B-1857 - the drying tint, on unless a region says otherwise.
    DRYING_TINT = REGIONS[region].get("drying_tint", True)
    DRYING_POST_SMOOTH_TOLERANCE_M = REGIONS[region].get(
        "drying_post_tol_m", POST_SMOOTH_TOLERANCE_M
    )
    LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = REGIONS[region]["bounds"]
    ERDDAP_URL = (
        "https://erddap.emodnet.eu/erddap/griddap/bathymetry_dtm_2024.csv"
        f"?elevation%5B({LAT_MIN}):({LAT_MAX})%5D%5B({LON_MIN}):({LON_MAX})%5D"
    )
    CACHE = HERE / "cache" / (
        f"dtm2024_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}.csv"
    )
    OUT = HERE.parent.parent / "assets" / "depth" / REGIONS[region]["out"]
    COASTLINE = (
        HERE.parent.parent / "assets" / "coastline" / f"{region}-coastline.json"
    )
    # Metres per degree at the box centre — good enough for
    # simplification tolerances; nobody navigates off a
    # Douglas-Peucker epsilon.
    lat_mid = (LAT_MIN + LAT_MAX) / 2.0
    M_PER_DEG_LON = 111_320.0 * math.cos(math.radians(lat_mid))


configure(DEFAULT_REGION)


# ── Any box: a chart area the helm cut himself (B-1845) ──────────
#
# **The bundled coasts are a table; a user cut is four numbers.** The
# relay snaps the helm's box to a 0.1 deg grid and the cut job in
# `jamesrolex/openracer-charts` runs this script with `--bbox`. The
# contour settings are the SHARED defaults every bundled coast starts
# from — 20 m shape, 4 m density, the 30,000 m2 ring floor — so a
# Cascais cut reads exactly like Abersoch. No per-box knobs: nobody is
# here to sweep them.

# **What EMODnet DTM 2024 actually covers**, from the dataset's own
# `actual_range` (ERDDAP `.das`, read 2026-10-07). European seas, the
# Med, the Black Sea and the North Atlantic shelf. A box is clipped to
# it first; a box wholly outside it is answered without a request.
EMODNET_LAT_RANGE = (15.0, 90.0)
EMODNET_LON_RANGE = (-36.0, 43.0)

# **The app's own size cap, mirrored** (`USER_DEPTH_MAX_BYTES` in
# `src/domain/userChartRegion.ts`). A bundled asset fails over
# 1.5 MB because it ships in the binary; a user depth file is a
# one-off download, so it may be larger. Over THIS cap the phone would
# refuse it anyway, so the script writes nothing instead.
USER_CUT_MAX_BYTES = 20 * 1024 * 1024

# **A parse guard, not a claim.** An arbitrary box can reach the
# abyssal plain off Portugal, so the British 150 m ceiling cannot
# apply. Deeper than the deepest trench is a parse bug.
BBOX_MAX_PLAUSIBLE_DEPTH_M = 11_000.0

# Set by `configure_bbox()`. The named-region mode never reads it.
BBOX_MODE = False


class NoDepthData(Exception):
    """EMODnet has nothing for this box. **Not an error** — the job
    writes no file, says why in one line, and exits 0."""


def parse_bbox(text: str) -> tuple[float, float, float, float]:
    """`W,S,E,N` in WGS84 degrees -> (west, south, east, north).

    Refuses a box the cutter would refuse: four finite numbers, west
    below east (no antimeridian wrap), south below north, on Earth."""
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 4:
        raise ValueError(f"--bbox wants W,S,E,N - got {text!r}")
    try:
        west, south, east, north = (float(p) for p in parts)
    except ValueError as err:
        raise ValueError(f"--bbox wants four numbers - got {text!r}") from err
    if not all(math.isfinite(v) for v in (west, south, east, north)):
        raise ValueError(f"--bbox wants four finite numbers - got {text!r}")
    if not (-180.0 <= west < east <= 180.0):
        raise ValueError(f"--bbox west must be below east, inside +-180 - got {text!r}")
    if not (-90.0 <= south < north <= 90.0):
        raise ValueError(f"--bbox south must be below north, inside +-90 - got {text!r}")
    return west, south, east, north


def clip_to_emodnet(
    west: float, south: float, east: float, north: float
) -> tuple[float, float, float, float] | None:
    """The part of the box EMODnet covers, or None when it covers none.

    **Clipped, not refused.** ERDDAP answers 404 to a box that pokes
    one cell past its axis, so a box straddling the edge would get no
    depth at all for the half it does cover."""
    s = max(south, EMODNET_LAT_RANGE[0])
    n = min(north, EMODNET_LAT_RANGE[1])
    w = max(west, EMODNET_LON_RANGE[0])
    e = min(east, EMODNET_LON_RANGE[1])
    if not (w < e and s < n):
        return None
    return w, s, e, n


def configure_bbox(
    west: float, south: float, east: float, north: float, out: Path
) -> None:
    """Point every box-derived global at an arbitrary box.

    **The same globals `configure()` binds**, so every helper below
    runs unchanged. Tolerances go back to the shared defaults — the
    bundled coasts' own starting point — whatever ran before."""
    global LAT_MIN, LAT_MAX, LON_MIN, LON_MAX
    global ERDDAP_URL, CACHE, OUT, REGION, M_PER_DEG_LON
    global BAND_SIMPLIFY_TOLERANCE_M, LINE_SIMPLIFY_TOLERANCE_M
    global MAX_PLAUSIBLE_DEPTH_M, POST_SMOOTH_TOLERANCE_M
    global SHOAL_MIN_RING_AREA_M2, BBOX_MODE
    global DRYING_TINT, DRYING_POST_SMOOTH_TOLERANCE_M, COASTLINE
    BBOX_MODE = True
    COASTLINE = None
    REGION = f"bbox {west},{south},{east},{north}"
    BAND_SIMPLIFY_TOLERANCE_M = SIMPLIFY_TOLERANCE_M
    LINE_SIMPLIFY_TOLERANCE_M = SIMPLIFY_TOLERANCE_M
    POST_SMOOTH_TOLERANCE_M = DEFAULT_POST_SMOOTH_TOLERANCE_M
    SHOAL_MIN_RING_AREA_M2 = MIN_RING_AREA_M2
    # B-1857 - a user's box always gets the drying tint: no budget,
    # and nobody here knows its tide.
    DRYING_TINT = True
    DRYING_POST_SMOOTH_TOLERANCE_M = DEFAULT_POST_SMOOTH_TOLERANCE_M
    MAX_PLAUSIBLE_DEPTH_M = BBOX_MAX_PLAUSIBLE_DEPTH_M
    LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = south, north, west, east
    ERDDAP_URL = (
        "https://erddap.emodnet.eu/erddap/griddap/bathymetry_dtm_2024.csv"
        f"?elevation%5B({LAT_MIN}):({LAT_MAX})%5D%5B({LON_MIN}):({LON_MAX})%5D"
    )
    CACHE = HERE / "cache" / (
        f"dtm2024_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}.csv"
    )
    OUT = out
    lat_mid = (LAT_MIN + LAT_MAX) / 2.0
    M_PER_DEG_LON = 111_320.0 * math.cos(math.radians(lat_mid))


# ── Download + parse ─────────────────────────────────────────────

def fetch_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (lats, lons, depth) — depth positive DOWN in metres."""
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {ERDDAP_URL}")
        try:
            with urllib.request.urlopen(ERDDAP_URL, timeout=180) as r:
                CACHE.write_bytes(r.read())
        except urllib.error.HTTPError as err:
            # **ERDDAP says 404 for "no matching results"** — a box
            # outside the grid. That is an answer, not a failure, and
            # it must never reach the curl fallback below (whose `-f`
            # would turn it into a crash).
            if err.code in (400, 404):
                raise NoDepthData(
                    f"EMODnet has no grid for this box (HTTP {err.code})"
                ) from err
            raise
        except urllib.error.URLError:
            # macOS framework Python often ships without CA certs.
            # curl uses the system trust store, so fall back to it.
            print("urllib SSL failed - falling back to curl")
            subprocess.run(
                ["curl", "-sSf", "--max-time", "180", "-o", str(CACHE),
                 ERDDAP_URL],
                check=True,
            )
        print(f"cached {CACHE.stat().st_size:,} bytes -> {CACHE}")
    else:
        print(f"using cache {CACHE} ({CACHE.stat().st_size:,} bytes)")

    # **Two passes, not one dictionary — B-1664.** The single-pass
    # form kept every row's lat, lon and elevation as Python objects
    # before the grid existed: about 250 bytes a cell. That is fine
    # for the Solent's 484,000 cells and it is 1 GB for the Ionian's
    # 3.8 million. Pass one learns the axes (a few thousand floats);
    # pass two fills the array straight from the file. **The grid is
    # identical either way** — same sorted unique axes, same
    # last-value-wins per cell — so the British cuts are untouched.
    def rows():
        with CACHE.open() as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header[:3] == ["latitude", "longitude", "elevation"], header
            next(reader)  # units row
            yield from reader

    lat_set: set[float] = set()
    lon_set: set[float] = set()
    for lat_s, lon_s, _elev_s in rows():
        lat_set.add(float(lat_s))
        lon_set.add(float(lon_s))

    ulats = np.array(sorted(lat_set))
    ulons = np.array(sorted(lon_set))
    depth = np.full((len(ulats), len(ulons)), NODATA_DEPTH_M)
    lat_idx = {v: i for i, v in enumerate(ulats)}
    lon_idx = {v: i for i, v in enumerate(ulons)}
    n_data = 0
    for lat_s, lon_s, elev_s in rows():
        if elev_s in ("", "NaN"):
            continue
        # Depth positive down; elevation positive up from LAT.
        depth[lat_idx[float(lat_s)], lon_idx[float(lon_s)]] = -float(elev_s)
        n_data += 1

    n_total = depth.size
    # **A user box can be all land, all blank, or one row thin
    # (B-1845).** None of those is a parse bug, and none has anything
    # to contour: say so, and let the caller write nothing. Named
    # coasts never get here — every one was cut with water in it.
    if BBOX_MODE:
        if len(ulats) < 2 or len(ulons) < 2:
            raise NoDepthData("the EMODnet grid for this box is too small to contour")
        if n_data == 0 or float(np.max(depth)) <= 0.0:
            raise NoDepthData("EMODnet has no sea depths inside this box")
    print(
        f"grid {len(ulats)} x {len(ulons)} = {n_total:,} cells, "
        f"{n_data:,} with data ({100 * n_data / n_total:.0f}%), "
        f"deepest {float(np.max(depth)):.1f} m"
    )
    if BBOX_MODE:
        # The ceiling only: a shallow box (the Wadden Sea, a lagoon)
        # is real water under the 20 m floor the named coasts assert.
        assert float(np.max(depth)) < MAX_PLAUSIBLE_DEPTH_M, (
            f"implausible max depth for {REGION}: {float(np.max(depth)):.1f} m"
        )
        return ulats, ulons, depth
    # Sanity: a wild max means a parse bug, not a discovery.
    # **The ceiling is per region** (B-1664): 150 m describes the
    # Channel and the Irish Sea and nothing else, and the Ionian
    # trench runs past 4,000 m inside its own box. Still a real guard
    # — a lost decimal or a flipped sign blows past any of them.
    assert 20.0 < float(np.max(depth)) < MAX_PLAUSIBLE_DEPTH_M, (
        f"implausible max depth for {REGION}: "
        f"{float(np.max(depth)):.1f} m against a "
        f"{MAX_PLAUSIBLE_DEPTH_M:.0f} m ceiling"
    )
    return ulats, ulons, depth


# ── Geometry helpers ─────────────────────────────────────────────

def to_metres(pts: np.ndarray) -> np.ndarray:
    """Lon/lat pairs -> local metre coordinates for tolerance maths."""
    out = np.empty_like(pts)
    out[:, 0] = (pts[:, 0] - LON_MIN) * M_PER_DEG_LON
    out[:, 1] = (pts[:, 1] - LAT_MIN) * M_PER_DEG_LAT
    return out


def simplify(pts: np.ndarray, tol_m: float) -> np.ndarray:
    """Iterative Douglas-Peucker on lon/lat points, tolerance in metres."""
    if len(pts) < 3:
        return pts
    m = to_metres(pts)
    keep = np.zeros(len(pts), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        seg = m[b] - m[a]
        seg_len = math.hypot(seg[0], seg[1])
        rel = m[a + 1 : b] - m[a]
        if seg_len < 1e-9:
            dist = np.hypot(rel[:, 0], rel[:, 1])
        else:
            dist = np.abs(rel[:, 0] * seg[1] - rel[:, 1] * seg[0]) / seg_len
        i = int(np.argmax(dist))
        if dist[i] > tol_m:
            split = a + 1 + i
            keep[split] = True
            stack.append((a, split))
            stack.append((split, b))
    return pts[keep]


def chaikin(pts: np.ndarray, iterations: int, closed: bool) -> np.ndarray:
    """Chaikin corner-cutting. Open lines keep their endpoints;
    closed rings cut every corner including the seam. Input and
    output are lon/lat point arrays (no re-closing point — callers
    close rings themselves)."""
    if len(pts) < 3:
        return pts
    out = pts
    for _ in range(iterations):
        if closed:
            # Treat the ring as cyclic. Drop a duplicated closing
            # point if present so the seam is cut like any corner.
            ring = out[:-1] if np.allclose(out[0], out[-1]) else out
            if len(ring) < 3:
                return out
            nxt = np.roll(ring, -1, axis=0)
            q = 0.75 * ring + 0.25 * nxt
            r = 0.25 * ring + 0.75 * nxt
            out = np.empty((len(ring) * 2, 2))
            out[0::2] = q
            out[1::2] = r
        else:
            seg_a = out[:-1]
            seg_b = out[1:]
            q = 0.75 * seg_a + 0.25 * seg_b
            r = 0.25 * seg_a + 0.75 * seg_b
            mid = np.empty((len(seg_a) * 2, 2))
            mid[0::2] = q
            mid[1::2] = r
            out = np.vstack([out[:1], mid, out[-1:]])
    return out


def smooth_line(pts: np.ndarray) -> np.ndarray:
    """Full smoothing pipeline for an open contour line: staircase
    kill (DP at the region's line tolerance) -> Chaikin -> light DP.

    `LINE_SIMPLIFY_TOLERANCE_M` is the shared 20 m for every British
    coast and always has been; only the Greek boxes raise it (B-1664),
    and they say so in the REGIONS table."""
    pts = simplify(pts, LINE_SIMPLIFY_TOLERANCE_M)
    pts = chaikin(pts, CHAIKIN_ITERATIONS, closed=False)
    return simplify(pts, POST_SMOOTH_TOLERANCE_M)


def smooth_ring(
    pts: np.ndarray,
    tol_m: float | None = None,
    post_m: float | None = None,
) -> np.ndarray:
    """Same pipeline for a closed polygon ring (not re-closed).

    `tol_m` overrides the staircase-kill tolerance — the band-only
    knob (B-1103). `post_m` overrides the post-smooth density — the
    drying-only knob (B-1857). Chaikin stays shared."""
    pts = simplify(
        pts, LINE_SIMPLIFY_TOLERANCE_M if tol_m is None else tol_m
    )
    pts = chaikin(pts, CHAIKIN_ITERATIONS, closed=True)
    return simplify(
        pts, POST_SMOOTH_TOLERANCE_M if post_m is None else post_m
    )


def line_length_m(pts: np.ndarray) -> float:
    m = to_metres(pts)
    return float(np.sum(np.hypot(np.diff(m[:, 0]), np.diff(m[:, 1]))))


def ring_area_m2(pts: np.ndarray) -> float:
    m = to_metres(pts)
    x, y = m[:, 0], m[:, 1]
    return abs(float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])) / 2.0)


def round_coords(pts: np.ndarray) -> list[list[float]]:
    out: list[list[float]] = []
    prev: list[float] | None = None
    for lon, lat in pts:
        p = [round(float(lon), COORD_DECIMALS), round(float(lat), COORD_DECIMALS)]
        if p != prev:
            out.append(p)
        prev = p
    return out


def label_points(pts: np.ndarray) -> list[list[float]]:
    """Anchor points along a contour line, LABEL_SPACING_M apart."""
    m = to_metres(pts)
    seg = np.hypot(np.diff(m[:, 0]), np.diff(m[:, 1]))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(cum[-1])
    if total < MIN_LABEL_LINE_M:
        return []
    anchors: list[list[float]] = []
    target = LABEL_SPACING_M / 2.0
    while target < total and len(anchors) < MAX_LABELS_PER_LINE:
        i = int(np.searchsorted(cum, target))
        i = min(max(i, 1), len(pts) - 1)
        f = (target - cum[i - 1]) / max(float(cum[i] - cum[i - 1]), 1e-9)
        lon = float(pts[i - 1, 0] + f * (pts[i, 0] - pts[i - 1, 0]))
        lat = float(pts[i - 1, 1] + f * (pts[i, 1] - pts[i - 1, 1]))
        anchors.append(
            [round(lon, COORD_DECIMALS), round(lat, COORD_DECIMALS)]
        )
        target += LABEL_SPACING_M
    return anchors


def filled_band(cg, lo: float, hi: float, kind: str) -> list[dict]:
    """Polygon features for the water between `lo` and `hi` metres.

    `FillType.OuterOffset`: one entry per outer boundary; the offsets
    mark where the outer ring ends and each hole ring begins.

    **Lifted out of the old inline shallow-tint loop unchanged**, so
    `kind="shallow"` still produces byte-identical geometry to the
    asset that shipped before the graduated bands existed. The
    regeneration check in `docs/plans/marine-cartography.md` asserts
    exactly that, hash and all.
    """
    out: list[dict] = []
    # Bands take the band-only tolerance (B-1103); the shallow tint
    # keeps the shared one, so `kind="shallow"` stays byte-identical
    # to the pre-band asset everywhere.
    # **The shoal tint takes the LINE tolerance, not the band one.**
    # It is the rung a helm reads as a warning, and both its edges —
    # the 0 m drying line and the 5 m contour — are drawn at the line
    # tolerance on the same locus, so matching them keeps the seam as
    # tight as the region allows. On every British coast
    # `LINE_SIMPLIFY_TOLERANCE_M` is the shared 20 m, so `shallow`
    # stays byte-identical to the pre-band asset exactly as B-1103
    # promised; only the Greek boxes move it (B-1664).
    tol_m = (
        BAND_SIMPLIFY_TOLERANCE_M if kind == "depth-band"
        else LINE_SIMPLIFY_TOLERANCE_M
    )
    # **The shoal band gets its own, lower floor (B-1743).** Every
    # other rung drops a ring under `MIN_RING_AREA_M2` as a speckle;
    # `kind="shallow"` (0-5 m) uses the region's `SHOAL_MIN_RING_AREA_M2`
    # instead, which defaults to the same 30,000 m2 everywhere and only
    # moves where a region says so. See the constant's own comment.
    area_floor = (
        SHOAL_MIN_RING_AREA_M2 if kind == "shallow" else MIN_RING_AREA_M2
    )
    # **B-1857 — the drying tint is the shoal rung carried up the
    # beach.** It takes the shoal band's floor and the line tolerance
    # (set above, since it is not a `depth-band`), plus its own
    # post-smooth density so a region can pay for it.
    post_m: float | None = None
    if kind == DRYING_KIND:
        area_floor = SHOAL_MIN_RING_AREA_M2
        post_m = DRYING_POST_SMOOTH_TOLERANCE_M
    points_list, offsets_list = cg.filled(lo, hi)
    for pts_arr, offsets in zip(points_list, offsets_list):
        rings: list[list[list[float]]] = []
        for a, b in zip(offsets[:-1], offsets[1:]):
            ring = np.asarray(pts_arr[a:b])
            if len(ring) < 4:
                continue
            ring = smooth_ring(ring, tol_m, post_m)
            if len(ring) < 4:
                continue
            # First ring is the outer boundary; holes below the area
            # floor just disappear into the tint.
            if len(rings) == 0 and ring_area_m2(ring) < area_floor:
                break
            if len(rings) > 0 and ring_area_m2(ring) < area_floor:
                continue
            coords = round_coords(ring)
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            if len(coords) >= 4:
                rings.append(coords)
        if rings:
            # A drying feature carries its kind and nothing else: see
            # `DRYING_KIND` for why an old app build needs exactly that.
            props: dict = (
                {"kind": kind}
                if kind == DRYING_KIND
                else {"kind": kind, "fromM": lo, "toM": hi}
            )
            out.append(
                {
                    "type": "Feature",
                    "properties": props,
                    "geometry": {"type": "Polygon", "coordinates": rings},
                }
            )
    return out


def _neighbour_sum(
    mask: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per cell: the sum of `values` over its 8 neighbours where `mask`
    holds, and how many such neighbours there are."""
    total = np.zeros(values.shape)
    count = np.zeros(values.shape)
    h, w = values.shape
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            # Target cell [i, j] reads neighbour [i + di, j + dj].
            ti = slice(max(-di, 0), h - max(di, 0))
            si = slice(max(di, 0), h - max(-di, 0))
            tj = slice(max(-dj, 0), w - max(dj, 0))
            sj = slice(max(dj, 0), w - max(-dj, 0))
            hit = mask[si, sj]
            total[ti, tj] += np.where(hit, values[si, sj], 0.0)
            count[ti, tj] += hit
    return total, count


def _inside_rings(
    xs: np.ndarray, ys: np.ndarray, rings: list[np.ndarray]
) -> np.ndarray:
    """Even-odd point-in-polygon for many points against many rings
    (outer rings and holes alike), in numpy. No new dependency."""
    inside = np.zeros(len(xs), dtype=bool)
    for ring in rings:
        x0, y0 = ring[:-1, 0], ring[:-1, 1]
        x1, y1 = ring[1:, 0], ring[1:, 1]
        # Edges that cannot reach any point are skipped up front.
        keep = (np.maximum(y0, y1) >= ys.min()) & (np.minimum(y0, y1) <= ys.max())
        x0, y0, x1, y1 = x0[keep], y0[keep], x1[keep], y1[keep]
        if len(x0) == 0:
            continue
        step = max(1, 4_000_000 // len(x0))
        for a in range(0, len(xs), step):
            px = xs[a : a + step, None]
            py = ys[a : a + step, None]
            spans = (y0 > py) != (y1 > py)
            with np.errstate(divide="ignore", invalid="ignore"):
                cross = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
            hits = spans & (px < cross)
            inside[a : a + step] ^= (hits.sum(axis=1) % 2).astype(bool)
    return inside


def load_water_rings(path: Path | None) -> list[np.ndarray] | None:
    """The coastline asset's WATER rings (B-1857), or None.

    `assets/coastline/<region>-coastline.json` is the OSM high-water
    line, cut by `tools/coastline-extract/`. Only named coasts have
    one; a user's `--bbox` cut never does."""
    if path is None or not path.exists():
        return None
    body = json.loads(path.read_text())
    rings: list[np.ndarray] = []
    for feat in body.get("features", []):
        if feat.get("properties", {}).get("kind") != "water":
            continue
        geom = feat["geometry"]
        polys = (
            [geom["coordinates"]]
            if geom["type"] == "Polygon"
            else geom["coordinates"] if geom["type"] == "MultiPolygon" else []
        )
        for poly in polys:
            rings.extend(np.asarray(r, dtype=float) for r in poly)
    return rings or None


def drying_grid(
    depth: np.ndarray,
    lats: np.ndarray | None = None,
    lons: np.ndarray | None = None,
    water: list[np.ndarray] | None = None,
) -> tuple[np.ndarray, float] | None:
    """The grid the drying tint is traced from, and its floor level.

    **A copy — the contour grid is never touched** (B-1857). Data
    cells keep their real values, so the tint's seaward edge is the
    same 0 m crossing the drying line and the 0-5 m band are drawn on.
    The tint is `floor <= z <= 0`, with `floor` strictly below every
    data value, so every above-LAT cell is in it whatever its height.

    **Where the coast has a coastline cut (`water`),** a no-data cell
    within `DRYING_SHORE_REACH_CELLS` of real data whose centre lies in
    the OSM WATER is foreshore EMODnet never measured: it keeps the
    contour grid's own `NODATA_DEPTH_M` and so is tinted. Its 0 m
    crossing against wet water is the very one the 0-5 m band is cut
    on, so the two meet without a seam of their own.

    Then, at the land edge: a no-data cell next to tinted ground (and
    to no wet cell) takes `2 * floor - mean(neighbours)`, which puts
    the boundary half way to it. Every other no-data cell sits far
    below the floor, so the tint never reaches it. See `DRYING_KIND`.

    None when there is nothing to tint at all."""
    nodata = depth == NODATA_DEPTH_M
    data = ~nodata
    if not data.any():
        return None
    wet = data & (depth > 0)
    tinted = data & (depth < 0)
    grid = depth.copy()
    # Strictly below every data value AND the shore cells' -5 m.
    floor = min(float(np.min(depth[data])), NODATA_DEPTH_M) - 1.0
    if water is not None and lats is not None and lons is not None:
        reach = data.copy()
        for _ in range(DRYING_SHORE_REACH_CELLS):
            reach |= _neighbour_sum(reach, depth)[1] > 0
        near = nodata & reach
        rows, cols = np.nonzero(near)
        shore = np.zeros_like(nodata)
        if len(rows):
            shore[rows, cols] = _inside_rings(
                np.asarray(lons, dtype=float)[cols],
                np.asarray(lats, dtype=float)[rows],
                water,
            )
        tinted |= shore
        # `grid` already holds NODATA_DEPTH_M there: kept, and tinted.
        nodata = nodata & ~shore
    if not tinted.any():
        return None
    grid[nodata] = floor - DRYING_FAR_BELOW_M
    t_sum, t_n = _neighbour_sum(tinted, grid)
    _wet_sum, wet_n = _neighbour_sum(wet, depth)
    edge = nodata & (t_n > 0) & (wet_n == 0)
    grid[edge] = 2.0 * floor - t_sum[edge] / t_n[edge]
    return grid, floor


def drying_band(
    lats: np.ndarray,
    lons: np.ndarray,
    depth: np.ndarray,
    water: list[np.ndarray] | None = None,
) -> list[dict]:
    """The drying-ground tint polygons (B-1857), or none."""
    traced = drying_grid(depth, lats, lons, water)
    if traced is None:
        return []
    grid, floor = traced
    cg = contour_generator(
        x=lons,
        y=lats,
        z=grid,
        line_type=LineType.Separate,
        fill_type=FillType.OuterOffset,
    )
    return filled_band(cg, floor, 0.0, DRYING_KIND)


# ── Extraction ───────────────────────────────────────────────────

def extract(
    lats: np.ndarray, lons: np.ndarray, depth: np.ndarray
) -> tuple[dict, dict[str, object]]:
    """The asset for one grid: (FeatureCollection, stats for the log).

    **One pipeline for every caller** — the named coasts, a user's box
    and the self-test all run exactly this, so a bbox file cannot drift
    from the bundled shape the app parses."""
    cg = contour_generator(
        x=lons,
        y=lats,
        z=depth,
        line_type=LineType.Separate,
        fill_type=FillType.OuterOffset,
    )

    features: list[dict] = []
    labels: list[dict] = []
    stats: dict[str, int] = {}

    # Contour lines, one level at a time.
    for level in DEPTH_LEVELS_M:
        n = 0
        for line in cg.lines(float(level)):
            pts = smooth_line(np.asarray(line))
            if len(pts) < 2 or line_length_m(pts) < MIN_LINE_LENGTH_M:
                continue
            features.append(
                {
                    "type": "Feature",
                    "properties": {"kind": "contour", "depthM": level},
                    "geometry": {
                        "type": "LineString",
                        "coordinates": round_coords(pts),
                    },
                }
            )
            n += 1
            for anchor in label_points(pts):
                labels.append(
                    {
                        "type": "Feature",
                        "properties": {"kind": "label", "depthM": level},
                        "geometry": {"type": "Point", "coordinates": anchor},
                    }
                )
        stats[f"{level} m"] = n

    # Shallow-water tint: filled band between 0 and 5 m.
    shallow = filled_band(cg, SHALLOW_BAND_M[0], SHALLOW_BAND_M[1], "shallow")
    features.extend(shallow)
    n_poly = len(shallow)

    # Graduated tints for the water outside it (B-1070). Deepest band
    # FIRST, so the feature order runs deep -> shoal and a reader of
    # the raw asset sees the ramp in the order the chart paints it.
    band_stats: dict[str, int] = {}
    for lo, hi in sorted(DEPTH_TINT_BANDS_M, reverse=True):
        band = filled_band(cg, lo, hi, "depth-band")
        features.extend(band)
        band_stats[f"{lo:g}-{hi:g} m"] = len(band)

    # The drying tint (B-1857), from its OWN grid, after every other
    # tint so the features above keep their bytes and their order.
    drying = (
        drying_band(lats, lons, depth, load_water_rings(COASTLINE))
        if DRYING_TINT
        else []
    )
    features.extend(drying)

    labels = labels[:MAX_LABELS_TOTAL]
    fc = {
        "type": "FeatureCollection",
        "features": features + labels,
    }
    return fc, {
        "lines": stats,
        "shallow": n_poly,
        "labels": len(labels),
        "bands": band_stats,
        "drying": len(drying) if DRYING_TINT else None,
    }


def encode(fc: dict) -> str:
    """The file body. Compact separators and a trailing newline: the
    bytes every bundled asset has always been written with."""
    return json.dumps(fc, separators=(",", ":")) + "\n"


def print_stats(stats: dict[str, object]) -> None:
    print(f"contour lines per level: {stats['lines']}")
    print(f"shallow polygons: {stats['shallow']}, labels: {stats['labels']}")
    print(f"depth-band polygons: {stats['bands']}")
    drying = stats["drying"]
    print(
        "drying polygons: off for this region"
        if drying is None
        else f"drying polygons: {drying}"
    )
    print(
        f"tolerances: lines {LINE_SIMPLIFY_TOLERANCE_M:g} m, "
        f"bands {BAND_SIMPLIFY_TOLERANCE_M:g} m"
    )


def main() -> None:
    """One named coast -> its bundled asset. Behaviour unchanged."""
    lats, lons, depth = fetch_grid()
    fc, stats = extract(lats, lons, depth)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(encode(fc))

    size = OUT.stat().st_size
    print_stats(stats)
    print(f"wrote {OUT} ({size:,} bytes)")
    if size > MAX_ASSET_BYTES:
        print(
            f"ASSET OVER BUDGET ({size:,} > {MAX_ASSET_BYTES:,}). "
            "Raise SIMPLIFY_TOLERANCE_M / MIN_LINE_LENGTH_M and re-run.",
            file=sys.stderr,
        )
        sys.exit(1)


def run_bbox(
    west: float,
    south: float,
    east: float,
    north: float,
    out: Path,
    cache: Path | None = None,
) -> bool:
    """One user box -> `out`, or nothing. True when a file was written.

    **Never fails for want of depth** (B-1845). Outside EMODnet, an
    all-land box, a box with nothing shallower than 50 m, or a file
    over the app's cap: each prints one plain line and returns False,
    and the caller exits 0. Only a real fault (a parse bug, a dead
    network) raises.

    `cache` points the grid read at a local CSV. The self-test uses
    it; the cut job does not."""
    clipped = clip_to_emodnet(west, south, east, north)
    if clipped is None:
        print(
            "No depth data for this area: EMODnet covers European waters "
            f"only (lat {EMODNET_LAT_RANGE[0]:g} to {EMODNET_LAT_RANGE[1]:g}, "
            f"lon {EMODNET_LON_RANGE[0]:g} to {EMODNET_LON_RANGE[1]:g}). "
            "Nothing written."
        )
        return False
    global CACHE
    configure_bbox(*clipped, out)
    if cache is not None:
        CACHE = cache
    print(f"-- bbox -- box {LAT_MIN},{LAT_MAX} {LON_MIN},{LON_MAX}")
    try:
        lats, lons, depth = fetch_grid()
    except NoDepthData as why:
        print(f"No depth data for this area: {why}. Nothing written.")
        return False
    fc, stats = extract(lats, lons, depth)
    print_stats(stats)
    if not fc["features"]:
        print(
            "No depth data for this area: the box holds no water "
            "shallower than 50 m to draw. Nothing written."
        )
        return False
    body = encode(fc)
    size = len(body.encode("utf-8"))
    if size > USER_CUT_MAX_BYTES:
        print(
            f"No depth data for this area: the file would be {size:,} bytes, "
            f"over the app's {USER_CUT_MAX_BYTES:,} byte cap. Nothing written."
        )
        return False
    out.parent.mkdir(parents=True, exist_ok=True)
    # **Write beside, then rename.** A job killed mid-write must never
    # leave a half file under the name the release upload looks for.
    tmp = out.with_name(out.name + ".part")
    tmp.write_text(body)
    tmp.replace(out)
    print(f"wrote {out} ({size:,} bytes)")
    return True


# ── Self-test: a tiny synthetic grid, no network ─────────────────

def _synthetic_csv(path: Path, lats: list[float], lons: list[float], fn) -> None:
    """An ERDDAP-shaped CSV: header row, units row, then one row per
    cell. `fn(lat, lon)` returns elevation (negative = sea), or None
    for a blank cell."""
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["latitude", "longitude", "elevation"])
        w.writerow(["degrees_north", "degrees_east", "m"])
        for lat in lats:
            for lon in lons:
                elev = fn(lat, lon)
                w.writerow(
                    [f"{lat:.6f}", f"{lon:.6f}", "NaN" if elev is None else f"{elev:.3f}"]
                )


def _check_shape(path: Path, box: tuple[float, float, float, float]) -> dict[str, int]:
    """**The app's contract, checked in Python.** Every feature is a
    kind `splitDepthAsset` routes, at a level or band it draws, and
    inside the box. Returns the count per kind."""
    west, south, east, north = box
    body = json.loads(path.read_text())
    assert body["type"] == "FeatureCollection", body.get("type")
    kinds: dict[str, int] = {}
    levels: set[float] = set()
    bands: set[tuple[float, float]] = set()
    for feat in body["features"]:
        assert feat["type"] == "Feature"
        props = feat["properties"]
        kind = props["kind"]
        kinds[kind] = kinds.get(kind, 0) + 1
        geom = feat["geometry"]
        if kind == "contour":
            assert geom["type"] == "LineString"
            levels.add(props["depthM"])
            pts = geom["coordinates"]
        elif kind in ("shallow", "depth-band"):
            assert geom["type"] == "Polygon"
            bands.add((props["fromM"], props["toM"]))
            for ring in geom["coordinates"]:
                assert ring[0] == ring[-1], "ring not closed"
            pts = [p for ring in geom["coordinates"] for p in ring]
        elif kind == DRYING_KIND:
            # B-1857: a polygon carrying its kind and nothing else, so
            # an older app build skips it rather than misreading it.
            assert geom["type"] == "Polygon"
            assert set(props) == {"kind"}, props
            for ring in geom["coordinates"]:
                assert ring[0] == ring[-1], "ring not closed"
            pts = [p for ring in geom["coordinates"] for p in ring]
        elif kind == "label":
            assert geom["type"] == "Point"
            assert props["depthM"] in DEPTH_LEVELS_M
            pts = [geom["coordinates"]]
        else:
            raise AssertionError(f"unknown kind {kind!r}")
        for lon, lat in pts:
            assert west - 1e-4 <= lon <= east + 1e-4, lon
            assert south - 1e-4 <= lat <= north + 1e-4, lat
    assert set(kinds) == {
        "contour", "shallow", "depth-band", DRYING_KIND, "label"
    }, kinds
    assert levels == set(DEPTH_LEVELS_M), levels
    assert bands == {SHALLOW_BAND_M, *DEPTH_TINT_BANDS_M}, bands
    return kinds


def _inside(lon: float, lat: float, feats: list[dict]) -> bool:
    """Even-odd point-in-polygon over every ring of `feats`. For the
    self-test only: no dependency beyond the standard library."""
    hit = False
    for feat in feats:
        for ring in feat["geometry"]["coordinates"]:
            for (x0, y0), (x1, y1) in zip(ring[:-1], ring[1:]):
                if (y0 > lat) != (y1 > lat):
                    if lon < x0 + (lat - y0) * (x1 - x0) / (y1 - y0):
                        hit = not hit
    return hit


def self_test() -> None:
    """Build contours from a 60 x 60 synthetic seabed and check the
    file is the shape the app parses. **No network.**

    The seabed: dry land along the west edge, then a FORESHORE that
    dries from 4.5 m above LAT down to the 0 m line (B-1857), shoaling
    out to 60 m in the east, with a round shoal in the middle that
    comes up to 1 m. Every contour level, every band, the drying tint
    and the label anchors appear."""
    south, north, west, east = 50.0, 50.3, -4.0, -3.6
    box = (west, south, east, north)
    n = 60
    lats = [south + (north - south) * i / (n - 1) for i in range(n)]
    lons = [west + (east - west) * j / (n - 1) for j in range(n)]
    cell = (east - west) / (n - 1)
    shore_fx = 0.08  # the land / data edge, as a fraction of the box
    beach_fx = 0.20  # where the foreshore meets the 0 m line

    def seabed(lat: float, lon: float) -> float | None:
        fx = (lon - west) / (east - west)  # 0 west .. 1 east
        if fx < shore_fx:
            return None  # land interior: blank, as ERDDAP sends it
        if fx < beach_fx:
            # The drying foreshore: 4.5 m above LAT at the top of the
            # beach, falling to LAT at the drying line.
            return 4.5 * (beach_fx - fx) / (beach_fx - shore_fx)
        depth = 60.0 * (fx - beach_fx) / (1.0 - beach_fx)
        # A round shoal at the centre that comes up to 1 m.
        d = math.hypot((lat - 50.15) / 0.05, (lon + 3.8) / 0.07)
        depth -= max(0.0, 1.0 - d) * (depth - 1.0)
        return -depth

    with tempfile.TemporaryDirectory() as tmp_s:
        tmp = Path(tmp_s)

        # 1. The box parser refuses what the relay refuses.
        assert parse_bbox("-9.5,38.6,-9.0,38.9") == (-9.5, 38.6, -9.0, 38.9)
        for bad in ("1,2,3", "a,b,c,d", "5,0,4,1", "0,5,1,4", "nan,0,1,1"):
            try:
                parse_bbox(bad)
            except ValueError:
                continue
            raise AssertionError(f"parse_bbox accepted {bad!r}")

        # 2. Outside EMODnet: nothing written, and no request made.
        assert clip_to_emodnet(-80.0, 25.0, -79.0, 26.0) is None
        assert clip_to_emodnet(-37.0, 38.0, -35.0, 39.0) == (-36.0, 38.0, -35.0, 39.0)
        far = tmp / "far.depth.json"
        assert run_bbox(-80.0, 25.0, -79.0, 26.0, far) is False
        assert not far.exists()

        # 3. A synthetic seabed, end to end: the same CSV parse, the
        #    same pipeline, the bundled shape.
        grid_csv = tmp / "grid.csv"
        _synthetic_csv(grid_csv, lats, lons, seabed)
        out = tmp / "cut_selftest.depth.json"
        assert run_bbox(*box, out, cache=grid_csv) is True
        assert not out.with_name(out.name + ".part").exists()
        kinds = _check_shape(out, box)

        # 3b. B-1857 — the shoal tint reaches the data edge. Before the
        #     fix the foreshore had no tint at all and read as bare
        #     water, the colour of the deepest sea.
        body = json.loads(out.read_text())
        drying = [f for f in body["features"] if f["properties"]["kind"] == DRYING_KIND]
        first_data = min(x for x in lons if (x - west) / (east - west) >= shore_fx)
        # The middle rows: Chaikin rounds the polygon's corners where
        # it meets the box edge, which is the box, not the beach.
        mid_lats = lats[20:40]
        for lat in mid_lats:
            # Up the beach, on the last data column, and mid-foreshore.
            for lon in (first_data + 0.25 * cell, first_data + 3 * cell):
                assert _inside(lon, lat, drying), ("not tinted", lon, lat)
            # Below the drying line it is the 0-5 m band's water, never
            # the drying tint's.
            assert not _inside(west + 0.30 * (east - west), lat, drying)
        reach = min(p[0] for f in drying for r in f["geometry"]["coordinates"] for p in r)
        # Half a cell past the data, never a whole one (the overspill
        # bound in `DRYING_KIND`'s header).
        assert first_data - 0.75 * cell < reach < first_data, (reach, first_data)

        # 3c. The drying tint is an ADDITION. Strip it and the file is
        #     byte-identical to a run with the tint switched off: the
        #     0 m line, every contour, label and band stay exactly put.
        global DRYING_TINT
        off = tmp / "no_drying.depth.json"
        lats_a, lons_a, depth_a = fetch_grid()
        DRYING_TINT = False
        try:
            fc_off, _ = extract(lats_a, lons_a, depth_a)
        finally:
            DRYING_TINT = True
        off.write_text(encode(fc_off))
        stripped = {
            "type": "FeatureCollection",
            "features": [
                f for f in body["features"] if f["properties"]["kind"] != DRYING_KIND
            ],
        }
        assert encode(stripped) == off.read_text(), "drying moved another feature"

        # 3d. Where a coastline cut says a no-data cell is WATER, it is
        #     unmeasured foreshore and is tinted, up to the reach; a
        #     no-data cell on LAND, or beyond the reach, never is.
        #     Columns: 0-4 blank, 5-6 drying, 7+ wet.
        g_lons = np.arange(10, dtype=float)
        g_lats = np.arange(5, dtype=float)
        g = np.full((5, 10), NODATA_DEPTH_M)
        g[:, 5:7] = -2.0
        g[:, 7:] = 3.0
        sea = [np.array([[1.5, -1.0], [20.0, -1.0], [20.0, 9.0], [1.5, 9.0], [1.5, -1.0]])]
        grid_w, floor_w = drying_grid(g, g_lats, g_lons, sea)
        tint = (grid_w >= floor_w) & (grid_w <= 0.0)
        # Cells 3 and 4 are within two cells of data and in the water.
        assert tint[:, 3:7].all(), grid_w
        # Cell 2 is in the water but three cells out: past the reach,
        # so it only takes the half-cell edge value, below the floor.
        assert not tint[:, 2].any() and (grid_w[:, 2] < floor_w).all()
        # Cells 0-1 are on land: never reached.
        assert not tint[:, :2].any()
        # Without the coastline cut, only the measured ground is tinted.
        grid_n, floor_n = drying_grid(g)
        tint_n = (grid_n >= floor_n) & (grid_n <= 0.0)
        assert tint_n[:, 5:7].all() and not tint_n[:, :5].any()

        # 4. Deterministic: the same grid writes the same bytes.
        again = tmp / "again.depth.json"
        assert run_bbox(*box, again, cache=grid_csv) is True
        assert again.read_text() == out.read_text(), "two runs differ"

        # 5. All blank: nothing written, and no failure.
        blank_csv = tmp / "blank.csv"
        _synthetic_csv(blank_csv, lats[:5], lons[:5], lambda _a, _b: None)
        none_out = tmp / "none.depth.json"
        assert run_bbox(*box, none_out, cache=blank_csv) is False
        assert not none_out.exists()

        # 6. All deep: water, but nothing shallower than 50 m to draw.
        deep_csv = tmp / "deep.csv"
        _synthetic_csv(deep_csv, lats[:10], lons[:10], lambda _a, _b: -900.0)
        deep_out = tmp / "deep.depth.json"
        assert run_bbox(*box, deep_out, cache=deep_csv) is False
        assert not deep_out.exists()

    print(
        "self-test passed: "
        + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items()))
    )


def cli(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Cut EMODnet depth contours for a named coast or any box.",
    )
    parser.add_argument(
        "regions",
        nargs="*",
        help=f"named coasts, or 'all' (default: {DEFAULT_REGION})",
    )
    parser.add_argument(
        "--bbox",
        help="W,S,E,N in degrees. Use the = form: --bbox=-9.5,38.6,-9.0,38.9",
    )
    parser.add_argument("--out", help="where --bbox writes its file")
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="build contours from a tiny synthetic grid, no network",
    )
    args = parser.parse_args(argv)

    if args.self_test:
        self_test()
        return 0

    if args.bbox is not None or args.out is not None:
        if args.bbox is None or args.out is None or args.regions:
            parser.error("--bbox and --out go together, without region names")
        try:
            box = parse_bbox(args.bbox)
        except ValueError as err:
            parser.error(str(err))
        run_bbox(*box, Path(args.out))
        # Written or not, this is success: no depth is an answer.
        return 0

    # `extract-depth-contours.py [region ...]` — no arguments re-cuts
    # the default region, exactly as before the table existed.
    # `all` cuts every region in one run.
    names = args.regions or [DEFAULT_REGION]
    if names == ["all"]:
        names = list(REGIONS)
    for name in names:
        configure(name)
        print(f"── {name} ── box {LAT_MIN},{LAT_MAX} {LON_MIN},{LON_MAX}")
        main()
    return 0


if __name__ == "__main__":
    sys.exit(cli(sys.argv[1:]))
