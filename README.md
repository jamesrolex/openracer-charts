# OpenYachtRacer — chart regions

Offline chart extracts for [OpenYachtRacer](https://github.com/jamesrolex/openracer),
an open-source sailing race navigation system.

**This repository holds map data only.** No application code. It exists
so the app can download its offline basemap without anyone needing a
GitHub account — the whole point of an offline-first chart is that it
works when nothing else does.

## Regions

| Region | File | Size | Covers |
|---|---|---|---|
| Cardigan Bay | `cardigan-bay.pmtiles` | 13.5 MB | Abersoch and Cardigan Bay, North Wales |

Download them from [Releases](https://github.com/jamesrolex/openracer-charts/releases).

## Format

[PMTiles](https://docs.protomaps.com/pmtiles/) — a single-file archive
of vector tiles, read directly by MapLibre with no tile server. One
file, one download, works in aeroplane mode afterwards.

## Licence and attribution

Map data © [OpenStreetMap](https://www.openstreetmap.org/copyright)
contributors, available under the
[Open Database Licence (ODbL)](https://opendatacommons.org/licenses/odbl/).

Extracted with [Protomaps](https://protomaps.com/) basemaps. Any
product using this data must carry the same attribution.

## Rebuilding

The extract is reproducible in about 22 seconds from a documented
command. See `tools/charts/EXTRACT-PMTILES.md` in the main repository.

## Areas cut on request

When a sailor draws a box on the chart in OpenYachtRacer and nobody has
asked for that water before, the app's relay starts the **Cut a chart
area** workflow (`.github/workflows/cut.yml`). It cuts the box from the
newest Protomaps build at max zoom 14, the same recipe as the regions
above, and publishes it to the
[`user-cuts`](https://github.com/jamesrolex/openracer-charts/releases/tag/user-cuts)
release. Boxes are snapped to a 0.1° grid and capped at about 60 nm a
side, so the same water is only ever cut once.
