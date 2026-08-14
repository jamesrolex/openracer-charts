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
