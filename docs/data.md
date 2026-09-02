# Data Setup

This project uses local data files for gravity models, Earth orientation, leap seconds, nutation, and TLE-based diagnostics. These files can be large or externally licensed, so they are not guaranteed to be committed to the repository.

Put downloaded files under the project-root `assets/` directory unless a scenario explicitly points somewhere else.

## Folder Layout

Recommended local folders:

```sh
mkdir -p assets/earth assets/moon assets/mars assets/mercury assets/output
```

Common expected files (diagnostic dependent):

```text
assets/earth/EGM2008.gfc.txt
assets/earth/GGM05S.gfc.txt
assets/egm2008_120.txt
assets/moon/JGL165P1.gfc.txt
assets/moon/gggrx_1200a_sha.tab.txt
assets/mars/ggm1025a.gfc.txt
assets/mars/jgm85f01.gfc.txt
assets/mercury/sha.Doppler.txt
assets/leap-seconds.list.txt
assets/EOP_20u24_C04_one_file_1962-now.txt
assets/EOP_C01_IAU1980_1846-now.txt
assets/EOP_C01_IAU2000_1846-now.txt
assets/nut_IAU1980.dat.txt
assets/nut_IERS1996.dat.txt
assets/tle_all.txt
assets/tle_iss.txt
```

Only a subset is needed for any one diagnostic or scenario.

## Scenario Paths

Scenario files can reference data with project-root-relative paths:

```json
"path": "/assets/earth/EGM2008.gfc.txt"
```

## Gravity Models

### Earth Spherical Harmonics

Common files:

```text
assets/earth/EGM2008.gfc.txt
assets/egm2008_120.txt
```

Sources:

- NGA / Earth Gravitational Model data: https://earth-info.nga.mil/
- NASA PDS Geosciences Node: https://pds-geosciences.wustl.edu/
- ICGEM gravity model archive: https://icgem.gfz.de/tom

Supported/experimental formats in the code include:

```text
gfc
icgem
egm
sha
nasa_sha
nasa
```

Scenario provider example:

```json
{
    "id": "EGM2008_gfc",
    "format": "gfc",
    "path": "assets/earth/EGM2008.gfc.txt",
    "lineskips": 0,
    "normalized": true
}
```

### Lunar and Planetary Gravity

Useful sources:

- NASA Planetary Data System: https://pds.nasa.gov/tools/about/
- PDS Geosciences Node: https://pds-geosciences.wustl.edu/
- PDS gravity models: https://pds-geosciences.wustl.edu/dataserv/gravity_models.htm

Example local files currently used by diagnostics/scenarios:

```text
assets/moon/JGL165P1.gfc.txt
assets/moon/gggrx_1200a_sha.tab.txt
assets/mars/ggm1025a.gfc.txt
assets/mars/jgm85f01.gfc.txt
assets/mercury/sha.Doppler.txt
```

## Leap Seconds

Expected file:

```text
assets/leap-seconds.list.txt
```

Source:

- IANA leap seconds list: https://data.iana.org/time-zones/data/leap-seconds.list

The current parser expects the downloaded text table and a caller-provided line-skip count.

## Earth Orientation Parameters

Expected files may include:

```text
assets/EOP_20u24_C04_one_file_1962-now.txt
assets/EOP_C01_IAU1980_1846-now.txt
assets/EOP_C01_IAU2000_1846-now.txt
assets/EOP2long.txt
```

Sources:

- IERS Earth orientation data: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop
- IERS data center: https://datacenter.iers.org/
- JPL EOP2: https://eop2-external.jpl.nasa.gov/

Common data products:

- IERS EOP C04, IAU 2000A: https://datacenter.iers.org/data/latestVersion/EOP_20u24_C04_one_file_1962-now.txt
- IERS EOP C01, IAU 1980: https://datacenter.iers.org/data/latestVersion/EOP_C01_IAU1980_1846-now.txt

Current parser usage:

```text
EOP_20u24_C04_one_file_1962-now.txt
  polar motion: xp, yp
  Earth rotation timing: UT1-UTC
  intended model path: IAU 2000A

EOP_C01_IAU1980_1846-now.txt
  polar motion: xp, yp
  Earth rotation timing: UT1-UTC, derived from the table value plus leap seconds
  intended model path: IAU 1980

EOP_C01_IAU2000_1846-now.txt
  polar motion: xp, yp
  Earth rotation timing: UT1-UTC, derived from the table value plus leap seconds
  intended model path: IAU 2000

EOP2long.txt
  polar motion: xp, yp
  Earth rotation timing: UT1-UTC, derived from TAI-UT1 plus leap seconds
  intended model path: JPL EOP2
```

Related quantities that are not read from these EOP files:

- Leap seconds / `TAI-UTC`: read from `assets/leap-seconds.list.txt`
- Nutation series coefficients: read from `assets/nut_IAU1980.dat.txt` or `assets/nut_IERS1996.dat.txt`
- Precession model coefficients: currently implemented from model formulas, not from an external EOP table

## Nutation Tables

Expected files:

```text
assets/nut_IAU1980.dat.txt
assets/nut_IERS1996.dat.txt
```

Sources:

- HPIERS nutation models: https://hpiers.obspm.fr/eop-pc/models/models.html
- HPIERS explanatory material: https://hpiers.obspm.fr/eoppc/bul/bulb/explanatory.html
- IAU 1980 table: https://hpiers.obspm.fr/eop-pc/models/nutations/nut_IAU1980.dat
- IERS 1996 table: https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3a.txt

## TLE Files

Expected files for current diagnostics include:

```text
assets/tle_all.txt
assets/tle_iss.txt
assets/tle_fake.txt
assets/tle_iss_fail.txt
assets/tle_trunc.txt
```

Useful sources:

- CelesTrak: https://celestrak.org/
- Space-Track: https://www.space-track.org/

TLE parsing supports raw TLE data and conversion to satellite initial states through the project ingest utilities.

## SPICE Toolkit and Kernels

Download the platform-specific CSPICE N0067 toolkit and the project's pinned
generic kernels:

```sh
./scripts/download_spice.sh
```

The script installs:

```text
external/cspice/             platform-specific CSPICE toolkit
assets/spice/de442s.bsp      compact DE442 planetary ephemeris
assets/spice/naif0012.tls    leap-second kernel
assets/spice/pck00011.tpc    planetary constants kernel
```

Existing validated files are reused. Pass `--force` to replace them, or use
`--toolkit-only` / `--kernels-only` to download one group. These directories are
ignored because CSPICE is platform-specific and the ephemeris is relatively
large.

The downloader selects the NAIF package for Apple Silicon macOS, Intel macOS,
x86-64 Linux, Windows/Visual C under Git Bash or MSYS, or Windows/Cygwin. Set
`CSPICE_PLATFORM` to an official NAIF platform directory when automatic
detection does not match the compiler used by CMake. Unsupported architectures
fail explicitly rather than installing an incompatible binary toolkit.

Enable the source adapter at configure time:

```sh
cmake -S . -B build -DASTROLIB_ENABLE_CSPICE=ON
```

References:

- SPK position API: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html
- NAIF target IDs: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
- NAIF frame IDs: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
- CSPICE toolkit downloads: https://naif.jpl.nasa.gov/naif/toolkit.html
- Generic kernels: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/

## Notes

- Data files are not committed to the repository. Check with data providers before redistributing downloaded assets.

## Recommended Future Data Sources

These sources are not all wired into the code yet, but they are useful candidates for future data providers, validation tools, or rendering/estimation features.

### Satellite Catalogs and Space-Object Data

- CelesTrak GP element sets: https://celestrak.org/NORAD/elements/
  - Use for GP/TLE, OMM XML/KVN, JSON, CSV, supplemental GP data, SATCAT metadata, and future public satellite ingestion.
- CelesTrak space data: https://celestrak.org/SpaceData/
  - Convenient EOP and space-weather CSV/legacy files.

### Solar-System Ephemerides

- JPL Horizons API: https://ssd-api.jpl.nasa.gov/doc/horizons.html
  - Sun, Moon, planet, asteroid, comet, and observer/vector/element ephemerides.
  - Useful for validating propagators and for generating table-backed ephemeris providers.
- NAIF SPICE generic kernels: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
  - Use later for higher-fidelity ephemeris and frame/orientation providers.
  - Relevant kernel families include `SPK` for trajectories, `PCK` for body orientation/shape, `LSK` for leap seconds, and `FK` for frames.

### GNSS and Space Geodesy

- NASA Earthdata / CDDIS GNSS: https://www.earthdata.nasa.gov/data/space-geodesy-techniques/gnss
  - GNSS precise orbit and clock products, RINEX observation/navigation data, and station products.

### Earth Orientation

- HPIERS EOP C04 directory: https://hpiers.obspm.fr/eoppc/eop/eopc04/
  - Direct C04 variants, 12-hour C04 files, and `dPsi/dEps` variants.

### Small Bodies

- Minor Planet Center data: https://minorplanetcenter.net/data
  - Asteroid/comet orbit catalogs, observation data, and small-body scenario generation.
  - Useful later for optical tracking examples or asteroid propagation/OD scenarios.

### Visualization and Star Backgrounds

- ESA Gaia Archive: https://gea.esac.esa.int/archive/
  - Star catalogs, sky backgrounds, inertial-frame visualization, and future optical sensor simulation.

### Space Weather and Atmosphere Inputs

- NOAA Space Weather Prediction Center products: https://www.spaceweather.gov/products/solar-and-geophysical-event-reports
  - Solar/geophysical activity inputs.
