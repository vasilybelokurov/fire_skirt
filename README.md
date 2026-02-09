# FIRE m12i -> SKIRT pipeline

This repository builds a reproducible SKIRT run for FIRE m12i snapshot 600 (z=0) with many viewing angles.

## One-command run

```bash
bash code/run_skirt.sh
```

That command will:
1. Build `run/stars.txt`, `run/gas.txt`, and `run/meta.json`.
2. Generate `run/views.json` (Fibonacci sphere).
3. Generate `run/m12i_600.ski`.
4. Run SKIRT and write outputs under `run/`.
5. Run quick checks before and after SKIRT.
6. Copy final images + metadata to `outputs/<timestamp>/`.

## Python environment

If you do not already have `numpy` and `h5py`, create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`code/run_skirt.sh` will `source ~/Work/venvs/.venv/bin/activate` if it exists. You can override this with:

```bash
VENV_ACTIVATE=/path/to/activate bash code/run_skirt.sh
```

You can also override data locations:

```bash
SNAPSHOT_GLOB=/path/to/snapshot_600.*.hdf5 HOST_COORDS=/path/to/host_coordinates.hdf5 bash code/run_skirt.sh
```

If your snapshot metallicity is stored as multiple columns (e.g., FIRE element abundances),
the pipeline uses column 0 by default. Override with:

```bash
METALLICITY_COLUMN=0 bash code/run_skirt.sh
```

If host coordinates are already in kpc while snapshot coordinates are in kpc/h, set:

```bash
HOST_POS_TO_KPC=1.0 POS_TO_KPC=1.0/0.702 bash code/run_skirt.sh
```

## Scripts

- `code/fire_to_skirt_tables.py`
  - Reads `data/output/snapshot_600.*.hdf5` and `data/track/host_coordinates.hdf5`.
  - Centers on the host at snapshot 600 and cuts to `R_CUT_KPC`.
  - Writes `run/stars.txt` and `run/gas.txt` with SKIRT column headers.
  - Writes `run/meta.json` with detected datasets, unit info, and assumptions.

- `code/make_views.py`
  - Generates `N` viewing directions on a Fibonacci sphere.
  - Writes `run/views.json` with `theta`/`phi` in degrees.

- `code/build_ski.py`
  - Builds `run/m12i_600.ski` from `run/meta.json` and `run/views.json`.
  - Uses `ParticleSource` + `BruzualCharlotSEDFamily` by default (no Cloud_* packs).
  - Uses `ParticleMedium` with `dust_to_metals` from meta and a built-in dust mix.
  - Creates one `FrameInstrument` per view with `instrumentName=view_###`.

- `code/run_skirt.sh`
  - Runs the full pipeline and logs SKIRT to `logs/skirt_run.log`.

- `code/quickcheck.py`
  - Sanity checks particle tables and (optionally) checks that images exist for all views.

- `code/convert_fits_to_png.py`
  - Converts SKIRT FITS output to PNG/JPG images.
  - Collapses the wavelength axis by summing, and applies an asinh stretch by default.

- `code/make_rgb.py`
  - Makes false-color RGB PNGs from SKIRT wavelength cubes using Lupton asinh stretch.
  - Includes presets for HST-like optical and JWST-like NIR false color.

- `code/camera_positions.py`
  - Computes camera positions (in Mpc and kpc) from `views.json`.
  - Writes `camera_positions.txt` next to `views.json` by default.

## Notes and knobs

- `fire_to_skirt_tables.py` uses `POS_TO_KPC` (default 1.0) and records header units in `run/meta.json`.
- If the metallicity dataset has multiple columns, `--metallicity-column` selects which column to use (default 0).
- Dust mass is derived in SKIRT with `Mdust = f_dtm * Z * Mgas` by setting `massFraction=f_dtm` and `importMetallicity=true`.
- The SKIRT simulation mode is `ExtinctionOnly` for a modest first run.

Grid refinement knobs (passed through `code/run_skirt.sh`):

```bash
MIN_LEVEL=3 MAX_LEVEL=9 MAX_DUST_FRACTION=1e-6 bash code/run_skirt.sh
```

### Bruzual-Charlot fallback

If the Bruzual-Charlot SSP resources are unavailable, you can switch to a built-in fallback SED:

```bash
python code/build_ski.py --sed-family BlackBodySEDFamily --bb-radius-km 6.96e5 --bb-temp-k 5000
```

This will write `run/stars_bb.txt` with constant blackbody parameters and use it in the generated `.ski`.
This is a crude fallback and should be treated as a last resort.

## FITS -> PNG/JPG

Convert all views in an output folder:

```bash
source ~/Work/venvs/.venv/bin/activate
python code/convert_fits_to_png.py --input-dir outputs/<timestamp>
```

This writes images to `outputs/<timestamp>/images/`.

## False-color RGB

Make HST-like optical RGB images:

```bash
source ~/Work/venvs/.venv/bin/activate
python code/make_rgb.py --input-dir outputs/<timestamp> --preset hst --out-dir outputs/<timestamp>/rgb_hst
```

Make JWST-like NIR RGB images:

```bash
source ~/Work/venvs/.venv/bin/activate
python code/make_rgb.py --input-dir outputs/<timestamp> --preset jwst --out-dir outputs/<timestamp>/rgb_jwst
```

You can also specify custom bands (micron):

```bash
python code/make_rgb.py --input-dir outputs/<timestamp> \
  --band-b 0.40 0.50 --band-g 0.55 0.70 --band-r 0.75 0.90 \
  --out-dir outputs/<timestamp>/rgb_custom
```

## Camera view docs

- `CUSTOM_VIEWS.md` explains how to define custom camera directions.
- `CAMERA_DISTANCE.md` explains where the instrument distance is set.
- `CAMERA_POSITIONS.md` explains how camera positions are computed.
- `RECIPE.md` records a working run configuration with visible dust lanes.
