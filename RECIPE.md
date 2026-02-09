# FIRE m12i SKIRT recipe (working run with visible dust lanes)

This document records the exact recipe that produced the clean, non‑stripy, blue‑band images with visible dust lanes.

## Environment

- Default venv: `source ~/Work/venvs/.venv/bin/activate`
- Repo root: `/Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/fire_skirt`

## Data

- Snapshots: `~/data/sims/fire/m12i_res7100/output/snapshot_600.*.hdf5`
- Host track: `~/data/sims/fire/m12i_res7100/track/host_coordinates.hdf5`

## Key fixes

1) **Metallicity column selection**
   - FIRE `Metallicity` is multi‑column (11 columns). Summing columns inflated Z.
   - Use **column 0 only** (default now) to get mass fraction metallicity.
   - Set `METALLICITY_COLUMN=0` explicitly in the run command.

2) **Photon starvation / striping fix**
   - Stripy/checkerboard artifacts were due to too few contributions per wavelength bin per pixel.
   - Fix by reducing sampling load:
     - `NPIX=512` (instead of 1024)
     - `NWAVES=10` within the blue band
   - Keep `NUM_PACKETS=5e6`.

3) **Grid refinement**
   - Use deeper refinement so dust grid does not imprint large slabs:
     - `MIN_LEVEL=3`
     - `MAX_LEVEL=9`
     - `MAX_DUST_FRACTION=1e-6`

## Exact run command (the one that worked)

```bash
MASS_TO_MSUN=1.4245e10 METALLICITY_SCALE=1.0 METALLICITY_COLUMN=0 \
NUM_PACKETS=5e6 NPIX=512 FOV_KPC=40 WMIN=0.35 WMAX=0.55 NWAVES=10 \
MIN_LEVEL=3 MAX_LEVEL=9 MAX_DUST_FRACTION=1e-6 \
bash code/run_skirt.sh
```

Notes:
- `MASS_TO_MSUN=1.4245e10` is the FIRE m12i conversion used in this project.
- `METALLICITY_SCALE=1.0` because snapshot Z is already a mass fraction.
- `POS_TO_KPC` is auto‑derived from `HubbleParam` in the snapshot header.
- `HOST_POS_TO_KPC=1.0` is correct for host positions in kpc.

## Output location (from this run)

- `outputs/20260209_214702/`
  - `m12i_600.ski`
  - `views.json`
  - `meta.json`
  - `m12i_600_view_###_total.fits`

## PNG conversion (blue band)

```bash
source ~/Work/venvs/.venv/bin/activate
python code/convert_fits_to_png.py \
  --input-dir outputs/20260209_214702 \
  --out-dir outputs/20260209_214702/images_blue_1000_smooth6 \
  --wmin 0.35 --wmax 0.55 \
  --fwhm-pix 6 --cmap gray --pmax 99.9 --Q 4 --stretch 0.6 \
  --size-px 1000 --resample bilinear
```

Example image:
- `outputs/20260209_214702/images_blue_1000_smooth6/m12i_600_view_015_total.png`

## Why this worked

- **Correct metallicity units** restored dust mass to realistic levels.
- **Lower NPIX + fewer wavelength bins** increased photons per pixel per bin.
- **Deeper grid refinement** reduced large‑cell dust artifacts.
- **PSF smoothing** (`fwhm=6`) reduced residual speckle without smearing lanes.

## Optional next steps (not required for the working recipe)

- Increase `NUM_PACKETS` to `2e7` for even cleaner blue‑band images.
- Add RGB false‑color using `code/make_rgb.py` once band‑separated images are produced.

