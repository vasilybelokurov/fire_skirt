#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p run outputs logs

VENV_ACTIVATE_DEFAULT="$HOME/Work/venvs/.venv/bin/activate"
VENV_ACTIVATE="${VENV_ACTIVATE:-$VENV_ACTIVATE_DEFAULT}"
if [[ -f "$VENV_ACTIVATE" ]]; then
  # shellcheck disable=SC1090
  source "$VENV_ACTIVATE"
fi

PYTHON_BIN="${PYTHON_BIN:-python}"

SNAPSHOT_GLOB_DEFAULT="$HOME/data/sims/fire/m12i_res7100/output/snapshot_600.*.hdf5"
HOST_COORDS_DEFAULT="$HOME/data/sims/fire/m12i_res7100/track/host_coordinates.hdf5"
SNAPSHOT_GLOB="${SNAPSHOT_GLOB:-$SNAPSHOT_GLOB_DEFAULT}"
HOST_COORDS="${HOST_COORDS:-$HOST_COORDS_DEFAULT}"
POS_TO_KPC="${POS_TO_KPC:-}"
HOST_POS_TO_KPC="${HOST_POS_TO_KPC:-1.0}"
MASS_TO_MSUN="${MASS_TO_MSUN:-}"
METALLICITY_SCALE="${METALLICITY_SCALE:-}"
METALLICITY_COLUMN="${METALLICITY_COLUMN:-}"
NUM_PACKETS="${NUM_PACKETS:-}"
NPIX="${NPIX:-}"
FOV_KPC="${FOV_KPC:-}"
WMIN="${WMIN:-}"
WMAX="${WMAX:-}"
NWAVES="${NWAVES:-}"
MIN_LEVEL="${MIN_LEVEL:-}"
MAX_LEVEL="${MAX_LEVEL:-}"
MAX_DUST_FRACTION="${MAX_DUST_FRACTION:-}"
export SNAPSHOT_GLOB

if [[ -z "$POS_TO_KPC" ]]; then
  POS_TO_KPC="$("$PYTHON_BIN" - <<'PY'
import glob
import os
import h5py

snap_glob = os.path.expanduser(os.environ.get("SNAPSHOT_GLOB", ""))
paths = sorted(glob.glob(snap_glob))
if not paths:
    raise SystemExit(1)
with h5py.File(paths[0], "r") as f:
    h = f["Header"].attrs.get("HubbleParam", None)
if h:
    print(1.0 / float(h))
PY
  )" || POS_TO_KPC=""
fi

"$PYTHON_BIN" code/fire_to_skirt_tables.py \
  --snapshot-glob "$SNAPSHOT_GLOB" \
  --host-coords "$HOST_COORDS" \
  ${POS_TO_KPC:+--pos-to-kpc "$POS_TO_KPC"} \
  --host-pos-to-kpc "$HOST_POS_TO_KPC" \
  ${MASS_TO_MSUN:+--mass-to-msun "$MASS_TO_MSUN"} \
  ${METALLICITY_SCALE:+--metallicity-scale "$METALLICITY_SCALE"} \
  ${METALLICITY_COLUMN:+--metallicity-column "$METALLICITY_COLUMN"}
"$PYTHON_BIN" code/make_views.py
"$PYTHON_BIN" code/build_ski.py \
  ${NUM_PACKETS:+--num-packets "$NUM_PACKETS"} \
  ${NPIX:+--pixels "$NPIX"} \
  ${FOV_KPC:+--fov-kpc "$FOV_KPC"} \
  ${WMIN:+--min-wavelength "$WMIN"} \
  ${WMAX:+--max-wavelength "$WMAX"} \
  ${NWAVES:+--num-wavelengths "$NWAVES"} \
  ${MIN_LEVEL:+--min-level "$MIN_LEVEL"} \
  ${MAX_LEVEL:+--max-level "$MAX_LEVEL"} \
  ${MAX_DUST_FRACTION:+--max-dust-fraction "$MAX_DUST_FRACTION"}

"$PYTHON_BIN" code/quickcheck.py --run-dir run --skip-images

skirt -i run -o run run/m12i_600.ski > logs/skirt_run.log 2>&1

"$PYTHON_BIN" code/quickcheck.py --run-dir run

TS="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="outputs/${TS}"
mkdir -p "$OUT_DIR"

cp run/m12i_600.ski "$OUT_DIR"/
cp run/views.json "$OUT_DIR"/
cp run/meta.json "$OUT_DIR"/

shopt -s nullglob
for f in run/*view_*_*.fits; do
  cp "$f" "$OUT_DIR"/
done

printf "Outputs copied to %s\n" "$OUT_DIR"
