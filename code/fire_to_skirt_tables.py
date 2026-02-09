#!/usr/bin/env python3
import argparse
import glob
import json
import os
import time

import h5py
import numpy as np

CM_PER_KPC = 3.085677581e21
G_PER_MSUN = 1.98847e33
SEC_PER_GYR = 3.15576e16

DEFAULT_SNAPSHOT_GLOB = "data/output/snapshot_600.*.hdf5"
DEFAULT_HOST_COORDS = "data/track/host_coordinates.hdf5"
DEFAULT_OUT_DIR = "run"
DEFAULT_R_CUT_KPC = 60.0
DEFAULT_DUST_TO_METALS = 0.4
DEFAULT_POS_TO_KPC = 1.0
DEFAULT_HSML_FALLBACK_KPC = 0.5

STAR_METALLICITY_KEYS = [
    "GFM_Metallicity",
    "Metallicity",
    "StellarMetallicity",
    "Z",
]
STAR_FORM_TIME_KEYS = [
    "GFM_StellarFormationTime",
    "StellarFormationTime",
    "FormationTime",
]
STAR_INIT_MASS_KEYS = [
    "GFM_InitialMass",
    "InitialMass",
    "BirthMass",
]
GAS_METALLICITY_KEYS = [
    "GFM_Metallicity",
    "Metallicity",
    "Z",
]
GAS_HSML_KEYS = [
    "SmoothingLength",
    "HSML",
    "Hsml",
]


def _jsonify(value):
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (bytes, bytearray)):
        try:
            return value.decode("utf-8")
        except Exception:
            return str(value)
    return value


def _read_header_attrs(snapshot_path):
    with h5py.File(snapshot_path, "r") as f:
        if "Header" not in f:
            return {}
        attrs = dict(f["Header"].attrs)
        return {k: _jsonify(v) for k, v in attrs.items()}


def _find_first_key(group, candidates):
    for key in candidates:
        if key in group:
            return key
    return None


def _flatten_metallicity(z, meta, column=0):
    if z is None:
        return None
    arr = np.asarray(z)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2 and arr.shape[1] == 1:
        return arr[:, 0]
    if arr.ndim == 2 and arr.shape[1] > 1:
        col = int(column)
        if col < 0 or col >= arr.shape[1]:
            col = 0
        meta["metallicity_note"] = (
            f"Metallicity array has {arr.shape[1]} columns; using column {col}."
        )
        return arr[:, col]
    raise ValueError(f"Unsupported metallicity array shape: {arr.shape}")


def _cosmic_time_gyr_from_scale_factor(a_vals, omega_m, omega_l, h0):
    a_vals = np.asarray(a_vals)
    a_max = float(np.nanmax(a_vals))
    a_max = max(a_max, 1e-4)
    a_grid = np.linspace(1e-4, a_max, 10000)
    omega_k = 1.0 - omega_m - omega_l
    e = np.sqrt(omega_m / a_grid**3 + omega_k / a_grid**2 + omega_l)
    integrand = 1.0 / (a_grid * e)
    da = np.diff(a_grid)
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * da)])

    mpc_km = 3.085677581e19
    h0_s = 100.0 * h0 / mpc_km
    t_h_gyr = (1.0 / h0_s) / SEC_PER_GYR
    t_grid = integ * t_h_gyr
    return np.interp(a_vals, a_grid, t_grid)


def _compute_ages_gyr(form_time, header, meta):
    if form_time is None:
        meta["age_note"] = "No formation time dataset found; setting ages to 0 Gyr."
        return None

    a_now = header.get("Time", None)
    redshift = header.get("Redshift", None)
    omega_m = header.get("Omega0", None)
    omega_l = header.get("OmegaLambda", None)
    h0 = header.get("HubbleParam", None)
    unit_time_s = header.get("UnitTimeInS", None)

    ft = np.asarray(form_time)
    if ft.size == 0:
        return ft.astype(np.float64)
    ft_max = float(np.nanmax(ft)) if ft.size else 0.0
    a_now_val = float(a_now) if a_now is not None else None

    if a_now_val is not None and a_now_val <= 1.1 and ft_max <= 1.1:
        if omega_m is not None and omega_l is not None and h0 is not None:
            t_now = _cosmic_time_gyr_from_scale_factor(a_now_val, omega_m, omega_l, h0)
            t_form = _cosmic_time_gyr_from_scale_factor(ft, omega_m, omega_l, h0)
            ages = t_now - t_form
            ages = np.clip(ages, 0.0, None)
            meta["age_note"] = (
                "Ages computed from scale factor using LCDM with Omega0, OmegaLambda, HubbleParam."
            )
            meta["age_now_gyr"] = _jsonify(t_now)
            return ages
        meta["age_note"] = (
            "Formation time looks like scale factor but cosmology missing; using delta(a) as age proxy."
        )
        ages = a_now_val - ft
        ages = np.clip(ages, 0.0, None)
        return ages

    if unit_time_s is not None and a_now_val is not None:
        t_unit_gyr = float(unit_time_s) / SEC_PER_GYR
        ages = (a_now_val - ft) * t_unit_gyr
        ages = np.clip(ages, 0.0, None)
        meta["age_note"] = "Ages computed from code time units using UnitTimeInS."
        return ages

    if a_now_val is not None:
        ages = a_now_val - ft
        ages = np.clip(ages, 0.0, None)
        meta["age_note"] = "Ages computed as Time - formation time; assuming both are in Gyr."
        return ages

    meta["age_note"] = "Could not interpret formation times; setting ages to 0 Gyr."
    return None


def _load_host_position(host_path, snapnum):
    if not os.path.exists(host_path):
        raise FileNotFoundError(f"Host coordinates file not found: {host_path}")

    with h5py.File(host_path, "r") as f:
        datasets = {}

        def _collect(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets[name] = obj

        f.visititems(_collect)

        snap_index = None
        snap_key = None
        for name, ds in datasets.items():
            lname = name.lower()
            if any(k in lname for k in ["snap", "snapshot", "snapnum", "snap_num", "index"]):
                arr = np.array(ds)
                if arr.ndim == 1 and arr.size:
                    if snapnum in arr:
                        snap_index = int(np.where(arr == snapnum)[0][0])
                        snap_key = name
                        break

        candidates = []
        for name, ds in datasets.items():
            if ds.ndim == 2:
                shape_ok = ds.shape[-1] == 3 or ds.shape[0] == 3
            elif ds.ndim == 3:
                shape_ok = 3 in ds.shape
            else:
                shape_ok = False

            if shape_ok:
                score = 0
                lname = name.lower()
                if "host" in lname:
                    score += 3
                if "coord" in lname or "pos" in lname:
                    score += 2
                candidates.append((score, name, ds))

        if not candidates:
            raise RuntimeError("Could not find a (N,3) or (3,N) dataset in host coordinates file.")

        candidates.sort(key=lambda x: (-x[0], x[1]))
        _, pos_key, pos_ds = candidates[0]

        if snap_index is None:
            # assume row index equals snapshot number
            if pos_ds.shape[0] > snapnum:
                snap_index = snapnum
                snap_key = None
            elif pos_ds.shape[-1] > snapnum and pos_ds.shape[0] == 3:
                snap_index = snapnum
                snap_key = None
            else:
                raise RuntimeError(
                    "Could not locate snapshot index in host file; add a snapshot index dataset or adjust snapnum."
                )

        if pos_ds.ndim == 2:
            if pos_ds.shape[-1] == 3:
                pos = np.array(pos_ds[snap_index])
            else:
                pos = np.array(pos_ds[:, snap_index])
        else:
            # handle shapes like (N,1,3) or (1,N,3) or (N,3,1)
            if pos_ds.shape[-1] == 3:
                if pos_ds.shape[1] == 1:
                    pos = np.array(pos_ds[snap_index, 0])
                elif pos_ds.shape[0] == 1:
                    pos = np.array(pos_ds[0, snap_index])
                else:
                    pos = np.array(pos_ds[snap_index, 0])
            elif pos_ds.shape[0] == 3:
                if pos_ds.shape[2] == 1:
                    pos = np.array(pos_ds[:, snap_index, 0])
                else:
                    pos = np.array(pos_ds[:, snap_index, 0])
            else:
                raise RuntimeError(f"Unhandled host position dataset shape: {pos_ds.shape}")

        return pos, {"host_position_dataset": pos_key, "host_index_dataset": snap_key}


def _write_table(path, header_lines, data):
    with open(path, "w", encoding="utf-8") as f:
        for line in header_lines:
            f.write(line + "\n")
        np.savetxt(f, data, fmt="%.8e")


def main():
    parser = argparse.ArgumentParser(description="Convert FIRE snapshot to SKIRT particle tables.")
    parser.add_argument("--snapshot-glob", default=DEFAULT_SNAPSHOT_GLOB)
    parser.add_argument("--host-coords", default=DEFAULT_HOST_COORDS)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--r-cut-kpc", type=float, default=DEFAULT_R_CUT_KPC)
    parser.add_argument("--dust-to-metals", type=float, default=DEFAULT_DUST_TO_METALS)
    parser.add_argument("--pos-to-kpc", type=float, default=DEFAULT_POS_TO_KPC)
    parser.add_argument("--host-pos-to-kpc", type=float, default=None)
    parser.add_argument("--mass-to-msun", type=float, default=None)
    parser.add_argument("--metallicity-scale", type=float, default=1.0)
    parser.add_argument("--metallicity-column", type=int, default=0)
    parser.add_argument("--hsml-fallback-kpc", type=float, default=DEFAULT_HSML_FALLBACK_KPC)
    parser.add_argument("--snapnum", type=int, default=600)
    args = parser.parse_args()

    snapshot_glob = os.path.expanduser(args.snapshot_glob)
    host_coords = os.path.expanduser(args.host_coords)

    snapshot_paths = sorted(glob.glob(snapshot_glob))
    if not snapshot_paths:
        raise FileNotFoundError(f"No snapshot files found for glob: {snapshot_glob}")

    os.makedirs(args.out_dir, exist_ok=True)

    header = _read_header_attrs(snapshot_paths[0])
    meta = {
        "snapshot_glob": snapshot_glob,
        "snapshot_files": snapshot_paths,
        "host_coords": host_coords,
        "snapnum": args.snapnum,
        "R_CUT_KPC": args.r_cut_kpc,
        "dust_to_metals": args.dust_to_metals,
        "pos_to_kpc": args.pos_to_kpc,
        "host_pos_to_kpc": args.host_pos_to_kpc if args.host_pos_to_kpc is not None else args.pos_to_kpc,
        "hsml_fallback_kpc": args.hsml_fallback_kpc,
        "mass_to_msun_arg": args.mass_to_msun,
        "metallicity_scale": args.metallicity_scale,
        "header_attrs": header,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    if "UnitLengthInCm" in header:
        meta["pos_to_kpc_from_header"] = float(header["UnitLengthInCm"]) / CM_PER_KPC
    if "UnitMassInG" in header:
        meta["mass_to_msun_from_header"] = float(header["UnitMassInG"]) / G_PER_MSUN
    if "UnitMassInG" not in header:
        meta["mass_to_msun_note"] = (
            "UnitMassInG missing; using --mass-to-msun if provided, else 1.0."
        )
    if "UnitLengthInCm" not in header:
        meta["pos_to_kpc_note"] = (
            "UnitLengthInCm missing; using POS_TO_KPC argument (default 1.0 unless overridden)."
        )
        if "HubbleParam" in header:
            meta["pos_to_kpc_hint_from_hubble"] = 1.0 / float(header["HubbleParam"])

    pos_to_kpc = args.pos_to_kpc
    host_pos_to_kpc = args.host_pos_to_kpc if args.host_pos_to_kpc is not None else args.pos_to_kpc
    if args.mass_to_msun is not None:
        mass_to_msun = float(args.mass_to_msun)
    else:
        mass_to_msun = float(meta.get("mass_to_msun_from_header", 1.0))
    meta["mass_to_msun_used"] = mass_to_msun
    meta["metallicity_column_used"] = int(args.metallicity_column)

    host_pos, host_meta = _load_host_position(host_coords, args.snapnum)
    meta.update(host_meta)
    meta["host_position_raw"] = _jsonify(host_pos)
    host_pos_kpc = host_pos * host_pos_to_kpc
    meta["host_position_kpc"] = _jsonify(host_pos_kpc)

    star_rows = []
    gas_rows = []

    star_keys_used = {}
    gas_keys_used = {}

    total_stars = 0
    total_gas = 0
    kept_stars = 0
    kept_gas = 0

    for path in snapshot_paths:
        with h5py.File(path, "r") as f:
            if "PartType4" in f:
                g = f["PartType4"]
                coords = np.array(g["Coordinates"], dtype=np.float64)
                masses = np.array(g["Masses"], dtype=np.float64)
                z_key = _find_first_key(g, STAR_METALLICITY_KEYS)
                ft_key = _find_first_key(g, STAR_FORM_TIME_KEYS)
                im_key = _find_first_key(g, STAR_INIT_MASS_KEYS)

                star_keys_used.setdefault("Coordinates", "Coordinates")
                star_keys_used.setdefault("Masses", "Masses")
                if z_key:
                    star_keys_used.setdefault("Metallicity", z_key)
                if ft_key:
                    star_keys_used.setdefault("FormationTime", ft_key)
                if im_key:
                    star_keys_used.setdefault("InitialMass", im_key)

                z = np.array(g[z_key], dtype=np.float64) if z_key else None
                z = _flatten_metallicity(z, meta, column=args.metallicity_column)
                ft = np.array(g[ft_key], dtype=np.float64) if ft_key else None
                im = np.array(g[im_key], dtype=np.float64) if im_key else None

                coords_kpc = coords * pos_to_kpc - host_pos_kpc
                r = np.sqrt((coords_kpc ** 2).sum(axis=1))
                mask = r <= args.r_cut_kpc

                total_stars += coords.shape[0]
                kept_stars += int(mask.sum())

                coords_kpc = coords_kpc[mask]
                masses = masses[mask]
                if z is not None:
                    z = z[mask]
                if ft is not None:
                    ft = ft[mask]
                if im is not None:
                    im = im[mask]

                hsml = np.full(coords_kpc.shape[0], args.hsml_fallback_kpc, dtype=np.float64)

                ages = _compute_ages_gyr(ft, header, meta) if ft is not None else None
                if ages is None:
                    ages = np.zeros(coords_kpc.shape[0], dtype=np.float64)

                if im is None:
                    meta["initial_mass_note"] = "Initial mass not found; using Masses as Minit."
                    im = masses

                minit = im * mass_to_msun
                z = np.zeros(coords_kpc.shape[0], dtype=np.float64) if z is None else z
                z = z * args.metallicity_scale
                ages = np.asarray(ages, dtype=np.float64)

                rows = np.column_stack([
                    coords_kpc[:, 0],
                    coords_kpc[:, 1],
                    coords_kpc[:, 2],
                    hsml,
                    minit,
                    z,
                    ages,
                ])
                star_rows.append(rows)

            if "PartType0" in f:
                g = f["PartType0"]
                coords = np.array(g["Coordinates"], dtype=np.float64)
                masses = np.array(g["Masses"], dtype=np.float64)
                z_key = _find_first_key(g, GAS_METALLICITY_KEYS)
                h_key = _find_first_key(g, GAS_HSML_KEYS)

                gas_keys_used.setdefault("Coordinates", "Coordinates")
                gas_keys_used.setdefault("Masses", "Masses")
                if z_key:
                    gas_keys_used.setdefault("Metallicity", z_key)
                if h_key:
                    gas_keys_used.setdefault("SmoothingLength", h_key)

                z = np.array(g[z_key], dtype=np.float64) if z_key else None
                z = _flatten_metallicity(z, meta, column=args.metallicity_column)
                hsml = np.array(g[h_key], dtype=np.float64) if h_key else None

                coords_kpc = coords * pos_to_kpc - host_pos_kpc
                r = np.sqrt((coords_kpc ** 2).sum(axis=1))
                mask = r <= args.r_cut_kpc

                total_gas += coords.shape[0]
                kept_gas += int(mask.sum())

                coords_kpc = coords_kpc[mask]
                masses = masses[mask]
                if z is not None:
                    z = z[mask]
                if hsml is not None:
                    hsml = hsml[mask]

                if hsml is None:
                    hsml = np.full(coords_kpc.shape[0], args.hsml_fallback_kpc, dtype=np.float64)
                    meta["hsml_note"] = "SmoothingLength missing; using constant fallback."
                else:
                    hsml = hsml * pos_to_kpc

                mgas = masses * mass_to_msun
                z = np.zeros(coords_kpc.shape[0], dtype=np.float64) if z is None else z
                z = z * args.metallicity_scale

                rows = np.column_stack([
                    coords_kpc[:, 0],
                    coords_kpc[:, 1],
                    coords_kpc[:, 2],
                    hsml,
                    mgas,
                    z,
                ])
                gas_rows.append(rows)

    meta["star_keys_used"] = star_keys_used
    meta["gas_keys_used"] = gas_keys_used
    meta["counts"] = {
        "stars_total": total_stars,
        "stars_kept": kept_stars,
        "gas_total": total_gas,
        "gas_kept": kept_gas,
    }

    stars = np.vstack(star_rows) if star_rows else np.zeros((0, 7))
    gas = np.vstack(gas_rows) if gas_rows else np.zeros((0, 6))

    star_header = [
        "# column 1: x (kpc)",
        "# column 2: y (kpc)",
        "# column 3: z (kpc)",
        "# column 4: hsml (kpc)",
        "# column 5: Minit (Msun)",
        "# column 6: Z (1)",
        "# column 7: age (Gyr)",
    ]
    gas_header = [
        "# column 1: x (kpc)",
        "# column 2: y (kpc)",
        "# column 3: z (kpc)",
        "# column 4: hsml (kpc)",
        "# column 5: Mgas (Msun)",
        "# column 6: Z (1)",
    ]

    stars_path = os.path.join(args.out_dir, "stars.txt")
    gas_path = os.path.join(args.out_dir, "gas.txt")
    meta_path = os.path.join(args.out_dir, "meta.json")

    _write_table(stars_path, star_header, stars)
    _write_table(gas_path, gas_header, gas)

    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote {stars.shape[0]} stars to {stars_path}")
    print(f"Wrote {gas.shape[0]} gas cells to {gas_path}")
    print(f"Wrote metadata to {meta_path}")


if __name__ == "__main__":
    main()
