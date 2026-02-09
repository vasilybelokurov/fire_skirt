#!/usr/bin/env python3
import argparse
import glob
import os

import numpy as np
from astropy.io import fits
import imageio.v2 as imageio


PRESETS = {
    "hst": {
        "B": (0.41, 0.49),
        "G": (0.54, 0.70),
        "R": (0.75, 0.90),
    },
    "jwst": {
        "B": (1.00, 1.30),
        "G": (1.80, 2.20),
        "R": (3.20, 3.90),
    },
}


def load_cube(path):
    with fits.open(path, memmap=False) as hdul:
        cube = np.asarray(hdul[0].data, dtype=np.float64)
        if cube.ndim != 3:
            raise RuntimeError(f"Expected 3D cube in {path}, got {cube.shape}")
        # Wavelength axis stored in extension 1 as a vector
        if len(hdul) > 1 and hdul[1].data is not None:
            waves = np.asarray(hdul[1].data, dtype=np.float64).reshape(-1)
        else:
            # fallback: try header WCS (not always present)
            hdr = hdul[0].header
            nlam = cube.shape[0]
            crval = hdr.get("CRVAL3", 0.0)
            cdelt = hdr.get("CDELT3", 1.0)
            crpix = hdr.get("CRPIX3", 1.0)
            waves = crval + (np.arange(nlam) + 1 - crpix) * cdelt
        return cube, waves


def band_image(cube, waves, wmin, wmax):
    mask = (waves >= wmin) & (waves <= wmax)
    if not np.any(mask):
        raise RuntimeError(f"No wavelengths in range {wmin}-{wmax} micron")
    return np.nansum(cube[mask, :, :], axis=0)


def lupton_rgb(r, g, b, Q=10.0, stretch=0.1):
    r = np.maximum(r, 0.0)
    g = np.maximum(g, 0.0)
    b = np.maximum(b, 0.0)
    I = (r + g + b) / 3.0
    I = np.maximum(I, 1e-12)
    f = np.arcsinh(Q * (I / stretch)) / Q
    R = f * (r / I)
    G = f * (g / I)
    B = f * (b / I)
    rgb = np.stack([R, G, B], axis=-1)
    rgb = np.clip(rgb, 0.0, None)
    return rgb


def normalize_rgb(rgb, p=99.5):
    vals = rgb[np.isfinite(rgb)]
    if vals.size == 0:
        return rgb
    vmax = np.percentile(vals, p)
    if vmax > 0:
        rgb = rgb / vmax
    return np.clip(rgb, 0.0, 1.0)


def main():
    parser = argparse.ArgumentParser(description="Make false-color RGB images from SKIRT cubes.")
    parser.add_argument("--input-dir", default="run")
    parser.add_argument("--glob", default="*view_*_total.fits")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--preset", choices=sorted(PRESETS.keys()), default="hst")
    parser.add_argument("--band-b", type=float, nargs=2, metavar=("BMIN", "BMAX"))
    parser.add_argument("--band-g", type=float, nargs=2, metavar=("GMIN", "GMAX"))
    parser.add_argument("--band-r", type=float, nargs=2, metavar=("RMIN", "RMAX"))
    parser.add_argument("--Q", type=float, default=10.0)
    parser.add_argument("--stretch", type=float, default=0.1)
    parser.add_argument("--percentile", type=float, default=99.5)
    parser.add_argument("--suffix", default=None)
    args = parser.parse_args()

    input_dir = os.path.expanduser(args.input_dir)
    pattern = os.path.join(input_dir, args.glob)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No FITS files found for pattern: {pattern}")

    out_dir = args.out_dir or os.path.join(input_dir, "rgb")
    os.makedirs(out_dir, exist_ok=True)

    bands = PRESETS[args.preset].copy()
    if args.band_b:
        bands["B"] = tuple(args.band_b)
    if args.band_g:
        bands["G"] = tuple(args.band_g)
    if args.band_r:
        bands["R"] = tuple(args.band_r)

    suffix = args.suffix or args.preset

    for path in paths:
        cube, waves = load_cube(path)
        b = band_image(cube, waves, *bands["B"])
        g = band_image(cube, waves, *bands["G"])
        r = band_image(cube, waves, *bands["R"])

        rgb = lupton_rgb(r, g, b, Q=args.Q, stretch=args.stretch)
        rgb = normalize_rgb(rgb, p=args.percentile)

        base = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(out_dir, f"{base}_rgb_{suffix}.png")
        imageio.imwrite(out_path, (rgb * 255).astype(np.uint8))
        print(out_path)


if __name__ == "__main__":
    main()
