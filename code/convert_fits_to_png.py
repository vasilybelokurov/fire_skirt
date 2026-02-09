#!/usr/bin/env python3
import argparse
import glob
import os

import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve_fft
import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from PIL import Image


DEFAULT_GLOB = "*view_*_total.fits"


def load_image(path):
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data

    if data is None:
        raise RuntimeError(f"No data in {path}")

    arr = np.asarray(data, dtype=np.float64)
    if arr.ndim < 2:
        raise RuntimeError(f"Unexpected data shape {arr.shape} in {path}")

    # Collapse all leading axes; assume last two are spatial
    if arr.ndim == 2:
        img = arr
    else:
        ny, nx = arr.shape[-2], arr.shape[-1]
        img = np.nansum(arr.reshape(-1, ny, nx), axis=0)

    return img


def _smooth_image(img, fwhm_pix):
    if fwhm_pix <= 0:
        return img
    sigma = fwhm_pix / 2.355
    kernel = Gaussian2DKernel(sigma)
    return convolve_fft(img, kernel, boundary="wrap", nan_treatment="fill", fill_value=0.0)


def _center_crop(img, frac):
    if frac >= 1.0:
        return img
    if frac <= 0.0:
        raise RuntimeError("crop fraction must be > 0")
    ny, nx = img.shape
    new_ny = max(1, int(round(ny * frac)))
    new_nx = max(1, int(round(nx * frac)))
    y0 = (ny - new_ny) // 2
    x0 = (nx - new_nx) // 2
    return img[y0 : y0 + new_ny, x0 : x0 + new_nx]


def _apply_stretch(img, mode, Q, stretch):
    if mode == "linear":
        return img
    if mode == "log":
        positive = img > 0
        if not np.any(positive):
            raise RuntimeError("No positive values for log scaling")
        min_pos = np.percentile(img[positive], 1.0)
        return np.log10(np.clip(img, min_pos, None))
    if mode == "asinh":
        return np.arcsinh(Q * (img / stretch)) / Q
    raise RuntimeError(f"Unknown stretch mode: {mode}")


def scale_image(img, stretch_mode, pmin, pmax, Q, stretch):
    img = np.array(img, dtype=np.float64)
    finite = np.isfinite(img)

    if not np.any(finite):
        raise RuntimeError("No finite values in image")

    img = _apply_stretch(img, stretch_mode, Q, stretch)
    finite = np.isfinite(img)

    vals = img[finite]
    vmin = np.percentile(vals, pmin)
    vmax = np.percentile(vals, pmax)

    return img, vmin, vmax


def main():
    parser = argparse.ArgumentParser(description="Convert SKIRT FITS images to PNG/JPG.")
    parser.add_argument("--input-dir", default="run")
    parser.add_argument("--glob", default=DEFAULT_GLOB)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--format", choices=["png", "jpg"], default="png")
    parser.add_argument("--stretch-mode", choices=["log", "linear", "asinh"], default="asinh")
    parser.add_argument("--Q", type=float, default=5.0)
    parser.add_argument("--stretch", type=float, default=0.4)
    parser.add_argument("--pmin", type=float, default=1.0)
    parser.add_argument("--pmax", type=float, default=99.0)
    parser.add_argument("--cmap", default="magma")
    parser.add_argument("--fwhm-pix", type=float, default=4.0)
    parser.add_argument("--crop-frac", type=float, default=1.0)
    parser.add_argument("--size-px", type=int, default=None)
    args = parser.parse_args()

    input_dir = os.path.expanduser(args.input_dir)
    pattern = os.path.join(input_dir, args.glob)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No FITS files found for pattern: {pattern}")

    out_dir = args.out_dir or os.path.join(input_dir, "images")
    os.makedirs(out_dir, exist_ok=True)

    for path in paths:
        img = load_image(path)
        img = _smooth_image(img, args.fwhm_pix)
        img = _center_crop(img, args.crop_frac)
        scaled, vmin, vmax = scale_image(img, args.stretch_mode, args.pmin, args.pmax, args.Q, args.stretch)

        base = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(out_dir, f"{base}.{args.format}")

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        cmap = cm.get_cmap(args.cmap)
        rgba = cmap(norm(scaled), bytes=True)
        img = Image.fromarray(rgba, mode="RGBA")
        if args.size_px:
            img = img.resize((args.size_px, args.size_px), resample=Image.LANCZOS)
        if args.format == "jpg":
            img = img.convert("RGB")
        img.save(out_path)
        print(out_path)


if __name__ == "__main__":
    main()
