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


def _read_wavelengths(hdul):
    if len(hdul) > 1 and hdul[1].data is not None:
        w = hdul[1].data
        if hasattr(w, "dtype") and w.dtype.fields:
            key = list(w.dtype.fields.keys())[0]
            return np.asarray(w[key], dtype=np.float64).reshape(-1)
        return np.asarray(w, dtype=np.float64).reshape(-1)
    hdr = hdul[0].header
    nlam = hdul[0].data.shape[0]
    crval = hdr.get("CRVAL3", 0.0)
    cdelt = hdr.get("CDELT3", 1.0)
    crpix = hdr.get("CRPIX3", 1.0)
    return crval + (np.arange(nlam) + 1 - crpix) * cdelt


def load_image(path, wmin=None, wmax=None):
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data
        if data is None:
            raise RuntimeError(f"No data in {path}")
        arr = np.asarray(data, dtype=np.float64)

        if arr.ndim < 2:
            raise RuntimeError(f"Unexpected data shape {arr.shape} in {path}")

        if arr.ndim == 2:
            return arr

        waves = _read_wavelengths(hdul)
        if wmin is not None or wmax is not None:
            if wmin is None:
                wmin = float(np.min(waves))
            if wmax is None:
                wmax = float(np.max(waves))
            mask = (waves >= wmin) & (waves <= wmax)
            if not np.any(mask):
                raise RuntimeError(f"No wavelengths in range {wmin}-{wmax} micron in {path}")
            return np.nansum(arr[mask, :, :], axis=0)

        ny, nx = arr.shape[-2], arr.shape[-1]
        return np.nansum(arr.reshape(-1, ny, nx), axis=0)


def _smooth_image(img, fwhm_pix):
    if fwhm_pix <= 0:
        return img
    sigma = fwhm_pix / 2.355
    kernel = Gaussian2DKernel(sigma)
    return convolve_fft(img, kernel, boundary="fill", nan_treatment="fill", fill_value=0.0)


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


def _replace_zeros(img, percentile):
    positive = img > 0
    if not np.any(positive):
        return img
    floor = np.percentile(img[positive], percentile)
    out = img.copy()
    out[out == 0] = floor
    return out


def scale_image(img, stretch_mode, pmin, pmax, Q, stretch, zero_fill, zero_percentile):
    img = np.array(img, dtype=np.float64)
    img = np.nan_to_num(img, nan=0.0, posinf=0.0, neginf=0.0)
    img = np.maximum(img, 0.0)
    if zero_fill:
        img = _replace_zeros(img, zero_percentile)
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
    parser.add_argument("--resample", choices=["lanczos", "bilinear", "bicubic", "nearest"], default="lanczos")
    parser.add_argument("--zero-fill", action="store_true", help="Replace zero pixels with low-percentile floor")
    parser.add_argument("--zero-percentile", type=float, default=1.0, help="Percentile for zero-fill floor")
    parser.add_argument("--wmin", type=float, default=None, help="Min wavelength (micron) for band integration")
    parser.add_argument("--wmax", type=float, default=None, help="Max wavelength (micron) for band integration")
    args = parser.parse_args()

    input_dir = os.path.expanduser(args.input_dir)
    pattern = os.path.join(input_dir, args.glob)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No FITS files found for pattern: {pattern}")

    out_dir = args.out_dir or os.path.join(input_dir, "images")
    os.makedirs(out_dir, exist_ok=True)

    for path in paths:
        img = load_image(path, wmin=args.wmin, wmax=args.wmax)
        img = _smooth_image(img, args.fwhm_pix)
        img = _center_crop(img, args.crop_frac)
        scaled, vmin, vmax = scale_image(
            img, args.stretch_mode, args.pmin, args.pmax, args.Q, args.stretch, args.zero_fill, args.zero_percentile
        )

        base = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(out_dir, f"{base}.{args.format}")

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        cmap = cm.get_cmap(args.cmap)
        rgba = cmap(norm(scaled), bytes=True)
        img = Image.fromarray(rgba, mode="RGBA")
        if args.size_px:
            resample_map = {
                "lanczos": Image.LANCZOS,
                "bilinear": Image.BILINEAR,
                "bicubic": Image.BICUBIC,
                "nearest": Image.NEAREST,
            }
            img = img.resize((args.size_px, args.size_px), resample=resample_map[args.resample])
        if args.format == "jpg":
            img = img.convert("RGB")
        img.save(out_path)
        print(out_path)


if __name__ == "__main__":
    main()
