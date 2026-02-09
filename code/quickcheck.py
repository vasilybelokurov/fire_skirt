#!/usr/bin/env python3
import argparse
import glob
import json
import os
import random

import numpy as np

DEFAULT_RUN_DIR = "run"


def _load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _stream_stats(path, usecols, sample_max=200000):
    mins = [float("inf")] * len(usecols)
    maxs = [float("-inf")] * len(usecols)
    sample = [[] for _ in usecols]
    count = 0
    rng = random.Random(0)

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            vals = np.fromstring(line, sep=" ")
            if vals.size <= max(usecols):
                continue
            count += 1
            for i, col in enumerate(usecols):
                v = float(vals[col])
                if v < mins[i]:
                    mins[i] = v
                if v > maxs[i]:
                    maxs[i] = v
                if len(sample[i]) < sample_max:
                    sample[i].append(v)
                else:
                    j = rng.randint(0, count - 1)
                    if j < sample_max:
                        sample[i][j] = v

    medians = [float(np.median(s)) if s else float("nan") for s in sample]
    return {
        "count": count,
        "mins": mins,
        "maxs": maxs,
        "medians": medians,
    }


def _check_images(output_dir, views):
    missing = []
    for view in views:
        idx = int(view["index"])
        pattern = os.path.join(output_dir, f"*view_{idx:03d}_*.fits")
        matches = glob.glob(pattern)
        if not matches:
            missing.append(idx)
    return missing


def main():
    parser = argparse.ArgumentParser(description="Quick sanity checks for SKIRT pipeline outputs.")
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--skip-images", action="store_true")
    parser.add_argument("--image-dir", default=None)
    args = parser.parse_args()

    run_dir = args.run_dir
    stars_path = os.path.join(run_dir, "stars.txt")
    gas_path = os.path.join(run_dir, "gas.txt")
    views_path = os.path.join(run_dir, "views.json")

    if not os.path.exists(stars_path):
        raise FileNotFoundError(f"Missing {stars_path}")
    if not os.path.exists(gas_path):
        raise FileNotFoundError(f"Missing {gas_path}")

    stars_stats = _stream_stats(stars_path, usecols=[0, 1, 2, 5])
    gas_stats = _stream_stats(gas_path, usecols=[0, 1, 2, 5])

    print("Stars:")
    print(f"  count: {stars_stats['count']}")
    print(f"  coord medians (kpc): {stars_stats['medians'][:3]}")
    print(f"  Z min/max: {stars_stats['mins'][3]} / {stars_stats['maxs'][3]}")

    print("Gas:")
    print(f"  count: {gas_stats['count']}")
    print(f"  coord medians (kpc): {gas_stats['medians'][:3]}")
    print(f"  Z min/max: {gas_stats['mins'][3]} / {gas_stats['maxs'][3]}")

    def _warn(cond, msg):
        if cond:
            print(f"WARNING: {msg}")

    _warn(stars_stats["count"] == 0, "stars.txt is empty")
    _warn(gas_stats["count"] == 0, "gas.txt is empty")

    _warn(stars_stats["mins"][3] < 0 or stars_stats["maxs"][3] > 0.1, "Star metallicity out of [0, 0.1]")
    _warn(gas_stats["mins"][3] < 0 or gas_stats["maxs"][3] > 0.1, "Gas metallicity out of [0, 0.1]")

    for axis, val in zip(["x", "y", "z"], stars_stats["medians"][:3]):
        _warn(abs(val) > 1.0, f"Star {axis}-median {val:.3f} kpc not centered")
    for axis, val in zip(["x", "y", "z"], gas_stats["medians"][:3]):
        _warn(abs(val) > 1.0, f"Gas {axis}-median {val:.3f} kpc not centered")

    if not args.skip_images:
        views = _load_json(views_path)
        if not views:
            raise FileNotFoundError(f"Missing {views_path}")
        output_dir = args.image_dir or run_dir
        missing = _check_images(output_dir, views.get("views", []))
        if missing:
            print(f"Missing images for {len(missing)} views: {missing[:10]}")
            raise SystemExit(1)
        print(f"Found images for all {views.get('num_views', len(views.get('views', [])))} views")


if __name__ == "__main__":
    main()
