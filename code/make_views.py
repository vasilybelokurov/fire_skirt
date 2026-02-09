#!/usr/bin/env python3
import argparse
import json
import math
import os
import time

import numpy as np

DEFAULT_OUT_DIR = "run"
DEFAULT_NUM_VIEWS = 32


def fibonacci_sphere(n):
    golden_ratio = (1.0 + 5.0 ** 0.5) / 2.0
    i = np.arange(n, dtype=np.float64)
    phi = 2.0 * math.pi * i / golden_ratio
    cos_theta = 1.0 - 2.0 * (i + 0.5) / n
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    return x, y, z


def main():
    parser = argparse.ArgumentParser(description="Generate Fibonacci sphere view directions.")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--num-views", type=int, default=DEFAULT_NUM_VIEWS)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    x, y, z = fibonacci_sphere(args.num_views)
    views = []
    for idx in range(args.num_views):
        theta = math.degrees(math.acos(max(-1.0, min(1.0, z[idx]))))
        phi = math.degrees(math.atan2(y[idx], x[idx]))
        views.append(
            {
                "index": idx,
                "theta_deg": theta,
                "phi_deg": phi,
                "dir": [float(x[idx]), float(y[idx]), float(z[idx])],
            }
        )

    payload = {
        "num_views": args.num_views,
        "method": "fibonacci_sphere",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "views": views,
    }

    out_path = os.path.join(args.out_dir, "views.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Wrote {args.num_views} views to {out_path}")


if __name__ == "__main__":
    main()
