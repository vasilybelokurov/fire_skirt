#!/usr/bin/env python3
import argparse
import json
import os


def _load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="Compute camera positions from views.json.")
    parser.add_argument("--views", default="run/views.json", help="Path to views.json")
    parser.add_argument("--distance-mpc", type=float, default=10.0, help="Camera distance in Mpc")
    parser.add_argument(
        "--out",
        default=None,
        help="Output text file. Defaults to <views_dir>/camera_positions.txt",
    )
    args = parser.parse_args()

    data = _load_json(args.views)
    views = data.get("views", [])
    if not views:
        raise SystemExit(f"No views found in {args.views}")

    out_path = args.out
    if out_path is None:
        out_path = os.path.join(os.path.dirname(args.views), "camera_positions.txt")

    d_mpc = float(args.distance_mpc)
    d_kpc = d_mpc * 1000.0

    header = (
        "# index theta_deg phi_deg  dx dy dz  cam_x_Mpc cam_y_Mpc cam_z_Mpc  "
        "cam_x_kpc cam_y_kpc cam_z_kpc"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for v in views:
            idx = int(v["index"])
            th = float(v["theta_deg"])
            ph = float(v["phi_deg"])
            dx, dy, dz = map(float, v["dir"])
            cam_mpc = (d_mpc * dx, d_mpc * dy, d_mpc * dz)
            cam_kpc = (d_kpc * dx, d_kpc * dy, d_kpc * dz)
            f.write(
                f"{idx:2d} {th:12.6f} {ph:12.6f} {dx: .8f} {dy: .8f} {dz: .8f} "
                f"{cam_mpc[0]: .8f} {cam_mpc[1]: .8f} {cam_mpc[2]: .8f} "
                f"{cam_kpc[0]: .8f} {cam_kpc[1]: .8f} {cam_kpc[2]: .8f}\n"
            )

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
