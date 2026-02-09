#!/usr/bin/env python3
import argparse
import json
import math
import os
import time

import numpy as np

DEFAULT_RUN_DIR = "run"
DEFAULT_SKI_NAME = "m12i_600.ski"
DEFAULT_NUM_PACKETS = 2e5
DEFAULT_MIN_WAVELENGTH = 0.09
DEFAULT_MAX_WAVELENGTH = 100.0
DEFAULT_NUM_WAVELENGTHS = 200
DEFAULT_DISTANCE_MPC = 10.0
DEFAULT_PIXELS = 256
DEFAULT_MIN_LEVEL = 2
DEFAULT_MAX_LEVEL = 6
DEFAULT_MAX_DUST_FRACTION = 1e-5
DEFAULT_SED_FAMILY = "BruzualCharlotSEDFamily"
DEFAULT_BB_RADIUS_KM = 6.96e5
DEFAULT_BB_TEMP_K = 5000.0


def _load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _write_blackbody_star_file(run_dir, radius_km, temp_k):
    src = os.path.join(run_dir, "stars.txt")
    dst = os.path.join(run_dir, "stars_bb.txt")
    if not os.path.exists(src):
        raise FileNotFoundError(f"Missing {src} required for BlackBodySEDFamily fallback.")

    data = []
    with open(src, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            x, y, z, hsml = map(float, parts[:4])
            data.append([x, y, z, hsml, radius_km, temp_k])

    header = [
        "# column 1: x (kpc)",
        "# column 2: y (kpc)",
        "# column 3: z (kpc)",
        "# column 4: hsml (kpc)",
        "# column 5: R (km)",
        "# column 6: T (K)",
    ]

    with open(dst, "w", encoding="utf-8") as f:
        for line in header:
            f.write(line + "\n")
        np.savetxt(f, np.array(data), fmt="%.8e")

    return "stars_bb.txt"


def main():
    parser = argparse.ArgumentParser(description="Build a SKIRT .ski file for m12i snapshot 600.")
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--ski-name", default=DEFAULT_SKI_NAME)
    parser.add_argument("--num-packets", type=float, default=DEFAULT_NUM_PACKETS)
    parser.add_argument("--min-wavelength", type=float, default=DEFAULT_MIN_WAVELENGTH)
    parser.add_argument("--max-wavelength", type=float, default=DEFAULT_MAX_WAVELENGTH)
    parser.add_argument("--num-wavelengths", type=int, default=DEFAULT_NUM_WAVELENGTHS)
    parser.add_argument("--distance-mpc", type=float, default=DEFAULT_DISTANCE_MPC)
    parser.add_argument("--pixels", type=int, default=DEFAULT_PIXELS)
    parser.add_argument("--min-level", type=int, default=DEFAULT_MIN_LEVEL)
    parser.add_argument("--max-level", type=int, default=DEFAULT_MAX_LEVEL)
    parser.add_argument("--max-dust-fraction", type=float, default=DEFAULT_MAX_DUST_FRACTION)
    parser.add_argument("--sed-family", default=DEFAULT_SED_FAMILY)
    parser.add_argument("--bb-radius-km", type=float, default=DEFAULT_BB_RADIUS_KM)
    parser.add_argument("--bb-temp-k", type=float, default=DEFAULT_BB_TEMP_K)
    args = parser.parse_args()

    os.makedirs(args.run_dir, exist_ok=True)

    meta = _load_json(os.path.join(args.run_dir, "meta.json")) or {}
    views = _load_json(os.path.join(args.run_dir, "views.json")) or {"views": []}

    r_cut = float(meta.get("R_CUT_KPC", 60.0))
    dust_to_metals = float(meta.get("dust_to_metals", 0.4))

    sed_family = args.sed_family
    star_filename = "stars.txt"

    if sed_family == "BlackBodySEDFamily":
        star_filename = _write_blackbody_star_file(args.run_dir, args.bb_radius_km, args.bb_temp_k)

    view_list = views.get("views", [])
    if not view_list:
        raise RuntimeError("No views found in run/views.json; run make_views.py first.")

    instrument_entries = []
    for view in view_list:
        idx = int(view["index"])
        theta = float(view["theta_deg"])
        phi = float(view["phi_deg"])
        instrument_entries.append(
            f"""                    <FrameInstrument instrumentName=\"view_{idx:03d}\" distance=\"{args.distance_mpc} Mpc\" inclination=\"{theta:.6f} deg\" azimuth=\"{phi:.6f} deg\" roll=\"0 deg\" fieldOfViewX=\"{2 * r_cut:.6f} kpc\" fieldOfViewY=\"{2 * r_cut:.6f} kpc\" numPixelsX=\"{args.pixels}\" numPixelsY=\"{args.pixels}\"/>"""
        )

    instrument_block = "\n".join(instrument_entries)

    sed_family_block = ""
    if sed_family == "BruzualCharlotSEDFamily":
        sed_family_block = (
            """                        <sedFamily type=\"SEDFamily\">
                            <BruzualCharlotSEDFamily imf=\"Chabrier\" resolution=\"Low\"/>
                        </sedFamily>"""
        )
    elif sed_family == "BlackBodySEDFamily":
        sed_family_block = (
            """                        <sedFamily type=\"SEDFamily\">
                            <BlackBodySEDFamily/>
                        </sedFamily>"""
        )
    else:
        raise ValueError(f"Unsupported sed family: {sed_family}")

    now = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())

    ski_text = f"""<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<!-- A SKIRT parameter file generated by build_ski.py -->
<skirt-simulation-hierarchy type=\"MonteCarloSimulation\" format=\"9\" producer=\"fire_skirt pipeline\" time=\"{now}\">
    <MonteCarloSimulation userLevel=\"Regular\" simulationMode=\"ExtinctionOnly\" iteratePrimaryEmission=\"false\" iterateSecondaryEmission=\"false\" numPackets=\"{args.num_packets}\">
        <random type=\"Random\">
            <Random seed=\"0\"/>
        </random>
        <units type=\"Units\">
            <ExtragalacticUnits wavelengthOutputStyle=\"Wavelength\" fluxOutputStyle=\"Frequency\"/>
        </units>
        <cosmology type=\"Cosmology\">
            <LocalUniverseCosmology/>
        </cosmology>
        <sourceSystem type=\"SourceSystem\">
            <SourceSystem minWavelength=\"{args.min_wavelength} micron\" maxWavelength=\"{args.max_wavelength} micron\" sourceBias=\"0.5\">
                <sources type=\"Source\">
                    <ParticleSource filename=\"{star_filename}\">
{sed_family_block}
                    </ParticleSource>
                </sources>
            </SourceSystem>
        </sourceSystem>
        <mediumSystem type=\"MediumSystem\">
            <MediumSystem>
                <media type=\"Medium\">
                    <ParticleMedium filename=\"gas.txt\" massFraction=\"{dust_to_metals}\" importMetallicity=\"true\">
                        <materialMix type=\"MaterialMix\">
                            <ThemisDustMix/>
                        </materialMix>
                    </ParticleMedium>
                </media>
                <grid type=\"SpatialGrid\">
                    <PolicyTreeSpatialGrid minX=\"{-r_cut} kpc\" maxX=\"{r_cut} kpc\" minY=\"{-r_cut} kpc\" maxY=\"{r_cut} kpc\" minZ=\"{-r_cut} kpc\" maxZ=\"{r_cut} kpc\">
                        <policy type=\"TreePolicy\">
                            <DensityTreePolicy minLevel=\"{args.min_level}\" maxLevel=\"{args.max_level}\" maxDustFraction=\"{args.max_dust_fraction}\"/>
                        </policy>
                    </PolicyTreeSpatialGrid>
                </grid>
            </MediumSystem>
        </mediumSystem>
        <instrumentSystem type=\"InstrumentSystem\">
            <InstrumentSystem>
                <defaultWavelengthGrid type=\"WavelengthGrid\">
                    <LogWavelengthGrid minWavelength=\"{args.min_wavelength} micron\" maxWavelength=\"{args.max_wavelength} micron\" numWavelengths=\"{args.num_wavelengths}\"/>
                </defaultWavelengthGrid>
                <instruments type=\"Instrument\">
{instrument_block}
                </instruments>
            </InstrumentSystem>
        </instrumentSystem>
        <probeSystem type=\"ProbeSystem\">
            <ProbeSystem/>
        </probeSystem>
    </MonteCarloSimulation>
</skirt-simulation-hierarchy>
"""

    out_path = os.path.join(args.run_dir, args.ski_name)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(ski_text)

    print(f"Wrote {out_path}")
    if sed_family == "BlackBodySEDFamily":
        print("BlackBodySEDFamily selected; wrote stars_bb.txt")


if __name__ == "__main__":
    main()
