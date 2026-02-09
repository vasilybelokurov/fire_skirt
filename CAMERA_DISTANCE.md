# Camera distance in this pipeline

This document explains where the camera distance is set and how to change it.

## Default distance

The SKIRT instruments are created in `code/build_ski.py`.

The default camera distance is:

```
DEFAULT_DISTANCE_MPC = 10.0
```

This value is used when generating the instrument entries, e.g.:

```xml
<FrameInstrument instrumentName="view_015" distance="10 Mpc" ... />
```

So unless you override it, the camera is placed at **10 Mpc** along each view direction.

## How to change the distance

You can override the default by passing `--distance-mpc` when building the `.ski`:

```bash
python code/build_ski.py --distance-mpc 20
```

This will write:

```xml
distance="20 Mpc"
```

into every instrument entry in `run/m12i_600.ski`.

## Pipeline note

If you run the full pipeline with `code/run_skirt.sh`, it currently uses the default distance (10 Mpc). If you want distance control in the one‑shot pipeline, add a passthrough to `code/run_skirt.sh` and re‑run `build_ski.py` with `--distance-mpc`.

## Camera positions file

If you use `code/camera_positions.py`, it also takes the camera distance as input and must match the value used in `build_ski.py`:

```bash
python code/camera_positions.py --views run/views.json --distance-mpc 10
```

