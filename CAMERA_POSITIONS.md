# Camera position computation

This document describes how camera positions are computed from the view directions in `views.json`.

## Inputs

- `views.json` contains one entry per view with:
  - `theta_deg`: polar angle in degrees (0 = +z, 90 = in the x–y plane)
  - `phi_deg`: azimuthal angle in degrees (in the x–y plane)
  - `dir`: unit direction vector `[dx, dy, dz]`

The pipeline generates these with a Fibonacci sphere so directions are approximately uniform on the sphere.

## Definition

We place each camera at a fixed distance `D` from the origin (galaxy center), along the unit direction vector.

Let:
- `D` be the distance in Mpc (default `10.0` in the SKIRT setup).
- `dir = (dx, dy, dz)` be the unit view direction.

Then the camera position in Mpc is:

```
cam = D * dir = (D*dx, D*dy, D*dz)
```

For convenience, the same position in kpc is:

```
cam_kpc = 1000 * cam
```

## Output format

The script writes a plain text file with one row per view:

```
# index theta_deg phi_deg  dx dy dz  cam_x_Mpc cam_y_Mpc cam_z_Mpc  cam_x_kpc cam_y_kpc cam_z_kpc
 0    14.361512     0.000000  0.24803919  0.00000000  0.96875000  2.48039185  0.00000000  9.68750000  2480.39185412  0.00000000  9687.50000000
 ...
```

## Script

Use the provided script:

```bash
python code/camera_positions.py --views run/views.json --distance-mpc 10
```

By default it writes `camera_positions.txt` next to `views.json`. You can override the output file with:

```bash
python code/camera_positions.py --views run/views.json --distance-mpc 10 --out outputs/camera_positions.txt
```

## Notes

- The camera distance `D` should match the `distance` used in the SKIRT instrument definition.
- The camera positions here are purely geometric; they do not include any roll angle or projection.
