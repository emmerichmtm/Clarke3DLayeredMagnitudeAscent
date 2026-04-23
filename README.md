# Layered Magnitude / Hypervolume 3D Runner

This repository contains a single self-contained Python script for running a Clarke-motivated projected ascent method on three-objective benchmark problems.

Current script:

- `layered_magnitude_3d_singlefile_bulged_hv_recovery.py`

## Features

- layered **magnitude** indicator gradients
- layered **hypervolume** indicator gradients
- exact small-front mode
- sweep-based large-front mode
- adjustable **bulge** parameter for the bulged three-peaks benchmark
- progress output during optimization
- step-size **recovery from stagnation**
- PNG and CSV outputs

## Requirements

Python 3 with:

- `numpy`
- `matplotlib`

Install with:

```bash
python -m pip install numpy matplotlib
```

## Quick start

Run the bulged three-peaks benchmark with magnitude gradients:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem bulged_three_peaks \
  --n-points 15 \
  --three-peaks-iters 200 \
  --bulge-gamma 0.25 \
  --indicator magnitude \
  --progress-every 10
```

Run the same benchmark with hypervolume gradients:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem bulged_three_peaks \
  --n-points 15 \
  --three-peaks-iters 200 \
  --bulge-gamma 0.25 \
  --indicator hypervolume \
  --progress-every 10
```

Run the original three-peaks benchmark:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem three_peaks \
  --n-points 15 \
  --three-peaks-iters 200 \
  --indicator magnitude
```

Run crashworthiness:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem vehicle_crashworthiness \
  --n-points 15 \
  --crash-iters 96 \
  --indicator magnitude
```

## Problems

### 1. `three_peaks`

Synthetic three-objective benchmark with peaks at the unit vectors.

### 2. `bulged_three_peaks`

A simplex-constrained benchmark with adjustable bulge parameter `gamma`:

```text
f_i(x) = 1 - (||x - e_i||_2^2 / 2)^gamma,   i = 1,2,3
```

where:

- `e1 = (1,0,0)`
- `e2 = (0,1,0)`
- `e3 = (0,0,1)`

and `x` is projected onto the simplex.

Smaller `gamma` gives a stronger bulge toward `(1,1,1)` in objective space.

### 3. `vehicle_crashworthiness`

Three-objective engineering benchmark with five decision variables and box constraints.

## Indicator option

Use:

- `--indicator magnitude`
- `--indicator hypervolume`

For `magnitude`, the 3D quantity includes:

- edge / axis terms
- face / projected-area terms
- volume / hypervolume term

For `hypervolume`, only the hypervolume part is optimized.

## Exact and switched modes

The code uses two regimes:

- **exact mode** for small first-front sizes
- **sweep-based mode** for larger first fronts

The switching threshold is controlled by:

```bash
--exact-front-threshold 10
```

## Progress output

The script prints intermediate progress such as:

- current iteration
- current objective value
- current step size
- nondominated count
- layer sizes
- mode used
- retry count
- whether the step was accepted or stalled

Control the print frequency with:

```bash
--progress-every 10
```

Silence progress output with:

```bash
--quiet
```

## Step-size recovery

The line search uses a reduction factor of `0.99`.

If the step size becomes too small and the run stagnates for a while, the code attempts a controlled increase of the step size again. This is intended to help the method recover from long stagnation phases near the alpha floor.

## Output files

The script writes files such as:

- `*_objective_space.png`
- `*_decision_space.png` (for three-peaks variants)
- `*_convergence.png`
- `*_initial_decisions.csv`
- `*_final_decisions.csv`
- `*_initial_objectives.csv`
- `*_final_objectives.csv`
- `*_reference_archive.csv`
- `*_history.csv`
- `benchmark_summary.json`

## Examples

Bulged front with stronger bend:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem bulged_three_peaks \
  --bulge-gamma 0.25 \
  --indicator magnitude
```

Hypervolume-driven run:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery.py \
  --problem bulged_three_peaks \
  --bulge-gamma 0.25 \
  --indicator hypervolume
```

## Clarke-motivated

The method is Clarke-motivated in the following practical sense:

- on chambers with fixed layer structure, the code uses exact or sweep-based indicator gradients;
- at layer switches, the hard-layer objective is nonsmooth;
- the implemented method is therefore a practical projected ascent scheme guided by generalized-gradient ideas, rather than a full Clarke-subgradient method with convergence guarantees.

## Notes

- exact 3D computations are intended for small or moderate front sizes
- larger fronts use a switched sweep-based regime
- longer runs can still become expensive, especially for denser fronts

## License

Add your preferred license here.
