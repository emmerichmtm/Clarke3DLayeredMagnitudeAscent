# Layered Magnitude / Hypervolume 3D Runner with Stochastic Moves

This repository contains a single self-contained Python script for running a Clarke-motivated projected ascent method on three-objective benchmark problems.

Current script:

- `layered_magnitude_3d_singlefile_bulged_hv_recovery_names_dd_stochastic.py`

## Features

- layered **magnitude** indicator gradients
- layered **hypervolume** indicator gradients
- **gradient** moves
- **stochastic hillclimbing** moves
- exact small-front mode
- sweep-based large-front mode
- adjustable **bulge** parameter for the bulged three-peaks benchmark
- progress output during optimization
- step-size **recovery from stagnation**
- exact **Das–Dennis simplex-grid initialization**
- compact output filenames with short suffixes for non-default settings
- PNG and CSV outputs

## Requirements

Python 3 with:

- `numpy`
- `matplotlib`

Install with:

```bash
python -m pip install numpy matplotlib
```

## Main script

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery_names_dd_stochastic.py
```

## Problems

Use one of:

- `--problem three_peaks`
- `--problem bulged_three_peaks`
- `--problem vehicle_crashworthiness`
- `--problem both`

## Indicator option

Choose:

- `--indicator magnitude`
- `--indicator hypervolume`

For `magnitude`, the 3D quantity includes:

- edge / axis terms
- face / projected-area terms
- volume / hypervolume term

For `hypervolume`, only the hypervolume part is optimized.

## Move option

Choose:

- `--move gradient`
- `--move stochastic`

### Gradient

Uses the projected gradient-like move based on the selected indicator.

### Stochastic hillclimbing

Perturbs **one point at a time** and accepts the move if the selected layered indicator improves.
The same adaptive step-size logic is used:

- shrink factor `0.99`
- alpha floor
- recovery from stagnation after repeated stalls

This makes stochastic hillclimbing a lightweight alternative to gradient moves.

## Bulged three-peaks benchmark

The bulged simplex benchmark uses

```text
f_i(x) = 1 - (||x - e_i||_2^2 / 2)^gamma,   i = 1,2,3
```

where:

- `e1 = (1,0,0)`
- `e2 = (0,1,0)`
- `e3 = (0,0,1)`

and `x` is projected onto the simplex.

Control the bulge with:

```bash
--bulge-gamma 0.25
```

Smaller `gamma` gives a stronger bend toward `(1,1,1)` in objective space.

## Initialization

Choose:

- `--initialization random`
- `--initialization dasdenis`

### Das–Dennis initialization

For simplex-based problems, the script supports an exact Das–Dennis grid with a small perturbation and reprojection to the simplex.

Control it with:

```bash
--initialization dasdenis
--dd-h 3
--dd-sigma 0.01
```

For 3 objectives, the exact Das–Dennis cardinality is

```text
(H + 1)(H + 2) / 2
```

Examples:

- `H=1` gives `3` points
- `H=2` gives `6` points
- `H=3` gives `10` points
- `H=4` gives `15` points

## Quick examples

Bulged three-peaks with magnitude gradients:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery_names_dd_stochastic.py \
  --problem bulged_three_peaks \
  --initialization dasdenis \
  --dd-h 3 \
  --dd-sigma 0.01 \
  --bulge-gamma 0.25 \
  --indicator magnitude \
  --move gradient \
  --three-peaks-iters 200 \
  --progress-every 10
```

Bulged three-peaks with hypervolume-based stochastic hillclimbing:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery_names_dd_stochastic.py \
  --problem bulged_three_peaks \
  --initialization dasdenis \
  --dd-h 3 \
  --dd-sigma 0.01 \
  --bulge-gamma 0.25 \
  --indicator hypervolume \
  --move stochastic \
  --three-peaks-iters 200 \
  --progress-every 10
```

Crashworthiness with stochastic magnitude hillclimbing:

```bash
python layered_magnitude_3d_singlefile_bulged_hv_recovery_names_dd_stochastic.py \
  --problem vehicle_crashworthiness \
  --n-points 15 \
  --crash-iters 96 \
  --indicator magnitude \
  --move stochastic
```

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
- whether the step was accepted, stalled, or recovered
- move type

Control the print frequency with:

```bash
--progress-every 10
```

Silence progress output with:

```bash
--quiet
```

## Step-size adaptation and recovery

The line search uses a reduction factor of `0.99`.

If the step size becomes too small and the run stagnates for a while, the code attempts a controlled increase of the step size again. This is intended to help the method recover from long stagnation phases near the alpha floor.

## Output filenames

Output filenames automatically include a compact suffix for **non-default settings only**.

Two-letter tags are used, for example:

- `se` = seed
- `np` = n-points
- `tp` = three-peaks-iters
- `cr` = crash-iters
- `ex` = exact-front-threshold
- `bu` = bulge-gamma
- `in` = indicator
- `ii` = initialization
- `dh` = dd-h
- `ds` = dd-sigma
- `mv` = move

So a run with non-default settings may produce names like:

```text
bulged_three_peaks_bu0p25_inhypervolume_iidasdenis_dh3_ds0p01_mvstochastic_objective_space.png
```

Default-valued settings are omitted from the suffix.

## Output files

The script writes files such as:

- `*_objective_space.png`
- `*_decision_space.png` (for simplex / three-peaks variants)
- `*_convergence.png`
- `*_initial_decisions.csv`
- `*_final_decisions.csv`
- `*_initial_objectives.csv`
- `*_final_objectives.csv`
- `*_reference_archive.csv`
- `*_history.csv`
- `benchmark_summary.json`

## Clarke-motivated

The method is Clarke-motivated in the following practical sense:

- on chambers with fixed layer structure, the code uses exact or sweep-based indicator gradients;
- at layer switches, the hard-layer objective is nonsmooth;
- the implemented method is therefore a practical projected ascent scheme guided by generalized-gradient ideas, rather than a full Clarke-subgradient method with convergence guarantees.

## Notes

- exact 3D computations are intended for small or moderate front sizes
- larger fronts use a switched sweep-based regime
- longer runs can still become expensive, especially for denser fronts
- stochastic hillclimbing can be slower to improve than gradient moves, but it is simple and robust

## License

Add your preferred license here.
