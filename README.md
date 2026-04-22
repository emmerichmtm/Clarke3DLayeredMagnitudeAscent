# Layered Magnitude 3D Single-File Runner

This repository contains a single self-contained Python script for running a Clarke-motivated layered-magnitude ascent method on two three-objective benchmark problems:

- `three_peaks` — a synthetic 3-objective warm-up problem
- `vehicle_crashworthiness` — a 3-objective engineering benchmark

The script uses an exact 3D magnitude / hypervolume engine and can write PNG and CSV files for the benchmark runs.

## File

- `layered_magnitude_3d_singlefile.py`

## Requirements

Python 3 with:

- `numpy`
- `matplotlib`

Install them with:

```bash
python -m pip install numpy matplotlib
```

## Quick start

Run both examples:

```bash
python layered_magnitude_3d_singlefile.py
```

Run only the three-peaks example:

```bash
python layered_magnitude_3d_singlefile.py --problem three_peaks
```

Run only the crashworthiness example:

```bash
python layered_magnitude_3d_singlefile.py --problem vehicle_crashworthiness
```

## Useful options

Set the output directory:

```bash
python layered_magnitude_3d_singlefile.py --outdir results
```

Change the number of points:

```bash
python layered_magnitude_3d_singlefile.py --n-points 15
```

Set the iteration budgets:

```bash
python layered_magnitude_3d_singlefile.py --three-peaks-iters 48 --crash-iters 24
```

Example:

```bash
python layered_magnitude_3d_singlefile.py --problem both --outdir results --three-peaks-iters 48 --crash-iters 24
```

## What the script prints

For each benchmark, the script reports:

- iteration count
- initial layered magnitude
- final layered magnitude
- accepted steps
- final step size
- nondominated count change
- IGD change

This makes it easy to monitor progress from the terminal without opening the output files first.

## Output files

The script writes files such as:

### Three-peaks
- `three_peaks_objective_space.png`
- `three_peaks_decision_space.png`
- `three_peaks_convergence.png`
- `three_peaks_initial_decisions.csv`
- `three_peaks_final_decisions.csv`
- `three_peaks_initial_objectives.csv`
- `three_peaks_final_objectives.csv`
- `three_peaks_reference_archive.csv`
- `three_peaks_history.csv`

### Vehicle crashworthiness
- `vehicle_crashworthiness_objective_space.png`
- `vehicle_crashworthiness_convergence.png`
- `vehicle_crashworthiness_initial_decisions.csv`
- `vehicle_crashworthiness_final_decisions.csv`
- `vehicle_crashworthiness_initial_objectives.csv`
- `vehicle_crashworthiness_final_objectives.csv`
- `vehicle_crashworthiness_initial_objectives_raw.csv`
- `vehicle_crashworthiness_final_objectives_raw.csv`
- `vehicle_crashworthiness_reference_archive.csv`
- `vehicle_crashworthiness_history.csv`

The script also writes:

- `benchmark_summary.json`

## Benchmarks

### 1. Three-peaks problem

The synthetic warm-up benchmark is defined by

    f_i(x) = 1 - ||x - e_i||_2,    i = 1,2,3

where

- `e1 = (1,0,0)`
- `e2 = (0,1,0)`
- `e3 = (0,0,1)`

and

- `x in [-2,2]^3`

The approximation points are initialized near the origin.

### 2. Vehicle crashworthiness

The crashworthiness benchmark is a standard 3-objective engineering problem with 5 decision variables and box constraints. The script transforms the original minimization objectives into a bounded maximization space before applying the layered-magnitude ascent method.

## Clarke-motivated: what it means here

The method is Clarke-motivated in a practical sense:

- on chambers with fixed layer structure, the exact magnitude engine provides an objective-space tie-shared subgradient;
- at layer switches, the hard-layer objective is nonsmooth;
- the implemented method should therefore be interpreted as a practical ascent scheme guided by generalized-gradient ideas, rather than as a full Clarke-subgradient method with convergence guarantees.

## Notes

- The exact 3D magnitude engine is intended for small to medium-sized point sets.
- Longer runs can become computationally expensive because exact 3D inclusion-exclusion is costly.
- The default settings are chosen to be useful while still remaining reasonably runnable.

## Suggested citation text

If you use this code in connection with the manuscript, cite the corresponding paper and mention that the repository contains a single-file reference implementation of the 3-objective layered-magnitude ascent method.

## License

Add your preferred license here, for example:

- MIT
- BSD-3-Clause
- GPL-3.0-only
