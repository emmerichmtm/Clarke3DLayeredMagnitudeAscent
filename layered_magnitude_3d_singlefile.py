#!/usr/bin/env python3
"""Standalone 3D layered-magnitude benchmark runner.

This single file combines:
- exact 2D/3D magnitude and hypervolume computation,
- the 3D layered ascent method,
- the three-peaks and vehicle crashworthiness examples,
- a simple CLI that prints iterations and achieved layered magnitude.
"""

from __future__ import annotations


"""
exact_magnitude_hypervolume.py

Exact 2D/3D hypervolume and magnitude for maximization with a given anchor point,
together with tie-shared inclusion-exclusion subgradients.

The implementation is intended for small and medium-sized point sets when one
wants exact values and exact piecewise-linear / piecewise-polynomial gradient
information for research experiments.

Definitions
-----------
For a finite approximation set P = {p^(1), ..., p^(n)} in R^d with anchor point a
and maximization, define the translated points

    q^(i) = p^(i) - a,

assuming q^(i) >= 0 componentwise.

The dominated region is the union of anchored boxes
    D(P) = union_i [0, q^(i)].

Hypervolume:
    HV_d(P; a) = vol_d(D(P))

computed exactly by inclusion-exclusion.

Magnitude:
For d = 2:
    Mag_2(P; a) = 1 + 1/2 (L_x + L_y) + 1/4 HV_2(P; a)

For d = 3:
    Mag_3(P; a) = 1
                  + 1/2 (L_x + L_y + L_z)
                  + 1/4 (A_xy + A_xz + A_yz)
                  + 1/8 HV_3(P; a)

where L_* are axis projection lengths, and A_xy etc. are 2D dominated areas of
the coordinate projections. These are the dominated-set formulas used in the
report in the l1 box setting.

Subgradients:
The exact inclusion-exclusion gradient is not unique at ties. We return a
tie-shared subgradient: when multiple points attain a coordinatewise minimum in a
subset term, that term's derivative is split equally among the tied points.
Likewise, maxima in the axis-projection terms are shared equally among ties.

Author: OpenAI assistant, prepared for Michael Emmerich.
License: MIT-style; feel free to adapt.
"""

from itertools import combinations
from typing import Callable, Iterable, Optional, Sequence, Tuple

import numpy as np


ArrayLike = Sequence[Sequence[float]]


def _as_array(points: ArrayLike, dim: Optional[int] = None) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2:
        raise ValueError("points must be a 2D array-like object of shape (n_points, dimension)")
    if dim is not None and arr.shape[1] != dim:
        raise ValueError(f"expected points in dimension {dim}, got {arr.shape[1]}")
    return arr


def _translate_and_validate(points: ArrayLike, anchor: Sequence[float], dim: Optional[int] = None) -> np.ndarray:
    pts = _as_array(points, dim=dim)
    a = np.asarray(anchor, dtype=float)
    if a.ndim != 1 or a.shape[0] != pts.shape[1]:
        raise ValueError("anchor must be a 1D array-like object with the same dimension as the points")
    q = pts - a
    if np.any(q < -1e-12):
        raise ValueError("all points must weakly dominate the anchor point componentwise")
    q[q < 0.0] = 0.0
    return q


def _prod_except(mins: np.ndarray, k: int) -> float:
    if mins.size == 1:
        return 1.0
    prod = 1.0
    for idx, val in enumerate(mins):
        if idx != k:
            prod *= float(val)
    return prod


def exact_hypervolume_max(points: ArrayLike, anchor: Sequence[float]) -> float:
    """
    Exact hypervolume for maximization with a given anchor point, in any dimension.
    Complexity is O(2^n * d) via inclusion-exclusion, intended for exact studies.
    """
    q = _translate_and_validate(points, anchor)
    n, d = q.shape
    total = 0.0
    for r in range(1, n + 1):
        sign = 1.0 if (r % 2 == 1) else -1.0
        for subset in combinations(range(n), r):
            mins = np.min(q[list(subset), :], axis=0)
            total += sign * float(np.prod(mins))
    return total


def exact_hypervolume_gradient_max(points: ArrayLike, anchor: Sequence[float]) -> np.ndarray:
    """
    Tie-shared exact inclusion-exclusion subgradient of hypervolume for maximization.

    Returns an array of shape (n_points, dimension).
    """
    q = _translate_and_validate(points, anchor)
    n, d = q.shape
    grad = np.zeros((n, d), dtype=float)

    for r in range(1, n + 1):
        sign = 1.0 if (r % 2 == 1) else -1.0
        for subset in combinations(range(n), r):
            sub = q[list(subset), :]
            mins = np.min(sub, axis=0)
            for k in range(d):
                min_val = mins[k]
                tied_positions = [pos for pos, idx in enumerate(subset) if abs(q[idx, k] - min_val) <= 1e-12]
                if not tied_positions:
                    continue
                coeff = sign * _prod_except(mins, k) / float(len(tied_positions))
                for pos in tied_positions:
                    idx = subset[pos]
                    grad[idx, k] += coeff
    return grad


def hypervolume_2d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0)) -> float:
    return exact_hypervolume_max(points, anchor)


def hypervolume_3d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0, 0.0)) -> float:
    return exact_hypervolume_max(points, anchor)


def hypervolume_gradient_2d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0)) -> np.ndarray:
    return exact_hypervolume_gradient_max(points, anchor)


def hypervolume_gradient_3d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0, 0.0)) -> np.ndarray:
    return exact_hypervolume_gradient_max(points, anchor)


def _axis_max_gradient(q: np.ndarray) -> np.ndarray:
    """
    Gradient of sum_k max_i q_{ik}, with ties shared equally.
    """
    n, d = q.shape
    grad = np.zeros((n, d), dtype=float)
    for k in range(d):
        m = float(np.max(q[:, k]))
        tied = np.where(np.abs(q[:, k] - m) <= 1e-12)[0]
        if tied.size:
            grad[tied, k] += 1.0 / float(tied.size)
    return grad


def magnitude_2d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0)) -> float:
    q = _translate_and_validate(points, anchor, dim=2)
    lx = float(np.max(q[:, 0])) if len(q) else 0.0
    ly = float(np.max(q[:, 1])) if len(q) else 0.0
    hv = exact_hypervolume_max(q, (0.0, 0.0))
    return 1.0 + 0.5 * (lx + ly) + 0.25 * hv


def magnitude_gradient_2d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0)) -> np.ndarray:
    q = _translate_and_validate(points, anchor, dim=2)
    axis_grad = _axis_max_gradient(q)
    hv_grad = exact_hypervolume_gradient_max(q, (0.0, 0.0))
    return 0.5 * axis_grad + 0.25 * hv_grad


def magnitude_3d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0, 0.0)) -> float:
    q = _translate_and_validate(points, anchor, dim=3)
    lx = float(np.max(q[:, 0])) if len(q) else 0.0
    ly = float(np.max(q[:, 1])) if len(q) else 0.0
    lz = float(np.max(q[:, 2])) if len(q) else 0.0

    area_xy = exact_hypervolume_max(q[:, [0, 1]], (0.0, 0.0))
    area_xz = exact_hypervolume_max(q[:, [0, 2]], (0.0, 0.0))
    area_yz = exact_hypervolume_max(q[:, [1, 2]], (0.0, 0.0))
    hv3 = exact_hypervolume_max(q, (0.0, 0.0, 0.0))

    return 1.0 + 0.5 * (lx + ly + lz) + 0.25 * (area_xy + area_xz + area_yz) + 0.125 * hv3


def magnitude_gradient_3d_max(points: ArrayLike, anchor: Sequence[float] = (0.0, 0.0, 0.0)) -> np.ndarray:
    q = _translate_and_validate(points, anchor, dim=3)
    n = q.shape[0]

    axis_grad = _axis_max_gradient(q)

    hv3_grad = exact_hypervolume_gradient_max(q, (0.0, 0.0, 0.0))

    hv_xy = exact_hypervolume_gradient_max(q[:, [0, 1]], (0.0, 0.0))
    hv_xz = exact_hypervolume_gradient_max(q[:, [0, 2]], (0.0, 0.0))
    hv_yz = exact_hypervolume_gradient_max(q[:, [1, 2]], (0.0, 0.0))

    proj_grad = np.zeros((n, 3), dtype=float)
    proj_grad[:, 0] += hv_xy[:, 0]
    proj_grad[:, 1] += hv_xy[:, 1]
    proj_grad[:, 0] += hv_xz[:, 0]
    proj_grad[:, 2] += hv_xz[:, 1]
    proj_grad[:, 1] += hv_yz[:, 0]
    proj_grad[:, 2] += hv_yz[:, 1]

    return 0.5 * axis_grad + 0.25 * proj_grad + 0.125 * hv3_grad


def normalize_rows(grad: np.ndarray, eps: float = 1e-15) -> np.ndarray:
    """
    Normalize each pointwise gradient vector to unit Euclidean norm when nonzero.
    Useful for the projected normalized-gradient methods described in the paper.
    """
    g = np.asarray(grad, dtype=float).copy()
    norms = np.linalg.norm(g, axis=1)
    mask = norms > eps
    if np.any(mask):
        g[mask] /= norms[mask][:, None]
    return g


def simplex_project_rows(points: np.ndarray, total: float = 1.0) -> np.ndarray:
    """
    Project each row onto the standard simplex {x >= 0, sum x = total}.
    """
    X = np.asarray(points, dtype=float)
    Y = np.zeros_like(X)
    for i, v in enumerate(X):
        u = np.sort(v)[::-1]
        cssv = np.cumsum(u) - total
        rho = np.nonzero(u - cssv / (np.arange(len(u)) + 1) > 0)[0]
        if len(rho) == 0:
            theta = 0.0
        else:
            rho = rho[-1]
            theta = cssv[rho] / (rho + 1.0)
        Y[i] = np.maximum(v - theta, 0.0)
    return Y


def projected_gradient_step(
    points: ArrayLike,
    gradient_fn: Callable[[ArrayLike, Sequence[float]], np.ndarray],
    anchor: Sequence[float],
    step_size: float,
    projector: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    normalize_pointwise: bool = True,
) -> np.ndarray:
    """
    One projected ascent step:
        P_new = Pi( P + alpha * G(P) )

    where G may be normalized row-wise if desired.
    """
    P = _as_array(points)
    G = np.asarray(gradient_fn(P, anchor), dtype=float)
    if normalize_pointwise:
        G = normalize_rows(G)
    P_new = P + step_size * G
    if projector is not None:
        P_new = projector(P_new)
    return P_new


def _demo() -> None:
    print("=== 2D example ===")
    P2 = np.array([(1.0, 3.0), (3.0, 2.0), (5.0, 1.0)])
    a2 = (0.0, 0.0)
    print("HV_2 =", hypervolume_2d_max(P2, a2))
    print("Mag_2 =", magnitude_2d_max(P2, a2))
    print("grad HV_2 =\n", hypervolume_gradient_2d_max(P2, a2))
    print("grad Mag_2 =\n", magnitude_gradient_2d_max(P2, a2))

    print("\n=== 3D example ===")
    P3 = np.array(
        [
            (1.0, 6.0, 4.0),
            (3.0, 5.0, 1.0),
            (4.0, 4.0, 6.0),
            (5.0, 2.0, 3.0),
            (6.0, 1.0, 5.0),
            (1.0, 3.0, 7.0),
            (2.0, 2.0, 8.0),
        ]
    )
    a3 = (0.0, 0.0, 0.0)
    print("HV_3 =", hypervolume_3d_max(P3, a3))
    print("Mag_3 =", magnitude_3d_max(P3, a3))
    print("grad HV_3 =\n", hypervolume_gradient_3d_max(P3, a3))
    print("grad Mag_3 =\n", magnitude_gradient_3d_max(P3, a3))



import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, List, Sequence

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


Array = np.ndarray

# Viewpoint approximating a look from the positive diagonal, e.g. eye=(3,3,3)
DEFAULT_ELEV = 28
DEFAULT_AZIM = 45


def normalize_rows_local(G: Array) -> Array:
    H = np.array(G, dtype=float, copy=True)
    norms = np.linalg.norm(H, axis=1)
    mask = norms > 1.0e-14
    H[mask] = H[mask] / norms[mask][:, None]
    return H



def nondomination_layers(Y: Array) -> List[List[int]]:
    remaining = list(range(len(Y)))
    layers: List[List[int]] = []
    while remaining:
        front: List[int] = []
        for i in remaining:
            yi = Y[i]
            dominated = False
            for j in remaining:
                if j != i and np.all(Y[j] >= yi) and np.any(Y[j] > yi):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        layers.append(front)
        front_set = set(front)
        remaining = [i for i in remaining if i not in front_set]
    return layers


def nondominated_subset(Y: Array) -> Array:
    return Y[nondomination_layers(Y)[0]].copy() if len(Y) else Y.copy()


def repulsion_value(Y: Array, sigma: float) -> float:
    val = 0.0
    for i in range(len(Y)):
        diff = Y[i + 1 :] - Y[i]
        if len(diff):
            val += float(np.sum(np.exp(-np.sum(diff * diff, axis=1) / (sigma * sigma))))
    return val


def repulsion_gradient(Y: Array, sigma: float) -> Array:
    G = np.zeros_like(Y)
    inv = 1.0 / (sigma * sigma)
    for i in range(len(Y)):
        for j in range(i + 1, len(Y)):
            diff = Y[i] - Y[j]
            e = np.exp(-inv * float(np.dot(diff, diff)))
            c = 2.0 * inv * e * diff
            G[i] += c
            G[j] -= c
    return G


def layered_value_and_gradient_obj(Y: Array, anchor: Sequence[float], eps_layer: float, tau: float, sigma: float):
    value = 0.0
    G = np.zeros_like(Y)
    for ell, front in enumerate(nondomination_layers(Y)):
        w = eps_layer**ell
        pts = Y[front]
        value += w * magnitude_3d_max(pts, anchor=anchor)
        gfront = magnitude_gradient_3d_max(pts, anchor=anchor)
        for loc, idx in enumerate(front):
            G[idx] += w * gfront[loc]
    if tau > 0:
        value -= tau * repulsion_value(Y, sigma)
        G -= tau * repulsion_gradient(Y, sigma)
    return float(value), G


def igd(approx: Array, reference: Array) -> float:
    if len(approx) == 0 or len(reference) == 0:
        return float("inf")
    return float(np.mean([np.sqrt(np.min(np.sum((approx - r) ** 2, axis=1))) for r in reference]))


def save_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[float]]):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(header))
        w.writerows(list(rows))


def plot_objective_space(path: Path, Y0: Array, Yf: Array, ref: Array, title: str, elev: float = DEFAULT_ELEV, azim: float = DEFAULT_AZIM):
    fig = plt.figure(figsize=(10, 4.6))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    for ax, Y, ttl in ((ax1, Y0, "Initial"), (ax2, Yf, "Final")):
        if len(ref):
            ax.scatter(ref[:, 0], ref[:, 1], ref[:, 2], s=4, alpha=0.14, depthshade=False)
        ax.scatter(Y[:, 0], Y[:, 1], Y[:, 2], s=34, depthshade=False)
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel("obj 1")
        ax.set_ylabel("obj 2")
        ax.set_zlabel("obj 3")
        ax.set_title(ttl)
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_three_peaks_decision_space(path: Path, X0: Array, Xf: Array, title: str, elev: float = DEFAULT_ELEV, azim: float = DEFAULT_AZIM):
    fig = plt.figure(figsize=(10, 4.6))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    E = np.eye(3)
    for ax, X, ttl in ((ax1, X0, "Initial"), (ax2, Xf, "Final")):
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], s=34, depthshade=False)
        ax.scatter(E[:, 0], E[:, 1], E[:, 2], s=46, marker="^", depthshade=False)
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(-2.0, 2.0)
        ax.set_zlim(-2.0, 2.0)
        ax.set_xlabel("x1")
        ax.set_ylabel("x2")
        ax.set_zlabel("x3")
        ax.set_title(ttl)
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_convergence(path: Path, values, alphas, title, stride=2):
    idx = np.arange(len(values))[:: max(1, stride)]
    v = np.asarray(values)[idx]
    a = np.asarray(alphas)[idx]
    fig = plt.figure(figsize=(7, 4.2))
    ax1 = fig.add_subplot(111)
    ax1.plot(idx, v, marker="o", markersize=2.2, linewidth=1.2)
    ax1.set_xlabel("iteration")
    ax1.set_ylabel("layered value")
    ax1.set_title(title)
    ax2 = ax1.twinx()
    ax2.plot(idx, a, "--", linewidth=1.0)
    ax2.set_ylabel("step size")
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


@dataclass
class RunResult:
    X0: Array
    Y0: Array
    Xf: Array
    Yf: Array
    values: list
    alphas: list
    accepted_steps: int


def run_projected_ascent(
    objective_fn,
    jacobian_fn,
    projector,
    X0,
    anchor,
    eps_layer=1e-3,
    tau=5e-4,
    sigma=0.03,
    alpha0=0.05,
    max_iter=48,
    normalize_per_point=True,
    shrink=0.96,
    max_retries=8,
):
    X = projector(np.asarray(X0, float).copy())
    Y = objective_fn(X)
    val, GY = layered_value_and_gradient_obj(Y, anchor, eps_layer, tau, sigma)
    values = [val]
    alphas = [alpha0]
    alpha = float(alpha0)
    accepted_steps = 0
    Xinit = X.copy()
    Yinit = Y.copy()

    for _ in range(max_iter):
        J = jacobian_fn(X)
        GX = np.einsum("nij,ni->nj", J, GY)
        if normalize_per_point:
            GX = normalize_rows(GX)
        trial_alpha = alpha
        accepted = False
        for _k in range(max_retries + 1):
            Xt = projector(X + trial_alpha * GX)
            Yt = objective_fn(Xt)
            vt, GYt = layered_value_and_gradient_obj(Yt, anchor, eps_layer, tau, sigma)
            if vt >= val - 1e-12:
                X, Y, val, GY = Xt, Yt, vt, GYt
                alpha = trial_alpha
                accepted_steps += 1
                accepted = True
                break
            trial_alpha *= shrink
        if not accepted:
            alpha = trial_alpha
        values.append(val)
        alphas.append(alpha)
    return RunResult(Xinit, Yinit, X, Y, values, alphas, accepted_steps)


# Three-peaks problem

def three_peaks_objective(X: Array) -> Array:
    E = np.eye(3)
    Y = np.zeros((len(X), 3))
    for i, e in enumerate(E):
        Y[:, i] = 1.0 - np.linalg.norm(X - e, axis=1)
    return Y


def three_peaks_jacobian(X: Array) -> Array:
    E = np.eye(3)
    J = np.zeros((len(X), 3, 3))
    for i, e in enumerate(E):
        diff = X - e
        norms = np.linalg.norm(diff, axis=1)
        safe = norms > 1e-14
        J[safe, i, :] = -diff[safe] / norms[safe][:, None]
    return J


# Vehicle crashworthiness benchmark

def crashworthiness_raw_objective(X: Array) -> Array:
    x1, x2, x3, x4, x5 = [X[:, i] for i in range(5)]
    f1 = 1640.2823 + 2.3573285 * x1 + 2.3220035 * x2 + 4.5688768 * x3 + 7.7213633 * x4 + 4.4559504 * x5
    f2 = 6.5856 + 1.15 * x1 - 1.0427 * x2 + 0.9738 * x3 + 0.8364 * x4 - 0.3695 * x1 * x4 + 0.0861 * x1 * x5 + 0.3628 * x2 * x4 - 0.1106 * x1**2 - 0.3437 * x3**2 + 0.1764 * x4**2
    f3 = -0.0551 + 0.0181 * x1 + 0.1024 * x2 + 0.0421 * x3 - 0.0073 * x1 * x2 + 0.0240 * x2 * x3 - 0.0118 * x2 * x4 - 0.0204 * x3 * x4 - 0.0080 * x3 * x5 - 0.0241 * x2**2 + 0.0109 * x4**2
    return np.column_stack([f1, f2, f3])


def crashworthiness_raw_jacobian(X: Array) -> Array:
    n = len(X)
    J = np.zeros((n, 3, 5))
    x1, x2, x3, x4, x5 = [X[:, i] for i in range(5)]
    J[:, 0, :] = np.array([2.3573285, 2.3220035, 4.5688768, 7.7213633, 4.4559504])
    J[:, 1, 0] = 1.15 - 0.3695 * x4 + 0.0861 * x5 - 0.2212 * x1
    J[:, 1, 1] = -1.0427 + 0.3628 * x4
    J[:, 1, 2] = 0.9738 - 0.6874 * x3
    J[:, 1, 3] = 0.8364 - 0.3695 * x1 + 0.3628 * x2 + 0.3528 * x4
    J[:, 1, 4] = 0.0861 * x1
    J[:, 2, 0] = 0.0181 - 0.0073 * x2
    J[:, 2, 1] = 0.1024 - 0.0073 * x1 + 0.0240 * x3 - 0.0118 * x4 - 0.0482 * x2
    J[:, 2, 2] = 0.0421 + 0.0240 * x2 - 0.0204 * x4 - 0.0080 * x5
    J[:, 2, 3] = -0.0118 * x2 - 0.0204 * x3 + 0.0218 * x4
    J[:, 2, 4] = -0.0080 * x3
    return J


def make_crash_transform(seed=0, calibration_samples=800):
    rng = np.random.default_rng(seed)
    X = rng.uniform(1.0, 3.0, size=(calibration_samples, 5))
    F = crashworthiness_raw_objective(X)
    ideal = np.min(F, axis=0)
    nadir = np.max(F, axis=0)
    scale = np.maximum(nadir - ideal, 1e-12)

    def obj(X):
        return (nadir - crashworthiness_raw_objective(X)) / scale

    def jac(X):
        return -crashworthiness_raw_jacobian(X) / scale[None, :, None]

    return obj, jac, {"ideal": ideal.tolist(), "nadir": nadir.tolist(), "scale": scale.tolist()}


# Utilities

def project_box(X: Array, lo: float, hi: float) -> Array:
    return np.clip(X, lo, hi)


def approx_reference(objective_fn, projector, n_dec, lo, hi, samples, rng):
    X = rng.uniform(lo, hi, size=(samples, n_dec))
    return nondominated_subset(objective_fn(projector(X)))


def run_three_peaks(prefix="three_peaks", outdir=".", seed=8, n_points=15, max_iter=48):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    X0 = np.clip(rng.normal(0.0, 0.06, size=(n_points, 3)), -2.0, 2.0)
    projector = lambda X: project_box(X, -2.0, 2.0)
    anchor = (-4.0, -4.0, -4.0)
    ref = approx_reference(three_peaks_objective, projector, 3, -2.0, 2.0, 1200, np.random.default_rng(seed + 100))
    res = run_projected_ascent(
        three_peaks_objective,
        three_peaks_jacobian,
        projector,
        X0,
        anchor,
        alpha0=0.05,
        max_iter=max_iter,
        sigma=0.08,
        tau=5e-4,
        shrink=0.96,
        max_retries=8,
    )
    save_csv(out / f"{prefix}_initial_decisions.csv", ["x1", "x2", "x3"], res.X0)
    save_csv(out / f"{prefix}_final_decisions.csv", ["x1", "x2", "x3"], res.Xf)
    save_csv(out / f"{prefix}_initial_objectives.csv", ["f1", "f2", "f3"], res.Y0)
    save_csv(out / f"{prefix}_final_objectives.csv", ["f1", "f2", "f3"], res.Yf)
    save_csv(out / f"{prefix}_reference_archive.csv", ["f1", "f2", "f3"], ref)
    save_csv(out / f"{prefix}_history.csv", ["iteration", "layered_value", "alpha"], zip(range(len(res.values)), res.values, res.alphas))
    plot_objective_space(out / f"{prefix}_objective_space.png", res.Y0, res.Yf, ref, "Three-peaks benchmark (objective space)")
    plot_three_peaks_decision_space(out / f"{prefix}_decision_space.png", res.X0, res.Xf, "Three-peaks benchmark (decision space)")
    plot_convergence(out / f"{prefix}_convergence.png", res.values, res.alphas, "Three-peaks convergence", stride=2)
    nd0 = nondominated_subset(res.Y0)
    ndf = nondominated_subset(res.Yf)
    return {
        "problem": "three_peaks",
        "population": n_points,
        "iterations": max_iter,
        "reference_archive_size": int(len(ref)),
        "layered_value_initial": res.values[0],
        "layered_value_final": res.values[-1],
        "igd_initial": igd(nd0, ref),
        "igd_final": igd(ndf, ref),
        "nd_count_initial": int(len(nd0)),
        "nd_count_final": int(len(ndf)),
        "accepted_steps": int(res.accepted_steps),
        "initial_alpha": float(res.alphas[0]),
        "final_alpha": float(res.alphas[-1]),
    }


def run_crashworthiness(prefix="vehicle_crashworthiness", outdir=".", seed=9, n_points=15, max_iter=24):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    obj, jac, meta = make_crash_transform(seed + 50, 800)
    projector = lambda X: project_box(X, 1.0, 3.0)
    X0 = rng.uniform(1.0, 3.0, size=(n_points, 5))
    ref = approx_reference(obj, projector, 5, 1.0, 3.0, 1500, np.random.default_rng(seed + 200))
    res = run_projected_ascent(
        obj,
        jac,
        projector,
        X0,
        (-0.2, -0.2, -0.2),
        alpha0=0.02,
        max_iter=max_iter,
        sigma=0.06,
        tau=5e-4,
        shrink=0.96,
        max_retries=8,
    )
    save_csv(out / f"{prefix}_initial_decisions.csv", [f"x{i}" for i in range(1, 6)], res.X0)
    save_csv(out / f"{prefix}_final_decisions.csv", [f"x{i}" for i in range(1, 6)], res.Xf)
    save_csv(out / f"{prefix}_initial_objectives.csv", ["g1", "g2", "g3"], res.Y0)
    save_csv(out / f"{prefix}_final_objectives.csv", ["g1", "g2", "g3"], res.Yf)
    save_csv(out / f"{prefix}_initial_objectives_raw.csv", ["f1", "f2", "f3"], crashworthiness_raw_objective(res.X0))
    save_csv(out / f"{prefix}_final_objectives_raw.csv", ["f1", "f2", "f3"], crashworthiness_raw_objective(res.Xf))
    save_csv(out / f"{prefix}_reference_archive.csv", ["g1", "g2", "g3"], ref)
    save_csv(out / f"{prefix}_history.csv", ["iteration", "layered_value", "alpha"], zip(range(len(res.values)), res.values, res.alphas))
    plot_objective_space(out / f"{prefix}_objective_space.png", res.Y0, res.Yf, ref, "Vehicle crashworthiness benchmark (objective space)")
    plot_convergence(out / f"{prefix}_convergence.png", res.values, res.alphas, "Vehicle crashworthiness convergence", stride=2)
    nd0 = nondominated_subset(res.Y0)
    ndf = nondominated_subset(res.Yf)
    return {
        "problem": "vehicle_crashworthiness",
        "population": n_points,
        "iterations": max_iter,
        "reference_archive_size": int(len(ref)),
        "layered_value_initial": res.values[0],
        "layered_value_final": res.values[-1],
        "igd_initial": igd(nd0, ref),
        "igd_final": igd(ndf, ref),
        "nd_count_initial": int(len(nd0)),
        "nd_count_final": int(len(ndf)),
        "accepted_steps": int(res.accepted_steps),
        "initial_alpha": float(res.alphas[0]),
        "final_alpha": float(res.alphas[-1]),
        "normalization": meta,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", choices=["three_peaks", "vehicle_crashworthiness", "both"], default="both")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--seed", type=int, default=8)
    ap.add_argument("--n-points", type=int, default=15)
    ap.add_argument("--three-peaks-iters", type=int, default=48)
    ap.add_argument("--crash-iters", type=int, default=24)
    args = ap.parse_args()
    summaries = {}
    if args.problem in ("three_peaks", "both"):
        summaries["three_peaks"] = run_three_peaks(outdir=args.outdir, seed=args.seed, n_points=args.n_points, max_iter=args.three_peaks_iters)
    if args.problem in ("vehicle_crashworthiness", "both"):
        summaries["vehicle_crashworthiness"] = run_crashworthiness(outdir=args.outdir, seed=args.seed + 1, n_points=args.n_points, max_iter=args.crash_iters)
    with open(Path(args.outdir) / "benchmark_summary.json", "w") as f:
        json.dump(summaries, f, indent=2)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()

import argparse
import json
from pathlib import Path


def report(name: str, s: dict) -> None:
    print(name)
    print(f"  iterations: {s['iterations']}")
    print(f"  final layered magnitude: {s['layered_value_final']:.12f}")
    print(f"  initial layered magnitude: {s['layered_value_initial']:.12f}")
    print(f"  accepted steps: {s['accepted_steps']}")
    print(f"  final step size: {s['final_alpha']:.12f}")
    print(f"  nondominated count: {s['nd_count_initial']} -> {s['nd_count_final']}")
    print(f"  IGD: {s['igd_initial']:.12f} -> {s['igd_final']:.12f}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Run the 3D layered-magnitude examples and print a concise summary.")
    ap.add_argument("--problem", choices=["three_peaks", "vehicle_crashworthiness", "both"], default="both")
    ap.add_argument("--outdir", default=".", help="Directory for PNG/CSV output files.")
    ap.add_argument("--three-peaks-iters", type=int, default=48)
    ap.add_argument("--crash-iters", type=int, default=24)
    ap.add_argument("--n-points", type=int, default=15)
    ap.add_argument("--seed", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summaries = {}

    if args.problem in ("three_peaks", "both"):
        s = run_three_peaks(outdir=str(outdir), max_iter=args.three_peaks_iters, n_points=args.n_points, seed=args.seed)
        summaries["three_peaks"] = s
        report("three_peaks", s)

    if args.problem in ("vehicle_crashworthiness", "both"):
        s = run_crashworthiness(outdir=str(outdir), max_iter=args.crash_iters, n_points=args.n_points, seed=args.seed + 1)
        summaries["vehicle_crashworthiness"] = s
        report("vehicle_crashworthiness", s)

    with open(outdir / "benchmark_summary.json", "w") as f:
        json.dump(summaries, f, indent=2)


if __name__ == "__main__":
    main()
