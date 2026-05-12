#!/usr/bin/env python3
"""
Standalone Python local-PSF builder for MORIA UVP products.

This script is intended to replace the old ``uvp2psf_simst*.xOg`` flow with a
pure-Python workflow that:

1. reads ``img2extract_wfc3uv_psflist_*.uvp(.gz)``
2. estimates a sky value per star/exposure from pixels with ``8.5 < r < 13.5``
3. applies a per-exposure recentring correction from the bright core pixels
4. iteratively stacks a 4x supersampled PSF
5. writes diagnostic FITS products analogous to ``show_psf`` / ``show_str``
   and a final ``psfout`` file

It intentionally does not call the Fortran PSF builders.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

import numpy as np
from astropy.io import fits


OVERSAMPLE = 4
GRID_SIZE = 101
GRID_CENTER = 50.0
FIT_RADIUS = 12.5
CORE_RADIUS = 2.5
ANNULUS_INNER = 8.5
ANNULUS_OUTER = 13.5


@dataclass
class GroupStats:
    star: int
    exposure: int
    sky: float
    flux: float
    pmax: float
    dx: float
    dy: float
    n_annulus: int
    n_core: int
    usable: bool


@dataclass
class IterationResult:
    psf_before: np.ndarray
    raw_scene: np.ndarray
    smooth_delta: np.ndarray
    phase_correction: np.ndarray
    psf_after: np.ndarray
    exposure_rms: Dict[Tuple[int, int], float]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def infer_base_name(workdir: Path, base_name: str | None) -> str:
    if base_name:
        return base_name
    if workdir.name.upper() == "F814W":
        return "simst"
    if workdir.name.upper() == "F606W":
        return "simstV"
    raise ValueError(
        "Could not infer base name from workdir. Pass --base-name explicitly."
    )


def infer_uvp_path(workdir: Path, base_name: str, uvp_path: Path | None) -> Path:
    if uvp_path is not None:
        return uvp_path
    candidates = [
        workdir / f"img2extract_wfc3uv_psflist_{base_name}.uvp",
        workdir / f"img2extract_wfc3uv_psflist_{base_name}.uvp.gz",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Could not find UVP file for base name '{base_name}' in {workdir}"
    )


def load_uvp(path: Path) -> Dict[str, np.ndarray]:
    with _open_text(path) as handle:
        data = np.loadtxt(handle, comments="#")

    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.shape[1] < 12:
        raise ValueError(f"Expected at least 12 UVP columns, got {data.shape[1]}")

    result = {
        "u": data[:, 0].astype(float),
        "v": data[:, 1].astype(float),
        "p": data[:, 2].astype(float),
        "nim": data[:, 3].astype(int),
        "i": data[:, 4].astype(int),
        "j": data[:, 5].astype(int),
        "xc": data[:, 6].astype(float),
        "yc": data[:, 7].astype(float),
        "du": data[:, 8].astype(float),
        "dv": data[:, 9].astype(float),
        "r": data[:, 10].astype(float),
        "nst": data[:, 11].astype(int),
    }
    if data.shape[1] >= 13:
        result["use_flag"] = data[:, 12].astype(int)
    else:
        result["use_flag"] = np.ones(data.shape[0], dtype=int)
    return result


def star_exposure_groups(nst: np.ndarray, nim: np.ndarray) -> Dict[Tuple[int, int], np.ndarray]:
    groups: Dict[Tuple[int, int], List[int]] = {}
    for idx, key in enumerate(zip(nst.tolist(), nim.tolist())):
        groups.setdefault(key, []).append(idx)
    return {key: np.asarray(indices, dtype=int) for key, indices in groups.items()}


def read_good_list(path: Path, nstars: int) -> Dict[int, bool]:
    flags = {star: True for star in range(1, nstars + 1)}
    if not path.exists():
        raise FileNotFoundError(f"Missing good-PSF list: {path}")

    with path.open("r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            star_s, flag_s = line.split()[:2]
            flags[int(star_s)] = bool(int(flag_s))
    return flags


def write_default_good_list(path: Path, nstars: int, target_star: int = 1) -> None:
    with path.open("w") as handle:
        for star in range(1, nstars + 1):
            flag = 0 if star == target_star else 1
            handle.write(f"{star:2d}   {flag}\n")


def deposit_bilinear(
    image: np.ndarray,
    weight: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray | None = None,
) -> None:
    if weights is None:
        weights = np.ones_like(values)

    ix0 = np.floor(x).astype(int)
    iy0 = np.floor(y).astype(int)
    fx = x - ix0
    fy = y - iy0

    contributions = (
        (ix0, iy0, (1.0 - fx) * (1.0 - fy)),
        (ix0 + 1, iy0, fx * (1.0 - fy)),
        (ix0, iy0 + 1, (1.0 - fx) * fy),
        (ix0 + 1, iy0 + 1, fx * fy),
    )

    ny, nx = image.shape
    for ix, iy, frac in contributions:
        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        if not np.any(valid):
            continue
        np.add.at(image, (iy[valid], ix[valid]), values[valid] * frac[valid] * weights[valid])
        np.add.at(weight, (iy[valid], ix[valid]), frac[valid] * weights[valid])


def sample_bilinear(image: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ix0 = np.floor(x).astype(int)
    iy0 = np.floor(y).astype(int)
    fx = x - ix0
    fy = y - iy0

    ny, nx = image.shape
    ix1 = ix0 + 1
    iy1 = iy0 + 1

    out = np.zeros_like(x, dtype=float)
    contributions = (
        (ix0, iy0, (1.0 - fx) * (1.0 - fy)),
        (ix1, iy0, fx * (1.0 - fy)),
        (ix0, iy1, (1.0 - fx) * fy),
        (ix1, iy1, fx * fy),
    )
    for ix, iy, frac in contributions:
        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        if np.any(valid):
            out[valid] += image[iy[valid], ix[valid]] * frac[valid]
    return out


def barsig(values: np.ndarray, niters: int = 30) -> Tuple[float, float, int]:
    if values.size == 0:
        return 0.0, 0.999, 0

    bar = 0.0
    sig = 9e9
    nuse = 0
    for _ in range(niters):
        keep = np.abs(values - bar) <= 2.25 * sig
        nuse = int(np.count_nonzero(keep))
        if nuse == 0:
            continue
        selected = values[keep]
        old_bar = bar
        bar = float(selected.mean())
        if nuse > 1:
            sig = float(np.abs(selected - old_bar).sum() / (nuse - 1))
    if nuse <= 1:
        sig = 0.999
    return bar, sig, nuse


def evaluate_local_psf(x: np.ndarray, y: np.ndarray, psf: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dd = np.hypot(x, y)
    out = np.zeros_like(x, dtype=float)

    mask_outer = (dd > 4.0) & (dd <= 12.0)
    if np.any(mask_outer):
        rx = 50.0 + x[mask_outer] * OVERSAMPLE
        ry = 50.0 + y[mask_outer] * OVERSAMPLE
        out[mask_outer] = sample_bilinear(psf, rx, ry)

    mask_core = dd <= 4.0
    if not np.any(mask_core):
        return out

    rx = 50.0 + x[mask_core] * OVERSAMPLE
    ry = 50.0 + y[mask_core] * OVERSAMPLE
    ix = np.floor(rx).astype(int)
    iy = np.floor(ry).astype(int)
    fx = rx - ix
    fy = ry - iy

    def pix(dx: int, dy: int) -> np.ndarray:
        xx = np.clip(ix + dx, 0, psf.shape[1] - 1)
        yy = np.clip(iy + dy, 0, psf.shape[0] - 1)
        return psf[yy, xx]

    a1 = pix(0, 0)
    b1 = (pix(1, 0) - pix(-1, 0)) / 2.0
    c1 = (pix(0, 1) - pix(0, -1)) / 2.0
    d1 = (pix(1, 0) + pix(-1, 0) - 2.0 * a1) / 2.0
    f1 = (pix(0, 1) + pix(0, -1) - 2.0 * a1) / 2.0
    e1 = pix(1, 1) - a1

    a2 = pix(1, 0)
    b2 = (pix(2, 0) - pix(0, 0)) / 2.0
    c2 = (pix(1, 1) - pix(1, -1)) / 2.0
    d2 = (pix(2, 0) + pix(0, 0) - 2.0 * a2) / 2.0
    f2 = (pix(1, 1) + pix(1, -1) - 2.0 * a2) / 2.0
    e2 = -(pix(0, 1) - a2)

    a3 = pix(0, 1)
    b3 = (pix(1, 1) - pix(-1, 1)) / 2.0
    c3 = (pix(0, 2) - pix(0, 0)) / 2.0
    d3 = (pix(1, 1) + pix(-1, 1) - 2.0 * a3) / 2.0
    f3 = (pix(0, 2) + pix(0, 0) - 2.0 * a3) / 2.0
    e3 = -(pix(1, 0) - a3)

    a4 = pix(1, 1)
    b4 = (pix(2, 1) - pix(0, 1)) / 2.0
    c4 = (pix(1, 2) - pix(1, 0)) / 2.0
    d4 = (pix(2, 1) + pix(0, 1) - 2.0 * a4) / 2.0
    f4 = (pix(1, 2) + pix(1, 0) - 2.0 * a4) / 2.0
    e4 = pix(0, 0) - a4

    v1 = a1 + b1 * fx + c1 * fy + d1 * fx**2 + e1 * fx * fy + f1 * fy**2
    v2 = a2 + b2 * (fx - 1.0) + c2 * fy + d2 * (fx - 1.0) ** 2 + e2 * (fx - 1.0) * fy + f2 * fy**2
    v3 = a3 + b3 * fx + c3 * (fy - 1.0) + d3 * fx**2 + e3 * fx * (fy - 1.0) + f3 * (fy - 1.0) ** 2
    v4 = a4 + b4 * (fx - 1.0) + c4 * (fy - 1.0) + d4 * (fx - 1.0) ** 2 + e4 * (fx - 1.0) * (fy - 1.0) + f4 * (fy - 1.0) ** 2

    out[mask_core] = (
        (1.0 - fx) * (1.0 - fy) * v1
        + fx * (1.0 - fy) * v2
        + (1.0 - fx) * fy * v3
        + fx * fy * v4
    )
    return out


def smooth_image(image: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array(
        [
            [1.0, 2.0, 1.0],
            [2.0, 4.0, 2.0],
            [1.0, 2.0, 1.0],
        ],
        dtype=float,
    )
    kernel /= kernel.sum()

    smoothed = image.copy()
    for _ in range(passes):
        padded = np.pad(smoothed, 1, mode="edge")
        out = np.zeros_like(smoothed)
        for dy in range(3):
            for dx in range(3):
                out += kernel[dy, dx] * padded[dy : dy + smoothed.shape[0], dx : dx + smoothed.shape[1]]
        smoothed = out
    return smoothed


def normalize_psf(psf: np.ndarray) -> np.ndarray:
    clipped = np.clip(psf, 0.0, None)
    total = clipped.sum()
    if total <= 0:
        raise ValueError("PSF normalization failed: non-positive total flux")
    return clipped / total


def shift_image_grid(image: np.ndarray, dx: float, dy: float) -> np.ndarray:
    yy, xx = np.indices(image.shape, dtype=float)
    return sample_bilinear(image, xx - dx, yy - dy)


def recenter_psf_image(image: np.ndarray) -> np.ndarray:
    positive = np.clip(image, 0.0, None)
    pmax = float(positive.max())
    if pmax <= 0.0:
        return image

    yy, xx = np.indices(image.shape, dtype=float)
    rr = np.hypot(xx - GRID_CENTER, yy - GRID_CENTER)
    focus = rr <= 24.0
    weights = np.where(focus & (positive >= 0.2 * pmax), positive, 0.0)
    if weights.sum() <= 0.0:
        weights = np.where(focus, positive, 0.0)
    if weights.sum() <= 0.0:
        return image

    cx = float((xx * weights).sum() / weights.sum())
    cy = float((yy * weights).sum() / weights.sum())
    return shift_image_grid(image, GRID_CENTER - cx, GRID_CENTER - cy)


def radialize_psf(image: np.ndarray) -> np.ndarray:
    positive = np.clip(image, 0.0, None)
    yy, xx = np.indices(image.shape, dtype=float)
    rr = np.rint(np.hypot(xx - GRID_CENTER, yy - GRID_CENTER)).astype(int)
    max_r = int(rr.max())
    profile = np.zeros(max_r + 1, dtype=float)
    counts = np.zeros(max_r + 1, dtype=float)
    np.add.at(profile, rr.ravel(), positive.ravel())
    np.add.at(counts, rr.ravel(), np.isfinite(positive).ravel().astype(float))
    profile = np.divide(profile, counts, out=np.zeros_like(profile), where=counts > 0)
    radial = profile[rr]
    return radial


def crop100(image101: np.ndarray) -> np.ndarray:
    return image101[:100, :100]


def _clamped_center(index: int, radius: int, size: int) -> int:
    return max(radius, min(index, size - radius - 1))


def sm7plan(image: np.ndarray) -> np.ndarray:
    out = np.zeros_like(image)
    ny, nx = image.shape
    for y in range(ny):
        for x in range(nx):
            xc = _clamped_center(x, 3, nx)
            yc = _clamped_center(y, 3, ny)
            patch = image[yc - 3 : yc + 4, xc - 3 : xc + 4]
            grid = np.arange(-3, 4, dtype=float)
            xx, yy = np.meshgrid(grid, grid)
            a = float(patch.sum() / xx.size)
            b = float((patch * xx).sum() / (xx * xx).sum())
            c = float((patch * yy).sum() / (yy * yy).sum())
            out[y, x] = a + b * (x - xc) + c * (y - yc)
    return out


def sm5plan(image: np.ndarray) -> np.ndarray:
    aa = np.array(
        [
            [0.2500, 0.5000, 0.5000, 0.5000, 0.2500],
            [0.5000, 1.0000, 1.0000, 1.0000, 0.5000],
            [0.5000, 1.0000, 1.0000, 1.0000, 0.5000],
            [0.5000, 1.0000, 1.0000, 1.0000, 0.5000],
            [0.2500, 0.5000, 0.5000, 0.5000, 0.2500],
        ],
        dtype=float,
    )
    bb = np.tile(np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float), (5, 1))
    cc = bb.T
    out = np.zeros_like(image)
    ny, nx = image.shape
    for y in range(ny):
        for x in range(nx):
            xc = _clamped_center(x, 2, nx)
            yc = _clamped_center(y, 2, ny)
            patch = image[yc - 2 : yc + 3, xc - 2 : xc + 3]
            a = float((aa * patch).sum() / 16.0)
            b = float((bb * patch).sum() / 25.0)
            c = float((cc * patch).sum() / 25.0)
            out[y, x] = a + b * (x - xc) + c * (y - yc)
    return out


def sm5quar(image: np.ndarray) -> np.ndarray:
    aa = np.array(
        [
            [1.0408, -2.0204, 1.9592, -2.0204, 1.0408],
            [-2.0204, -0.4898, 5.0204, -0.4898, -2.0204],
            [1.9592, 5.0204, 11.0408, 5.0204, 1.9592],
            [-2.0204, -0.4898, 5.0204, -0.4898, -2.0204],
            [1.0408, -2.0204, 1.9592, -2.0204, 1.0408],
        ],
        dtype=float,
    )
    out = np.zeros_like(image)
    ny, nx = image.shape
    for y in range(ny):
        for x in range(nx):
            xc = _clamped_center(x, 2, nx)
            yc = _clamped_center(y, 2, ny)
            patch = image[yc - 2 : yc + 3, xc - 2 : xc + 3]
            out[y, x] = float((aa * patch).sum() / 25.0)
    return out


def sm5quad(image: np.ndarray) -> np.ndarray:
    aa = np.array(
        [
            [-1.8571, 0.2857, 1.0000, 0.2857, -1.8571],
            [0.2857, 2.4286, 3.1429, 2.4286, 0.2857],
            [1.0000, 3.1429, 3.8571, 3.1429, 1.0000],
            [0.2857, 2.4286, 3.1429, 2.4286, 0.2857],
            [-1.8571, 0.2857, 1.0000, 0.2857, -1.8571],
        ],
        dtype=float,
    )
    bb = np.tile(np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float), (5, 1))
    cc = bb.T
    dd = np.tile(np.array([0.7143, -0.3571, -0.7143, -0.3571, 0.7143], dtype=float), (5, 1))
    ee = np.array(
        [
            [1.0000, 0.5000, 0.0000, -0.5000, -1.0000],
            [0.5000, 0.2500, 0.0000, -0.2500, -0.5000],
            [0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
            [-0.5000, -0.2500, 0.0000, 0.2500, 0.5000],
            [-1.0000, -0.5000, 0.0000, 0.5000, 1.0000],
        ],
        dtype=float,
    )
    ff = np.tile(np.array([0.7143, -0.3571, -0.7143, -0.3571, 0.7143], dtype=float)[:, None], (1, 5))
    out = np.zeros_like(image)
    ny, nx = image.shape
    for y in range(ny):
        for x in range(nx):
            xc = _clamped_center(x, 2, nx)
            yc = _clamped_center(y, 2, ny)
            patch = image[yc - 2 : yc + 3, xc - 2 : xc + 3]
            a = float((aa * patch).sum() / 25.0)
            b = float((bb * patch).sum() / 25.0)
            c = float((cc * patch).sum() / 25.0)
            d = float((dd * patch).sum() / 25.0)
            e = float((ee * patch).sum() / 25.0)
            f = float((ff * patch).sum() / 25.0)
            dx = x - xc
            dy = y - yc
            out[y, x] = a + b * dx + c * dy + d * dx**2 + e * dx * dy + f * dy**2
    return out


def psf_smooth_like_fortran(image: np.ndarray) -> np.ndarray:
    ma = sm5quar(image)
    mb = sm5quad(image)
    mc = sm5plan(image)
    md = sm7plan(image)
    yy, xx = np.indices(image.shape, dtype=float)
    r = np.hypot(xx - GRID_CENTER, yy - GRID_CENTER) / OVERSAMPLE
    out = md.copy()
    out[r < 8.0] = (mc[r < 8.0] + md[r < 8.0]) / 2.0
    out[r < 7.0] = mc[r < 7.0]
    out[r < 6.0] = (mb[r < 6.0] + mc[r < 6.0]) / 2.0
    out[r < 5.0] = mb[r < 5.0]
    out[r < 3.0] = (ma[r < 3.0] + mb[r < 3.0]) / 2.0
    out[r < 2.0] = ma[r < 2.0]
    return out


def psf_phase_normalize_like_fortran(image: np.ndarray) -> np.ndarray:
    psfc = np.zeros_like(image)
    ny, nx = image.shape
    for y in range(ny):
        for x in range(nx):
            acc = 0.0
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    ff = 1.0 / 16.0
                    if abs(dx) == 2:
                        ff *= 0.5
                    if abs(dy) == 2:
                        ff *= 0.5
                    xx = np.clip(x + dx, 0, nx - 1)
                    yy = np.clip(y + dy, 0, ny - 1)
                    acc += ff * image[yy, xx]
            psfc[y, x] = acc

    strobe = np.zeros((OVERSAMPLE, OVERSAMPLE), dtype=float)
    strobc = np.zeros((OVERSAMPLE, OVERSAMPLE), dtype=float)
    yy, xx = np.indices(image.shape, dtype=float)
    rr = np.hypot(xx - GRID_CENTER, yy - GRID_CENTER) / OVERSAMPLE
    for y in range(ny):
        for x in range(nx):
            if rr[y, x] >= 5.0:
                continue
            strobe[x % OVERSAMPLE, y % OVERSAMPLE] += image[y, x]
            strobc[x % OVERSAMPLE, y % OVERSAMPLE] += psfc[y, x]

    ratio = np.ones_like(strobe)
    valid = strobc > 0
    ratio[valid] = strobe[valid] / strobc[valid]
    out = image.copy()
    for y in range(ny):
        for x in range(nx):
            phase = ratio[x % OVERSAMPLE, y % OVERSAMPLE]
            if phase != 0.0:
                out[y, x] = image[y, x] / phase
    return out


def build_residual_map(
    groups: Dict[Tuple[int, int], np.ndarray],
    used_groups: Dict[Tuple[int, int], bool],
    star_fluxes: Dict[int, float],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    current_psf: np.ndarray,
) -> np.ndarray:
    buckets: List[List[float]] = [[[] for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
    for key, indices in groups.items():
        if not used_groups.get(key, False):
            continue
        flux = star_fluxes.get(key[0], 1.0)
        if flux <= 0.0:
            continue
        du = du_corr[indices]
        dv = dv_corr[indices]
        obs = psub[indices] / flux
        pred = evaluate_local_psf(du, dv, current_psf)
        resid = obs - pred
        for xval, yval, rval in zip(du.tolist(), dv.tolist(), resid.tolist()):
            if abs(xval) > FIT_RADIUS or abs(yval) > FIT_RADIUS:
                continue
            xbase = xval * OVERSAMPLE + GRID_CENTER
            ybase = yval * OVERSAMPLE + GRID_CENTER
            ix0 = max(0, int(np.floor(xbase - 2.0)))
            ix1 = min(GRID_SIZE - 1, int(np.ceil(xbase + 2.0)))
            iy0 = max(0, int(np.floor(ybase - 2.0)))
            iy1 = min(GRID_SIZE - 1, int(np.ceil(ybase + 2.0)))
            for ix in range(ix0, ix1 + 1):
                if abs(xval - (ix - GRID_CENTER) / OVERSAMPLE) >= 0.5:
                    continue
                for iy in range(iy0, iy1 + 1):
                    if abs(yval - (iy - GRID_CENTER) / OVERSAMPLE) < 0.5:
                        buckets[iy][ix].append(rval)

    psfa = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)
    for iy in range(GRID_SIZE):
        for ix in range(GRID_SIZE):
            if buckets[iy][ix]:
                psfa[iy, ix], _, _ = barsig(np.asarray(buckets[iy][ix], dtype=float))
    return psfa


def analyze_star_exposures(
    star: int,
    exposures: List[int],
    groups: Dict[Tuple[int, int], np.ndarray],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    psf: np.ndarray,
) -> Tuple[float, np.ndarray]:
    z_values = np.full(len(exposures), np.nan, dtype=float)
    for idx, exposure in enumerate(exposures):
        indices = groups.get((star, exposure))
        if indices is None:
            continue
        mask = np.hypot(du_corr[indices], dv_corr[indices]) < CORE_RADIUS
        if not np.any(mask):
            continue
        core_idx = indices[mask]
        ptot = float(psub[core_idx].sum())
        ftot = float(evaluate_local_psf(du_corr[core_idx], dv_corr[core_idx], psf).sum())
        if ftot > 0.0:
            z_values[idx] = ptot / ftot

    use = np.isfinite(z_values).astype(int)
    if use.sum() == 0:
        return 1.0, use

    while True:
        candidate = -1
        candidate_sig = -np.inf
        for idx in range(len(exposures)):
            if use[idx] != 1:
                continue
            others = [j for j in range(len(exposures)) if j != idx and use[j] == 1]
            if len(others) <= 2:
                mean = float(np.nanmean(z_values[use == 1]))
                break
            other_vals = z_values[others]
            mean = float(np.nanmean(other_vals))
            sig = float(np.nanstd(other_vals, ddof=1))
            if sig <= 0.0:
                continue
            score = abs(z_values[idx] - mean) / sig
            if score > candidate_sig:
                candidate_sig = score
                candidate = idx
        else:
            if candidate_sig > 4.0:
                use[candidate] = 0
                if use.sum() >= 4:
                    continue
            mean = float(np.nanmean(z_values[use == 1]))
        break

    if not np.isfinite(mean) or mean <= 0.0:
        mean = 1.0
    return mean, use
def build_group_stats(
    uvp: Dict[str, np.ndarray],
    groups: Dict[Tuple[int, int], np.ndarray],
    recenter_frac: float,
    min_flux: float,
) -> Tuple[
    Dict[Tuple[int, int], GroupStats],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    psub = np.zeros_like(uvp["p"], dtype=float)
    du_corr = uvp["du"].astype(float).copy()
    dv_corr = uvp["dv"].astype(float).copy()
    r_corr = uvp["r"].astype(float).copy()

    stats: Dict[Tuple[int, int], GroupStats] = {}
    for (star, exposure), indices in groups.items():
        radii = uvp["r"][indices]
        pixels = uvp["p"][indices]

        annulus = indices[(radii > ANNULUS_INNER) & (radii < ANNULUS_OUTER)]
        if annulus.size:
            sky, _, _ = barsig(uvp["p"][annulus].astype(float))
        else:
            sky, _, _ = barsig(pixels.astype(float))

        psub_group = pixels - sky
        psub[indices] = psub_group

        core_mask = radii <= CORE_RADIUS
        core_indices = indices[core_mask]
        core_values = psub_group[core_mask]
        pmax = float(np.max(core_values)) if core_values.size else 0.0

        dx = 0.0
        dy = 0.0
        flux = float(core_values.sum())
        usable = bool(pmax > 0.0 and flux > min_flux)

        stats[(star, exposure)] = GroupStats(
            star=star,
            exposure=exposure,
            sky=sky,
            flux=flux,
            pmax=pmax,
            dx=dx,
            dy=dy,
            n_annulus=int(annulus.size),
            n_core=int(core_indices.size),
            usable=usable,
        )

    return stats, psub, du_corr, dv_corr, r_corr


def compute_star_fluxes(
    stats: Dict[Tuple[int, int], GroupStats],
    nstars: int,
    num_exposures: int,
) -> Dict[int, float]:
    fluxes: Dict[int, float] = {}
    for star in range(1, nstars + 1):
        values = [
            group.flux
            for key, group in stats.items()
            if key[0] == star and group.flux > 0.0 and group.usable
        ]
        if values:
            fluxes[star] = float(np.sum(values) / max(num_exposures, 1))
        else:
            fluxes[star] = 1.0
    return fluxes


def get_group_observations(
    indices: np.ndarray,
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    flux: float,
    shift: Tuple[float, float] = (0.0, 0.0),
    fit_radius: float = FIT_RADIUS,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if flux <= 0.0:
        return np.array([]), np.array([]), np.array([])

    du = du_corr[indices] - shift[0]
    dv = dv_corr[indices] - shift[1]
    obs = psub[indices] / flux
    mask = np.abs(du) <= fit_radius
    mask &= np.abs(dv) <= fit_radius
    if not np.any(mask):
        return np.array([]), np.array([]), np.array([])
    return du[mask], dv[mask], obs[mask]


def build_group_stamp(
    indices: np.ndarray,
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    flux: float,
    shift: Tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    du, dv, obs = get_group_observations(indices, psub, du_corr, dv_corr, flux, shift=shift)
    stamp = np.full((GRID_SIZE, GRID_SIZE), np.nan, dtype=float)
    if obs.size == 0:
        return stamp

    accum = np.zeros_like(stamp)
    weight = np.zeros_like(stamp)
    x = du * OVERSAMPLE + GRID_CENTER
    y = dv * OVERSAMPLE + GRID_CENTER
    deposit_bilinear(accum, weight, x, y, obs)
    valid = weight > 0
    stamp[valid] = accum[valid] / weight[valid]
    return stamp


def estimate_group_shift(
    current_psf: np.ndarray,
    indices: np.ndarray,
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    flux: float,
) -> Tuple[float, float, float]:
    if flux <= 0.0:
        return 0.0, 0.0, np.inf

    du = du_corr[indices]
    dv = dv_corr[indices]
    obs = psub[indices] / flux
    mask = np.hypot(du, dv) <= 4.0
    if np.count_nonzero(mask) < 4:
        return 0.0, 0.0, np.inf

    du = du[mask]
    dv = dv[mask]
    obs = obs[mask]
    weights = np.clip(obs, 0.0, None)
    if weights.sum() <= 0.0:
        return 0.0, 0.0, np.inf

    if current_psf is None:
        bright = weights >= 0.6 * float(weights.max())
        if np.count_nonzero(bright) < 1:
            bright = weights > 0.0
        dx = float(np.sum(du[bright] * weights[bright]) / np.sum(weights[bright]))
        dy = float(np.sum(dv[bright] * weights[bright]) / np.sum(weights[bright]))
        return dx, dy, np.inf

    best_e = np.inf
    best_dx = 0.0
    best_dy = 0.0
    for dx in np.arange(-4.0, 4.0001, 0.25):
        x = (du - dx) * OVERSAMPLE + GRID_CENTER
        for dy in np.arange(-4.0, 4.0001, 0.25):
            y = (dv - dy) * OVERSAMPLE + GRID_CENTER
            pred = sample_bilinear(current_psf, x, y)
            resid = obs - pred
            e = float(np.sum(weights * resid**2) / np.sum(weights))
            if e < best_e:
                best_e = e
                best_dx = float(dx)
                best_dy = float(dy)
    return best_dx, best_dy, float(np.sqrt(best_e))


def stack_group_stamps(
    groups: Dict[Tuple[int, int], np.ndarray],
    used_groups: Dict[Tuple[int, int], bool],
    star_fluxes: Dict[int, float],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    current_psf: np.ndarray | None,
    exposure_rms_prev: Dict[Tuple[int, int], float] | None = None,
) -> Tuple[np.ndarray, Dict[Tuple[int, int], Tuple[float, float]], Dict[Tuple[int, int], float]]:
    stamps: List[np.ndarray] = []
    align_shifts: Dict[Tuple[int, int], Tuple[float, float]] = {}
    seed_rms: Dict[Tuple[int, int], float] = {}
    for key, indices in groups.items():
        if not used_groups.get(key, False):
            continue

        star = key[0]
        flux = star_fluxes[star]
        if flux <= 0.0:
            continue

        dx, dy, rms = estimate_group_shift(current_psf, indices, psub, du_corr, dv_corr, flux)
        align_shifts[key] = (dx, dy)
        seed_rms[key] = rms
        stamp = build_group_stamp(indices, psub, du_corr, dv_corr, flux, shift=(dx, dy))
        if not np.isfinite(stamp).any():
            continue
        if exposure_rms_prev is not None:
            prev = exposure_rms_prev.get(key, np.inf)
            if np.isfinite(prev) and prev > 0.5:
                continue
        stamps.append(stamp)

    if not stamps:
        raise ValueError("No usable group stamps available for PSF construction")

    cube = np.stack(stamps, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        raw = np.nanmedian(cube, axis=0)
    raw = np.nan_to_num(raw, nan=0.0)
    if current_psf is not None:
        raw = np.where(np.isfinite(raw) & (raw != 0.0), raw, current_psf)
    return raw, align_shifts, seed_rms


def compute_exposure_rms(
    psf: np.ndarray,
    groups: Dict[Tuple[int, int], np.ndarray],
    used_groups: Dict[Tuple[int, int], bool],
    star_fluxes: Dict[int, float],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    align_shifts: Dict[Tuple[int, int], Tuple[float, float]],
    fit_radius: float = 4.0,
) -> Dict[Tuple[int, int], float]:
    exposure_rms: Dict[Tuple[int, int], float] = {}
    for key, indices in groups.items():
        star = key[0]
        flux = star_fluxes.get(star, 1.0)
        if flux <= 0.0:
            exposure_rms[key] = np.inf
            continue

        mask = np.hypot(du_corr[indices], dv_corr[indices]) <= fit_radius
        if not np.any(mask):
            exposure_rms[key] = np.inf
            continue

        idx = indices[mask]
        obs = psub[idx] / flux
        pred = evaluate_local_psf(du_corr[idx], dv_corr[idx], psf)
        resid = obs - pred
        exposure_rms[key] = float(np.sqrt(np.mean(resid**2)))
    return exposure_rms


def iterate_psf(
    groups: Dict[Tuple[int, int], np.ndarray],
    stats: Dict[Tuple[int, int], GroupStats],
    good_flags: Dict[int, bool],
    star_fluxes: Dict[int, float],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    iterations: int,
    reject_sigma: float,
) -> Tuple[
    np.ndarray,
    List[IterationResult],
    Dict[Tuple[int, int], bool],
    Dict[Tuple[int, int], float],
    Dict[Tuple[int, int], Tuple[float, float]],
]:
    exposures = sorted({key[1] for key in groups})
    nstars = max(key[0] for key in groups)
    used_groups = {key: (good_flags.get(key[0], True) and group.usable) for key, group in stats.items()}
    all_iterations: List[IterationResult] = []
    exposure_rms_prev: Dict[Tuple[int, int], float] = {key: np.inf for key in groups}
    align_shifts: Dict[Tuple[int, int], Tuple[float, float]] = {key: (0.0, 0.0) for key in groups}
    current_psf = np.zeros((GRID_SIZE, GRID_SIZE), dtype=float)
    cumulative_psf = np.zeros_like(current_psf)
    yy, xx = np.indices(current_psf.shape, dtype=float)
    support = np.hypot(xx - GRID_CENTER, yy - GRID_CENTER) / OVERSAMPLE <= 12.0

    for _ in range(iterations):
        psf_before = current_psf.copy()
        psfa = build_residual_map(
            groups,
            used_groups,
            star_fluxes,
            psub,
            du_corr,
            dv_corr,
            current_psf,
        )
        cumulative_psf = cumulative_psf + psfa
        psfs = psf_smooth_like_fortran(cumulative_psf)
        current_psf = psf_phase_normalize_like_fortran(psfs)
        current_psf = np.where(support, current_psf, 0.0)
        smooth_delta = psfs - psf_before
        phase_correction = current_psf - psfs

        next_fluxes: Dict[int, float] = {}
        next_used = used_groups.copy()
        for star in range(1, nstars + 1):
            flux, use_mask = analyze_star_exposures(
                star,
                exposures,
                groups,
                psub,
                du_corr,
                dv_corr,
                current_psf,
            )
            next_fluxes[star] = flux
            for idx, exposure in enumerate(exposures):
                key = (star, exposure)
                if key in next_used:
                    next_used[key] = bool(good_flags.get(star, True) and stats[key].usable and use_mask[idx] == 1)

        star_fluxes = next_fluxes
        used_groups = next_used

        exposure_rms = compute_exposure_rms(
            current_psf,
            groups,
            used_groups,
            star_fluxes,
            psub,
            du_corr,
            dv_corr,
            align_shifts,
        )
        all_iterations.append(
            IterationResult(
                psf_before=psf_before,
                raw_scene=psf_before + psfa,
                smooth_delta=smooth_delta,
                phase_correction=phase_correction,
                psf_after=current_psf.copy(),
                exposure_rms=exposure_rms,
            )
        )

        exposure_rms_prev = exposure_rms

    return current_psf, all_iterations, used_groups, exposure_rms_prev, align_shifts


def render_stamp(
    du: np.ndarray,
    dv: np.ndarray,
    values: np.ndarray,
    panel_size: int = 100,
) -> np.ndarray:
    grid = np.zeros((panel_size, panel_size), dtype=float)
    weight = np.zeros((panel_size, panel_size), dtype=float)
    center = (panel_size - 1) / 2.0
    x = du * OVERSAMPLE + center
    y = dv * OVERSAMPLE + center
    deposit_bilinear(grid, weight, x, y, values)
    return np.divide(grid, weight, out=np.zeros_like(grid), where=weight > 0)


def draw_box(panel: np.ndarray, value: float = 0.1) -> None:
    panel[0, :] = value
    panel[-1, :] = value
    panel[:, 0] = value
    panel[:, -1] = value


def build_show_psf(iterations: List[IterationResult]) -> np.ndarray:
    canvas = np.zeros((1001, 501), dtype=np.float32)
    for it_index, iteration in enumerate(iterations[:6]):
        row0 = (it_index) * 100
        panels = [
            crop100(iteration.psf_before),
            crop100(iteration.raw_scene - iteration.psf_before),
            crop100(iteration.smooth_delta),
            crop100(iteration.phase_correction),
            crop100(iteration.psf_after - iteration.psf_before),
        ]
        for col, panel in enumerate(panels):
            x0 = col * 100
            canvas[row0 : row0 + 100, x0 : x0 + 100] = panel.astype(np.float32)
    return canvas


def build_show_str(
    uvp: Dict[str, np.ndarray],
    groups: Dict[Tuple[int, int], np.ndarray],
    stats: Dict[Tuple[int, int], GroupStats],
    good_flags: Dict[int, bool],
    used_groups: Dict[Tuple[int, int], bool],
    star_fluxes: Dict[int, float],
    psub: np.ndarray,
    du_corr: np.ndarray,
    dv_corr: np.ndarray,
    final_psf: np.ndarray,
    align_shifts: Dict[Tuple[int, int], Tuple[float, float]],
) -> np.ndarray:
    nstars = int(uvp["nst"].max())
    canvas = np.zeros((2001, 4001), dtype=np.float32)
    exposures = sorted({int(x) for x in uvp["nim"]})

    for star in range(1, nstars + 1):
        x0 = (star - 1) * 100
        star_panels_scene: List[np.ndarray] = []
        star_panels_resid: List[np.ndarray] = []

        for exp_idx, exposure in enumerate(exposures[:16], start=3):
            key = (star, exposure)
            if key not in groups:
                continue
            indices = groups[key]
            flux = star_fluxes.get(star, 1.0)
            du = du_corr[indices]
            dv = dv_corr[indices]
            scene = psub[indices] / flux
            pred = evaluate_local_psf(du, dv, final_psf)
            resid = scene - pred

            scene_panel = render_stamp(du, dv, scene)
            resid_panel = render_stamp(du, dv, resid)
            if not used_groups.get(key, False):
                draw_box(resid_panel)

            y0 = exp_idx * 100
            canvas[y0 : y0 + 100, x0 : x0 + 100] = resid_panel.astype(np.float32)
            star_panels_scene.append(scene_panel)
            star_panels_resid.append(resid_panel)

        if star_panels_scene:
            stacked_scene = np.mean(star_panels_scene, axis=0)
            stacked_resid = np.mean(star_panels_resid, axis=0)
        else:
            stacked_scene = np.zeros((100, 100), dtype=float)
            stacked_resid = np.zeros((100, 100), dtype=float)

        psf_panel = crop100(final_psf).astype(np.float32)
        canvas[200:300, x0 : x0 + 100] = psf_panel
        canvas[100:200, x0 : x0 + 100] = stacked_scene.astype(np.float32)
        canvas[0:100, x0 : x0 + 100] = stacked_resid.astype(np.float32)

        if not good_flags.get(star, True):
            draw_box(canvas[0:100, x0 : x0 + 100])

    return canvas


def write_diagnostics_csv(
    path: Path,
    nstars: int,
    exposures: Iterable[int],
    good_flags: Dict[int, bool],
    stats: Dict[Tuple[int, int], GroupStats],
    used_groups: Dict[Tuple[int, int], bool],
    exposure_rms: Dict[Tuple[int, int], float],
    star_fluxes: Dict[int, float],
    align_shifts: Dict[Tuple[int, int], Tuple[float, float]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "star",
                "good_flag",
                "star_flux",
                "exposure",
                "used_in_psf",
                "usable_group",
                "sky",
                "flux",
                "pmax",
                "core_dx",
                "core_dy",
                "align_dx",
                "align_dy",
                "rms",
            ]
        )
        for star in range(1, nstars + 1):
            for exposure in exposures:
                key = (star, exposure)
                group = stats.get(key)
                writer.writerow(
                    [
                        star,
                        int(good_flags.get(star, True)),
                        star_fluxes.get(star, 1.0),
                        exposure,
                        int(used_groups.get(key, False)),
                        int(group.usable) if group else 0,
                        group.sky if group else "",
                        group.flux if group else "",
                        group.pmax if group else "",
                        group.dx if group else "",
                        group.dy if group else "",
                        align_shifts.get(key, ("", ""))[0],
                        align_shifts.get(key, ("", ""))[1],
                        exposure_rms.get(key, ""),
                    ]
                )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a local PSF from MORIA UVP products without Fortran."
    )
    parser.add_argument(
        "workdir",
        type=Path,
        help="Filter directory, e.g. .../04.EXTRACT_PSF/F814W",
    )
    parser.add_argument(
        "--uvp",
        type=Path,
        default=None,
        help="Explicit UVP/UVP.GZ path. Defaults to img2extract_wfc3uv_psflist_<base>.uvp(.gz).",
    )
    parser.add_argument(
        "--base-name",
        default=None,
        help="Base tag such as simst or simstV. Inferred from workdir if omitted.",
    )
    parser.add_argument(
        "--good-list",
        type=Path,
        default=None,
        help="Path to IN.good_psf_list-style file. Defaults to IN.good_psf_list in the workdir.",
    )
    parser.add_argument(
        "--pass-index",
        type=int,
        default=1,
        help="Diagnostic pass index used in numbered output copies.",
    )
    parser.add_argument(
        "--output-prefix",
        default="",
        help="Optional prefix for all output files, e.g. py_.",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=6,
        help="Number of PSF refinement iterations.",
    )
    parser.add_argument(
        "--recenter-frac",
        type=float,
        default=0.5,
        help="Use core pixels above this fraction of the group peak when recentring.",
    )
    parser.add_argument(
        "--min-flux",
        type=float,
        default=1.0,
        help="Minimum positive core flux required to use a star/exposure in the PSF build.",
    )
    parser.add_argument(
        "--reject-sigma",
        type=float,
        default=3.0,
        help="Sigma threshold for exposure-level residual rejection between iterations.",
    )
    parser.add_argument(
        "--write-default-good-list",
        action="store_true",
        help="Write IN.good_psf_list.1-style defaults (star 1 off, all others on) and exit.",
    )
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    base_name = infer_base_name(workdir, args.base_name)
    uvp_path = infer_uvp_path(workdir, base_name, args.uvp.resolve() if args.uvp else None)
    uvp = load_uvp(uvp_path)
    nstars = int(uvp["nst"].max())
    exposures = sorted({int(x) for x in uvp["nim"]})

    if args.write_default_good_list:
        target = workdir / "IN.good_psf_list.1"
        write_default_good_list(target, nstars)
        print(f"Wrote default good-list to {target}")
        return

    good_list_path = args.good_list.resolve() if args.good_list else workdir / "IN.good_psf_list"
    good_flags = read_good_list(good_list_path, nstars)
    groups = star_exposure_groups(uvp["nst"], uvp["nim"])
    stats, psub, du_corr, dv_corr, r_corr = build_group_stats(
        uvp,
        groups,
        recenter_frac=args.recenter_frac,
        min_flux=args.min_flux,
    )
    star_fluxes = compute_star_fluxes(stats, nstars, len(exposures))

    final_psf, iterations, used_groups, exposure_rms, align_shifts = iterate_psf(
        groups,
        stats,
        good_flags,
        star_fluxes,
        psub,
        du_corr,
        dv_corr,
        iterations=args.iterations,
        reject_sigma=args.reject_sigma,
    )

    show_psf = build_show_psf(iterations)
    show_str = build_show_str(
        uvp,
        groups,
        stats,
        good_flags,
        used_groups,
        star_fluxes,
        psub,
        du_corr,
        dv_corr,
        final_psf,
        align_shifts,
    )

    prefix = args.output_prefix
    pass_suffix = str(args.pass_index) if args.pass_index is not None else ""

    show_psf_base = workdir / f"{prefix}show_psf_{base_name}.fits"
    show_psf_pass = workdir / f"{prefix}show_psf_{base_name}{pass_suffix}.fits"
    show_str_base = workdir / f"{prefix}show_str_{base_name}.fits"
    show_str_pass = workdir / f"{prefix}show_str_{base_name}{pass_suffix}.fits"
    psfout_base = workdir / f"{prefix}psfout_{base_name}.fits"
    psfout_pass = workdir / f"{prefix}psfout_{base_name}{pass_suffix}.fits"
    diag_csv = workdir / f"{prefix}psf_diagnostics_{base_name}.csv"

    fits.writeto(show_psf_base, show_psf.astype(np.float32), overwrite=True)
    fits.writeto(show_psf_pass, show_psf.astype(np.float32), overwrite=True)
    fits.writeto(show_str_base, show_str.astype(np.float32), overwrite=True)
    fits.writeto(show_str_pass, show_str.astype(np.float32), overwrite=True)
    fits.writeto(psfout_base, final_psf.astype(np.float32), overwrite=True)
    fits.writeto(psfout_pass, final_psf.astype(np.float32), overwrite=True)

    write_diagnostics_csv(
        diag_csv,
        nstars,
        exposures,
        good_flags,
        stats,
        used_groups,
        exposure_rms,
        star_fluxes,
        align_shifts,
    )

    print(f"UVP file:         {uvp_path}")
    print(f"Good list:        {good_list_path}")
    print(f"Base name:        {base_name}")
    print(f"Show PSF:         {show_psf_base.name}")
    print(f"Show STR:         {show_str_base.name}")
    print(f"Final PSF:        {psfout_base.name}")
    print(f"Diagnostics CSV:  {diag_csv.name}")
    print(f"Stars:            {nstars}")
    print(f"Exposures:        {len(exposures)}")
    print(f"Used groups:      {sum(used_groups.values())}")


if __name__ == "__main__":
    main()
