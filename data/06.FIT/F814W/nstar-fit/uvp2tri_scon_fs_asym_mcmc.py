#!/usr/bin/env python3
"""
N-star PSF fit with optional separation constraint between stars 1–2 and MCMC.
"""

from __future__ import annotations
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

import numpy as np

try:
    from astropy.io import fits
except ImportError as e:
    raise SystemExit("Install astropy.") from e


def lenc(s: str) -> int:
    s = s.rstrip()
    return len(s)


DEFAULT_INFILE = "IN.uvp2pri_NOscon_cs_asym_mcmc"
DEFAULT_PSF_INFILE = "psfout_simst.fits"
MAX_STARS = 20


@dataclass
class InfileParams:
    """Batch / IN-file inputs matching Fortran interactive reads.

    After optional line 1 (``N`` or ``NSTAR N``), the block is the same order as
    ``read(5,*)`` in the Fortran file:

    * PSF path, ``tagin``, ``tag`` (output tag for ``uvp2tri_scon_fsky_<tag>.*``)
    * ``sep12x, sig12x, sep12y, sig12y`` — constraint on stars 1–2 separation
    * ``x1 y1 … xN yN f1 … f_{N-1}`` (see ``_parse_star_geometry``)
    * ``dr1 … drN, df1 … df_{N-1}``
    * fudge, ``Nmcmc``, ``dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut``
    """

    psffile: str
    tagin: str
    tag: str
    sep12x: float
    sig12x: float
    sep12y: float
    sig12y: float
    nstar: int
    xy_init: tuple[float, ...]
    f_prefix_init: tuple[float, ...]
    dr: tuple[float, ...]
    df: tuple[float, ...]
    fudge: float
    nmcmc: int
    dufitmn: float
    dufitmx: float
    dvfitmn: float
    dvfitmx: float
    chi2cut: float


def _noncomment_lines(stream: Iterable[str]) -> list[str]:
    out: list[str] = []
    for line in stream:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        out.append(s)
    return out


def _parse_star_geometry(vals: list[float]) -> tuple[int, tuple[float, ...], tuple[float, ...]]:
    """Return N, flat xy (2N values), flux prefix (N-1 values).

    Generalizes Fortran::

        enter initial x1, y1, x2, y2, x3, y3, f1, f2 values:
        If abs(f1+f2-1.d0) < 1.d-6, source 3 flux f3 == 0

    to N stars: ``3N-1`` numbers with the last flux implied as ``1 - sum(f_prefix)``.
    """
    L = len(vals)
    if L < 2:
        raise ValueError("Line 5: need at least x1 y1 for one star.")
    if (L + 1) % 3 != 0:
        raise ValueError(
            "Line 5: expected 3N−1 numbers "
            f"got {L} numbers."
        )
    nstar = (L + 1) // 3
    if nstar > MAX_STARS:
        raise ValueError(f"N={nstar} exceeds MAX_STARS={MAX_STARS}")
    if nstar < 1:
        raise ValueError("Invalid star count.")
    n_xy = 2 * nstar
    xy_init = tuple(vals[:n_xy])
    f_prefix_init = tuple(vals[n_xy:])
    assert len(f_prefix_init) == nstar - 1
    return nstar, xy_init, f_prefix_init


def _leading_integer_star_count(line: str) -> int | None:
    """
    If ``line`` opens with an integer N (e.g. ``3`` or ``3 stars``), return N.

    Returns None when line 1 should be treated as legacy PSF path (e.g. ends with
    ``.fits`` or contains ``/``).
    """
    s = line.strip()
    if not s:
        return None
    low = s.lower()
    if "/" in s or "\\" in low:
        return None
    if low.endswith(".fits") or low.endswith(".fit"):
        return None
    parts = s.split()
    if not parts:
        return None
    try:
        n = int(float(parts[0]))
    except ValueError:
        return None
    if 1 <= n <= MAX_STARS:
        return n
    return None


def load_infile_params(stream: TextIO) -> InfileParams:
    """Parse stdin / ``IN.*``: 9-line block; optional line 1 = N or ``NSTAR N`` (PSF follows).

    Fortran equivalents (after PSF / tags)::

        enter constraint parameters sep12x, sig12x, sep12y, sig12y:
        enter initial x1, y1, … f1, f2 values:
        enter dr1, dr2, … df1, df2:
        enter error bar fudge factor:
        enter Nmcmc:
        enter dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut:
    """
    rows = _noncomment_lines(stream)
    row0_parts = rows[0].split()
    nstar_declared: int | None = None
    base = 0
    if row0_parts and row0_parts[0].upper() == "NSTAR":
        if len(row0_parts) < 2:
            raise ValueError('First line starts with NSTAR')
        nstar_declared = int(float(row0_parts[1]))
        base = 1
        if len(rows) < 9 + base:
            raise ValueError(f"IN file needs at least {9 + base} non-comment lines after NSTAR")
    else:
        n_lead = _leading_integer_star_count(rows[0])
        if n_lead is not None:
            nstar_declared = n_lead
            base = 1
            if len(rows) < 10:
                raise ValueError()
        elif len(rows) < 9:
            raise ValueError()

    psffile = rows[base + 0].strip()
    tagin = rows[base + 1].strip()
    tag = rows[base + 2].strip()
    s4 = rows[base + 3].split()
    if len(s4) != 4:
        raise ValueError(f"Separation line needs 4 numbers, got: {rows[base + 3]!r}")
    sep12x, sig12x, sep12y, sig12y = map(float, s4)
    g5 = list(map(float, rows[base + 4].split()))
    nstar, xy_init, f_prefix_init = _parse_star_geometry(g5)
    if nstar_declared is not None and nstar_declared != nstar:
        raise ValueError()
    g6 = list(map(float, rows[base + 5].split()))
    need6 = 2 * nstar - 1
    if len(g6) != need6:
        raise ValueError(f"Line 6: expected {need6} numbers (N dr + N−1 df), got {len(g6)}")
    dr = tuple(g6[:nstar])
    df = tuple(g6[nstar:])
    fudge = float(rows[base + 6].split()[0])
    nmcmc = int(float(rows[base + 7].split()[0]))
    s9 = rows[base + 8].split()
    if len(s9) != 5:
        raise ValueError(f"Window/chi2cut line needs 5 numbers, got: {rows[base + 8]!r}")
    dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut = map(float, s9)

    return InfileParams(
        psffile=psffile,
        tagin=tagin,
        tag=tag,
        sep12x=sep12x,
        sig12x=sig12x,
        sep12y=sep12y,
        sig12y=sig12y,
        nstar=nstar,
        xy_init=xy_init,
        f_prefix_init=f_prefix_init,
        dr=dr,
        df=df,
        fudge=fudge,
        nmcmc=nmcmc,
        dufitmn=dufitmn,
        dufitmx=dufitmx,
        dvfitmn=dvfitmn,
        dvfitmx=dvfitmx,
        chi2cut=chi2cut,
    )


def freeze_last_star(f_prefix: np.ndarray) -> bool:
    """True if implied last flux is negligible (freeze last star position).

    Matches Fortran ``ifit_f3`` when ``abs(f1+f2-1.d0) <= 1.d-6`` for three stars:
    the last source flux is zero and its position is not updated in MCMC
    (``x3 = x3o``, ``y3 = y3o`` in the 3-star code).
    """
    tail = 1.0 - float(np.sum(f_prefix))
    return tail <= 1e-6


def n_flux_work(nstar: int, freeze_last: bool) -> int:
    """Count of independent flux parameters in the MCMC state (``df`` slice length)."""
    if nstar == 1:
        return 0
    return (nstar - 2) if freeze_last else (nstar - 1)


def f_work_from_prefix(f_prefix: np.ndarray, nstar: int, freeze_last: bool) -> np.ndarray:
    """Reduced flux vector used in MCMC proposals."""
    if nstar == 1:
        return np.zeros(0, dtype=np.float64)
    if freeze_last:
        return np.asarray(f_prefix[: nstar - 2], dtype=np.float64).copy()
    return np.asarray(f_prefix, dtype=np.float64).copy()


def assemble_full_flux(f_work: np.ndarray, nstar: int, freeze_last: bool) -> np.ndarray:
    """Length-N fractions summing to 1."""
    f = np.zeros(nstar, dtype=np.float64)
    if nstar == 1:
        f[0] = 1.0
        return f
    if freeze_last:
        f[: max(0, nstar - 2)] = f_work
        if nstar >= 2:
            f[nstar - 2] = 1.0 - float(np.sum(f_work))
        f[nstar - 1] = 0.0
        return f
    f[: nstar - 1] = f_work
    f[nstar - 1] = 1.0 - float(np.sum(f_work))
    return f


def barsig(xlist: np.ndarray) -> tuple[float, float, int]:
    """Iterative robust mean / scatter (Fortran ``subroutine barsig``).
    """
    xlist = np.asarray(xlist, dtype=np.float64)
    ntot = int(xlist.shape[0])
    bar = 0.0
    sig = 9e9
    nsum = 0
    for _nit in range(30):
        bsum = 0.0
        ssum = 0.0
        nsum = 0
        for n in range(ntot):
            if abs(float(xlist[n]) - bar) <= 2.25 * sig:
                bsum += float(xlist[n])
                ssum += abs(float(xlist[n]) - bar)
                nsum += 1
        if nsum > 0:
            bar = bsum / nsum
        if nsum > 1:
            sig = ssum / (nsum - 1)
    if nsum <= 1:
        sig = 0.999
    return float(bar), float(sig), int(nsum)


def rpsf_phot(x: float, y: float, psf: np.ndarray) -> float:
    """PSF fraction at offset (x, y) from a star center (Fortran ``rpsf_phot``).

    This is the general function that evaluates the PSF for a given offset from
    the center.  I
    
    f you have a pixel at (ix, iy) that is (dx, dy) from the
    center of a star, this routine returns what fraction of the light should
    fall in that pixel.

    ``psf`` is shape (101, 101); sampling uses ``rx = 51 + x*4``, ``ry = 51 + y*4``.
    Returns 0 beyond ``sqrt(x^2+y^2) > 12``.
    """
    rx = 51.0 + x * 4.0
    ry = 51.0 + y * 4.0
    ix = int(rx)
    iy = int(ry)
    fx = rx - ix
    fy = ry - iy

    dd = np.sqrt(x * x + y * y)
    if dd > 12.0:
        return 0.0

    # Fortran 1-based indexing into psf(ix, iy); Python 0-based.
    # ix, iy from int(rx), int(ry) — use ix0 = ix - 1 for corner cell (ix,iy).
    ix0 = ix - 1
    iy0 = iy - 1
    if dd > 4.0:
        return float(
            (1 - fx) * (1 - fy) * psf[ix0, iy0]
            + fx * (1 - fy) * psf[ix0 + 1, iy0]
            + (1 - fx) * fy * psf[ix0, iy0 + 1]
            + fx * fy * psf[ix0 + 1, iy0 + 1]
        )

    def pf(ii: int, jj: int) -> float:
        return float(psf[ii, jj])

    A1 = pf(ix0, iy0)
    B1 = (pf(ix0 + 1, iy0) - pf(ix0 - 1, iy0)) / 2
    C1 = (pf(ix0, iy0 + 1) - pf(ix0, iy0 - 1)) / 2
    D1 = (pf(ix0 + 1, iy0) + pf(ix0 - 1, iy0) - 2 * A1) / 2
    F1 = (pf(ix0, iy0 + 1) + pf(ix0, iy0 - 1) - 2 * A1) / 2
    E1 = pf(ix0 + 1, iy0 + 1) - A1

    A2 = pf(ix0 + 1, iy0)
    B2 = (pf(ix0 + 2, iy0) - pf(ix0, iy0)) / 2
    C2 = (pf(ix0 + 1, iy0 + 1) - pf(ix0 + 1, iy0 - 1)) / 2
    D2 = (pf(ix0 + 2, iy0) + pf(ix0, iy0) - 2 * A2) / 2
    F2 = (pf(ix0 + 1, iy0 + 1) + pf(ix0 + 1, iy0 - 1) - 2 * A2) / 2
    E2 = -(pf(ix0, iy0 + 1) - A2)

    A3 = pf(ix0, iy0 + 1)
    B3 = (pf(ix0 + 1, iy0 + 1) - pf(ix0 - 1, iy0 + 1)) / 2
    C3 = (pf(ix0, iy0 + 2) - pf(ix0, iy0)) / 2
    D3 = (pf(ix0 + 1, iy0 + 1) + pf(ix0 - 1, iy0 + 1) - 2 * A3) / 2
    F3 = (pf(ix0, iy0 + 2) + pf(ix0, iy0) - 2 * A3) / 2
    E3 = -(pf(ix0 + 1, iy0) - A3)

    A4 = pf(ix0 + 1, iy0 + 1)
    B4 = (pf(ix0 + 2, iy0 + 1) - pf(ix0, iy0 + 1)) / 2
    C4 = (pf(ix0 + 1, iy0 + 2) - pf(ix0 + 1, iy0)) / 2
    D4 = (pf(ix0 + 2, iy0 + 1) + pf(ix0, iy0 + 1) - 2 * A4) / 2
    F4 = (pf(ix0 + 1, iy0 + 2) + pf(ix0 + 1, iy0) - 2 * A4) / 2
    E4 = pf(ix0, iy0) - A4

    V1 = A1 + B1 * fx + C1 * fy + D1 * fx**2 + E1 * fx * fy + F1 * fy**2
    V2 = (
        A2
        + B2 * (fx - 1)
        + C2 * fy
        + D2 * (fx - 1) ** 2
        + E2 * (fx - 1) * fy
        + F2 * fy**2
    )
    V3 = (
        A3
        + B3 * fx
        + C3 * (fy - 1)
        + D3 * fx**2
        + E3 * fx * (fy - 1)
        + F3 * (fy - 1) ** 2
    )
    V4 = (
        A4
        + B4 * (fx - 1)
        + C4 * (fy - 1)
        + D4 * (fx - 1) ** 2
        + E4 * (fx - 1) * (fy - 1)
        + F4 * (fy - 1) ** 2
    )

    return float((1 - fx) * (1 - fy) * V1 + fx * (1 - fy) * V2 + (1 - fx) * fy * V3 + fx * fy * V4)


def read_psf_fits(path: Path) -> np.ndarray:
    """Read in the PSF model file (Fortran ``readfits_r4`` on ``psffile`` → ``psfu(101,101)``)."""
    with fits.open(path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(f"PSF FITS must be 2-D, got shape {data.shape}")
    if data.shape != (101, 101):
        # Allow (101,101) regardless of axis order warning
        if data.shape[0] == 101 and data.shape[1] == 101:
            pass
        else:
            raise ValueError(f"Expected 101×101 PSF, got {data.shape}")
    return data


def read_uvp_star(path: Path, star_number: int = 1) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read ``img2extract_wfc3uv_psflist_<tag>.uvp`` for one NEARBY_SIM star.

    First 12 fields in each row (Fortran ``read(63,*)`` on the ``.uvp`` file):

    (1) u position in the reference frame (u,v are ref-frame coords)
    (2) v position in the reference frame
    (3) p, the pixel value for the pixel in question
    (4) exposure number
    (5) i coordinate of pixel
    (6) j coordinate of pixel
    (7) x-coordinate of distortion-corrected location of the pixel
    (8) y-coordinate of distortion-corrected location of the pixel
    (9)  delta-u (offset from the center of this object)
    (10) delta-v (offset from the center of this object)
    (11) sqrt[(delta-u)^2+(delta-v)^2]
    (12) star number in the NEARBY_SIM list (ignore last column if present)

    This selects the target star: for photometry of other stars pick a different
    ``star_number`` (Fortran keeps rows with ``nstl(Ls+1).ne.1`` → ``goto 1``).
    """
    ul, vl, pl = [], [], []
    niml, dul, dvl = [], [], []
    nstl, rl = [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            u = float(parts[0])
            v = float(parts[1])
            p = float(parts[2])
            nim = int(parts[3])
            du = float(parts[8])
            dv = float(parts[9])
            r = float(parts[10])
            nst = int(parts[11])
            if nst != star_number:
                continue
            ul.append(u)
            vl.append(v)
            pl.append(p)
            niml.append(nim)
            dul.append(du)
            dvl.append(dv)
            rl.append(r)
            nstl.append(nst)
    if not ul:
        raise FileNotFoundError(f"No pixels with star number {star_number} in {path}")
    return (
        np.asarray(ul),
        np.asarray(vl),
        np.asarray(pl, dtype=np.float64),
        np.asarray(niml, dtype=np.int32),
        np.asarray(dul, dtype=np.float64),
        np.asarray(dvl, dtype=np.float64),
        np.asarray(rl, dtype=np.float64),
        np.asarray(nstl, dtype=np.int32),
    )


def solve_flux_sky(
    fu: np.ndarray,
    pp_pskyu: np.ndarray,
    sgu: np.ndarray,
    nimu: np.ndarray,
    nims: int,
) -> np.ndarray:
    """Simultaneous fit sky method (Fortran ``SVDCMP`` / ``SVBKSB`` block).
    """
    us = fu.shape[0]
    nim_maxp1 = nims + 1
    sg = sgu.astype(np.float64)
    X = np.zeros((us, nim_maxp1), dtype=np.float64)
    X[:, 0] = fu.astype(np.float64) / sg
    for u in range(us):
        nim = int(nimu[u])  # 1-based exposure index as in Fortran
        col = nim  # sky column index matching Fortran Aa(nimu+1) -> Python Aa[nim]
        if col < 1 or col > nims:
            raise ValueError(f"Bad exposure index nimu={nim} for NIMs={nims}")
        X[u, col] += 1.0 / sg[u]
    bb = pp_pskyu.astype(np.float64) / sg
    coef, *_ = np.linalg.lstsq(X, bb, rcond=None)
    return coef.astype(np.float64)


def chi2_tot(
    coef: np.ndarray,
    fu: np.ndarray,
    pp_pskyu: np.ndarray,
    sgu: np.ndarray,
    nimu: np.ndarray,
    chi2norm: float,
    sep12x: float,
    sig12x: float,
    sep12y: float,
    sig12y: float,
    xys: np.ndarray,
) -> tuple[float, float]:
    """Calculate chi^2 (data term / ``chi2norm`` + optional separation constraint).
    """
    sg = sgu.astype(np.float64)
    c0 = coef[0]
    s = 0.0
    for u in range(fu.shape[0]):
        nim = int(nimu[u])
        sky = coef[nim]
        sm = c0 * float(fu[u]) + sky
        d = (float(pp_pskyu[u]) - sm) / sg[u]
        s += d * d
    s = s / chi2norm
    chi2scon = 0.0
    if xys.shape[0] >= 2:
        chi2scon = ((sep12x - (xys[0, 0] - xys[1, 0])) / sig12x) ** 2 + (
            (sep12y - (xys[0, 1] - xys[1, 1])) / sig12y
        ) ** 2
    return s + chi2scon, chi2scon


def model_fu(
    duu: np.ndarray,
    dvu: np.ndarray,
    psf: np.ndarray,
    xys: np.ndarray,
    f_full: np.ndarray,
) -> np.ndarray:
    """Normalized PSF model at each used pixel (Fortran ``fu(U)`` loop before sky fit).
    """
    nstar = xys.shape[0]
    out = np.zeros(duu.shape[0], dtype=np.float64)
    for u in range(duu.shape[0]):
        du = float(duu[u])
        dv = float(dvu[u])
        acc = 0.0
        for k in range(nstar):
            acc += float(f_full[k]) * rpsf_phot(
                du - float(xys[k, 0]), dv - float(xys[k, 1]), psf
            )
        out[u] = acc
    return out


def propose_flux(
    f_work_o: np.ndarray,
    df: np.ndarray,
    rng: np.random.Generator,
    nstar: int,
    freeze_last: bool,
    max_try: int = 400,
) -> np.ndarray | None:
    """
    Reject when any fraction is negative or the implied sum exceeds 1
    (Fortran: ``if((f1+f2.gt.1.0).or.(f1.lt.0.).or.(f2.lt.0.)) go to 22``).
    Returns ``None`` if no valid draw after ``max_try``.
    """
    nw = f_work_o.shape[0]
    if nw == 0:
        return f_work_o.copy()
    dfv = np.asarray(df[:nw], dtype=np.float64)
    for _ in range(max_try):
        cand = f_work_o + dfv * (2.0 * rng.random(nw) - 1.0)
        if freeze_last:
            penult = 1.0 - float(np.sum(cand))
            if np.any(cand < 0.0) or penult < 0.0:
                continue
            return cand
        tail = 1.0 - float(np.sum(cand))
        if np.any(cand < 0.0) or tail < 0.0:
            continue
        return cand
    return None


def mcmc_ascii_header(nstar: int, nims: int) -> str:
    """Header lines for ``.07.mcmc`` / ``.04.probe_fit`` (read by ``mcmc_expand_average``).
    Documents ``NSTAR`` and ``NIMS`` for variable-N chains; Fortran format 993 had no header.
    """
    return f"# NSTAR {nstar}\n# NIMS {nims}\n"


def format_mcmc_row(
    nrep: int,
    xys: np.ndarray,
    f_full: np.ndarray,
    chi2tot_: float,
    aa1: float,
    skies_: list[float],
) -> str:
    """One MCMC/probe line: **i4** + **(3N−1)×f10.5** + **2f13.4** + skies (f10.4).
    """
    parts_xy = "".join(
        f"{float(xys[k, 0]):10.5f}{float(xys[k, 1]):10.5f}" for k in range(xys.shape[0])
    )
    nstar = int(xys.shape[0])
    n_flux_prefix = max(0, nstar - 1)
    parts_f = "".join(
        f"{float(f_full[k]):10.5f}" for k in range(n_flux_prefix)
    )
    tail = "".join(f"{s:10.4f}" for s in skies_)
    return (
        f"{nrep:4d}"
        f"{parts_xy}"
        f"{parts_f}"
        f"{chi2tot_:13.4f}{aa1:13.4f}"
        f"{tail}\n"
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="N-star MCMC PSF fit.")
    p.add_argument(
        "infile",
        nargs="?",
        default=None,
        help=(
            f"Fortran-format parameter file (default: ./{DEFAULT_INFILE} if present). "
            "Use '-' to read stdin."
        ),
    )
    p.add_argument("--psf", help="PSF FITS (101×101 float)")
    p.add_argument("--tagin", help="Input tag for img2extract_wfc3uv_psflist_<tagin>.uvp")
    p.add_argument("--tag", help="Output tag for product filenames")
    p.add_argument("--star", type=int, default=1, help="Star index in .uvp last column (default 1)")
    p.add_argument("--sep12x", type=float)
    p.add_argument("--sig12x", type=float)
    p.add_argument("--sep12y", type=float)
    p.add_argument("--sig12y", type=float)
    p.add_argument(
        "--xyf",
        nargs="*",
        type=float,
        metavar="XYF",
        help="Line 5 numbers: x1 y1 … xN yN f1 … f_{N−1}",
    )
    p.add_argument(
        "--drdf",
        nargs="*",
        type=float,
        metavar="DRDF",
        help="Line 6 numbers: dr1…drN df1…df_{N−1}",
    )
    p.add_argument("--fudge", type=float, default=None, help="Override IN-file fudge")
    p.add_argument("--Nmcmc", type=int)
    p.add_argument("--duv", nargs=5, metavar=("dufitmn", "dufitmx", "dvfitmn", "dvfitmx", "chi2cut"), type=float)
    p.add_argument("--seed", type=int, default=45780)
    p.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="Prompt for all inputs on a TTY; do not read the default IN file.",
    )
    return p


def prompt(msg: str) -> str:
    sys.stdout.write(msg)
    sys.stdout.flush()
    return sys.stdin.readline().rstrip("\n")


def main(argv: list[str] | None = None) -> None:
    args_ns = build_parser().parse_args(argv)

    loaded: InfileParams | None = None
    if args_ns.interactive:
        if not sys.stdin.isatty():
            raise SystemExit("--interactive requires an interactive terminal (stdin TTY).")
        loaded = None
    elif args_ns.infile == "-":
        loaded = load_infile_params(sys.stdin)
    elif args_ns.infile is not None:
        p = Path(args_ns.infile)
        with open(p, encoding="utf-8", errors="replace") as fh:
            loaded = load_infile_params(fh)
    elif not sys.stdin.isatty():
        loaded = load_infile_params(sys.stdin)
    else:
        default_path = Path(DEFAULT_INFILE)
        if default_path.is_file():
            with default_path.open(encoding="utf-8", errors="replace") as fh:
                loaded = load_infile_params(fh)

    full_cli = (
        args_ns.psf
        and args_ns.tagin
        and args_ns.tag
        and args_ns.sep12x is not None
        and args_ns.sig12x is not None
        and args_ns.sep12y is not None
        and args_ns.sig12y is not None
        and args_ns.xyf is not None
        and len(args_ns.xyf) >= 2
        and args_ns.drdf is not None
        and len(args_ns.drdf) >= 1
        and args_ns.Nmcmc is not None
        and args_ns.duv is not None
    )

    interactive = loaded is None and not full_cli

    star = args_ns.star
    seed = args_ns.seed

    geom_vals: list[float]
    drdf_vals: list[float]

    if interactive:
        sys.stderr.write(
            "This builds the same columns as an IN file whose first line is N "
            f"(10 non-comment lines total; geometry is line 6).\n\n"
        )
        n_in = int(
            float(prompt(f"Line 1 — N (number of stars) [1–{MAX_STARS}]: ").strip())
        )
        if n_in < 1 or n_in > MAX_STARS:
            raise SystemExit(f"N must be 1…{MAX_STARS}")
        psffile = (
            prompt(f"Line 2 — PSF [{DEFAULT_PSF_INFILE}]: ").strip() or DEFAULT_PSF_INFILE
        )
        tagin = prompt(
            "Line 3 — tag for img2extract_wfc3uv_psflist_<tag>.uvp: ",
        ).strip()
        tag = prompt("Line 4 — output tag for uvp2tri_scon_fsky_<tag>.*: ").strip()
        # enter constraint parameters sep12x, sig12x, sep12y, sig12y
        sys.stderr.write(
            "\nLine 5 — four numbers before the star positions: a soft prior on the "
            "separation of stars 1 and 2 in du and dv (not your star x,y coordinates). "
            "Order: target Δu(2−1), σ for that in u, target Δv(2−1), σ for that in v. "
            "To weaken the prior, use huge σ, e.g.: 0 99999 0 99999\n"
        )
        sep12x, sig12x, sep12y, sig12y = map(
            float,
            prompt("Line 5 — four numbers (space-separated): ").split(),
        )

        need_geom = 3 * n_in - 1
        if n_in == 1:
            sys.stderr.write(
                "\nLine 6 — one line with 2 numbers: x1 y1 (same layout as IN geometry row).\n"
            )
        else:
            # enter initial x1, y1, … f1, f2 values; last flux = 1 - sum(prefix)
            sys.stderr.write(
                f"\nLine 6 — one line with {need_geom} numbers, same as in your IN file: "
                f"x1 y1 … x{n_in} y{n_in}  f1 … f{n_in - 1} "
                "(last flux omitted; it is 1 − sum of these). "
                "If abs(sum(prefix)-1) < 1e-6, last star flux is 0 (Fortran ifit_f3).\n"
            )
        geom_vals = list(
            map(float, prompt("Line 6 — paste geometry line: ").split())
        )

        need_dr = 2 * n_in - 1
        # enter dr1, dr2, … df1, df2
        sys.stderr.write(
            f"\nLine 7 — one line with {need_dr} numbers: dr1 … dr{n_in}, "
            f"then  df1 … df_{max(0, n_in - 1)}  (positional step sizes, then flux "
            "step sizes).\n"
        )
        drdf_vals = list(
            map(float, prompt("Line 7 — paste dr/df line: ").split())
        )
        if len(geom_vals) != need_geom:
            raise SystemExit(
                f"Line 6: expected {need_geom} numbers for N={n_in}; got {len(geom_vals)}"
            )
        if len(drdf_vals) != need_dr:
            raise SystemExit(
                f"Line 7: expected {need_dr} numbers for N={n_in}; got {len(drdf_vals)}"
            )

        fudge = float(prompt("Line 8 — fudge: ").strip())
        Nmcmc = int(float(prompt("Line 9 — Nmcmc: ").strip()))
        dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut = map(
            float,
            prompt(
                "Line 10 — dufitmn dufitmx dvfitmn dvfitmx chi2cut: ",
            ).split(),
        )

        outp = prompt(
            f"Save IN file path (Enter = {DEFAULT_INFILE}): ",
        ).strip()
        in_path = Path(outp or DEFAULT_INFILE)
        lines = (
            [str(n_in), psffile, tagin, tag]
            + [f"{sep12x} {sig12x} {sep12y} {sig12y}"]
            + [" ".join(f"{v:g}" for v in geom_vals)]
            + [" ".join(f"{v:g}" for v in drdf_vals)]
            + [f"{fudge:g}", str(Nmcmc)]
            + [f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"]
        )
        in_path.parent.mkdir(parents=True, exist_ok=True)
        in_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        sys.stderr.write(f"Saved IN parameter file → {in_path.resolve()}\n")
    elif loaded is not None:
        lp = loaded
        psffile = args_ns.psf or lp.psffile
        tagin = args_ns.tagin or lp.tagin
        tag = args_ns.tag or lp.tag
        sep12x = lp.sep12x if args_ns.sep12x is None else args_ns.sep12x
        sig12x = lp.sig12x if args_ns.sig12x is None else args_ns.sig12x
        sep12y = lp.sep12y if args_ns.sep12y is None else args_ns.sep12y
        sig12y = lp.sig12y if args_ns.sig12y is None else args_ns.sig12y
        if args_ns.xyf:
            geom_vals = list(args_ns.xyf)
        else:
            geom_vals = list(lp.xy_init) + list(lp.f_prefix_init)
        if args_ns.drdf:
            drdf_vals = list(args_ns.drdf)
        else:
            drdf_vals = list(lp.dr) + list(lp.df)
        fudge = args_ns.fudge if args_ns.fudge is not None else lp.fudge
        Nmcmc = lp.nmcmc if args_ns.Nmcmc is None else args_ns.Nmcmc
        if args_ns.duv is None:
            dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut = (
                lp.dufitmn,
                lp.dufitmx,
                lp.dvfitmn,
                lp.dvfitmx,
                lp.chi2cut,
            )
        else:
            dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut = args_ns.duv
    else:
        psffile = args_ns.psf
        tagin = args_ns.tagin
        tag = args_ns.tag
        sep12x = args_ns.sep12x
        sig12x = args_ns.sig12x
        sep12y = args_ns.sep12y
        sig12y = args_ns.sig12y
        geom_vals = list(args_ns.xyf)
        drdf_vals = list(args_ns.drdf)
        fudge = args_ns.fudge if args_ns.fudge is not None else 1.0
        Nmcmc = args_ns.Nmcmc
        dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut = args_ns.duv

    nstar, _xy_flat, _f_pref_t = _parse_star_geometry(geom_vals)
    need_drdf = 2 * nstar - 1
    if len(drdf_vals) != need_drdf:
        raise SystemExit(
            f"Line 6: expected {need_drdf} numbers for N={nstar}, got {len(drdf_vals)}"
        )
    dr_arr = np.asarray(drdf_vals[:nstar], dtype=np.float64)
    df_arr = np.asarray(drdf_vals[nstar:], dtype=np.float64)
    xy_in = np.asarray(_xy_flat, dtype=np.float64).reshape(nstar, 2)
    f_prefix_in = np.asarray(_f_pref_t, dtype=np.float64)
    freeze_last = freeze_last_star(f_prefix_in)
    f_work_in = f_work_from_prefix(f_prefix_in, nstar, freeze_last)
    nw = n_flux_work(nstar, freeze_last)
    if df_arr.shape[0] < nstar - 1:
        raise SystemExit("Internal error: df array too short.")
    print(f"NSTAR fit: N={nstar} (freeze_last_star={freeze_last}, n_flux_work={nw})")

    lpsff = lenc(psffile)
    ltag = lenc(tag)
    ltagin = lenc(tagin)
    rng = np.random.default_rng(seed)

    pixoutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.01.pix_all"
    rmpixfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.08.rm_pix"
    infile = f"img2extract_wfc3uv_psflist_{tagin[:ltagin]}.uvp"

    ul, vl, pl, niml, dul, dvl, rl, nstl = read_uvp_star(Path(infile), star_number=star)
    ls = ul.shape[0]
    nims = int(np.max(niml))
    nsts = int(np.max(nstl))

    with open(pixoutfile, "w", encoding="utf-8") as pf:
        for ell in range(ls):
            pf.write(
                f"{ell + 1:6d} {dul[ell]:7.3f} {dvl[ell]:7.3f} {pl[ell]:9.1f} "
                f"{niml[ell]:4d} {nstl[ell]:4d}\n"
            )

    print(f"NSTs: {nsts}")
    print(f"NIMs: {nims}")
    print(f"  Ls: {ls}")

    psf = read_psf_fits(Path(psffile))
    print(f"READIN {psffile[:lpsff]}")
    print(f"psfu(51,51): {float(psf[50, 50])}")

    pp_psky = pl.copy()

    puseoutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.03.pix_use"
    du_arr = []
    dv_arr = []
    nim_use = []
    nst_use = []
    pp_use = []
    sg_use = []

    # Pixels in fitting window: dufitmn < du < dufitmx, dvfitmn < dv < dvfitmx
    with open(puseoutfile, "w", encoding="utf-8") as uf:
        us_count = 0
        for ell in range(ls):
            du = float(dul[ell])
            dv = float(dvl[ell])
            if du > dufitmn and du < dufitmx and dv > dvfitmn and dv < dvfitmx:
                us_count += 1
                nim_use.append(int(niml[ell]))
                nst_use.append(int(nstl[ell]))
                du_arr.append(du)
                dv_arr.append(dv)
                # Fortran: sgu = sqrt(16 + max(pl,0) + (0.01*max(pl,0))^2)
                sgu = np.sqrt(
                    16.0 + max(pl[ell], 0.0) + (0.01 * max(pl[ell], 0.0)) ** 2
                )
                sg_use.append(sgu)
                pp_use.append(float(pp_psky[ell]))
                uf.write(
                    f"{ell + 1:6d}{us_count:4d}"
                    f"{du:10.4f}{dv:10.4f}{pp_use[-1]:10.1f}"
                    f"{nim_use[-1]:4d}{nst_use[-1]:4d}\n"
                )

    nimu = np.asarray(nim_use, dtype=np.int32)
    duu = np.asarray(du_arr, dtype=np.float64)
    dvu = np.asarray(dv_arr, dtype=np.float64)
    pp_pskyu = np.asarray(pp_use, dtype=np.float64)
    sgu = np.asarray(sg_use, dtype=np.float64)
    us0 = us = int(duu.shape[0])
    print(f"Ls: {ls}")
    print(f"Us: {us}")

    if us == 0:
        raise SystemExit("No pixels in fitting window; check dufitmn/mx and dvfitmn/mx.")

    probeoutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.04.probe_fit"
    mcmcoutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.07.mcmc"

    chi2norm = 1.0
    irenorm = 0

    # Store full-pixel lists for pix_show residual map (Fortran ``pix_show`` block)
    dul_full = dul.astype(np.float64)
    dvl_full = dvl.astype(np.float64)

    n_fit_params = (2 * nstar - (2 if freeze_last else 0)) + nw + 1

    # Label 180: restart after pixel removal or chi2norm renormalization
    while True:
        xys_o = xy_in.copy()
        f_work_o = f_work_in.copy()

        # MCMC starting from this point; open .04.probe_fit, .07.mcmc, .08.rm_pix
        with open(probeoutfile, "w", encoding="utf-8") as f24, open(
            mcmcoutfile, "w", encoding="utf-8"
        ) as f25, open(rmpixfile, "w", encoding="utf-8") as f26:
            hdr_mcmc = mcmc_ascii_header(nstar, nims)
            f24.write(hdr_mcmc)
            f25.write(hdr_mcmc)
            print(f"chi2norm = {chi2norm}")
            f_full = assemble_full_flux(f_work_o, nstar, freeze_last)
            fu = model_fu(duu, dvu, psf, xys_o, f_full)
            Aa = solve_flux_sky(fu, pp_pskyu, sgu, nimu, nims)
            chi2tot, chi2scon = chi2_tot(
                Aa,
                fu,
                pp_pskyu,
                sgu,
                nimu,
                chi2norm,
                sep12x,
                sig12x,
                sep12y,
                sig12y,
                xys_o,
            )
            chi2o = chi2tot

            chi2min = chi2tot
            chi2scmin = chi2scon
            fu_min = fu.copy()
            pp_pskyu_min = pp_pskyu.copy()
            Aamin = Aa.copy()
            xys_min = xys_o.copy()
            f_full_min = f_full.copy()

            skies = [float(Aa[j]) for j in range(1, nims + 1)]
            # write initial state to MCMC file (Fortran write(24,993) / write(25,993))
            row0 = format_mcmc_row(0, xys_o, f_full, chi2tot, float(Aa[0]), skies)
            f24.write(row0)
            f25.write(row0)

            naccept = nreject = 0
            nrepeat = 0
            idof = max(1, us - n_fit_params - nims)

            for i in range(1, Nmcmc + 1):
                xys = xys_o.copy()
                for k in range(nstar):
                    if freeze_last and k == nstar - 1:
                        continue
                    xys[k, 0] = xys_o[k, 0] + dr_arr[k] * (2.0 * rng.random() - 1.0)
                    xys[k, 1] = xys_o[k, 1] + dr_arr[k] * (2.0 * rng.random() - 1.0)

                f_work_cand = propose_flux(f_work_o, df_arr, rng, nstar, freeze_last)
                if f_work_cand is None:
                    nreject += 1
                    nrepeat += 1
                    f25.write(f"{nrepeat:8d}\n")
                    f25.flush()
                    if (i // 10000) * 10000 == i and i > 0:
                        print(f"{i} steps completed")
                    continue

                f_full = assemble_full_flux(f_work_cand, nstar, freeze_last)
                fu = model_fu(duu, dvu, psf, xys, f_full)
                Aa = solve_flux_sky(fu, pp_pskyu, sgu, nimu, nims)
                chi2tot, chi2scon = chi2_tot(
                    Aa,
                    fu,
                    pp_pskyu,
                    sgu,
                    nimu,
                    chi2norm,
                    sep12x,
                    sig12x,
                    sep12y,
                    sig12y,
                    xys,
                )

                if chi2tot < chi2min:
                    xys_min[:] = xys
                    f_full_min[:] = f_full
                    chi2min = chi2tot
                    chi2scmin = chi2scon
                    fu_min[:] = fu
                    pp_pskyu_min[:] = pp_pskyu
                    Aamin[:] = Aa

                if chi2tot <= chi2o:
                    naccept += 1
                    xys_o[:] = xys
                    f_work_o[:] = f_work_cand
                    chi2o = chi2tot
                    nrepeat = 0
                    skies = [float(Aa[j]) for j in range(1, nims + 1)]
                    f25.write(
                        format_mcmc_row(
                            nrepeat, xys, f_full, chi2tot, float(Aa[0]), skies
                        )
                    )
                    f25.flush()
                else:
                    prob_rand = rng.random()
                    prob = np.exp(-(chi2tot - chi2o) / 2.0)
                    if prob > prob_rand:
                        naccept += 1
                        xys_o[:] = xys
                        f_work_o[:] = f_work_cand
                        chi2o = chi2tot
                        nrepeat = 0
                        skies = [float(Aa[j]) for j in range(1, nims + 1)]
                        f25.write(
                            format_mcmc_row(
                                nrepeat,
                                xys,
                                f_full,
                                chi2tot,
                                float(Aa[0]),
                                skies,
                            )
                        )
                        f25.flush()
                    else:
                        nreject += 1
                        nrepeat += 1
                        f25.write(f"{nrepeat:8d}\n")
                        f25.flush()
                if (i // 10000) * 10000 == i and i > 0:
                    print(f"{i} steps completed")

            print(f"accepted, rejected MCMC steps: {naccept} {nreject}")

            # check if a pixel should be removed
            chi2max = 0.0
            ichi2 = 0
            for u in range(us):
                sm = Aamin[0] * fu_min[u] + Aamin[int(nimu[u])]
                chi2p = ((pp_pskyu_min[u] - sm) / sgu[u]) ** 2 / chi2norm
                if chi2p > chi2max:
                    ichi2 = u + 1
                    chi2max = chi2p

            # possibly remove pixel; else renormalize chi2norm once (Fortran label 180)
            restart = False
            if ichi2 > 0:
                if chi2max >= chi2cut:
                    udrop = ichi2 - 1
                    f26.write(
                        f"removing pix at img, du, dv, val, chi2 ="
                        f"{nimu[udrop]:4d}{duu[udrop]:10.4f}{dvu[udrop]:10.4f}"
                        f"{pp_pskyu[udrop]:10.2f}{chi2max:13.5g}\n"
                    )
                    keep = np.ones(us, dtype=bool)
                    keep[udrop] = False
                    nimu = nimu[keep]
                    duu = duu[keep]
                    dvu = dvu[keep]
                    pp_pskyu = pp_pskyu[keep]
                    sgu = sgu[keep]
                    us = int(duu.shape[0])
                    idof = max(1, us - n_fit_params - nims)
                    chi2norm = (chi2min - chi2max) * chi2norm / idof
                    irenorm = 1
                    restart = True
                elif irenorm == 0:
                    chi2norm = chi2min * chi2norm / idof
                    irenorm = 1
                    restart = True

            if restart:
                continue
        break

    # Last MCMC iteration model (Fortran ``Aa`` / ``fu`` at loop exit)
    fu_report = fu.copy()
    Aa_report = Aa.copy()

    hdr_xy_in = "".join(
        f"{float(xy_in[k, 0]):9.5f}{float(xy_in[k, 1]):9.5f}" for k in range(nstar)
    )
    f_full_input = assemble_full_flux(f_work_in, nstar, freeze_last)
    hdr_f_in = "".join(f"{float(f_full_input[k]):9.5f}" for k in range(nstar))
    hdr_dr = "".join(f"{float(dr_arr[k]):9.5f}" for k in range(nstar))
    hdr_df = "".join(f"{float(df_arr[k]):9.5f}" for k in range(nstar - 1))

    # Final report — same chi2 assembly as Fortran lines 653–661 (uses exit ``Aa``, ``fu``)
    finaloutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.05.final_fit"
    with open(finaloutfile, "w", encoding="utf-8") as ff:
        ff.write(
            f"# NSTAR = {nstar}  (freeze_last_star={freeze_last})\n"
            f"#input " + " ".join(f"x{k+1}   y{k+1}      " for k in range(nstar))
            + " ".join(f"f{k+1}       " for k in range(nstar))
            + "\n"
            f"#{hdr_xy_in}{hdr_f_in}\n"
            f"# " + " ".join(f"dr{k+1}    " for k in range(nstar))
            + " ".join(f"df{k+1}    " for k in range(nstar - 1))
            + "Nmcmc\n"
            f"#{hdr_dr}{hdr_df}{Nmcmc:8d}\n"
            f"#    sep12x         sig12x         sep12y         sig12y      fudge\n"
            f"#{sep12x:15.7g}{sig12x:15.7g}{sep12y:15.7g}{sig12y:15.7g}{fudge:15.7g}\n"
            f"#   dufitmn   dufitmx   dvfitmn  dvfitmx  chi2cut\n"
            f"#{dufitmn:10.5f}{dufitmx:10.5f}{dvfitmn:10.5f}{dvfitmx:10.5f}{chi2cut:10.5f}\n"
        )
        ff.write("# \n")
        ff.write("# input file:\n")
        ff.write(f"#  PSF  file ={psffile[:lpsff]}\n")
        ff.write(f"# pixel output file ={pixoutfile}\n")
        ff.write("# \n")
        ff.write("# BEST-FIT MODEL \n")
        ff.write(f"# pix0  = {us0:10d}\n")
        ff.write(f"# pix     {us:10d}\n")
        ff.write(f"# n_fit_params (for idof) = {n_fit_params}\n")
        for k in range(nstar):
            ff.write(
                f"#   X{k + 1}MIN = {float(xys_min[k, 0]):10.5f}\n"
                f"#   Y{k + 1}MIN = {float(xys_min[k, 1]):10.5f}\n"
            )
        for k in range(nstar):
            ff.write(f"#   F{k + 1}MIN = {float(f_full_min[k]):10.5f}\n")
        ff.write(f"# chi2min = {chi2min:10.4f}\n")
        ff.write(f"#chi2scmin= {chi2scmin:10.4f}\n")
        ff.write(f"#  A1_min = {Aamin[0]:15.7g}\n")
        for img in range(1, nims + 1):
            ff.write(f"# Sky({img:2d})  ={Aamin[img]:10.4f}\n")
        ff.write("# \n")
        ff.write("# INDIVIDUAL PIXELS:\n")
        ff.write("# \n")

        if nstar >= 2:
            chi2scon_f = ((sep12x - (xys_min[0, 0] - xys_min[1, 0])) / sig12x) ** 2 + (
                (sep12y - (xys_min[0, 1] - xys_min[1, 1])) / sig12y
            ) ** 2
        else:
            chi2scon_f = 0.0
        chi2tot_f = chi2scon_f
        for u in range(us):
            sm = Aa_report[0] * fu_report[u] + Aa_report[int(nimu[u])]
            chi2p = ((pp_pskyu[u] - sm) / sgu[u]) ** 2 / chi2norm
            chi2tot_f += chi2p
            ff.write(
                f"{u + 1:4d}{duu[u]:9.4f}{dvu[u]:9.4f}{pp_pskyu_min[u]:9.1f}"
                f"{fu_min[u]:9.5f}{Aa_report[0] * fu_report[u]:9.1f}{sgu[u]:9.2f}"
                f"{rpsf_phot(float(duu[u]), float(dvu[u]), psf):9.5f}{chi2p:13.5g}\n"
            )

    print("GENERATE PIX_SHOW...")
    pix_show = np.zeros((2001, 4001), dtype=np.float32)
    nbuf = max(ls, 1)
    po = np.zeros(nbuf, dtype=np.float32)
    ps_arr = np.zeros(nbuf, dtype=np.float32)
    for ii in range(1, 2002):
        if ii % 100 == 0:
            print(f"---> i: {ii}")
        for jj in range(1, 2002):
            x = 0.01 * (ii - 1001)
            y = 0.01 * (jj - 1001)
            ss = 0
            for ell in range(ls):
                if abs(float(dul_full[ell]) - x) <= 0.5 and abs(float(dvl_full[ell]) - y) <= 0.5:
                    f_mod = 0.0
                    for k in range(nstar):
                        f_mod += float(f_full_min[k]) * rpsf_phot(
                            float(dul_full[ell]) - float(xys_min[k, 0]),
                            float(dvl_full[ell]) - float(xys_min[k, 1]),
                            psf,
                        )
                    nim_ix = int(niml[ell])
                    po[ss] = float(pp_psky[ell]) - float(Aamin[nim_ix])
                    ps_arr[ss] = po[ss] - f_mod * float(Aamin[0])
                    ss += 1
            if ss > 0:
                pbar, _, _ = barsig(ps_arr[:ss])
                pbar2, _, _ = barsig(po[:ss])
            else:
                pbar = 0.0
                pbar2 = 0.0
            pix_show[jj - 1, ii - 1] = float(pbar)
            pix_show[jj - 1, ii + 2000 - 1] = float(pbar2)

    showoutfile = f"uvp2tri_scon_fsky_{tag[:ltag]}.06.pix_show.fits"
    fits.writeto(showoutfile, pix_show.astype(np.float32), overwrite=True)
    print(f"Wrote {showoutfile}")


if __name__ == "__main__":
    main()
