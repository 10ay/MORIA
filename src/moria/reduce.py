import numpy as np
import os
import glob
from pathlib import Path
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from astropy.table import Table
import pdb
import subprocess
from datetime import datetime as dt
import pandas as pd
import shutil
import sys
import matplotlib.pyplot as plt
import argparse
from importlib import resources
import urllib.request
import urllib.error
import bz2
import urllib.parse
import re


def _get_cal_matches4_field_specs():
    return (
        [("i6", 6)]
        + [("F9.3", 9)] * 2
        + [("i5", 5)] * 2
        + [("f8.3", 8)] * 5
        + [("f9.4", 9)]
        + [("f8.3", 8)] * 5
        + [("f9.4", 9)]
        + [("f8.3", 8)] * 6
    )


def _format_cal_matches4_field(spec: str, value) -> str:
    if spec == "i6":
        return f"{int(value):6d}"
    if spec == "i5":
        return f"{int(value):5d}"
    if spec == "F9.3":
        return f"{float(value):9.3f}"
    if spec == "f8.3":
        return f"{float(value):8.3f}"
    if spec == "f9.4":
        return f"{float(value):9.4f}"
    raise ValueError(f"unknown Cal_matches4 field spec: {spec}")


def _format_cal_matches4_line(values) -> str:
    specs = _get_cal_matches4_field_specs()
    return "".join(
        _format_cal_matches4_field(spec, value)
        for (spec, _width), value in zip(specs, values)
    )


def _parse_cal_matches4_field_token(spec: str, token: str):
    token = token.strip()
    if token and set(token) == {"*"}:
        return 99.999
    if spec in {"i6", "i5"}:
        return int(token)
    return float(token)


def _parse_cal_matches4_values(line: str):
    """Parse one format-999 row from fixed-width or repaired free-form text."""
    specs = _get_cal_matches4_field_specs()
    repaired = line.rstrip()
    repaired = re.sub(r"(\.\d{4})(\*{6,})", r"\1 \2", repaired)
    repaired = re.sub(r"(\.\d{4})(-)", r"\1 \2", repaired)

    expected_len = sum(width for _, width in specs)
    if "*" not in repaired and len(repaired) == expected_len:
        values = []
        pos = 0
        for spec, width in specs:
            chunk = repaired[pos : pos + width]
            values.append(_parse_cal_matches4_field_token(spec, chunk))
            pos += width
        return values

    tokens = repaired.split()
    if len(tokens) < len(specs):
        raise ValueError(f"expected {len(specs)} fields, got {len(tokens)}")
    return [
        _parse_cal_matches4_field_token(spec, token)
        for (spec, _width), token in zip(specs, tokens)
    ]


def _repair_vi_hst_cal_matches4_data_line(line: str) -> str:
    """
    Repair one VI_HST_ogle_Cal_matches4 data row while preserving Fortran
    format-999 fixed-width columns.
    """
    specs = _get_cal_matches4_field_specs()
    stripped = line.rstrip()
    if not stripped or stripped.startswith("#"):
        return stripped
    if not re.match(r"^\s*\d", stripped):
        return stripped

    expected_len = sum(width for _, width in specs)
    if "*" not in stripped and len(stripped) == expected_len:
        try:
            values = _parse_cal_matches4_values(stripped)
            if not (
                (values[9] == 0.0 and values[8] != 0.0)
                or (values[15] == 0.0 and values[14] != 0.0)
            ):
                return stripped
        except ValueError:
            return stripped

    try:
        values = _parse_cal_matches4_values(stripped)
    except ValueError:
        return stripped
    # psf_star_mags with sky_model=0 writes ZMIN instead of A1_min, leaving I_hfs=0.
    if values[9] == 0.0 and values[8] != 0.0:
        values[9] = values[8]
    if values[15] == 0.0 and values[14] != 0.0:
        values[15] = values[14]
    return _format_cal_matches4_line(values)


def fix_vi_hst_ogle_cal_matches4_spacing(directory):
    """
    VI_HST_ogle_man_match4 writes Cal_matches4 with glued f9.4/f8.3 fields and
    Fortran overflow masks (********). Normalize spacing and replace overflow
    tokens so fit_HST_IV_ogle_col can parse the file.
    """
    file_path = Path(directory).resolve() / "07.CALIBRATION" / "VI_HST_ogle_Cal_matches4.dat"
    if not file_path.exists():
        return

    out_lines = [
        _repair_vi_hst_cal_matches4_data_line(line)
        for line in file_path.read_text(encoding="utf-8").splitlines()
    ]
    file_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


_OGLE_AS_PIX = 0.26


def _parse_psf_star_fsky_table(psf_path):
    """Return PSF-star fsky mags keyed by star number."""
    psf_path = Path(psf_path)
    if not psf_path.is_file():
        return {}
    band = "Imag" if "_I." in psf_path.name or "Cal_I" in psf_path.name else "Vmag"
    lg_key = "lg_c2_I" if band == "Imag" else "lg_c2_V"
    stars = {}
    kstar = None
    chi2 = -1.0
    for line in psf_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# BEST-FIT MODEL"):
            kstar = int(line.split(":")[-1])
            chi2 = -1.0
        elif line.startswith("# chi2MIN =") and kstar is not None:
            tail = line.split("=", 1)[1].strip()
            chi2 = -1.0 if "*" in tail else float(tail)
        elif line.startswith("#  A1_min =") and kstar is not None:
            a1 = float(line.split("=", 1)[1])
            mag = -2.5 * np.log10(a1)
            lg = 99.999 if chi2 <= 0 else np.log10(chi2) + 0.25 * mag
            stars.setdefault(kstar, {})[band] = mag
            stars[kstar][lg_key] = lg
    return stars


def _fit_bar_to_ogle_from_vi_hst_in(in_path):
    """Affine bar(MATCHUP) -> OGLE map fit from manual IN anchors."""
    pairs = []
    for line in Path(in_path).read_text(encoding="utf-8").splitlines()[7:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        bar = (float(parts[0]), float(parts[1]))
        ogle = (float(parts[2]), float(parts[3]))
        if bar[0] == 0.0 and bar[1] == 0.0:
            continue
        pairs.append((bar, ogle))
    if len(pairs) < 3:
        raise ValueError(f"need at least 3 manual anchors in {in_path}")
    design = []
    target = []
    for (xb, yb), (xo, yo) in pairs:
        design.append([xb, yb, 1.0, 0.0, 0.0, 0.0])
        design.append([0.0, 0.0, 0.0, xb, yb, 1.0])
        target.extend([xo, yo])
    coeff, _, _, _ = np.linalg.lstsq(np.array(design), np.array(target), rcond=None)

    def bar_to_ogle(xb, yb):
        return (
            coeff[0] * xb + coeff[1] * yb + coeff[2],
            coeff[3] * xb + coeff[4] * yb + coeff[5],
        )

    return bar_to_ogle, pairs


def _parse_all_matches4_row(line: str):
    parts = line.split()
    if len(parts) < 17:
        return None
    try:
        return {
            "cal": int(parts[0]),
            "x": float(parts[1]),
            "y": float(parts[2]),
            "nmat": int(parts[3]),
            "nbmat": int(parts[4]),
            "I_ogle": float(parts[5]),
            "V_ogle": float(parts[6]),
            "I_hst1": float(parts[7]),
            "I_hst": float(parts[8]),
            "I_hstB": float(parts[9]),
            "V_hst1": float(parts[10]),
            "V_hst": float(parts[11]),
            "V_hstB": float(parts[12]),
            "Io_Ih1": float(parts[13]),
            "Io_Ih": float(parts[14]),
            "Vo_Vh1": float(parts[15]),
            "Vo_Vh": float(parts[16]),
        }
    except ValueError:
        return None


def _cal_matches4_header_lines(cal_dir):
    headers = []
    path = cal_dir / "VI_HST_ogle_Cal_matches4.dat"
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped and (stripped[0].isdigit() or (line.startswith(" ") and stripped[0].isdigit())):
                break
            headers.append(line)
    if len(headers) < 6:
        headers = [
            "# I_fsky file =",
            "psf_star_mags_Cal_I.05.final_fits",
            "# V_fsky file =",
            "psf_star_mags_Cal_V.05.final_fits",
            "# match radius in arc sec =   50.00 HWHMa =   0.40 HWHMb =   1.00  A_lg_chi2_I, A_lg_chi2_V =   0.250   0.250",
            "# Cal#   x        y     nmat nbmat I_ogle  V_ogle  I_hst1  I_hst   I_hfs   Io-Ihfs lg_c2Imx I_hstB  V_hst1  V_hst   V_hfs   Vo-Vhfs lg_c2Vmx V_hstB  Io-Ih1  Io-Ih   Vo-Vh1  Vo-Vh",
        ]
    return headers


def _assign_psf_star_index(row, psf_i, psf_v):
    best = None
    best_score = None
    for kstar, data in psf_i.items():
        imag = data.get("Imag")
        if imag is None:
            continue
        vmag = psf_v.get(kstar, {}).get("Vmag", row["V_hst"])
        score = abs(row["I_hst"] - imag) + abs(row["V_hst"] - vmag)
        if best_score is None or score < best_score:
            best_score = score
            best = kstar
    return best or 1


def _build_cal_matches4_values(row, psf_num, psf_i, psf_v):
    i_fsky = psf_i.get(psf_num, {}).get("Imag", row["I_hst"])
    v_fsky = psf_v.get(psf_num, {}).get("Vmag", row["V_hst"])
    if i_fsky == 0.0:
        i_fsky = row["I_hst"]
    if v_fsky == 0.0:
        v_fsky = row["V_hst"]
    lg_c2_i = psf_i.get(psf_num, {}).get("lg_c2_I", 0.0)
    lg_c2_v = psf_v.get(psf_num, {}).get("lg_c2_V", 0.0)
    return [
        psf_num,
        row["x"],
        row["y"],
        row["nmat"],
        row["nbmat"],
        row["I_ogle"],
        row["V_ogle"],
        row["I_hst1"],
        row["I_hst"],
        i_fsky,
        row["I_ogle"] - i_fsky,
        lg_c2_i,
        row["I_hstB"],
        row["V_hst1"],
        row["V_hst"],
        v_fsky,
        row["V_ogle"] - v_fsky,
        lg_c2_v,
        row["V_hstB"],
        row["Io_Ih1"],
        row["Io_Ih"],
        row["Vo_Vh1"],
        row["Vo_Vh"],
    ]


def expand_vi_hst_cal_matches_for_fit(
    directory, *, min_stars=12, match_arcsec=2.0, max_stars=20
):
    """
    VI_HST only tags 11 bar cal stars, but few land on OGLE sources at HWHMa.
    Rebuild Cal_matches4 from all_matches near manual anchors + NOTFAR positions.
    """
    cal_dir = Path(directory).resolve() / "07.CALIBRATION"
    all_path = cal_dir / "VI_HST_ogle_all_matches4.dat"
    out_path = cal_dir / "VI_HST_ogle_Cal_matches4.dat"
    in_path = cal_dir / "IN.VI_HST_ogle_man_match4"
    if not all_path.is_file() or not in_path.is_file():
        return

    bar_to_ogle, anchor_pairs = _fit_bar_to_ogle_from_vi_hst_in(in_path)
    ogle_targets = [ogle for _bar, ogle in anchor_pairs]
    notfar = cal_dir / "NOTFAR_CAL_STARS.XYIVB_targ"
    if notfar.is_file():
        for row in np.atleast_2d(np.loadtxt(notfar, skiprows=2)):
            ogle_targets.append(bar_to_ogle(float(row[0]), float(row[1])))

    match_rad = match_arcsec / _OGLE_AS_PIX
    candidates = []
    for line in all_path.read_text(encoding="utf-8").splitlines():
        row = _parse_all_matches4_row(line)
        if row is None:
            continue
        if row["nmat"] <= 0 or row["nbmat"] <= 0:
            continue
        if row["I_hst"] > 90 or row["V_hst"] > 90 or row["I_ogle"] > 30:
            continue
        nearest = min(
            ((row["x"] - tx) ** 2 + (row["y"] - ty) ** 2) ** 0.5
            for tx, ty in ogle_targets
        )
        if nearest <= match_rad:
            candidates.append((nearest, row))

    if len(candidates) < min_stars:
        extra = []
        for line in all_path.read_text(encoding="utf-8").splitlines():
            row = _parse_all_matches4_row(line)
            if row is None:
                continue
            if row["nmat"] <= 0 or row["nbmat"] <= 0:
                continue
            if row["I_hst"] > 90 or row["V_hst"] > 90:
                continue
            nearest = min(
                ((row["x"] - tx) ** 2 + (row["y"] - ty) ** 2) ** 0.5
                for tx, ty in ogle_targets
            )
            extra.append((nearest, row))
        extra.sort(key=lambda item: item[0])
        seen = {(round(r["x"], 3), round(r["y"], 3)) for _, r in candidates}
        for dist, row in extra:
            key = (round(row["x"], 3), round(row["y"], 3))
            if key in seen:
                continue
            candidates.append((dist, row))
            seen.add(key)
            if len(candidates) >= min_stars:
                break

    candidates.sort(key=lambda item: item[0])
    candidates = candidates[:max_stars]
    if not candidates:
        return

    psf_i = _parse_psf_star_fsky_table(cal_dir / "psf_star_mags_Cal_I.05.final_fits")
    psf_v = _parse_psf_star_fsky_table(cal_dir / "psf_star_mags_Cal_V.05.final_fits")
    used_psf = set()
    data_lines = []
    for _dist, row in candidates:
        psf_num = _assign_psf_star_index(row, psf_i, psf_v)
        while psf_num in used_psf and psf_num < 99:
            psf_num += 1
        used_psf.add(psf_num)
        data_lines.append(_format_cal_matches4_line(_build_cal_matches4_values(row, psf_num, psf_i, psf_v)))

    headers = _cal_matches4_header_lines(cal_dir)
    out_path.write_text("\n".join(headers + data_lines) + "\n", encoding="utf-8")
    print(f"Expanded {out_path.name} to {len(data_lines)} calibration stars")


parser = argparse.ArgumentParser()

def get_fortran_dir():
    """
    Return the path to the compiled Fortran executables bundled with MORIA.
    """
    return resources.files("moria") / "fortran_compile"

def copy_files(source, destination, extensions=[".fits"]):
    """
    Copy files from one folder to another with the same extension.
    """
    if not os.path.exists(source):
        raise FileNotFoundError()

    for f in os.listdir(source):
        path = os.path.join(source, f)
        if os.path.isfile(path) and (any(f.endswith(e) for e in extensions)):
            shutil.copy2(source / f, os.path.join(destination, f))

def copy_entire_files(source, destination, filename):
    """
    Copy entire files from one folder to another. 
    """
    if not os.path.exists(source):
        raise FileNotFoundError()

    path = os.path.join(source, filename)
    if os.path.isfile(path):
        shutil.copy2(source / filename, os.path.join(destination, filename))


def xympquvwrd_to_physical_xym(xym_dir):
    """
    Read each ``*xympquvwrd`` in ``xym_dir`` and write a sibling ``*.xym``.

    hst1pass ``OUT=xympquvwrd`` rows are: ``x_raw y_raw m p q u v w r d ...``.
    DS9 **Physical** circles in ``REG=uv`` use ``u,v`` (cols 6–7).  The output
    ``.xym`` repeats the same numeric layout but replaces the first two columns
    with ``u,v`` so list **x,y** match DS9 Physical while ``m p q u v w r d`` stay
    aligned with the source (``u,v`` appear in both positions 1–2 and 6–7).

    Example: ``visit_flc.F814W_xympquvwrd`` → ``visit_flc.F814W.xym``.

    Parameters
    ----------
    xym_dir : str or pathlib.Path
        Directory containing the hst1pass lists (often ``01.XYM``).

    Returns
    -------
    list[pathlib.Path]
        Paths of written ``.xym`` files.
    """
    d = Path(xym_dir).resolve()
    if not d.is_dir():
        raise FileNotFoundError(str(d))

    annot = (
        "# MORIA: cols 01-02 = u,v (DS9 Physical, hst1pass REG=uv); "
        "raw xy were former 01-02 of this xympquvwrd.\n"
    )
    written: list[Path] = []

    for src in sorted(d.glob("*xympquvwrd")):
        if not src.is_file():
            continue
        name = src.name
        if not name.endswith("xympquvwrd"):
            continue
        stem = name[: -len("_xympquvwrd")]
        out_path = src.with_name(stem + ".xym")

        out_lines: list[str] = []
        inserted_annot = False
        with src.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if (not inserted_annot) and "#  FILEOUT:" in line:
                    out_lines.append(line)
                    out_lines.append(annot)
                    inserted_annot = True
                    continue
                s = line.strip()
                if not s or s.startswith("#"):
                    out_lines.append(line)
                    continue
                parts = s.split()
                if len(parts) < 10:
                    out_lines.append(line)
                    continue
                _xr, _yr = parts[0], parts[1]
                m, p, q = parts[2], parts[3], parts[4]
                u, v, w, r, d = parts[5], parts[6], parts[7], parts[8], parts[9]
                tail = parts[10:]
                try:
                    uf, vf = float(u), float(v)
                    mf = float(m)
                    pf = float(p)
                    qf = float(q)
                    wf = float(w)
                    rf = float(r)
                    ddf = float(d)
                except ValueError:
                    out_lines.append(line)
                    continue
                # Match ~hst1pass column spacing (see xympquvwrd examples).
                body = (
                    f"{uf:11.4f}{vf:11.4f}"
                    f"{mf:8.3f}{pf:9.2f}{qf:10.3f}"
                    f"{uf:11.4f}{vf:11.4f}"
                    f"{wf:8.3f}{rf:14.8f}{ddf:14.8f}"
                )
                if tail:
                    body += " " + " ".join(tail)
                out_lines.append(body + "\n")

        if not inserted_annot:
            out_lines.insert(0, annot)

        out_path.write_text("".join(out_lines), encoding="utf-8")
        written.append(out_path)

    return written


_WFC3UV_PSF_STDPSF_BASE = (
    "https://www.stsci.edu/~jayander/WFC3/WFC3UV_PSFs/STDPSF/"
)
_WFC3UV_PSF_BASE = "https://www.stsci.edu/~jayander/WFC3/WFC3UV_PSFs/"
_WFC3UV_STDGDC_BASE = (
    "https://www.stsci.edu/~jayander/HST1PASS/LIB/GDCs/STDGDCs/WFC3UV/"
)

# (base URL, filename, per-file timeout seconds) — Jay Anderson WFC3/UV reference FITS
_WFC3UV_REFERENCE_DOWNLOADS = (
    (_WFC3UV_PSF_STDPSF_BASE, "PSFSTD_WFC3UV_F814W.fits", 120.0),
    (_WFC3UV_PSF_STDPSF_BASE, "PSFSTD_WFC3UV_F606W.fits", 120.0),
    (_WFC3UV_PSF_BASE, "PSFEFF_WFC3UV_F814W_C0.fits", 120.0),
    (_WFC3UV_PSF_BASE, "PSFEFF_WFC3UV_F606W_C0.fits", 120.0),
    (_WFC3UV_STDGDC_BASE, "STDGDC_WFC3UV_F814W.fits", 7200.0),
    (_WFC3UV_STDGDC_BASE, "STDGDC_WFC3UV_F606W.fits", 7200.0),
)

_ACS_PSF_STDPSF_BASE = ("https://www.stsci.edu/~jayander/HST1PASS/LIB/PSFs/STDPSFs/ACSWFC")
_ACS_STDGDC_BASE = ("https://www.stsci.edu/~jayander/HST1PASS/LIB/GDCs/STDGDCs/ACSWFC")
_ACS_REFERENCE_DOWNLOADS = (
    (_ACS_PSF_STDPSF_BASE, "STDPSF_ACSWFC_F814W_SM4.fits", 120.0),
    (_ACS_PSF_STDPSF_BASE, "STDPSF_ACSWFC_F606W_SM4.fits", 120.0),
    (_ACS_STDGDC_BASE, "STDGDC_OFFICIAL_JFRAME_ACSWFC_F606W.fits", 7200.0),
    (_ACS_STDGDC_BASE, "STDGDC_OFFICIAL_JFRAME_ACSWFC_F814W.fits", 7200.0),
)



MATCHUP_JAY_HEADER = (
    "#.......... ........... ........... ........... ........... "
    "........... ........... ................ ....... ..... ..... ............\n"
    "#    xbar        ybar        mbar        xsig        ysig        msig       "
    " qbar    Nf  Ng  Nm Nmin Nstar     pki   pkj  pkp pkn pku\n"
    "#    (01)        (02)        (03)        (04)        (05)        (06)      "
    "  (07)   (08)(09)(10)(11)   (12)   (13)  (14) (15)(16)(17)\n"
    "#.......... ........... ........... ........... ........... "
    "........... ........... ................ ....... ..... ..... ............\n"
)

MATCHUP_CAL_ONLY_HEADER = (
    "#.......... ........... .........\n"
    "#    xbar        ybar        mbar\n"
    "#    (01)        (02)        (03)\n"
    "#.......... ........... .........\n"
)

UVP2TRI_FSKY_OUTPUTS_F814W = (
    "uvp2tri_scon_fsky_I_KeckNOcon.01.pix_all",
    "uvp2tri_scon_fsky_I_KeckNOcon.03.pix_use",
    "uvp2tri_scon_fsky_I_KeckNOcon.04.probe_fit",
    "uvp2tri_scon_fsky_I_KeckNOcon.05.final_fit",
    "uvp2tri_scon_fsky_I_KeckNOcon.06.pix_show.fits",
    "uvp2tri_scon_fsky_I_KeckNOcon.07.mcmc",
    "uvp2tri_scon_fsky_I_KeckNOcon.08.rm_pix",
)
UVP2TRI_FSKY_OUTPUTS_F606W = (
    "uvp2tri_scon_fsky_V_KeckNOcon.01.pix_all",
    "uvp2tri_scon_fsky_V_KeckNOcon.03.pix_use",
    "uvp2tri_scon_fsky_V_KeckNOcon.04.probe_fit",
    "uvp2tri_scon_fsky_V_KeckNOcon.05.final_fit",
    "uvp2tri_scon_fsky_V_KeckNOcon.06.pix_show.fits",
    "uvp2tri_scon_fsky_V_KeckNOcon.07.mcmc",
    "uvp2tri_scon_fsky_V_KeckNOcon.08.rm_pix",
)
_UVP2TRI_MCMC_SCRIPT = "run_uvp2tri_NOscon_fs_asym_mcmc.src"
_UVP2TRI_MCMC_ALT_SCRIPT = "run_uvp2tri_NOscon_fs_asym_mcmc_alt.src"


def uvp2tri_fsky_outputs_missing(fit_dir, output_names):
    return [name for name in output_names if not (fit_dir / name).is_file()]


def remove_uvp2tri_fsky_outputs(fit_dir, output_names):
    for name in output_names:
        (fit_dir / name).unlink(missing_ok=True)



def run_uvp2tri_mcmc_csh(directory, filter_name, fit_folder, script, log_basename):
    base_dir = Path(directory).resolve() / "06.FIT" / filter_name
    subdir = base_dir / fit_folder
    log_file = subdir / "log_files" / f"{log_basename}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    script_path = subdir / script
    with open(log_file, "w") as logf:
        subprocess.run(
            ["csh", str(script_path)],
            cwd=subdir,
            stdout=logf,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )


def run_uvp2tri_mcmc(directory, filter_name, fit_folders, output_names):
    """Run uvp2tri MCMC; retry with alt script for fit folders missing expected outputs."""
    base_dir = Path(directory).resolve() / "06.FIT" / filter_name
    for fit_folder in fit_folders:
        run_uvp2tri_mcmc_csh(
            directory, filter_name, fit_folder, _UVP2TRI_MCMC_SCRIPT, "uvp2tri_scon_fs_asym_mcmc"
        )
    for fit_folder in fit_folders:
        fit_dir = base_dir / fit_folder
        if uvp2tri_fsky_outputs_missing(fit_dir, output_names):
            remove_uvp2tri_fsky_outputs(fit_dir, output_names)
            run_uvp2tri_mcmc_csh(
                directory,
                filter_name,
                fit_folder,
                _UVP2TRI_MCMC_ALT_SCRIPT,
                "uvp2tri_scon_fs_asym_mcmc_alt",
            )


_MCMC_ACCEPTANCE_RE = re.compile(
    r"accepted,\s*rejected\s+MCMC\s+steps:\s*(\d+)\s+(\d+)",
    re.IGNORECASE,
)



def print_uvp2tri_mcmc_acceptance_rate(directory, filter_name, fit_folders):
    """Read uvp2tri_scon_fs_asym_mcmc.log and print the MCMC acceptance rate."""
    base_dir = Path(directory).resolve() / "06.FIT" / filter_name
    fit_folders = ['1star-fit', '2star-fit']
    for fit_folder in fit_folders:
        log_file = base_dir / fit_folder / "log_files" / "uvp2tri_scon_fs_asym_mcmc.log"
        if not log_file.is_file():
            print(f"{filter_name}/{fit_folder}: log file not found: {log_file}")
            continue
        matches = _MCMC_ACCEPTANCE_RE.findall(log_file.read_text())
        naccept, nreject = map(int, matches[-1])
        total = naccept + nreject
        rate = naccept / total if total else 0.0
        if fit_folder == '1star-fit':
            string_1star = str(f"{filter_name}/{fit_folder}: MCMC acceptance rate = {rate:.4f} "f"({naccept} accepted, {nreject} rejected)")
        elif fit_folder == '2star-fit':
            string_2star = str(f"{filter_name}/{fit_folder}: MCMC acceptance rate = {rate:.4f} "f"({naccept} accepted, {nreject} rejected)")
    return string_1star, string_2star

def print_uvp2tri_mcmc_acceptance_rate_3star(directory, filter_name, fit_folders):
    """Read uvp2tri_scon_fs_asym_mcmc.log and print the MCMC acceptance rate."""
    base_dir = Path(directory).resolve() / "06.FIT" / filter_name
    fit_folders = ['3star-fit']
    for fit_folder in fit_folders:
        log_file = base_dir / fit_folder / "log_files" / "uvp2tri_scon_fs_asym_mcmc.log"
        if not log_file.is_file():
            print(f"{filter_name}/{fit_folder}: log file not found: {log_file}")
            continue
        matches = _MCMC_ACCEPTANCE_RE.findall(log_file.read_text())
        naccept, nreject = map(int, matches[-1])
        total = naccept + nreject
        rate = naccept / total if total else 0.0
        string_3star = str(f"{filter_name}/{fit_folder}: MCMC acceptance rate = {rate:.4f} "f"({naccept} accepted, {nreject} rejected)")
    return string_3star

def _download_url_to_path(url, out_path, timeout_s, *, verbose):
    """Download ``url`` to ``out_path`` (write via a ``.part`` temp file)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    if verbose:
        print(f"Downloading {out_path.name} …", flush=True)
    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as resp:
            with open(tmp_path, "wb") as f:
                shutil.copyfileobj(resp, f)
        tmp_path.replace(out_path)
        if verbose:
            mb = out_path.stat().st_size / (1024 * 1024)
            print(f"  → {out_path} ({mb:.1f} MiB)", flush=True)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise


def download_wfc3uv_psf_libraries(
    camera: str,
    directory: str | Path,
    dest_subdirs: tuple[str, ...] = ("01.XYM", "03.LOC_TRANS"),
    overwrite: bool = False,
    timeout_s: float | None = None,
    verbose: bool = True,
) -> dict[str, list[str]]:
    """
    Fetch WFC3/UV PSF and GDC reference FITS from Jay Anderson's STScI pages.

    Writes each file into every ``dest_subdirs`` folder under ``directory``
    (default ``01.XYM`` and ``03.LOC_TRANS``).
    """
    root = Path(directory).resolve()
    results: dict[str, list[str]] = {"downloaded": [], "skipped": []}

    if verbose:
        print(f"MORIA reference FITS → {root}", flush=True)

    if camera == "WFC3UV":
        reference_downloads = _WFC3UV_REFERENCE_DOWNLOADS
    elif camera == "ACS":
        reference_downloads = _ACS_REFERENCE_DOWNLOADS
    else:
        raise ValueError(f"Invalid camera: {camera}")

    for sub in dest_subdirs:
        dest_dir = root / sub
        dest_dir.mkdir(parents=True, exist_ok=True)
        if verbose:
            print(f"  [{sub}]", flush=True)
        for base_url, name, file_timeout in reference_downloads:
            url = base_url.rstrip("/") + "/" + name
            out_path = dest_dir / name
            if out_path.exists() and not overwrite:
                results["skipped"].append(str(out_path))
                if verbose:
                    print(f"  skip (exists): {out_path.name}", flush=True)
                continue
            tmo = file_timeout if timeout_s is None else timeout_s
            _download_url_to_path(url, out_path, timeout_s=tmo, verbose=verbose)
            results["downloaded"].append(str(out_path))

    if verbose:
        print(
            f"Done: {len(results['downloaded'])} downloaded, "
            f"{len(results['skipped'])} skipped.",
            flush=True,
        )

    return results


def data_prep_early(camera, destination):
    download_wfc3uv_psf_libraries(camera, destination)

    fortran_src = get_fortran_dir()
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "00.DATA" / "F814W", extensions=[".xOg"])
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "00.DATA" / "F606W", extensions=[".xOg"])
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM", filename="dex_no_gaia.e")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM", filename="1exp_no_gaia.e")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM", filename="hst1pass.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM", filename="hst1pass.F")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS", filename="img2extract_wfc3uv_psflist.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS", filename="xym2mat_new.e")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS", filename="1exp_xym2mat_new.e")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS", filename="img2extract_wfc3uv_psflist.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "04.EXTRACT_PSF" / "F814W", filename="uvp2psf_simst.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "04.EXTRACT_PSF" / "F606W", filename="uvp2psf_simstV.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F814W" / "1star-fit", filename="mcmc_expand_average.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F814W" / "2star-fit", filename="mcmc_expand_average.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F606W" / "1star-fit", filename="mcmc_expand_average.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F606W" / "2star-fit", filename="mcmc_expand_average.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F814W" / "1star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F814W" / "2star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F606W" / "1star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "06.FIT" / "F606W" / "2star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "07.CALIBRATION", filename="VI_HST_ogle_man_match4.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "07.CALIBRATION", filename="fit_HST_IV_ogle_col.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "07.CALIBRATION", filename="psf_star_mags_mcmc.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "07.CALIBRATION", filename="cal_star_num_2_MATCHUP.xOg")

    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "05.KECK_TRANS", filename="matched_HST_Keck_stars.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "05.KECK_TRANS", filename="HST_Keck_coord_trans.xOg")


def data_prep_module_five(destination):
    fortran_src = get_fortran_dir()
    copy_entire_files(source=Path(destination).resolve() / "01.XYM" / "F814W", destination=Path(destination).resolve() / "05.KECK_TRANS", filename="outputq_F814W.fits")
    copy_entire_files(source=Path(destination).resolve() / "01.XYM" / "F606W", destination=Path(destination).resolve() / "05.KECK_TRANS", filename="outputq_F606W.fits")

    copy_entire_files(source=Path(destination).resolve() / "02.CMD", destination=Path(destination).resolve() / "05.KECK_TRANS", filename="MATCHUP.F814W.XYM.02")
    copy_entire_files(source=Path(destination).resolve() / "02.CMD", destination=Path(destination).resolve() / "05.KECK_TRANS", filename="MATCHUP.F606W.XYM")



def run_xgf_conversion(directory, script='run_convert_C1K1C.src'):
    """
    Run the conversion script on _flc files in F814W and F606W subdirectories.
    Parameters
    ----------
    directory - The root directory to operate on. There is a specific directory
    structure that is expected within <dir>:
        00.DATA/
        01.XYM/

    Returns
    -------
    Convert from _flc files to _WJ2 files
    """
    base_dir = Path(directory).resolve() / '00.DATA'
    filters = ['F814W', 'F606W']
    
    for f in filters:
        subdir = base_dir / f
        script_path = subdir / script
        
        if not subdir.exists():
            print(f'Missing subdirectory.')
            continue

        flc_files = list(subdir.glob('*_flc.fits'))

        if not flc_files:
            print(f'Missing flc files')
            continue

        subprocess.run(
            ['csh', str(script_path)],  
            cwd=subdir,
            check=False
        )

    return

def data_prep(directory):
    """
    Preparae IN.* files in respective files using _WJ2.fits files in 00.DATA
    """
    
    def _estimate_nstars_from_xym(path: Path, max_lines: int = 200000) -> int:
        """
        Rough estimate of number of sources in an .xym file.
        Ignores blank lines and comment lines (starting with '#').
        """
        try:
            n = 0
            with path.open("r", encoding="utf-8", errors="ignore") as fh:
                for i, line in enumerate(fh):
                    if i >= max_lines:
                        break
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue
                    n += 1
            return n
        except OSError:
            return 0

    def _choose_mag_clip_for_xym(xym_path: Path, default_clip: str) -> str:
        """
        Pick a more inclusive magnitude clip for sparse lists.
        Instrumental magnitudes are negative; "mA,B" keeps A < m < B.
        """
        nstars = _estimate_nstars_from_xym(xym_path)
        if nstars and nstars < 200:
            if default_clip == "m-13.75,-8.5":
                return "m-13.75,-6.0"
            if default_clip == "m-14.75,-5.5":
                return "m-14.75,-3.5"
        return default_clip
        
        
    def data_prep_F814W(directory, f= 'F814W'):
        base_dir = Path(directory).resolve()
        subdir = base_dir / f
        in_img2sam_wfc3uv = 'IN.img2sam_wfc3uv'
        in_xym2bar_1 = 'IN.xym2bar.1'
        in_xym2bar_2 = 'IN.xym2bar.2'
        in_xym2mat = 'IN.xym2mat'
        in_xym2bar = 'IN.xym2bar'
        in_xym2mat_1 = 'IN.xym2mat.1'
        in_xym2mat_2 = 'IN.xym2mat.2'

        base_dir_one = base_dir / '01.XYM'/ f
        files = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.xym')])
        files_two = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.fits')])

        default_clip = "m-13.75,-8.5"
        clip = default_clip
        if files:
            clip = _choose_mag_clip_for_xym(base_dir_one / files[0], default_clip)


        
        output_file_dir = base_dir / '01.XYM' / f
        output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
        output_file_xym2mat = os.path.join(output_file_dir, in_xym2mat)
        output_file_xym2bar = os.path.join(output_file_dir, in_xym2bar)
        output_file_xym2mat1 = os.path.join(output_file_dir, in_xym2mat_1)
        output_file_xym2mat2 = os.path.join(output_file_dir, in_xym2mat_2)
        output_file_xym2bar1 = os.path.join(output_file_dir, in_xym2bar_1)
        output_file_xym2bar2 = os.path.join(output_file_dir, in_xym2bar_2)

        

        with open(output_file_xym2mat1, "w") as f:
            f.write("#00 MATCHUP.XYM.01 c0\n")
            for i, filename in enumerate(files, start=1):
                if i == 1:
                    f.write(f"{0:02d} {filename} c8 f8 \"{clip}\" \n")
                f.write(f"{i:02d} {filename} c8 f8 \"{clip}\" \n")


        with open(output_file_xym2bar1, "w") as f:
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 z0\n")
    
        with open(output_file_img2sam, "w") as t:
            for i, filename in enumerate(files_two, start=1):
                t.write(f"{i:02d} \"{filename}\" 8 0\n")

        with open(output_file_xym2bar2, "w") as f:
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 z0\n")

        with open(output_file_xym2mat2, "w") as f:
            f.write("00 MATCHUP.XYM.01 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 \"{clip}\" \n")
        return
    
    
    def data_prep_F606W(directory, f = "F606W"):
        base_dir = Path(directory).resolve()
        subdir = base_dir / f
        in_img2sam_wfc3uv = 'IN.img2sam_wfc3uv'
        in_xym2mat = 'IN.xym2mat'
        in_xym2bar = 'IN.xym2bar'

        base_dir_one = base_dir / '01.XYM'/ f
        files = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.xym')])
        files_two = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.fits')])

        default_clip = "m-14.75,-5.5"
        clip = default_clip
        if files:
            clip = _choose_mag_clip_for_xym(base_dir_one / files[0], default_clip)

        
        output_file_dir = base_dir / '01.XYM' / f
        output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
        output_file_xym2mat = os.path.join(output_file_dir, in_xym2mat)
        output_file_xym2bar = os.path.join(output_file_dir, in_xym2bar)
        
        with open(output_file_xym2mat, "w") as f:
            f.write("00 MATCHUP.F814W.XYM.02 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f6 \"{clip}\" \n")


        with open(output_file_xym2bar, "w") as f:
            f.write("00 MATCHUP.F814W.XYM.02 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f6\n")


        with open(output_file_img2sam, "w") as t:
            for i, filename in enumerate(files_two, start=1):
                t.write(f"{i:02d} \"{filename}\" 6 0\n")
        
                
        return
    
    #data_prep_F814W(directory)
    #data_prep_F606W(directory)
    
def _count_wj2_xym(field_dir: Path) -> int:
    """Number of WFC3 *WJ2.xym lists in 01.XYM/Filter (one per exposure/stack)."""
    if not field_dir.is_dir():
        return 0
    return len(sorted(field_dir.glob("*WJ2.xym")))


def _write_xym2bar_run_scripts(field_root: Path) -> None:
    """
    xym2bar's first argument is NIMMIN
    """
    # F814W scripts (produce MATCHUP.XYM.01/02 and MATCHUP.F814W.XYM.02)
    fdir_I = (field_root / "01.XYM" / "F814W").resolve()
    nexp_I = _count_wj2_xym(fdir_I)
    if nexp_I >= 1:
        nim_I = str(nexp_I)
        (fdir_I / "run_xym2bar_1.src").write_text(f"cp -p IN.xym2bar.1 IN.xym2bar\n./xym2bar.xOg {nim_I}\ncp -p MATCHUP.XYMEEE MATCHUP.XYM.01\n", 
        encoding="utf-8")
        (fdir_I / "run_xym2bar_2.src").write_text(
            f"cp -p IN.xym2bar.2 IN.xym2bar\n./xym2bar.xOg {nim_I}\n"
            f"cp -p MATCHUP.XYMEEE MATCHUP.XYM.02\n"
            f"cp -p MATCHUP.XYMEEE MATCHUP.F814W.XYM.02\n",
            encoding="utf-8")

    # F606W script (produces MATCHUP.F606W.XYM using I-band MATCHUP as input list)
    fdir_V = (field_root / "01.XYM" / "F606W").resolve()
    nexp_V = _count_wj2_xym(fdir_V)
    if nexp_V >= 1:
        nim_V = str(nexp_V)
        # Keep "I" so xym2bar reads the NIM=0 input list from IN.xym2bar.
        (fdir_V / "run_xym2bar.src").write_text(f"./xym2bar.xOg {nim_V} I DMATCH=8\ncp -p MATCHUP.XYMEEE MATCHUP.F606W.XYM\n", encoding="utf-8")


def _f606_flc_count(directory) -> int:
    """Number of F606W *_flc.fits exposures in 00.DATA."""
    data_dir = Path(directory).resolve() / "00.DATA" / "F606W"
    return len(list(data_dir.glob("*_flc.fits")))


def _hstpass_catalog_counts(directory):
    """Return sorted F606W and F814W hst1pass list paths under 01.XYM."""
    base_dir = Path(directory).resolve() / "01.XYM"
    f606 = sorted(base_dir.glob("*.F606W_xympquvwrd"))
    f814 = sorted(base_dir.glob("*.F814W_xympquvwrd"))
    return f606, f814


def _use_single_f606w_mode(directory) -> bool:
    """True when the field has one F606W exposure (1 F606W + 2 F814W)."""
    return _f606_flc_count(directory) == 1

    
def matchup_files(camera,directory):
    """
    Run the scripts to create MATCHUP Files on _WJ2 files in F814W and F606W subdirectories.
    """

    def reduce_wfc3(directory, script='reduce_wfc3.src'):
        """
        Run hst1pass on exposures in the F814W and F606W filters.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / "reduce_wfc3.log"
        base_dir = Path(directory).resolve() / "01.XYM"
        script_path = base_dir / script
        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )

    def reduce_acs(directory, script='reduce_acs.src'):
        """
        Run hst1pass on exposures in the F814W and F606W filters.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / "reduce_acs.log"
        base_dir = Path(directory).resolve() / "01.XYM"
        script_path = base_dir / script
        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )
    def reduce_wfc3_1exp(directory, script='reduce_wfc3_1exp.src'):
        """
        Run hst1pass on exposures in the F814W and F606W filters.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / "reduce_wfc3_1exp.log"
        base_dir = Path(directory).resolve() / "01.XYM"
        script_path = base_dir / script
        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )

    def reduce_acs_1exp(directory, script='reduce_acs_1exp.src'):
        """
        Run hst1pass on exposures in the F814W and F606W filters.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / "reduce_acs_1exp.log"
        base_dir = Path(directory).resolve() / "01.XYM"
        script_path = base_dir / script
        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )

    def no_gaia_matchup(directory, script='no_gaia_match.src'):
        """
        Write ``no_gaia_match.src`` and run ``dex_no_gaia.e`` on the four hst1pass lists
        produced by ``reduce_wfc3`` (two ``*.F606W_xympquvwrd``, two ``*.F814W_xympquvwrd``
        in ``01.XYM``).
        """
        base_dir = Path(directory).resolve() / "01.XYM"
        log_dir = base_dir / "log_files"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "matchup.log"

        f606 = sorted(base_dir.glob("*.F606W_xympquvwrd"))
        f814 = sorted(base_dir.glob("*.F814W_xympquvwrd"))
        if len(f606) != 2 or len(f814) != 2:
            raise ValueError(
                f"Need exactly 2 '*.F606W_xympquvwrd' and 2 '*.F814W_xympquvwrd'"
                            )
        names606 = [p.name for p in f606]
        names814 = [p.name for p in f814]
        field_root = Path(directory).resolve()
        ref_candidates = [
            (base_dir / "NEARBY_REF_STARS.XYIVB_targ", "NEARBY_REF_STARS.XYIVB_targ"),
            (base_dir / "NEARBY_SIM_STARS.XYIVB_targ", "NEARBY_SIM_STARS.XYIVB_targ"),
            (field_root / "02.CMD" / "NEARBY_REF_STARS.XYIVB_targ", "../02.CMD/NEARBY_REF_STARS.XYIVB_targ"),
            (field_root / "02.CMD" / "NEARBY_SIM_STARS.XYIVB_targ", "../02.CMD/NEARBY_SIM_STARS.XYIVB_targ"),
        ]
        ref_prefix = ""
        for path, argv_token in ref_candidates:
            if path.is_file():
                ref_prefix = argv_token + " "
                break
        line = "./dex_no_gaia.e " + ref_prefix + " ".join(names606 + names814) + "\n"
        script_path = base_dir / script
        script_path.write_text(line, encoding="utf-8")

        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )

    def no_gaia_matchup_1exp(directory, script='1exp_no_gaia_match.src'):
        """
        Write ``1exp_no_gaia_match.src`` and run ``1exp_no_gaia.e`` on one
        ``*.F606W_xympquvwrd`` and two ``*.F814W_xympquvwrd`` files in ``01.XYM``.
        """
        base_dir = Path(directory).resolve() / "01.XYM"
        log_dir = base_dir / "log_files"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "matchup.log"

        f606, f814 = _hstpass_catalog_counts(directory)
        if len(f606) != 1 or len(f814) != 2:
            raise ValueError(
                "Need exactly 1 '*.F606W_xympquvwrd' and 2 '*.F814W_xympquvwrd' "
                f"in {base_dir}, found F606W={len(f606)} F814W={len(f814)}"
            )
        names606 = [p.name for p in f606]
        names814 = [p.name for p in f814]
        field_root = Path(directory).resolve()
        ref_candidates = [
            (base_dir / "NEARBY_REF_STARS.XYIVB_targ", "NEARBY_REF_STARS.XYIVB_targ"),
            (base_dir / "NEARBY_SIM_STARS.XYIVB_targ", "NEARBY_SIM_STARS.XYIVB_targ"),
            (field_root / "02.CMD" / "NEARBY_REF_STARS.XYIVB_targ", "../02.CMD/NEARBY_REF_STARS.XYIVB_targ"),
            (field_root / "02.CMD" / "NEARBY_SIM_STARS.XYIVB_targ", "../02.CMD/NEARBY_SIM_STARS.XYIVB_targ"),
        ]
        ref_prefix = ""
        for path, argv_token in ref_candidates:
            if path.is_file():
                ref_prefix = argv_token + " "
                break
        line = "./1exp_no_gaia.e " + ref_prefix + " ".join(names606 + names814) + "\n"
        script_path = base_dir / script
        script_path.write_text(line, encoding="utf-8")

        with open(log_file, "w") as logf:
            subprocess.run(
                ["bash", str(script_path)],
                cwd=base_dir,
                stdout=logf,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=os.environ,
            )

    use_1exp = _use_single_f606w_mode(directory)

    if camera == "WFC3UV":
        reduce_fn = reduce_wfc3_1exp if use_1exp else reduce_wfc3
    elif camera == "ACS":
        reduce_fn = reduce_acs_1exp if use_1exp else reduce_acs
    else:
        raise ValueError(f"Invalid camera: {camera}")

    reduce_fn(directory)
    if use_1exp:
        no_gaia_matchup_1exp(directory)
    else:
        no_gaia_matchup(directory)

def run_output_stack(directory, script='run_img2sam_wfc3uv_379.src'):
    """
    Create a stack of the scene in the reference frame.
    """
    log_file = Path(directory).resolve() /  "01.XYM" / "log_files" / f"run_img2sam_wfc3uv_379.log"
    with open(log_file, "w") as logf:
        sys.stdout = sys.stderr = logf
        try:
            base_dir = Path(directory).resolve() / "01.XYM"
            filters = ['F814W', 'F606W']
            for f in filters:
                subdir = base_dir / f
                script_path = subdir / script if f == "F814W" else subdir / script
                subprocess.run(
                    ["csh", str(script_path)],  # assumes csh script
                    cwd=subdir,                        
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
        finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__


def write_run_xym2mat_src_loc_trans(
    output_path: Path,
    xym_dir: Path,
    *,
    use_1exp: bool | None = None,
) -> None:
    """
    Write ``run_xym2mat.src`` for LOC_TRANS.

    Invokes ``xym2mat_new.e`` (2 F606W + 2 F814W) or ``1exp_xym2mat_new.e``
    (1 F606W + 2 F814W) with catalog filenames sorted within each filter.
    """
    f606 = sorted(
        p.name
        for p in xym_dir.iterdir()
        if p.is_file() and p.name.endswith(".F606W_xympquvwrd")
    )
    f814 = sorted(
        p.name
        for p in xym_dir.iterdir()
        if p.is_file() and p.name.endswith(".F814W_xympquvwrd")
    )
    if use_1exp is None:
        use_1exp = len(f606) == 1
    if use_1exp:
        if len(f606) != 1 or len(f814) != 2:
            raise FileNotFoundError(
                "write_run_xym2mat_src_loc_trans: need 1 *.F606W_xympquvwrd and "
                f"2 *.F814W_xympquvwrd under {xym_dir}, found "
                f"F606W={len(f606)} F814W={len(f814)}"
            )
        executable = "./1exp_xym2mat_new.e"
    else:
        if len(f606) != 2 or len(f814) != 2:
            raise FileNotFoundError(
                "write_run_xym2mat_src_loc_trans: need 2 *.F606W_xympquvwrd and "
                f"2 *.F814W_xympquvwrd under {xym_dir}, found "
                f"F606W={len(f606)} F814W={len(f814)}"
            )
        executable = "./xym2mat_new.e"
    files = f606 + f814
    output_path.parent.mkdir(parents=True, exist_ok=True)
    line = executable + " " + " ".join(files) + "\n"
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(line)


def write_in_img2sam_wfc3uv_loc_trans(
    camera: str,
    output_path: Path,
    data_dir: Path,
) -> None:
    """
    Write ``IN.img2sam_wfc3uv`` for LOC_TRANS.

    Lists all F814W WJC/WJ2 files first, then all F606W files from ``00.DATA``.
    Supports one or two exposures per filter.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        idx = 1
        for filter_name, chip in (("F814W", 8), ("F606W", 6)):
            fits_dir = data_dir / filter_name
            if camera == "WFC3UV":
                names = sorted(
                    p.name
                    for p in fits_dir.iterdir()
                    if p.is_file() and p.name.endswith("_WJC.fits")
                )
            else:
                names = sorted(
                    p.name
                    for p in fits_dir.iterdir()
                    if p.is_file() and p.name.endswith("_WJ2.fits")
                )
            if not names:
                raise FileNotFoundError(
                    f"write_in_img2sam_wfc3uv_loc_trans: no WJC/WJ2 files under {fits_dir}"
                )
            for filename in names:
                fh.write(f'{idx:02d} "{filename}" {chip} 0\n')
                idx += 1


def data_prep_loc_trans(camera,directory, filters = 'F814W'):
    """
    Prepare LOC_TRANS inputs: copy star lists and FITS into ``03.LOC_TRANS``,
    write ``IN.img2sam_wfc3uv``, and write ``run_xym2mat.src``.

    Uses ``1exp_xym2mat_new.e`` when only one F606W catalog is present
    (1 F606W + 2 F814W); otherwise ``xym2mat_new.e`` (2+2).

    Parameters
    ----------
    directory - The root directory to operate on. There is a specific directory
    structure that is expected within <dir>:
        00.DATA/
        01.XYM/

    Returns
    -------
    several IN.* files and ``run_xym2mat.src``.
    """

    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "03.LOC_TRANS")
    copy_files(source=Path(directory).resolve() / "01.XYM", extensions=[".F606W_xympquvwrd"], destination=Path(directory).resolve() / "03.LOC_TRANS")
    copy_files(source=Path(directory).resolve() / "01.XYM", extensions=[".F814W_xympquvwrd"], destination=Path(directory).resolve() / "03.LOC_TRANS")
    copy_files(source=Path(directory).resolve() / "00.DATA" / "F606W", extensions=[".fits"], destination=Path(directory).resolve() / "03.LOC_TRANS")
    copy_files(source=Path(directory).resolve() / "00.DATA" / "F814W", extensions=[".fits"], destination=Path(directory).resolve() / "03.LOC_TRANS")

    base_dir = Path(directory).resolve()

    in_img2sam_wfc3uv = 'IN.img2sam_wfc3uv'

    base_dir_one = base_dir / '01.XYM'
    f606, f814 = _hstpass_catalog_counts(directory)
    use_1exp = len(f606) == 1
    output_file_dir = base_dir / '03.LOC_TRANS' 

    output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
    output_run_xym2mat = os.path.join(output_file_dir, "run_xym2mat.src")
    write_run_xym2mat_src_loc_trans(
        Path(output_run_xym2mat), base_dir_one, use_1exp=use_1exp
    )
    write_in_img2sam_wfc3uv_loc_trans(camera,Path(output_file_img2sam), base_dir / "00.DATA")

    return

def loc_trans(camera,directory):
    """
    Run the local transformation scripts in F814W and F606W subdirectories to extract 
    the pixels from each exposure and accurately transform their locations into the 
    reference frame so that we can use them to solve for a PSF and then use this PSF to model the target star.

    Automatically selects ``1exp_xym2mat_new.e`` when ``01.XYM`` has one F606W and
    two F814W hst1pass catalogs (single-F606W mode).

    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Produces pixel data for PSF generation.
    """

    def run_xym2mat(directory, script='run_xym2mat.src'):
        """Produces TRANS.xym2mat, as well as 16 MAT.0 files."""
        log_file = Path(directory).resolve() / "03.LOC_TRANS" / "log_files" / f"loc_trans_{script.replace('.src','')}.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                subdir = base_dir 
                script_path = subdir / script
                subprocess.run(
                    ["bash", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    def run_img2extract_wfc3uv_psflist(directory):
        """Run extraction for PSF list generation (simulation)."""
        log_file = Path(directory).resolve() /  "03.LOC_TRANS" / "log_files" / f"loc_trans_run_img2extract_wfc3uv_psflist.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                subdir = base_dir 
                script = 'run_img2extract_wfc3uv_psflist_simst.src' 
                script_path = base_dir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    def run_img2extract_wfc3uv_psflist_Cal(directory):
        """Run extraction for PSF list generation (calibration)."""

        log_file = Path(directory).resolve() /   "03.LOC_TRANS" / "log_files" / f"loc_trans_run_img2extract_wfc3uv_psflist_Cal.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                subdir = base_dir 
                script = 'run_img2extract_wfc3uv_psflist_Cal.src' 
                script_path = base_dir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                    )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                
    data_prep_loc_trans(camera, directory)
    run_xym2mat(directory)
    run_img2extract_wfc3uv_psflist(directory)
    run_img2extract_wfc3uv_psflist_Cal(directory)


def keck_trans(directory):
    """
    Run the Keck transformation scripts in F814W and F606W subdirectories.
    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Produces pixel data for PSF generation.
    """

    def run_matched_HST_Keck_stars_F814W(directory, script='run_matched_HST_Keck_stars_F814W.src'):
        """Produces matched HST keck stars."""
        log_file = Path(directory).resolve() / "05.KECK_TRANS" / "log_files" / f"run_matched_HST_Keck_stars_F814W.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "05.KECK_TRANS"
                subdir = base_dir
                script_path = subdir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    def run_matched_HST_Keck_stars_F606W(directory, script='run_matched_HST_Keck_stars_F606W.src'):
        """Produces matched HST keck stars."""
        log_file = Path(directory).resolve() / "05.KECK_TRANS" / "log_files" / f"run_matched_HST_Keck_stars_F606W.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "05.KECK_TRANS"
                subdir = base_dir
                script_path = subdir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    def run_HST_Keck_coord_trans_F814W(directory, script='run_HST_Keck_coord_trans_F814W.src'):
        """Produces matched HST keck stars."""
        log_file = Path(directory).resolve() / "05.KECK_TRANS" / "log_files" / f"run_HST_Keck_coord_trans_F814W.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "05.KECK_TRANS"
                subdir = base_dir
                script_path = subdir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    def run_HST_Keck_coord_trans_F606W(directory, script='run_HST_Keck_coord_trans_F606W.src'):
        """Produces matched HST keck stars."""
        log_file = Path(directory).resolve() / "05.KECK_TRANS" / "log_files" / f"run_HST_Keck_coord_trans_F606W.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "05.KECK_TRANS"
                subdir = base_dir
                script_path = subdir / script
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    run_matched_HST_Keck_stars_F814W(directory)
    run_matched_HST_Keck_stars_F606W(directory)
    run_HST_Keck_coord_trans_F814W(directory)
    run_HST_Keck_coord_trans_F606W(directory)



def cmd_partition_matchup_raw_lines(path: Path):
    """Split MATCHUP into preamble (# / blank) and data rows; omit xym2bar echo lines (fixed-width data unchanged)."""

    def _is_xym2bar_echo(line: str) -> bool:
        """Echo of xym2bar's command/input deck — drop; keep real column # headers."""
        s = line.strip()
        if not s.startswith("#"):
            return False
        inner = s.lstrip("#").strip()
        if inner and re.fullmatch(r"[-]+", inner):
            return True
        up = inner.upper()
        return up.startswith("ARG") or up.startswith("INP")

    preamble, data = [], []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.rstrip("\r\n")
            if not s.strip():
                preamble.append(s)
            elif s.lstrip().startswith("#"):
                if _is_xym2bar_echo(s):
                    continue
                preamble.append(s)
            else:
                data.append(s)
    return preamble, data


def cmd_write_matchup_raw_lines(path: Path, preamble: list[str], data: list[str]) -> None:
    path.write_text("\n".join(preamble + data) + "\n", encoding="utf-8")


def cmd_rewrite_matchup_drop_xym2bar_echo(path: Path) -> None:
    """Rewrite file without xym2bar echoes (same data line bytes)."""
    if not path.is_file():
        return
    pre, dat = cmd_partition_matchup_raw_lines(path)
    cmd_write_matchup_raw_lines(path, pre, dat)



def resolve_dex_xyvieeee_path(directory):
    """Return dex_no_gaia_STEP08_A.xyvieeee from 07.CALIBRATION, 02.CMD, or 01.XYM."""
    root = Path(directory).resolve()
    for sub in ("07.CALIBRATION", "02.CMD", "01.XYM"):
        candidate = root / sub / "dex_no_gaia_STEP08_A.xyvieeee"
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"dex_no_gaia_STEP08_A.xyvieeee not found under "
        f"{root}/07.CALIBRATION, {root}/02.CMD, or {root}/01.XYM"
    )


_MATCHUP_DATA_LINE_LEN = 133


def _format_matchup_data_line(
    xbar, ybar, mbar, xsig, ysig, msig, detection_count, star_row, qbar=9.9990
):
    """
    One xym2bar MATCHUP row (format 125, 133 chars).

    pkp/pkn/pku occupy the final three i4 fields; pku is column (17).  Cal-star
    tags are appended after column 133 so VI_HST_ogle_man_match4 reads them from
    cols 135-140 without confusing pkp with the calibration number.
    """
    nim = max(1, int(round(detection_count)))
    pki = int(round(xbar))
    pkj = int(round(ybar))
    line = (
        f"{xbar:11.4f} {ybar:11.4f} {mbar:11.4f} {xsig:11.4f} {ysig:11.4f} "
        f"{msig:11.4f} {qbar:11.4f} {nim:4d}{nim:4d}{nim:4d}{nim:4d} "
        f"N{star_row:06d} {pki:5d} {pkj:5d} {nim:4d}{nim:4d}{nim:4d}"
    )
    if len(line) != _MATCHUP_DATA_LINE_LEN:
        raise ValueError(
            f"MATCHUP row length {len(line)} != {_MATCHUP_DATA_LINE_LEN}: {line!r}"
        )
    return line


def _parse_matchup_data_line(line: str):
    """Parse one Jay Anderson MATCHUP data row; drop trailing cal annotation."""
    base = line.rstrip()
    if len(base) > _MATCHUP_DATA_LINE_LEN:
        tail = base[_MATCHUP_DATA_LINE_LEN:].strip()
        if tail.isdigit():
            base = base[:_MATCHUP_DATA_LINE_LEN]
    parts = base.split()
    if len(parts) < 17:
        return None
    star_token = parts[11]
    if not star_token.startswith("N"):
        return None
    return {
        "xbar": float(parts[0]),
        "ybar": float(parts[1]),
        "mbar": float(parts[2]),
        "xsig": float(parts[3]),
        "ysig": float(parts[4]),
        "msig": float(parts[5]),
        "qbar": float(parts[6]),
        "detection_count": int(parts[7]),
        "star_row": int(star_token[1:]),
    }


def _matchup_line_from_fields(fields, *, xbar=None, ybar=None) -> str:
    """Rebuild a fixed-width MATCHUP row from parsed fields."""
    return _format_matchup_data_line(
        xbar if xbar is not None else fields["xbar"],
        ybar if ybar is not None else fields["ybar"],
        fields["mbar"],
        fields["xsig"],
        fields["ysig"],
        fields["msig"],
        fields["detection_count"],
        fields["star_row"],
        fields["qbar"],
    )


def _annotate_matchup_cal_line(base_line: str, cal_num: int) -> str:
    """Append cal-star index after the 133-char xym2bar row."""
    base = base_line.rstrip()
    if len(base) > _MATCHUP_DATA_LINE_LEN and base[_MATCHUP_DATA_LINE_LEN:].strip().isdigit():
        base = base[:_MATCHUP_DATA_LINE_LEN]
    elif len(base) < _MATCHUP_DATA_LINE_LEN:
        fields = _parse_matchup_data_line(base)
        if fields is None:
            return base_line
        base = _matchup_line_from_fields(fields)
    if len(base) < _MATCHUP_DATA_LINE_LEN:
        base = base + " " * (_MATCHUP_DATA_LINE_LEN - len(base))
    return base[:_MATCHUP_DATA_LINE_LEN] + f"{cal_num:7d}"


def _notfar_is_bar_space(notfar_path):
    """NOTFAR bar catalogs use x/y > 1000; uuref lists stay near the image origin."""
    data = np.atleast_2d(np.loadtxt(notfar_path, skiprows=2))
    return float(data[0, 0]) > 1000.0


def _resolve_calibration_notfar(directory, cal_dir, *, f814_matchup=None):
    """Use bar-space NOTFAR for bar MATCHUP annotation."""
    root = Path(directory).resolve()
    cal_dir = Path(cal_dir)
    dest = cal_dir / "NOTFAR_CAL_STARS.XYIVB_targ"
    if dest.is_file() and not _notfar_is_bar_space(dest):
        dest.unlink()
    candidates = [
        dest,
        root / "04.EXTRACT_PSF" / "F814W" / "NOTFAR_CAL_STARS.XYIVB_targ",
        root / "02.CMD" / "NOTFAR_CAL_STARS.XYIVB_targ",
        *sorted(root.glob("07.CALIBRATION*/NOTFAR_CAL_STARS.XYIVB_targ")),
    ]
    seen = set()
    bar_candidates = []
    for path in candidates:
        path = path.resolve()
        if path in seen or not path.is_file():
            continue
        seen.add(path)
        if _notfar_is_bar_space(path):
            bar_candidates.append(path)
    if not bar_candidates:
        raise FileNotFoundError(
            "Bar-space NOTFAR_CAL_STARS.XYIVB_targ not found. Place the bar-frame "
            f"catalog in {cal_dir} (coords > 1000) before calibration_new_matchup."
        )
    chosen = bar_candidates[0]
    if chosen.resolve() != dest.resolve():
        shutil.copy2(chosen, dest)
    return dest


def _ensure_psf_mcmc_fit_sky(directory):
    """psf_star_mags_mcmc needs sky_model=1 (line 7) to write A1_min for I_hfs."""
    cal_dir = Path(directory).resolve() / "07.CALIBRATION"
    for name in ("IN.psf_star_mags_mcmc_I", "IN.psf_star_mags_mcmc_V"):
        path = cal_dir / name
        if not path.is_file():
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) < 7:
            continue
        if lines[6].strip() != "1":
            lines[6] = "1"
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _matchup_has_bar_anchor(matchup_path):
    """True when the catalog already carries the IN-file bar anchor near 6252/4738."""
    for line in Path(matchup_path).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = _parse_matchup_data_line(line)
        if fields is None:
            continue
        return fields["xbar"] > 5000.0
    return False


def _reformat_matchup_exposure_inplace(matchup_path, detection_count):
    """Rewrite MATCHUP rows in place, preserving the 133-char xym2bar layout."""
    matchup_path = Path(matchup_path)
    preamble, data_lines = cmd_partition_matchup_raw_lines(matchup_path)
    if not preamble:
        preamble = MATCHUP_JAY_HEADER.rstrip("\n").split("\n")
    out_lines = []
    for line in data_lines:
        fields = _parse_matchup_data_line(line)
        if fields is None:
            out_lines.append(line)
            continue
        fields["detection_count"] = detection_count
        rebuilt = _matchup_line_from_fields(fields)
        tail = line.rstrip()
        if len(tail) > _MATCHUP_DATA_LINE_LEN and tail[_MATCHUP_DATA_LINE_LEN:].strip().isdigit():
            rebuilt += tail[_MATCHUP_DATA_LINE_LEN:]
        out_lines.append(rebuilt)
    cmd_write_matchup_raw_lines(matchup_path, preamble, out_lines)
    return matchup_path


def _matchup_detection_count(lv, li, band):
    """Use dex_no_gaia exposure counts directly for synthetic MATCHUP rows."""
    if str(band).upper() == "F606W":
        return max(1, int(round(lv)))
    return max(1, int(round(li)))


def _find_matchup_catalog_candidates(directory, band):
    """Return existing MATCHUP catalogs for a field (most specific first)."""
    root = Path(directory).resolve()
    name = "MATCHUP.F814W.XYM.02" if band == "F814W" else "MATCHUP.F606W.XYM"
    filt = "F814W" if band == "F814W" else "F606W"
    return [
        p
        for p in (
            root / "07.CALIBRATION" / name,
            root / "01.XYM" / filt / name,
            root / "02.CMD" / name,
            root / name,
        )
        if p.is_file()
    ]


def _notfar_position_match_fraction(matchup_path, notfar_path, *, match_tol=0.1):
    """Fraction of NOTFAR cal stars that match MATCHUP x/y (uuref frame)."""
    matchup_path = Path(matchup_path)
    notfar_path = Path(notfar_path)
    if not matchup_path.is_file() or not notfar_path.is_file():
        return 0.0
    try:
        notfar = np.atleast_2d(np.loadtxt(notfar_path, skiprows=2))
        if notfar.size == 0:
            return 0.0
        stars = []
        for line in cmd_partition_matchup_raw_lines(matchup_path)[1]:
            fields = _parse_matchup_data_line(line)
            if fields is not None:
                stars.append((fields["xbar"], fields["ybar"]))
        if not stars:
            return 0.0
        tol2 = match_tol ** 2
        hits = 0
        for row in notfar:
            x_ref, y_ref = float(row[0]), float(row[1])
            if min((x - x_ref) ** 2 + (y - y_ref) ** 2 for x, y in stars) < tol2:
                hits += 1
        return hits / notfar.shape[0]
    except Exception:
        return 0.0


def _resolve_matchup_catalog_source(directory, band, notfar_path=None):
    """
    Pick the best MATCHUP catalog for calibration.

    Prefer catalogs already in 07.CALIBRATION for this field. Fall back to
    xym2bar outputs that agree with NOTFAR positions, then 02.CMD / field root.
    Dex-built bar-coordinate catalogs are only used when no suitable catalog exists.
    """
    root = Path(directory).resolve()
    if notfar_path is None:
        for candidate in (
            root / "07.CALIBRATION" / "NOTFAR_CAL_STARS.XYIVB_targ",
            root / "04.EXTRACT_PSF" / "F814W" / "NOTFAR_CAL_STARS.XYIVB_targ",
            root / "02.CMD" / "NOTFAR_CAL_STARS.XYIVB_targ",
        ):
            if candidate.is_file():
                notfar_path = candidate
                break

    cal_name = "MATCHUP.F814W.XYM.02" if band == "F814W" else "MATCHUP.F606W.XYM"
    cal_matchup = root / "07.CALIBRATION" / cal_name
    if cal_matchup.is_file():
        return cal_matchup, 1.0

    best = None
    best_score = -1.0
    for candidate in _find_matchup_catalog_candidates(directory, band):
        if candidate.resolve() == cal_matchup.resolve():
            continue
        score = _notfar_position_match_fraction(candidate, notfar_path)
        if score > best_score:
            best_score = score
            best = candidate
    return best, best_score


def _matchup_needs_rebuild_from_dex(matchup_path, dex_path, band, *, notfar_path=None):
    """
    Rebuild obviously synthetic/bad MATCHUP files.

    For 1- or 2-exposure datasets the dex_no_gaia counts should be small.
    Older placeholder catalogs used hard-coded values like 10/8/10/8, which
    break downstream assumptions in calibration.

    Do not replace xym2bar catalogs that already match NOTFAR cal-star positions.
    """
    matchup_path = Path(matchup_path)
    dex_path = Path(dex_path)
    if not matchup_path.is_file():
        return True

    if _notfar_position_match_fraction(matchup_path, notfar_path) >= 0.8:
        return False

    try:
        dex = np.loadtxt(dex_path)
        if dex.ndim == 1:
            dex = dex.reshape(1, -1)
        max_expected = 1
        for row in dex:
            max_expected = max(
                max_expected, _matchup_detection_count(float(row[12]), float(row[13]), band)
            )

        with open(matchup_path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 17:
                    return True
                counts = [int(parts[7]), int(parts[8]), int(parts[9]), int(parts[10])]
                return any(count > max_expected for count in counts)
    except Exception:
        return True

    return False


def write_matchup_from_dex_xyvieeee(dex_path, out_path, band="F814W"):
    """
    dex_no_gaia_STEP08_A.xyvieeee
    -----------------------------
    Each line from dex_no_gaia STEP08 looks like::

        1559.795 4389.189   -7.284   -8.555      0.015  0.039  9.990  9.990   268.6038096  -28.7303315  2 2 1 1 000001

    Column mapping to MATCHUP:

    ====  ============  ==========================================
    Col   Field         Used for MATCHUP
    ====  ============  ==========================================
    0     uu_bar        x position (xbar)
    1     vv_bar        y position (ybar)
    2     mv_bar        V magnitude (F606W)
    3     mi_bar        I magnitude (F814W)
    4-5   uu_sig, vv_sig  xsig, ysig
    6-7   mv_sig, mi_sig  msig (V or I, depending on band)
    12-13 Lv, Li        exposure-count fields -> Nf, Ng, Nm, Nmin
    14    star id (Q)   not used; row number i becomes N000001, etc.
    ====  ============  ==========================================

    Output:
    --------------
    - Header: 4-line Jay Anderson header (MATCHUP_JAY_HEADER)
    - One data row per dex line, formatted like xym2bar format 125::

        1559.7950   4389.1890     -8.5550      0.0150      0.0390      0.0390      9.9990   10    8   10    8 N000001  1560  4389   10   10    8

    Mapping logic:

    - band="F814W" (default): mbar = col 3 (I), msig = col 7
    - band="F606W": mbar = col 2 (V), msig = col 6
    - Nstar: N{i:06d} where i is the 1-based row in the dex file (row 1 -> N000001)
    - Nf/Ng/Nm/Nmin are set from the dex exposure counts for the requested band
      (Lv for F606W, Li for F814W), so 1- and 2-exposure fields stay consistent
      with the actual dataset.

    Parameters
    ----------
    dex_path : path to dex_no_gaia_STEP08_A.xyvieeee
    out_path : output MATCHUP file (e.g. MATCHUP.F814W.XYM.02 or MATCHUP.F606W.XYM)
    band : ``F814W`` uses I magnitudes (col 3); ``F606W`` uses V magnitudes (col 2)
    """
    dex_path = Path(dex_path)
    out_path = Path(out_path)
    data = np.loadtxt(dex_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    use_v = str(band).upper() == "F606W"
    out_lines = MATCHUP_JAY_HEADER.rstrip("\n").split("\n")
    for i, row in enumerate(data, start=1):
        xbar, ybar = float(row[0]), float(row[1])
        mbar = float(row[2] if use_v else row[3])
        xsig, ysig = float(row[4]), float(row[5])
        msig = float(row[6] if use_v else row[7])
        lv, li = float(row[12]), float(row[13])
        detection_count = _matchup_detection_count(lv, li, band)
        out_lines.append(
            _format_matchup_data_line(
                xbar, ybar, mbar, xsig, ysig, msig, detection_count, i
            )
        )
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return out_path


def write_cal_matchup_files(matchup_in, notfar_path, out_cal, out_cal_only, *, match_tol=0.1):
    """
  Annotate a MATCHUP catalog with calibration-star numbers from NOTFAR_CAL_STARS.

  Matches cal stars by position (same criterion as cal_star_num_2_MATCHUP).
  Falls back to the iMuref row index when the catalog was reordered after
  NOTFAR was written. MATCHUP bar coordinates are always preserved; NOTFAR
  positions live in a different frame and must not overwrite them.
  """
    matchup_in = Path(matchup_in)
    notfar_path = Path(notfar_path)
    out_cal = Path(out_cal)
    out_cal_only = Path(out_cal_only)

    preamble, data_lines = cmd_partition_matchup_raw_lines(matchup_in)
    if not preamble:
        preamble = MATCHUP_JAY_HEADER.rstrip("\n").split("\n")

    notfar = np.atleast_2d(np.loadtxt(notfar_path, skiprows=2))
    parsed = []
    for line in data_lines:
        fields = _parse_matchup_data_line(line)
        if fields is None:
            parsed.append(None)
            continue
        formatted = _matchup_line_from_fields(fields)
        parsed.append({**fields, "line": formatted})

    cal_only_rows = []
    tol2 = match_tol ** 2
    for ic in range(notfar.shape[0]):
        x_ref, y_ref = float(notfar[ic, 0]), float(notfar[ic, 1])
        i_muref = int(notfar[ic, 5]) if notfar.shape[1] > 5 else 0
        hit = None
        for j, rec in enumerate(parsed):
            if rec is None:
                continue
            dr2 = (rec["xbar"] - x_ref) ** 2 + (rec["ybar"] - y_ref) ** 2
            if dr2 < tol2:
                hit = j
                break
        if hit is None:
            continue
        rec = parsed[hit]
        if rec is None:
            continue
        annotated = _annotate_matchup_cal_line(rec["line"], ic + 1)
        parsed[hit]["line"] = annotated
        mbar = rec["mbar"]
        cal_only_rows.append(
            f"{rec['xbar']:11.4f} {rec['ybar']:11.4f} {float(mbar):11.4f}{ic + 1:7d}"
        )

    out_data = [rec["line"] if rec else "" for rec in parsed]
    cmd_write_matchup_raw_lines(out_cal, preamble, out_data)
    cal_only_body = MATCHUP_CAL_ONLY_HEADER + "".join(
        row + "\n" for row in cal_only_rows
    )
    out_cal_only.write_text(cal_only_body, encoding="utf-8")
    return out_cal, out_cal_only


def ensure_matchup_catalogs_from_dex(directory):
    """
  Ensure 02.CMD contains MATCHUP.F814W.XYM.02 and MATCHUP.F606W.XYM.

  Prefer xym2bar catalogs from 01.XYM or the field root when they agree with
  NOTFAR cal-star positions. Only build from dex_no_gaia_STEP08_A.xyvieeee as
  a last resort (dex bar coordinates/magnitudes break VI_HST_ogle_man_match4).
  """
    root = Path(directory).resolve()
    cmd_dir = root / "02.CMD"
    cmd_dir.mkdir(parents=True, exist_ok=True)
    dex_path = resolve_dex_xyvieeee_path(directory)

    notfar = None
    for candidate in (
        root / "07.CALIBRATION" / "NOTFAR_CAL_STARS.XYIVB_targ",
        root / "04.EXTRACT_PSF" / "F814W" / "NOTFAR_CAL_STARS.XYIVB_targ",
        root / "02.CMD" / "NOTFAR_CAL_STARS.XYIVB_targ",
    ):
        if candidate.is_file():
            notfar = candidate
            break

    f814 = cmd_dir / "MATCHUP.F814W.XYM.02"
    f606 = cmd_dir / "MATCHUP.F606W.XYM"
    for band, dest in (("F814W", f814), ("F606W", f606)):
        source, score = _resolve_matchup_catalog_source(directory, band, notfar)
        if source is not None and source.resolve() != dest.resolve():
            shutil.copy2(source, dest)
            print(f"Using {source} for {dest.name} (NOTFAR match {score:.0%})")
        if _matchup_needs_rebuild_from_dex(dest, dex_path, band, notfar_path=notfar):
            write_matchup_from_dex_xyvieeee(dex_path, dest, band=band)
            print(f"Wrote {dest} from {dex_path}")
    return f814, f606, dex_path


def cmd_diagram(directory):

#    fortran_src = get_fortran_dir()
    copy_entire_files(source=Path(directory).resolve() / "01.XYM", destination=Path(directory).resolve() / "02.CMD", filename="dex_no_gaia_STEP08_A.xyvieeee")

    subdir = Path(directory).resolve() / "02.CMD"
    path_match_I = subdir / "dex_no_gaia_STEP08_A.xyvieeee"
    #cmd_rewrite_matchup_drop_xym2bar_echo(path_match_I)
    xv, yv, mv, mi = np.loadtxt(path_match_I, unpack=True, usecols=(0, 1, 2, 3))
    xi = xv
    yi = yv
    #Establish target parameters
    response = str(input("Do you have a target? Enter 'Yes' if you do."))
    if response == 'yes' or response == 'Yes':
        xtarg = float(input("Enter x-coord of your target (from dex_no_gaia_STEP08_A.xyvieeee)"))
        ytarg = float(input("Enter y-coord of your target (from dex_no_gaia_STEP08_A.xyvieeee)"))
        Vtarg = float(input("Enter V magnitude of your target (from dex_no_gaia_STEP08_A.xyvieeee)"))
        Itarg = float(input("Enter I magnitude of your target (from dex_no_gaia_STEP08_A.xyvieeee)"))
    else:
#        pdb.set_trace()
        xtarg, ytarg = xi[0], yi[0]
        Vtarg, Itarg = mv[0], mi[0]
#        xtarg, ytarg = 535.2320, 623.8950
#        Vtarg, Itarg = -10.2651, -10.3052
    VmItarg = Vtarg - Itarg

    # If the user supplied a target, put that star first in the MATCHUP files (preserve fixed-width lines).
    if response == "yes" or response == "Yes":
        path_I = path_match_I
        #path_V = path_match_V
        preamble_I, lines_I = cmd_partition_matchup_raw_lines(path_I)
        #preamble_V, lines_V = cmd_partition_matchup_raw_lines(path_V)
        n = len(xi)
        dist2 = (xi - xtarg) ** 2 + (yi - ytarg) ** 2
        idx = int(np.argmin(dist2))
        order_list = [idx] + [i for i in range(n) if i != idx]
        if idx != 0:
            reordered_I = [lines_I[i] for i in order_list]
            #reordered_V = [lines_V[i] for i in order_list]
            cmd_write_matchup_raw_lines(path_I, preamble_I, reordered_I)
            #cmd_write_matchup_raw_lines(path_V, preamble_V, reordered_V)
        xv, yv, mv, mi = np.loadtxt(path_match_I, unpack=True, usecols=(0, 1, 2, 3))
        xi = xv
        yi = yv

    #Function to find the CMD of the target and get our list of Sim+Ref stars.

    def show_cmd_targ(directory, response, xi, yi, xv, yv, mi, mv):
    
        #Plotting parameters 
        if response == 'yes' or response == 'Yes':
    
            rad_max = float(input("Enter maximum plotting radius for PSF selection (recommended: 150)"))
            box_max = float(input("Enter maximum box radius for PSF selection (recommended: 150)"))
            mag_range = float(input("Enter magnitude range for plotting of PSF selection (recommended: 0.60)"))
            col_range = float(input("Enter color range for plotting of PSF selection (recommended: 0.60)"))
            ref_st_Imx = float(input("Enter reference star input I max (recommended: F814W mag + 2)"))
            ref_st_Imn = float(input("Enter reference star input I min (recommended: F814W mag - 2)"))
            ref_st_Vmx = float(input("Enter reference star input V max (recommended: F606W mag + 2)"))
            ref_st_Vmn = float(input("Enter reference star input V min (recommended: F606W mag - 2)"))
        else:
            rad_max, box_max = 300, 300
            mag_range, col_range = 0.50, 0.30
            ref_st_Imx, ref_st_Imn = mi + 4, mi -4 
            ref_st_Vmx, ref_st_Vmn = mv + 4. , mv-4
        
       # rad_max, box_max = 300, 300
        #mag_range, col_range = 0.50, 0.30
        #ref_st_Imx, ref_st_Imn = -8.8, -12.75
        #ref_st_Vmx, ref_st_Vmn = -8.7, -12.75
    
        d = np.sqrt((xi-xtarg)**2 + (yi-ytarg)**2)
        n = np.arange(1, len(mi)+1)
    
        #Plotting parameter tlo decide how many stars around the target should be selected 
    
        vprox = np.abs(mi - Itarg) < mag_range
        cprox = np.abs(mv - mi - VmItarg) < col_range
    
        #vprox and cprox are boolean arrays. Stars close to our target within a given window 
    
        u = vprox & cprox #These are the stars close to the target in both magnitude and color space.
    
    
        uref = (d < rad_max) & cprox & (mv < ref_st_Vmx) & (mv > ref_st_Vmn) & (mi < ref_st_Imx) & (mi > ref_st_Imn) & (n > 1) #More selective than u
    
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))
        ax_cmd = ax[0,0]
        ax_xy_I = ax[1,1]
        ax_xy_zoom = ax[1,0]
        ax_xy_V = ax[0, 1]
    
        # CMD 
        ax_cmd.scatter(mv-mi, mi, s=5, c='k', label='All Stars')
        ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='purple', label='Selected Stars')
        ax_cmd.scatter([VmItarg], [Itarg], marker='x', lw = 5, s=100, c='grey', label='Target')
        ax_cmd.set_xlim(-0.75, 1.75)
        ax_cmd.set_ylim(-15, -7)   
        ax_cmd.set_xlabel('F606W - F814W')
        ax_cmd.set_ylabel('F814W')
        ax_cmd.legend(loc = 'lower right')
        ax_cmd.set_title('CMD')
    
        # XY for I filter
        ax_xy_I.scatter(xi, yi, s=5, c='k', label = 'All Stars') # All stars
        #ax_xy_I.scatter(xi[u], yi[u], s=15, c='r') # I avoid this from SM and prefer uref instead
        ax_xy_I.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual selected stars')
        ax_xy_I.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target') #Target
        circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2) 
        ax_xy_I.add_patch(circle)
        ax_xy_I.set_xlim(xtarg-box_max, xtarg+box_max)
        ax_xy_I.set_ylim(ytarg-box_max, ytarg+box_max)
        ax_xy_I.set_xlabel('x coord')
        ax_xy_I.set_ylabel('y coord')
        ax_xy_I.legend(loc = 'lower right')
        ax_xy_I.set_title('XY Coord I-Filter')
    
        # XY for V filter
        ax_xy_V.scatter(xv, yv, s=5, c='k', label = 'All Stars') # All stars
        ax_xy_V.scatter(xv[uref], yv[uref], s=30, c='purple', label = 'Actual selected stars')
        ax_xy_V.scatter([xtarg], [ytarg], marker='x',lw = 5, s=100, c='grey', label = 'Target') #Target
        circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2) 
        ax_xy_V.add_patch(circle)
        ax_xy_V.set_xlim(xtarg-box_max, xtarg+box_max)
        ax_xy_V.set_ylim(ytarg-box_max, ytarg+box_max)
        ax_xy_V.set_xlabel('x coord')
        ax_xy_V.set_ylabel('y coord')
        ax_xy_V.legend(loc = 'lower right')
        ax_xy_V.set_title('XY Coord V-Filter')
    
    
        # XY zoom for I-filter
        mask_zoom = u & (d < rad_max)
        ax_xy_zoom.scatter(xi, yi, s=5, c='k', label = 'All Stars')
        #ax_xy_zoom.scatter(xi[mask_zoom], yi[mask_zoom], s=15, c='r')
        ax_xy_zoom.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual selected stars')
        ax_xy_zoom.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target')
        circle2 = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2)
        ax_xy_zoom.add_patch(circle2)
        ax_xy_zoom.set_xlim(xtarg-box_max, xtarg+box_max)
        ax_xy_zoom.set_ylim(ytarg-box_max, ytarg+box_max)
        ax_xy_zoom.set_xlabel('x coord')
        ax_xy_zoom.set_ylabel('y coord')
        ax_xy_zoom.legend(loc='lower right')
        ax_xy_zoom.set_title('XY Coord I-Filter Zoom')
    
        plt.tight_layout()
        plt.savefig(Path(directory).resolve() / "02.CMD" / "show_cmd_targ.pdf")
        plt.close(fig)
    
        xu = xi[u & (d < rad_max)]
        yu =  yi[u & (d < rad_max)]
        miu = mi[u & (d < rad_max)]
        mvu = mv[u & (d < rad_max)]
    
        upsf = np.ones_like(xu, dtype=int)
    
        assert len(upsf) > 0
    #    if len(upsf) > 0 could be another condition
    
        upsf[0] = 0  # first star (the target) can't be used for PSF 
     
    
        np.savetxt(Path(directory).resolve() / "02.CMD" / 'NEARBY_SIM_STARS.XYIVB_targ', np.column_stack([xu, yu, miu, mvu, upsf]), fmt='%10.3f %10.3f %8.4f %8.4f %1d', header="xu         yu      miu      mvu   upsf \n")
    
        np.savetxt(Path(directory).resolve() / "02.CMD" / 'NEARBY_REF_STARS.XYIVB_targ', np.column_stack([xi[uref], yi[uref], mi[uref], mv[uref]]), fmt='%10.3f %10.3f %8.4f %8.4f', header = "xuref      yuref   miuref      mvuref \n")
    
    #Function to give calibration stars
    def show_cmd_Cal(directory, response, xi, yi, xv, yv, mi, mv):
        #Plotting parameters 
        if response == 'yes' or response == 'Yes':

            rad_max = float(input("Enter maximum plotting radius for calibration star selection (recommended: 300)"))
            box_max = float(input("Enter maximum box radius for calibration star selection (recommended: 300)"))
            Vcalc = float(input("Enter F606W magnitude for calibration (recommended: -13)"))
            Icalc = float(input("Enter I magnitude for calibration (recommended: -13)"))
            mag_range = float(input("Enter magnitude range for plotting (recommended: 1.0)"))
            col_range = float(input("Enter color range for plotting (recommended: 1.0)"))
        else:
            rad_max, box_max = 300, 300
            mag_range, col_range = 0.50, 0.30
            Vcalc, Icalc = -10, -10
        ref_st_Imx, ref_st_Imn = Icalc + mag_range, Icalc - mag_range
        ref_st_Vmx, ref_st_Vmn = Vcalc + mag_range, Vcalc - mag_range

        d = np.sqrt((xi-xtarg)**2 + (yi-ytarg)**2)
        vprox = np.abs(mi - Icalc) < mag_range
        cprox = np.abs(mv - mi - VmItarg) < col_range
    
        u = vprox & cprox
    
        uref = (d < rad_max) & cprox & (mv < ref_st_Vmx) & (mv > ref_st_Vmn) & (mi < ref_st_Imx) & (mi > ref_st_Imn)
    
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    
        # CMD
        ax_cmd = ax[0]
        ax_cmd.scatter(mv-mi, mi, s=5, c='k', label='All Stars')
    #    ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='r', label='Selected Stars')
        ax_cmd.scatter((mv-mi)[u], mi[u], s=15, c='purple', label='Selected Stars')
        ax_cmd.scatter([VmItarg], [Itarg], marker='x', lw =5, s=100, c='grey', label='Target')
        ax_cmd.set_xlim(-0.75, 1.25)
        ax_cmd.set_ylim(-15, -7)
        ax_cmd.set_xlabel('F606W - F814W')
        ax_cmd.set_ylabel('F814W')
        ax_cmd.set_title('CMD')
        ax_cmd.legend(loc = 'lower right')
    
        # XY
        ax_xy = ax[1]
        ax_xy.scatter(xi, yi, s=5, c='k', label = 'All Stars')
        #ax_xy.scatter(xi[u], yi[u], s=15, c='purple')
        ax_xy.scatter(xi[uref], yi[uref], s=30, c='purple', label = 'Actual Selected Stars')
        ax_xy.scatter([xtarg], [ytarg], marker='x', lw = 5, s=100, c='grey', label = 'Target')
        circle = plt.Circle((xtarg, ytarg), rad_max, color='hotpink', fill=False, lw=2)
        ax_xy.add_patch(circle)
        ax_xy.set_xlim(xtarg-box_max, xtarg+box_max)
        ax_xy.set_ylim(ytarg-box_max, ytarg+box_max)
        ax_xy.set_xlabel('x coord')
        ax_xy.set_ylabel('y coord')
        ax_xy.set_title('Calibration Stars')
        ax_xy.legend(loc = 'lower right')
        plt.tight_layout()
        plt.savefig(Path(directory).resolve() / "02.CMD" / 'show_cmd_Cal.pdf')
        plt.close(fig)
    
        xuref = xi[uref]
        yuref =  yi[uref]
        miuref = mi[uref]
        mvuref = mv[uref]
    
        iMuref = np.arange(1, len(xi) + 1)[uref]
        upsf = np.ones_like(xuref, dtype=int)
    
        np.savetxt(Path(directory).resolve() / "02.CMD" / 'NOTFAR_CAL_STARS.XYIVB_targ',  np.column_stack([xuref, yuref, miuref, mvuref, upsf, iMuref]), fmt='%10.3f %10.3f %8.4f %8.4f %1d %6d', header=" #  xuref      yuref    miuref   mvuref  upsf iMuref  \n")

    show_cmd_targ(directory, response, xi, yi, xv, yv, mi, mv)

    show_cmd_Cal(directory, response, xi, yi, xv, yv, mi, mv)
    
    
        
    



def extract_psf_1(directory):
    """
    Generate a local PSF for each filter, using the stars that are similar in magnitude and color to the target star. The magnitude similarity is thought to be important because the CTE losses are expected to make PSF shapes magnitude dependent. 
    The initial selection of PSF stars in directory 02.CMD_final did not take into account any blending. 

    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def run_uvp2psf_simst(directory, iteration = 1, script='run_uvp2psf_simst_1.src'):
        """Finds a sky value for each star in each exposure using the pixels between 8.5 and 13.5 pixels of the center."""
        log_file = Path(directory).resolve() / "04.EXTRACT_PSF" / "log_files" / f"run_uvp2psf_simst_1.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "04.EXTRACT_PSF"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_uvp2psf_simst_1.src' if f == "F814W" else 'run_uvp2psf_simstV_1.src'
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W")
    sim_stars_file = Path(directory).resolve() / "02.CMD" / "NEARBY_SIM_STARS.XYIVB_targ"
    sim_stars = np.atleast_2d(np.loadtxt(sim_stars_file, skiprows=1))
    num_simstars = sim_stars.shape[0]
    def prepare_data(images, directory, f= 'F814W'):
        base_dir = Path(directory).resolve()
        subdir = base_dir / f
        in_good_psf_list = 'IN.good_psf_list.1'
        output_file_dir = base_dir / '04.EXTRACT_PSF' / f
        output_file_img = os.path.join(output_file_dir, in_good_psf_list)
        with open(output_file_img, "w") as f:
            for i in range(1, images + 1):
                value = 0 if i == 1 else 1
                f.write(f"{i:2d}   {value}\n")
                
    prepare_data(num_simstars, directory)
    prepare_data(num_simstars, directory, f = 'F606W')
    run_uvp2psf_simst(directory)
        
def extract_psf_2(good_psf, directory):
    """
    Generate a local PSF for each filter, using the stars that are similar in magnitude and color to the target star. The magnitude similarity is thought to be important because the CTE losses are expected to make PSF shapes magnitude dependent. 

    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def prepare_data(good_psf, directory, f= 'F814W'):
        base_dir = Path(directory).resolve()
        subdir = base_dir / f
        in_good_psf_list = 'IN.good_psf_list.2'
        output_file_dir = base_dir / '04.EXTRACT_PSF' / f
        output_file_img = os.path.join(output_file_dir, in_good_psf_list)


        with open(output_file_img, "w") as f:
            for i in range(1, len(good_psf) + 1):
                value = good_psf[i-1]
                f.write(f"{i:2d}   {value}\n")
                
    def run_uvp2psf_simst(directory, script='run_uvp2psf_simst_2.src'):
        """Finds a sky value for each star in each exposure using the pixels between 8.5 and 13.5 pixels of the center."""
        log_file = Path(directory).resolve() / "04.EXTRACT_PSF" / "log_files" / f"run_uvp2psf_simst_2.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "04.EXTRACT_PSF"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_uvp2psf_simst_2.src' if f == "F814W" else 'run_uvp2psf_simstV_2.src'
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
    
    sim_stars_file = Path(directory).resolve() / "02.CMD" / "NEARBY_SIM_STARS.XYIVB_targ"
    sim_stars = np.atleast_2d(np.loadtxt(sim_stars_file, skiprows=1))
    num_simstars = sim_stars.shape[0]

    prepare_data(good_psf, directory)
    prepare_data(good_psf, directory, f='F606W')
    run_uvp2psf_simst(directory)



def hst_fit_dataprep_twostar(directory, f = 'F814W'):

    x1 = float(input("Initial x position for object 1 in pixels (example: -0.5)"))
    y1 = float(input("Initial y position for object 1 in pixels (example: 0.25)"))
    x2 = float(input("Initial x position for object 2 in pixels (example: 0.5)"))
    y2 = float(input("Initial y position for object 2 in pixels (example: -0.17)"))
    x3 = float(0)
    y3 = float(0)
    f1 = float(input("Initial flux for object 1 (example: 0.55)"))
    f2 = float(input("Initial flux for object 2 (example: 0.45)"))
    mcmc_dr1 = float(input("MCMC step size for object 1's position in pixels. (example: 0.01)"))
    mcmc_dr2 = float(input("MCMC step size for object 2's position in pixels. (example: 0.01)"))
    mcmc_dr3 = float(0)
    mcmc_df1 = float(input("MCMC step size for flux of object 1 (example: 0.01)"))
    mcmc_df2 = float(input("MCMC step size for flux of object 2 (example: 0.01)"))

    nmcmc = int(input("Total MCMC steps (recommended > 50000)"))
    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))
    dufitmn = float(input("Minimum du cut (minimum distance in x-direction) in pixels (example: -4.5)"))
    dufitmx = float(input("Maximum du cut (maximum distance in x-direction) in pixels (example: 4.5)"))
    dvfitmn =float(input("Minimum dv cut (minimum distance in y-direction) in pixels (example: -4.5)"))
    dvfitmx = float(input("Maximum dv cut (maximum distance in y-direction) in pixels (example: 4.5)"))
    chi2cut = float(input("Chi-square cut for pixel selection (example: 45)"))
    

    
    if f == 'F814W':
        content = [
            "psfout_simst.fits",
            "simst",
            "I_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]
    else:
        content = [
            "psfout_simstV.fits",
            "simstV",
            "V_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]

    filename = 'IN.uvp2tri_NOscon_fs_asym_mcmc'
    base_dir = Path(directory).resolve()
    output_file_dir = base_dir / '06.FIT' / f / '2star-fit'
    output_file = os.path.join(output_file_dir, filename)
    
    with open(output_file, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{filename}' successfully created.")



def hst_fit_dataprep_threestar(directory, f = 'F814W'):

    x1 = float(input("Initial x position for object 1 (example: -0.5)"))
    y1 = float(input("Initial y position for object 1 (example: 0.01"))

    x2 = float(input("Initial x position for object 2 (example: 0.4)"))
    y2 = float(input("Initial y position for object 2 (example: -0.01)"))
 
    x3 = float(input("Initial x position for object 3 (example: -0.001)"))
    y3 = float(input("Initial y position for object 3 (example: 0.004)"))

    f1 = float(input("Initial flux for object 1 (example: 0.2)"))
    f2 = float(input("Initial flux for object 2 (example: 0.4)"))

    mcmc_dr1 = float(input("MCMC step size for object 1's position in pixels. (example: 0.01)"))
    mcmc_dr2 = float(input("MCMC step size for object 2's position in pixels. (example: 0.01)"))
    mcmc_dr3 = float(input("MCMC step size for object 3's position in pixels. (example: 0.01)"))
    mcmc_df1 = float(input("MCMC step size for flux of object 1 (example: 0.01)"))
    mcmc_df2 = float(input("MCMC step size for flux of object 2 (example: 0.01)"))

    nmcmc = int(input("Total MCMC steps (recommended > 50000)"))
    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))
    dufitmn = float(input("Minimum du cut (minimum distance in x-direction) in pixels (example: -4.5)"))
    dufitmx = float(input("Maximum du cut (maximum distance in x-direction) in pixels (example: 4.5)"))
    dvfitmn =float(input("Minimum dv cut (minimum distance in y-direction) in pixels (example: -4.5)"))
    dvfitmx = float(input("Maximum dv cut (maximum distance in y-direction) in pixels (example: 4.5)"))
    chi2cut = float(input("Chi-square cut for pixel selection (example: 45)"))

    
    if f == 'F814W':
        content = [
            "psfout_simst.fits",
            "simst",
            "I_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]
    else:
        content = [
            "psfout_simstV.fits",
            "simstV",
            "V_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]

    filename = 'IN.uvp2tri_NOscon_fs_asym_mcmc'
    base_dir = Path(directory).resolve()
    output_file_dir = base_dir / '06.FIT' / f / '3star-fit'
    output_file = os.path.join(output_file_dir, filename)
    
    with open(output_file, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{filename}' successfully created.")


def hst_fit_dataprep_onestar(directory, f = 'F814W'):

    x1 = float(input("Initial x position for object 1 (example: 0.0)"))
    y1 = float(input("Initial y position for object 1 (example: 0.0)"))

    x2 = float(0)
    y2 = float(0)
 
    x3 = float(0)
    y3 = float(0)

    f1 = float(1)
    f2 = float(0)

    mcmc_dr1 = float(input("MCMC step size for object 1's position in pixels. (example: 0.01)"))
    mcmc_dr2 = float(0)
    mcmc_dr3 = float(0)
    mcmc_df1 = float(0)
    mcmc_df2 = float(0)

    nmcmc = int(input("Total MCMC steps (recommended > 50000)"))
    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))
    dufitmn = float(input("Minimum du cut (minimum distance in x-direction) in pixels (example: -4.5)"))
    dufitmx = float(input("Maximum du cut (maximum distance in x-direction) in pixels (example: 4.5)"))
    dvfitmn =float(input("Minimum dv cut (minimum distance in y-direction) in pixels (example: -4.5)"))
    dvfitmx = float(input("Maximum dv cut (maximum distance in y-direction) in pixels (example: 4.5)"))
    chi2cut = float(input("Chi-square cut for pixel selection (example: 45)"))

    if f == 'F814W':
        content = [
            "psfout_simst.fits",
            "simst",
            "I_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]
    else:
        content = [
            "psfout_simstV.fits",
            "simstV",
            "V_KeckNOcon",
            "0.0  99999.  0.0  99999.",
            f"{x1} {y1} {x2} {y2} {x3} {y3} {f1} {f2}",
            f"{mcmc_dr1} {mcmc_dr2} {mcmc_dr3} {mcmc_df1} {mcmc_df2}",
            "1.0",
            f"{nmcmc}",
            f"{dufitmn} {dufitmx} {dvfitmn} {dvfitmx} {chi2cut}"
        ]

    filename = 'IN.uvp2tri_NOscon_fs_asym_mcmc'
    base_dir = Path(directory).resolve()
    output_file_dir = base_dir / '06.FIT' / f / '1star-fit'
    output_file = os.path.join(output_file_dir, filename)
    
    with open(output_file, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{filename}' successfully created.")


def tri_fit_F814W_opt(directory):
    """
    Fit the pixels of the target star with the PSF to determine the best-fit 2 or 3-star model in the F814W filter. 
    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def strip_star_lines_from_uvp2tri_mcmc_814W(directory):
        mcmc_path = (
            Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit" / "uvp2tri_scon_fsky_I_KeckNOcon.07.mcmc")
        if not mcmc_path.is_file():
            raise FileNotFoundError(f"Expected MCMC file not found: {mcmc_path}")
        lines = mcmc_path.read_text().splitlines(keepends=True)
        kept = [ln for ln in lines if not ln.lstrip().startswith("***")]
        mcmc_path.write_text("".join(kept))

    def run_mcmc_expand_average_814W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / "06.FIT" /  "F814W" / "3star-fit" / "log_files" /  f"run_mcmc_expand_average.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

    _fit_3star = Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit"
    for name in UVP2TRI_FSKY_OUTPUTS_F814W:
        (_fit_3star / name).unlink(missing_ok=True)

    fortran_src = get_fortran_dir()
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename="mcmc_expand_average.xOg")

    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc_alt.xOg")
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", extensions=[".fits"])

    run_uvp2tri_mcmc(directory, "F814W", ["3star-fit"], UVP2TRI_FSKY_OUTPUTS_F814W)
    #run_uvp2psf_simst_2(directory)
    strip_star_lines_from_uvp2tri_mcmc_814W(directory)
    run_mcmc_expand_average_814W(directory)
    string_3star =  print_uvp2tri_mcmc_acceptance_rate_3star(directory, "F814W", fit_folders = ['3star-fit'])
    return string_3star
    #run_mcmc_expand_average_606W(directory)



def hst_fit_final_F814W(directory):
    """
    Fit the pixels of the target star with the PSF to determine the best-fit 2 or 3-star model in the F814W filter. 
    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def strip_star_lines_from_uvp2tri_mcmc_814W(directory):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            mcmc_path = (Path(directory).resolve() / "06.FIT" / "F814W" / f / "uvp2tri_scon_fsky_I_KeckNOcon.07.mcmc")
            if not mcmc_path.is_file():
                raise FileNotFoundError(f"Expected MCMC file not found: {mcmc_path}")
            lines = mcmc_path.read_text().splitlines(keepends=True)
            kept = [ln for ln in lines if not ln.lstrip().startswith("***")]
            mcmc_path.write_text("".join(kept))

    def run_mcmc_expand_average_814W(directory, script='run_mcmc_expand_average.src'):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            log_file = Path(directory).resolve() / "06.FIT" / "F814W" / f / "log_files" / f"run_mcmc_expand_average.log"
            with open(log_file, "w") as logf:
                try:
                    base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                    subdir = base_dir / f
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
                finally:
                    sys.stdout = sys.__stdout__
                    sys.stderr = sys.__stderr__

    folders = ['1star-fit', '2star-fit']
    for f in folders:
        _fit_2star = Path(directory).resolve() / "06.FIT" / "F814W" / f
        for name in UVP2TRI_FSKY_OUTPUTS_F814W:
            (_fit_2star / name).unlink(missing_ok=True)

    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit")

    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit")

    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit", extensions=[".fits"])
  

    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit", extensions=[".fits"])

    run_uvp2tri_mcmc(directory, "F814W", folders, UVP2TRI_FSKY_OUTPUTS_F814W)
    #pdb.set_trace()
    string_1star, string_2star =  print_uvp2tri_mcmc_acceptance_rate(directory, "F814W", folders)
    strip_star_lines_from_uvp2tri_mcmc_814W(directory)
    run_mcmc_expand_average_814W(directory)
    return string_1star, string_2star




def hst_fit_final_F606W(directory):
    """
    Fit the pixels of the target star with the PSF to determine the best-fit 2 or 3-star model in the 606W filter. 
    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def strip_star_lines_from_uvp2tri_mcmc_606W(directory):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            mcmc_path = (Path(directory).resolve() / "06.FIT" / "F606W" / f / "uvp2tri_scon_fsky_V_KeckNOcon.07.mcmc")
            if not mcmc_path.is_file():
                raise FileNotFoundError(f"Expected MCMC file not found: {mcmc_path}")
            lines = mcmc_path.read_text().splitlines(keepends=True)
            kept = [ln for ln in lines if not ln.lstrip().startswith("***")]
            mcmc_path.write_text("".join(kept))

        
    def run_mcmc_expand_average_606W(directory, script='run_mcmc_expand_average.src'):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            log_file = Path(directory).resolve() / "06.FIT" / "F606W" / f / "log_files" / f"run_mcmc_expand_average.log"
            with open(log_file, "w") as logf:
                try:
                    base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                    subdir = base_dir / f
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
                finally:
                    sys.stdout = sys.__stdout__
                    sys.stderr = sys.__stderr__

    folders = ['1star-fit', '2star-fit']
    for f in folders:
        _fit_2star = Path(directory).resolve() / "06.FIT" / "F606W" / f
        for name in UVP2TRI_FSKY_OUTPUTS_F606W:
            (_fit_2star / name).unlink(missing_ok=True)

    run_uvp2tri_mcmc(directory, "F606W", folders, UVP2TRI_FSKY_OUTPUTS_F606W)
    star_1fit, star_2fit = print_uvp2tri_mcmc_acceptance_rate(directory, "F606W", folders)
    strip_star_lines_from_uvp2tri_mcmc_606W(directory)
    run_mcmc_expand_average_606W(directory)
    return star_1fit, star_2fit

def tri_fit_F606W_opt(directory):
    """
    Fit the pixels of the target star with the PSF to determine the best-fit 2 or 3-star model in the 606W filter. 
    Parameters
    ----------
    directory : str or Path
        Root directory containing 00.DATA/ and 01.XYM/ folders.

    Returns
    -------
    Local PSF for each filter
    """

    def strip_star_lines_from_uvp2tri_mcmc_606W(directory):
        mcmc_path = (Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit" / "uvp2tri_scon_fsky_V_KeckNOcon.07.mcmc")
        if not mcmc_path.is_file():
            raise FileNotFoundError(f"Expected MCMC file not found: {mcmc_path}")
        lines = mcmc_path.read_text().splitlines(keepends=True)
        kept = [ln for ln in lines if not ln.lstrip().startswith("***")]
        mcmc_path.write_text("".join(kept))
        
    def run_mcmc_expand_average_606W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit" / "log_files" / f"run_mcmc_expand_average.log" 
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
                    subprocess.run(
                        ["csh", str(script_path)],
                        cwd=subdir,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
    
        
    _fit_3star = Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit"
    for name in UVP2TRI_FSKY_OUTPUTS_F606W:
        (_fit_3star / name).unlink(missing_ok=True)

    fortran_src = get_fortran_dir()
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc_alt.xOg")
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", filename="mcmc_expand_average.xOg")
    run_uvp2tri_mcmc(directory, "F606W", ["3star-fit"], UVP2TRI_FSKY_OUTPUTS_F606W)
    strip_star_lines_from_uvp2tri_mcmc_606W(directory)
    run_mcmc_expand_average_606W(directory)
    string_3star =  print_uvp2tri_mcmc_acceptance_rate_3star(directory, "F606W", fit_folders = ['3star-fit'])
    return string_3star




def notebook_complete(message):
    """
    Call at the end of a MORIA notebook when the final pipeline step is done.

    Prints the Doors of Durin (ASCII art), plus an optional completion message.

    Parameters
    ----------
    message : str, optional
        Short note to print under the portrait (e.g. notebook name or step).
    """
    portrait = (Path(__file__).resolve().parent / "doors_of_durin.txt").read_text(
        encoding="utf-8"
    )
    print(portrait)
    if message:
        print(f"  >> {message}")
    print("  >> MORIA notebook complete. The doors are open.\n")




def notebook_complete_fit():
    """
    Call at the end of fitting_psfs.ipynb when PSF fitting is done.

    Prints gollum.
    """
    portrait = (
        Path(__file__).resolve().parent / "gollum.txt"
    ).read_text(encoding="utf-8")
    print(portrait)
    print("  >> Gollum Calls You His Precious")


def calibration_input_file_one(directory):
    
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "psfout_simst.fits")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "psfout_simstV.fits")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "img2extract_wfc3uv_psflist_Cal.uvp.gz")

    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "img2extract_wfc3uv_psflist_simst.uvp.gz")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "img2extract_wfc3uv_psflist_CalV.uvp.gz")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "img2extract_wfc3uv_psflist_simstV.uvp.gz")


    base_dir = Path(directory).resolve() 
    subdir = base_dir / "07.CALIBRATION"
    in_psf_star_mags_mcmc_I = 'IN.psf_star_mags_mcmc_I'
    output_file_psf_star_mags_mcmc_I = os.path.join(subdir, in_psf_star_mags_mcmc_I)
    markov_chain_steps = int(input("Enter the number of Markov chain steps"))
    maximum_size_mcmc = float(input("Enter the maximum size of MCMC coordinate steps in pixels"))
    fudge = float(input("Enter the error bar fudge factor. Let fuge be 1 by default"))
    maximum_distance_x = float(input("Enter maximum distance in x from the star center for pixels to be included in the fit"))
    maximum_distance_y = float(input("Enter maximum distance in y from the star center for pixels to be included in the fit"))
    chi2cut = float(input("Enter the χ2 threshold to define an outlier pixel."))
    sky_model = int(input("The sky model to use"))
    star_numbers = int(input("A list of star numbers to produce PIX_SHOW files for"))
    content = [
            "psfout_simst.fits",
            "Cal",
            "Cal_I",
            f"{markov_chain_steps}",
            f"{maximum_size_mcmc} {fudge}",
            f"{maximum_distance_x} {maximum_distance_y} {chi2cut}",
            f"{sky_model}",
            f"{star_numbers}",
            "0"
        ]

    with open(output_file_psf_star_mags_mcmc_I, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{in_psf_star_mags_mcmc_I}' successfully created.")


    
    return



    
def calibration_input_file_two(directory):
    
    base_dir = Path(directory).resolve() 
    subdir = base_dir / "07.CALIBRATION"
    in_psf_star_mags_mcmc_V = 'IN.psf_star_mags_mcmc_V'
    output_file_psf_star_mags_mcmc_V = os.path.join(subdir, in_psf_star_mags_mcmc_V)
    markov_chain_steps = int(input("Enter the number of Markov chain steps"))
    maximum_size_mcmc = float(input("Enter the maximum size of MCMC coordinate steps in pixels"))
    fudge = float(input("Enter the error bar fudge factor. Let fuge be 1 by default"))
    maximum_distance_x = float(input("Enter maximum distance in x from the star center for pixels to be included in the fit"))
    maximum_distance_y = float(input("Enter maximum distance in y from the star center for pixels to be included in the fit"))
    chi2cut = float(input("Enter the χ2 threshold to define an outlier pixel."))
    sky_model = int(input("The sky model to use"))
    star_numbers = int(input("A list of star numbers to produce PIX_SHOW files for"))
    content = [
            "psfout_simstV.fits",
            "CalV",
            "Cal_V",
            f"{markov_chain_steps}",
            f"{maximum_size_mcmc} {fudge}",
            f"{maximum_distance_x} {maximum_distance_y} {chi2cut}",
            f"{sky_model}",
            f"{star_numbers}",
            "0"
        ]

    with open(output_file_psf_star_mags_mcmc_V, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{in_psf_star_mags_mcmc_V}' successfully created.")


    
    return

    
def calibration_new_matchup(directory):
    """
    Prepare 07.CALIBRATION matchup files for OGLE calibration.

    Keeps an existing bar-space MATCHUP catalog when present (6252/4738 anchor),
    otherwise builds from 07.CALIBRATION/dex_no_gaia_STEP08_A.xyvieeee. Rows are
    written in the 133-char xym2bar layout with exposure counts 1 (F814W) or 2
    (F606W), then MATCHUP.F814W_cal.XYM is annotated from bar-space NOTFAR.
    """
    root = Path(directory).resolve()
    cal_dir = root / "07.CALIBRATION"
    cal_dir.mkdir(parents=True, exist_ok=True)

    for name in (
        "NEARBY_SIM_STARS.XYIVB_targ",
        "NEARBY_REF_STARS.XYIVB_targ",
    ):
        src = root / "04.EXTRACT_PSF" / "F814W" / name
        if not src.is_file():
            src = root / "02.CMD" / name
        if src.is_file():
            shutil.copy2(src, cal_dir / name)

    dex_path = cal_dir / "dex_no_gaia_STEP08_A.xyvieeee"
    if not dex_path.is_file():
        dex_path = resolve_dex_xyvieeee_path(directory)
        shutil.copy2(dex_path, cal_dir / "dex_no_gaia_STEP08_A.xyvieeee")
        dex_path = cal_dir / "dex_no_gaia_STEP08_A.xyvieeee"

    f814 = cal_dir / "MATCHUP.F814W.XYM.02"
    f606 = cal_dir / "MATCHUP.F606W.XYM"
    if f814.is_file() and _matchup_has_bar_anchor(f814):
        _reformat_matchup_exposure_inplace(f814, 1)
        _reformat_matchup_exposure_inplace(f606, 2)
        print(f"Reformatted exposure counts on existing {f814} and {f606}")
    else:
        write_matchup_from_dex_xyvieeee(dex_path, f814, band="F814W")
        write_matchup_from_dex_xyvieeee(dex_path, f606, band="F606W")
        print(f"Wrote {f814} and {f606} from {dex_path}")

    loc_trans = root / "03.LOC_TRANS" / "outputq.fits"
    if loc_trans.is_file():
        shutil.copy2(loc_trans, cal_dir / "outputq.fits")
    else:
        for alt in (
            root / "02.CMD" / "outputq_F814W.fits",
            root / "01.XYM" / "F814W" / "outputq_F814W.fits",
        ):
            if alt.is_file():
                shutil.copy2(alt, cal_dir / "outputq.fits")
                break

    notfar = _resolve_calibration_notfar(
        directory, cal_dir, f814_matchup=f814
    )
    write_cal_matchup_files(
        f814,
        notfar,
        cal_dir / "MATCHUP.F814W_cal.XYM",
        cal_dir / "MATCHUP.F814W_cal_only.XYM",
    )
    print("Wrote MATCHUP.F814W_cal.XYM and MATCHUP.F814W_cal_only.XYM")

    _ensure_psf_mcmc_fit_sky(directory)

    def psf_star_mags_mcmc(directory, script='run_psf_star_Imags_mcmc.src'):
        
        log_file = Path(directory).resolve() / "07.CALIBRATION" / "log_files" / f"run_psf_star_Imags_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "07.CALIBRATION"
                subdir = base_dir 
                script = 'run_psf_star_Imags_mcmc.src'
                script_path = base_dir / script
                print(base_dir)
                print(script_path)
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        print(f"Finished run_psf_star_Imags_mcmc")

        log_file = Path(directory).resolve() / "07.CALIBRATION" / "log_files" / f"run_psf_star_Vmags_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "07.CALIBRATION"
                subdir = base_dir 
                script = 'run_psf_star_Vmags_mcmc.src'
                script_path = base_dir / script
                print(base_dir)
                print(script_path)
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        print(f"Finished run_psf_star_Vmags_mcmc")
        return

    psf_star_mags_mcmc(directory)


def calibration_hst_ogle_match(directory):
    """
    The goal here is to calibrate the HST photometry to the OGLE-III database.
    """

    def VI_HST_ogle_man_match4(directory, script='run_VI_HST_ogle_man_match4.src'):
        
        log_file = Path(directory).resolve() / "07.CALIBRATION" / "log_files" / f"run_VI_HST_ogle_man_match4.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "07.CALIBRATION"
                subdir = base_dir 
                script = 'run_VI_HST_ogle_man_match4.src'
                script_path = base_dir / script
                print(base_dir)
                print(script_path)
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        print(f"Finished run_VI_HST_ogle_man_match4")
        return

    VI_HST_ogle_man_match4(directory)
    fix_vi_hst_ogle_cal_matches4_spacing(directory)
    expand_vi_hst_cal_matches_for_fit(directory, min_stars=15, max_stars=20)


def fit_calibration(directory):
    """
    The goal here is to calibrate the HST photometry to the OGLE-III database.
    """
    def fit_HST_IV_ogle_col_1(directory, script='run_fit_HST_IV_ogle_col_1.src'):
        
        log_file = Path(directory).resolve() / "07.CALIBRATION" / "log_files" / f"run_fit_VI_HST_ogle_man_match4.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "07.CALIBRATION"
                subdir = base_dir 
                script = 'run_fit_HST_IV_ogle_col_1.src'
                script_path = base_dir / script
                print(base_dir)
                print(script_path)
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
                script = 'run_fit_HST_IV_ogle_col_2.src'
                script_path = base_dir / script
                print(base_dir)
                print(script_path)
                subprocess.run(
                    ["csh", str(script_path)],
                    cwd=subdir,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
                
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        print(f"Finished fitting calibration")
        return
    fix_vi_hst_ogle_cal_matches4_spacing(directory)
    fit_HST_IV_ogle_col_1(directory)

def get_chip_number(ogle_ra_deg, ogle_dec_deg):
    """
    Get the chip number from the OGLE-III field finder for you target
    """
    if ogle_ra_deg is not None and ogle_dec_deg is not None:
        candidates = ogle_field_chip_candidates_from_coords(
            ogle_ra_deg,
            ogle_dec_deg,
            phase=3,
            epoch=2000.0,
        )
        # Most events land in a single chip; if multiple are returned, we take the first.
        ogle_field_number = candidates[0]["field_number"]
        ogle_chip_number = candidates[0]["chip_number"]
        print("OGLE field candidates:", candidates)
    else:
        ogle_field_number = 0
        ogle_chip_number = 0
    return ogle_field_number, ogle_chip_number



def ogle_map_and_reference_filenames(ogle_field_number, ogle_chip_number, ogle_band="I",
 map_prefix_case="upper", ref_prefix_case="lower"):
    """
    Construct the OGLE-III blg filenames from a field + chip identifier.

    Example (as described in `demo/in_dev.ipynb`):
    - catalog map: blg226.7.map
    - reference image (I band): blg226.I.7.fts
    """
    field = int(ogle_field_number)
    chip = int(ogle_chip_number)
    band = str(ogle_band).strip()

    map_prefix = "BLG" if map_prefix_case.lower() == "upper" else "blg"
    ref_prefix = "BLG" if ref_prefix_case.lower() == "upper" else "blg"

    # Filenames on the remote server (lowercase `blg`).
    map_filename_remote = f"blg{field}.{chip}.map"
    ref_filename_remote = f"blg{field}.{band}.{chip}.fts"

    # Filenames saved locally (sometimes Fortran/xgf pipelines are case-sensitive).
    map_filename_local = f"{map_prefix}{field}.{chip}.map"
    ref_filename_local = f"{ref_prefix}{field}.{band}.{chip}.fts"

    return {
        "map_filename_remote": map_filename_remote,
        "ref_filename_remote": ref_filename_remote,
        "map_filename_local": map_filename_local,
        "ref_filename_local": ref_filename_local,
    }


def download_ogle_map_and_reference(directory, ogle_field_number, ogle_chip_number, ogle_band="I", 
                                    destination_subdir="07.CALIBRATION", overwrite=False, map_prefix_case="upper", ref_prefix_case="lower", 
                                    maps_base_url="http://www.astrouw.edu.pl/ogle/ogle3/maps/blg/maps/", ref_images_base_url="http://www.astrouw.edu.pl/ogle/ogle3/maps/blg/ref_images/", timeout_s=120.0):
                                        
    """
    Download the OGLE-III photometry map and reference image needed for calibration.

    The OGLE server typically hosts these as compressed files:
    - `blg{field}.{chip}.map.bz2`  -> decompressed to `blg{field}.{chip}.map`
    - `blg{field}.{band}.{chip}.fts.bz2` -> decompressed to `blg{field}.{band}.{chip}.fts`
    """
    directory = Path(directory).resolve()
    dest_dir = directory / destination_subdir
    dest_dir.mkdir(parents=True, exist_ok=True)

    fn = ogle_map_and_reference_filenames(
        ogle_field_number=ogle_field_number,
        ogle_chip_number=ogle_chip_number,
        ogle_band=ogle_band,
        map_prefix_case=map_prefix_case,
        ref_prefix_case=ref_prefix_case,
    )

    map_path = dest_dir / fn["map_filename_local"]
    ref_path = dest_dir / fn["ref_filename_local"]

    map_candidates = [
        maps_base_url.rstrip("/") + "/" + fn["map_filename_remote"],
        maps_base_url.rstrip("/") + "/" + fn["map_filename_remote"] + ".bz2",
    ]
    ref_candidates = [
        ref_images_base_url.rstrip("/") + "/" + fn["ref_filename_remote"],
        ref_images_base_url.rstrip("/") + "/" + fn["ref_filename_remote"] + ".bz2",
    ]

    def download_to_path(url: str, out_path: Path):
        with urllib.request.urlopen(url, timeout=timeout_s) as resp:
            with open(out_path, "wb") as f:
                shutil.copyfileobj(resp, f)

    def try_download_candidates(candidates: list[str], out_path: Path) -> str:
        if out_path.exists() and not overwrite:
            return "already-present"

        last_err = None
        for url in candidates:
            try:
                if url.endswith(".bz2"):
                    tmp_path = out_path.with_name(out_path.name + ".bz2")
                    download_to_path(url, tmp_path)
                    with bz2.open(tmp_path, "rb") as fin, open(out_path, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
                    # Keep or remove the compressed copy depending on your preference;
                    # default to cleanup so the destination stays tidy.
                    try:
                        tmp_path.unlink()
                    except OSError:
                        pass
                    return url
                else:
                    download_to_path(url, out_path)
                    return url
            except urllib.error.HTTPError as e:
                last_err = e
            except Exception as e:
                last_err = e

        tried = "\n".join(f"- {c}" for c in candidates)
        raise RuntimeError(
            f"Failed to download OGLE files for {out_path.name}. Tried:\n{tried}\nLast error: {last_err}"
        )

    used_map_url = try_download_candidates(map_candidates, map_path)
    used_ref_url = try_download_candidates(ref_candidates, ref_path)

    subset_map_path = dest_dir / f"New_{fn['map_filename_local']}"
    columns = ["ogle_x", "ogle_y", "V", "I"]
    df = pd.read_csv( map_path, header=None, sep=r"\s+", usecols=[3, 4, 5, 7], names=columns)
    with open(subset_map_path, "w", encoding="utf-8") as f:
        f.write("ogle_x ogle_y V I\n")
        for row in df.itertuples(index=False):
            f.write(
                f"{row.ogle_x:7.2f} {row.ogle_y:7.2f} {row.V:6.3f} {row.I:6.3f}\n"
            )

    return {
        "map_url": used_map_url,
        "ref_url": used_ref_url,
        "map_path": str(map_path),
        "subset_map_path": str(subset_map_path),
        "ref_path": str(ref_path),
    }


def ogle_field_chip_candidates_from_coords(ra_deg, dec_deg, phase=3, epoch=2000.0,
 assume_ra_is_hours_if_lt_24=True, base_url="https://ogle.astrouw.edu.pl/cgi-ogle/uncgi.cgi/radec2field", 
 timeout_s=30.0):
    """
    Query the OGLE Field Finder to get candidate OGLE-III fields + chip numbers.

    Parameters
    ----------
    ra_deg, dec_deg
        Sky position in degrees 
    phase
        OGLE phase to query
    epoch
        Epoch passed through to OGLE Field Finder

    Returns
    -------
    List of dicts like:
      {"field_name": "BLG226.7", "field_number": 226, "chip_number": 7, "x": ..., "y": ...}
    """
    ra_val = float(ra_deg)
    dec_val = float(dec_deg)
    if assume_ra_is_hours_if_lt_24 and 0.0 <= ra_val <= 24.0:
     #"Interpreting input RA as hours because ra_deg={ra_val} is in [0,24]. (Will send RA={ra_hours} hours to OGLE.)"
        ra_hours = ra_val
    else:
        ra_hours = None

    payload = {
        "phase": str(phase),
        # OGLE Field Finder expects RA in hours (hh.hhhh or hh:mm:ss)
        "ra": f"{ra_hours:.6f}",
        "dec": f"{dec_val:.6f}",
        "epoch": f"{float(epoch):.1f}",
    }

    data = urllib.parse.urlencode(payload).encode("utf-8")

    req = urllib.request.Request(base_url, data=data, method="POST")

    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        html = resp.read().decode("iso-8859-1", errors="replace")

    m = re.search(r"<PRE>(.*?)</PRE>", html, flags=re.S | re.I)

    if not m:
        raise RuntimeError("Unexpected OGLE Field Finder response: missing <PRE> section.")

    pre = m.group(1)
    rows = [ln.strip() for ln in pre.splitlines() if ln.strip()]

    # First row is usually: "field phase x y"
    candidates: list[dict] = []
    for ln in rows:
        parts = ln.split()
        if not parts:
            continue
        if parts[0].lower() == "field":
            continue

        # Expected: field phase x y
        field_token = parts[0]
        # field_token like "BLG126.6"
        mf = re.match(r"([A-Za-z]+)(\d+)\.(\d+)", field_token)
        if not mf:
            continue

        field_prefix = mf.group(1)
        field_number = int(mf.group(2))
        chip_number = int(mf.group(3))

        x = float(parts[2]) if len(parts) > 2 else None
        y = float(parts[3]) if len(parts) > 3 else None

        candidates.append(
            {
                "field_name": field_token,
                "field_prefix": field_prefix,
                "field_number": field_number,
                "chip_number": chip_number,
                "x": x,
                "y": y,
            }
        )

    if not candidates:
        raise ValueError(f"No OGLE field candidates found for ra_deg={ra_deg}, dec_deg={dec_deg}.")

    return candidates

