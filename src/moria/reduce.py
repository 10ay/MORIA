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

def data_prep_early(destination):
    fortran_src = get_fortran_dir()
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "00.DATA" / "F814W", extensions=[".xOg"])
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "00.DATA" / "F606W", extensions=[".xOg"])
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM" / "F814W", extensions=[".xOg"])
    copy_files(source=fortran_src, destination=Path(destination).resolve() / "01.XYM" / "F606W", extensions=[".xOg"])
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS" / "F814W", filename="xym2mat.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS" / "F814W", filename="img2extract_wfc3uv_psflist.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS" / "F606W", filename="xym2mat.xOg")
    copy_entire_files(source=fortran_src, destination=Path(destination).resolve() / "03.LOC_TRANS" / "F606W", filename="img2extract_wfc3uv_psflist.xOg")
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
                    f.write(f"{0:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")
                f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")


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
                f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")
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

        
        output_file_dir = base_dir / '01.XYM' / f
        output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
        output_file_xym2mat = os.path.join(output_file_dir, in_xym2mat)
        output_file_xym2bar = os.path.join(output_file_dir, in_xym2bar)
        
        with open(output_file_xym2mat, "w") as f:
            f.write("00 MATCHUP.F814W.XYM.02 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f6 \"m-14.75,-5.5\" \n")


        with open(output_file_xym2bar, "w") as f:
            f.write("00 MATCHUP.F814W.XYM.02 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f6\n")


        with open(output_file_img2sam, "w") as t:
            for i, filename in enumerate(files_two, start=1):
                t.write(f"{i:02d} \"{filename}\" 6 0\n")
        
                
        return
    
    data_prep_F814W(directory)
    data_prep_F606W(directory)
    
def matchup_files(directory):
    """
    Run the scripts to create MATCHUP Files on _WJ2 files in F814W and F606W subdirectories.
    """

    def run_img2xym(directory, script='run_img2xym_wfc3uv.src'):
        """
        Run img2xym_wfc3uv on exposures in the F814W subdirectory.
        Produces .XYM files for each exposure using the PSFEFF_WFC3UV_F814W_C0.fits library PSF.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / "run_img2xym_wfc3uv.log"
        with open(log_file, "w") as logf:
            sys.stdout = sys.stderr = logf
            try:
                base_dir = Path(directory).resolve() / "01.XYM"
                filters = ['F814W', 'F606W']

                for f in filters:
                    subdir = base_dir / f
                    script_path = subdir / script if f == "F814W" else subdir / script
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

    def run_xym2mat(directory, filters = ['F814W'], script='run_xym2mat_1.src'):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / f"{script.replace('.src','')}.log"

        with open(log_file, "w") as logf:
            sys.stdout = sys.stderr = logf
            try:
                base_dir = Path(directory).resolve() / "01.XYM"

                for f in filters:
                    subdir = base_dir / f
                    script_to_use = script if f == 'F814W' else 'run_xym2mat_VI.src'
                    script_path = subdir / script_to_use
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
        

    def run_xym2bar(directory, filters = ['F814W'], script='run_xym2bar_1.src'):
        """
        Get a final list of photometry that allowed small zeropoint shifts for each exposure
        """
        log_file = Path(directory).resolve() / "01.XYM" / "log_files" / f"{script.replace('.src','')}.log"

        with open(log_file, "w") as logf:
            sys.stdout = sys.stderr = logf
            try:
                base_dir = Path(directory).resolve() / "01.XYM"
                for f in filters:
                    subdir = base_dir / f
                    script_to_use = script if f == 'F814W' else 'run_xym2bar.src'
                    script_path = subdir / script_to_use
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

    copy_files(source=Path(directory).resolve() / "00.DATA" / "F814W", destination=Path(directory).resolve() / "01.XYM" / "F814W", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "00.DATA" / "F606W", destination=Path(directory).resolve() / "01.XYM" / "F606W", extensions=[".fits"])
    run_img2xym(directory)
    data_prep(directory)
    run_xym2mat(directory)
    run_xym2bar(directory)
    run_xym2mat(directory, script='run_xym2mat_2.src')
    run_xym2bar(directory, script='run_xym2bar_2.src')
    copy_files(source=Path(directory).resolve() / "01.XYM" / "F814W", destination=Path(directory).resolve() / "01.XYM" / "F606W", extensions=[".02"])
    run_xym2mat(directory, filters = ['F606W'])
    run_xym2bar(directory, filters = ['F606W'])


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

def data_prep_loc_trans(directory, filters = 'F814W'):
    """
    Preparae IN.* files in respective files using _WJ2.fits files in 00.DATA
    Parameters
    ----------
    directory - The root directory to operate on. There is a specific directory
    structure that is expected within <dir>:
        00.DATA/
        01.XYM/

    Returns
    -------
    several IN.* files.
    """

    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "03.LOC_TRANS" / filters)
    copy_files(source=Path(directory).resolve() / "01.XYM" / filters, extensions=[".xym"], destination=Path(directory).resolve() / "03.LOC_TRANS" / filters)
    copy_files(source=Path(directory).resolve() / "01.XYM" / filters, extensions=[".fits"], destination=Path(directory).resolve() / "03.LOC_TRANS" / filters)

    base_dir = Path(directory).resolve()

    subdir = base_dir / filters
    in_img2sam_wfc3uv = 'IN.img2sam_wfc3uv'
    in_xym2mat = 'IN.xym2mat'

    base_dir_one = base_dir / '01.XYM'/ filters
    f=filters
    files = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.xym')])    
    files_two = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.fits')])    
    output_file_dir = base_dir / '03.LOC_TRANS' / f

    output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
    output_file_xym2mat = os.path.join(output_file_dir, in_xym2mat)
    if filters == 'F814W':
        with open(output_file_xym2mat, "w") as f:
            f.write("00 NEARBY_REF_STARS.XYIVB_targ c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")
    
        with open(output_file_img2sam, "w") as t:
            for i, filename in enumerate(files_two, start=1):
                t.write(f"{i:02d} \"{filename}\" 8 0\n")
    else:
        with open(output_file_xym2mat, "w") as f:
            f.write("00 NEARBY_REF_STARS.XYIVB_targ c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 \"m-14.75,-5.5\" \n")
    
        with open(output_file_img2sam, "w") as t:
            for i, filename in enumerate(files_two, start=1):
                t.write(f"{i:02d} \"{filename}\" 6 0\n")
            
    return

def loc_trans(directory):
    """
    Run the local transformation scripts in F814W and F606W subdirectories to extract 
    the pixels from each exposure and accurately transform their locations into the 
    reference frame so that we can use them to solve for a PSF and then use this PSF to model the target star.

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
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script_path = subdir / script if f == "F814W" else base_dir / "F606W" / script
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

    def run_img2extract_wfc3uv_psflist(directory):
        """Run extraction for PSF list generation (simulation)."""
        log_file = Path(directory).resolve() /  "03.LOC_TRANS" / "log_files" / f"loc_trans_run_img2extract_wfc3uv_psflist.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_img2extract_wfc3uv_psflist_simst.src' if f == "F814W" else 'run_img2extract_wfc3uv_psflist_simstV.src'
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

    def run_img2extract_wfc3uv_psflist_Cal(directory):
        """Run extraction for PSF list generation (calibration)."""

        log_file = Path(directory).resolve() /   "03.LOC_TRANS" / "log_files" / f"loc_trans_run_img2extract_wfc3uv_psflist_Cal.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_img2extract_wfc3uv_psflist_Cal.src' if f == "F814W" else 'run_img2extract_wfc3uv_psflist_CalV.src'
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
                
    data_prep_loc_trans(directory)
    data_prep_loc_trans(directory, filters = 'F606W')
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
    pre, dat = _cmd_partition_matchup_raw_lines(path)
    _cmd_write_matchup_raw_lines(path, pre, dat)


def cmd_diagram(directory):

    fortran_src = get_fortran_dir()
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "01.XYM" / "F814W", filename="MATCHUP.F814W.XYM.02")
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "01.XYM" / "F606W", filename="MATCHUP.F606W.XYM")

    subdir = Path(directory).resolve() / "02.CMD"
    path_match_I = subdir / "MATCHUP.F814W.XYM.02"
    path_match_V = subdir / "MATCHUP.F606W.XYM"

    cmd_rewrite_matchup_drop_xym2bar_echo(path_match_I)
    cmd_rewrite_matchup_drop_xym2bar_echo(path_match_V)

    xv, yv, mv = np.loadtxt(path_match_V, unpack=True, usecols=(0, 1, 2))
    xi, yi, mi = np.loadtxt(path_match_I, unpack=True, usecols=(0, 1, 2))
    
    #Establish target parameters
    response = str(input("Do you have a target? Enter 'Yes' if you do."))
    if response == 'yes' or response == 'Yes':
        xtarg = float(input("Enter x-coord of your target"))
        ytarg = float(input("Enter y-coord of your target"))
        Vtarg = float(input("Enter V magnitude of your target"))
        Itarg = float(input("Enter I magnitude of your target"))
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
        path_V = path_match_V
        preamble_I, lines_I = cmd_partition_matchup_raw_lines(path_I)
        preamble_V, lines_V = cmd_partition_matchup_raw_lines(path_V)
        n = len(xi)
        dist2 = (xi - xtarg) ** 2 + (yi - ytarg) ** 2
        idx = int(np.argmin(dist2))
        order_list = [idx] + [i for i in range(n) if i != idx]
        if idx != 0:
            reordered_I = [lines_I[i] for i in order_list]
            reordered_V = [lines_V[i] for i in order_list]
            cmd_write_matchup_raw_lines(path_I, preamble_I, reordered_I)
            cmd_write_matchup_raw_lines(path_V, preamble_V, reordered_V)
        xv, yv, mv = np.loadtxt(path_V, unpack=True, usecols=(0, 1, 2))
        xi, yi, mi = np.loadtxt(path_I, unpack=True, usecols=(0, 1, 2))

    #Function to find the CMD of the target and get our list of Sim+Ref stars.

    def show_cmd_targ(directory, response, xi, yi, xv, yv, mi, mv):
    
        #Plotting parameters 
        if response == 'yes' or response == 'Yes':
    
            rad_max = float(input("Enter maximum plotting radius"))
            box_max = float(input("Enter maximum box radius"))
            mag_range = float(input("Enter magnitude range for plotting"))
            col_range = float(input("Enter color range for plotting"))
            ref_st_Imx = float(input("Enter reference star input I max"))
            ref_st_Imn = float(input("Enter reference star input I min"))
            ref_st_Vmx = float(input("Enter reference star input V max"))
            ref_st_Vmn = float(input("Enter reference star input V min"))
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

            rad_max = float(input("Enter maximum plotting radius"))
            box_max = float(input("Enter maximum box radius"))
            Vcalc = float(input("Enter V magnitude for calibration"))
            Icalc = float(input("Enter I magnitude for calibration"))
            mag_range = float(input("Enter magnitude range for plotting"))
            col_range = float(input("Enter color range for plotting"))
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
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F606W", destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W")
    f814_images   = int(input("Enter the number of images for the F814W filter: "))
    f606_images   = int(input("Enter the number of images for the F606W filter: "))
    
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
                
    prepare_data(f814_images, directory)
    prepare_data(f606_images, directory, f = 'F606W')
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
        
    prepare_data(good_psf, directory)
    prepare_data(good_psf, directory, f='F606W')
    run_uvp2psf_simst(directory)


def hst_fit_dataprep_twostar(directory, f = 'F814W'):

    x1 = float(input("Initial x position for object 1"))
    y1 = float(input("Initial y position for object 1"))

    x2 = float(input("Initial x position for object 2"))
    y2 = float(input("Initial y position for object 2"))
 
    x3 = float(input("Initial x position for object 3"))
    y3 = float(input("Initial y position for object 3"))

    f1 = float(input("Initial flux for object 1"))
    f2 = float(input("Initial flux for object 2"))

    mcmc_dr1 = float(input("Maximum MCMC jump size for object 1"))
    mcmc_dr2 = float(input("Maximum MCMC jump size for object 2"))
    mcmc_dr3 = float(input("Maximum MCMC jump size for object 3"))
    mcmc_df1 = float(input("Maximum MCMC jump size for flux of object 1"))
    mcmc_df2 = float(input("Maximum MCMC jump size for flux of object 2"))

    nmcmc = int(input("MCMC step sizes (Recommended: >50000)"))

    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))

    dufitmn = float(input("Minimum du cut"))
    dufitmx = float(input("Maximum du cut"))
    dvfitmn =float(input("Minimum dv cut"))
    dvfitmx = float(input("Maximum dv cut"))
    chi2cut = float(input("Chi-squared cut"))

    
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

    x1 = float(input("Initial x position for object 1"))
    y1 = float(input("Initial y position for object 1"))

    x2 = float(input("Initial x position for object 2"))
    y2 = float(input("Initial y position for object 2"))
 
    x3 = float(input("Initial x position for object 3"))
    y3 = float(input("Initial y position for object 3"))

    f1 = float(input("Initial flux for object 1"))
    f2 = float(input("Initial flux for object 2"))

    mcmc_dr1 = float(input("Maximum MCMC jump size for object 1"))
    mcmc_dr2 = float(input("Maximum MCMC jump size for object 2"))
    mcmc_dr3 = float(input("Maximum MCMC jump size for object 3"))
    mcmc_df1 = float(input("Maximum MCMC jump size for flux of object 1"))
    mcmc_df2 = float(input("Maximum MCMC jump size for flux of object 2"))

    nmcmc = int(input("MCMC step sizes"))

    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))

    dufitmn = float(input("Minimum du cut"))
    dufitmx = float(input("Maximum du cut"))
    dvfitmn =float(input("Minimum dv cut"))
    dvfitmx = float(input("Maximum dv cut"))
    chi2cut = float(input("Chi-squared cut"))

    
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

    x1 = float(input("Initial x position for object 1"))
    y1 = float(input("Initial y position for object 1"))

    x2 = float(0)
    y2 = float(0)
 
    x3 = float(0)
    y3 = float(0)

    f1 = float(1)
    f2 = float(0)

    mcmc_dr1 = float(input("Maximum MCMC jump size for object 1"))
    mcmc_dr2 = float(0)
    mcmc_dr3 = float(0)
    mcmc_df1 = float(input("Maximum MCMC jump size for flux of object 1"))
    mcmc_df2 = float(0)

    nmcmc = int(input("MCMC step sizes"))

    fudge = float(input("Input fudge factor. Input 1.0 if you don't know what this is"))

    dufitmn = float(input("Minimum du cut"))
    dufitmx = float(input("Maximum du cut"))
    dvfitmn =float(input("Minimum dv cut"))
    dvfitmx = float(input("Maximum dv cut"))
    chi2cut = float(input("Chi-squared cut"))

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

    def run_uvp2psf_simst_1(directory, script='run_uvp2tri_NOscon_fs_asym_mcmc.src'):
        log_file = Path(directory).resolve() / "06.FIT" /  "F814W" / "3star-fit" / "log_files" / f"uvp2tri_scon_fs_asym_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" /  "F814W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script = script
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
    _uvp2tri_fsky_outputs = (
        "uvp2tri_scon_fsky_I_KeckNOcon.01.pix_all",
        "uvp2tri_scon_fsky_I_KeckNOcon.03.pix_use",
        "uvp2tri_scon_fsky_I_KeckNOcon.04.probe_fit",
        "uvp2tri_scon_fsky_I_KeckNOcon.05.final_fit",
        "uvp2tri_scon_fsky_I_KeckNOcon.06.pix_show.fits",
        "uvp2tri_scon_fsky_I_KeckNOcon.07.mcmc",
        "uvp2tri_scon_fsky_I_KeckNOcon.08.rm_pix",
    )
    for name in _uvp2tri_fsky_outputs:
        (_fit_3star / name).unlink(missing_ok=True)

    fortran_src = get_fortran_dir()
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename="mcmc_expand_average.xOg")

    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", extensions=[".fits"])

    
    run_uvp2psf_simst_1(directory)
    #run_uvp2psf_simst_2(directory)
    strip_star_lines_from_uvp2tri_mcmc_814W(directory)
    run_mcmc_expand_average_814W(directory)
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

    def run_uvp2psf_simst_1(directory, script='run_uvp2tri_NOscon_fs_asym_mcmc.src'):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            log_file = Path(directory).resolve() / "06.FIT" / "F814W" / f /"log_files" / f"uvp2tri_scon_fs_asym_mcmc.log"
            with open(log_file, "w") as logf:
                try:
                    base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                    subdir = base_dir / f
                    script = script
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
        _uvp2tri_fsky_outputs = (
            "uvp2tri_scon_fsky_I_KeckNOcon.01.pix_all",
            "uvp2tri_scon_fsky_I_KeckNOcon.03.pix_use",
            "uvp2tri_scon_fsky_I_KeckNOcon.04.probe_fit",
            "uvp2tri_scon_fsky_I_KeckNOcon.05.final_fit",
            "uvp2tri_scon_fsky_I_KeckNOcon.06.pix_show.fits",
            "uvp2tri_scon_fsky_I_KeckNOcon.07.mcmc",
            "uvp2tri_scon_fsky_I_KeckNOcon.08.rm_pix",
        )
        for name in _uvp2tri_fsky_outputs:
            (_fit_2star / name).unlink(missing_ok=True)

    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit")

    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit")

    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "1star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "2star-fit", extensions=[".fits"])
  

    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "1star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "2star-fit", extensions=[".fits"])

    
    run_uvp2psf_simst_1(directory)
    strip_star_lines_from_uvp2tri_mcmc_814W(directory)
    run_mcmc_expand_average_814W(directory)




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

    def run_uvp2psf_simst_2(directory, script='run_uvp2tri_NOscon_fs_asym_mcmc.src'):
        folders = ['1star-fit', '2star-fit']
        for f in folders:
            log_file = Path(directory).resolve() / "06.FIT" / "F606W" / f / "log_files" /f"uvp2tri_scon_fs_asym_mcmc.log"
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
        _uvp2tri_fsky_outputs = (
            "uvp2tri_scon_fsky_V_KeckNOcon.01.pix_all",
            "uvp2tri_scon_fsky_V_KeckNOcon.03.pix_use",
            "uvp2tri_scon_fsky_V_KeckNOcon.04.probe_fit",
            "uvp2tri_scon_fsky_V_KeckNOcon.05.final_fit",
            "uvp2tri_scon_fsky_V_KeckNOcon.06.pix_show.fits",
            "uvp2tri_scon_fsky_V_KeckNOcon.07.mcmc",
            "uvp2tri_scon_fsky_V_KeckNOcon.08.rm_pix",
        )
        for name in _uvp2tri_fsky_outputs:
            (_fit_2star / name).unlink(missing_ok=True)
        
    run_uvp2psf_simst_2(directory)
    strip_star_lines_from_uvp2tri_mcmc_606W(directory)
    run_mcmc_expand_average_606W(directory)


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

    def run_uvp2psf_simst_2(directory, script='run_uvp2tri_NOscon_fs_asym_mcmc.src'):
        log_file = Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit" / "log_files" / f"uvp2tri_scon_fs_asym_mcmc.log"
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
    
        
    _fit_2star = Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit"
    _uvp2tri_fsky_outputs = (
        "uvp2tri_scon_fsky_V_KeckNOcon.01.pix_all",
        "uvp2tri_scon_fsky_V_KeckNOcon.03.pix_use",
        "uvp2tri_scon_fsky_V_KeckNOcon.04.probe_fit",
        "uvp2tri_scon_fsky_V_KeckNOcon.05.final_fit",
        "uvp2tri_scon_fsky_V_KeckNOcon.06.pix_show.fits",
        "uvp2tri_scon_fsky_V_KeckNOcon.07.mcmc",
        "uvp2tri_scon_fsky_V_KeckNOcon.08.rm_pix",
    )
    for name in _uvp2tri_fsky_outputs:
        (_fit_2star / name).unlink(missing_ok=True)
    
    fortran_src = get_fortran_dir()
    copy_entire_files(source=fortran_src, destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", filename="uvp2tri_scon_fs_asym_mcmc.xOg")
    run_uvp2psf_simst_2(directory)
    strip_star_lines_from_uvp2tri_mcmc_606W(directory)
    run_mcmc_expand_average_606W(directory)


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
    The goal here is to calibrate the HST photometry to the OGLE-III database.
    """

    #def prepare_data(good_psf, directory, f= 'F814W'):
    #    base_dir = Path(directory).resolve()
     #   subdir = base_dir / f
      #  in_good_psf_list = 'IN.good_psf_list.2'
      #  output_file_dir = base_dir / '04.EXTRACT_PSF' / f
       # output_file_img = os.path.join(output_file_dir, in_good_psf_list)


#        with open(output_file_img, "w") as f:
 #           for i in range(1, len(good_psf) + 1):
  #              value = good_psf[i-1]
   #             f.write(f"{i:2d}   {value}\n")

    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "NOTFAR_CAL_STARS.XYIVB_targ")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "NEARBY_SIM_STARS.XYIVB_targ")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "NEARBY_REF_STARS.XYIVB_targ")

    copy_entire_files(source=Path(directory).resolve() / "02.CMD", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "MATCHUP.F606W.XYM")
    copy_entire_files(source=Path(directory).resolve() / "02.CMD", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "MATCHUP.F814W.XYM.02")

    copy_entire_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "outputq.fits")


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
    

    def cal_star_num(directory, script='run_cal_star_num_2_MATCHUP.src'):
        
        log_file = Path(directory).resolve() / "07.CALIBRATION" / "log_files" / f"run_cal_star_num_2_MATCHUP.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "07.CALIBRATION"
                subdir = base_dir 
                script = 'run_cal_star_num_2_MATCHUP.src'
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
        print(f"Finished cal_star_num_2_MATCHUP")
        return

    psf_star_mags_mcmc(directory)
    cal_star_num(directory)
    #VI_HST_ogle_man_match4(directory)
    #fit_HST_IV_ogle_col_1(directory)
    
    

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



def ogle_map_and_reference_filenames(
    ogle_field_number: int,
    ogle_chip_number: int,
    ogle_band: str = "I",
    *,
    map_prefix_case: str = "upper",
    ref_prefix_case: str = "lower",
):
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


def download_ogle_map_and_reference(
    directory: str | Path,
    ogle_field_number: int,
    ogle_chip_number: int,
    ogle_band: str = "I",
    destination_subdir: str = "07.CALIBRATION",
    overwrite: bool = False,
    map_prefix_case: str = "upper",
    ref_prefix_case: str = "lower",
    maps_base_url: str = "http://www.astrouw.edu.pl/ogle/ogle3/maps/blg/maps/",
    ref_images_base_url: str = "http://www.astrouw.edu.pl/ogle/ogle3/maps/blg/ref_images/",
    timeout_s: float = 120.0,
):
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

    return {
        "map_url": used_map_url,
        "ref_url": used_ref_url,
        "map_path": str(map_path),
        "ref_path": str(ref_path),
    }


def ogle_field_chip_candidates_from_coords(
    ra_deg: float,
    dec_deg: float,
    phase: str | int = 3,
    epoch: str | float = 2000.0,
    assume_ra_is_hours_if_lt_24: bool = True,
    base_url: str = "https://ogle.astrouw.edu.pl/cgi-ogle/uncgi.cgi/radec2field",
    timeout_s: float = 30.0,
) -> list[dict]:
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

