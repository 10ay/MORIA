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
parser = argparse.ArgumentParser()



"""
The HST fly star reduction pipeline depends on a certain directory structure.

<root_dir>/
    00.DATA/
    01.XYM/
"""

fcode_dir = '/Users/tmbhadra/Documents/Work/NASA/moira'

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

def data_prep_early(source, destination):
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "00.DATA" / "F814W", extensions=[".xOg"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "00.DATA" / "F606W", extensions=[".xOg"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "01.XYM" / "F814W", extensions=[".xOg"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "01.XYM" / "F606W", extensions=[".xOg"])
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "03.LOC_TRANS" / "F814W", filename = "xym2mat.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "03.LOC_TRANS" / "F814W", filename = "img2extract_wfc3uv_psflist.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "03.LOC_TRANS" / "F606W", filename = "xym2mat.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "03.LOC_TRANS" / "F606W", filename = "img2extract_wfc3uv_psflist.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "04.EXTRACT_PSF" / "F814W", filename = "uvp2psf_simst.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "04.EXTRACT_PSF" / "F606W", filename = "uvp2psf_simstV.xOg")

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F814W" / "1star-fit", filename = "mcmc_expand_average.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F814W" / "2star-fit", filename = "mcmc_expand_average.xOg")

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F606W" / "1star-fit", filename = "mcmc_expand_average.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F606W" / "2star-fit", filename = "mcmc_expand_average.xOg")

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F814W" / "1star-fit", filename = "uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F814W" / "2star-fit", filename = "uvp2tri_scon_fs_asym_mcmc.xOg")

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F606W" / "1star-fit", filename = "uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "06.FIT" / "F606W" / "2star-fit", filename = "uvp2tri_scon_fs_asym_mcmc.xOg")


    
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
        log_file = Path(directory).resolve() / f"loc_trans_run_xym2mat_{script.replace('.src','')}.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script_path = subdir / script if f == "F814W" else base_dir / "F606W" / script
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
        print(f"Finished run_xym2mat.src")

    def run_img2extract_wfc3uv_psflist(directory):
        """Run extraction for PSF list generation (simulation)."""
        log_file = Path(directory).resolve() / f"loc_trans_run_img2extract_wfc3uv_psflist.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_img2extract_wfc3uv_psflist_simst.src' if f == "F814W" else 'run_img2extract_wfc3uv_psflist_simstV.src'
                    script_path = base_dir / f / script
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
        print(f"Finished run_img2extract_wfc3uv_psflist")

    def run_img2extract_wfc3uv_psflist_Cal(directory):
        """Run extraction for PSF list generation (calibration)."""

        log_file = Path(directory).resolve() / f"loc_trans_run_img2extract_wfc3uv_psflist_Cal.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "03.LOC_TRANS"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_img2extract_wfc3uv_psflist_Cal.src' if f == "F814W" else 'run_img2extract_wfc3uv_psflist_CalV.src'
                    script_path = base_dir / f / script
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
        print(f"Finished run_img2extract_wfc3uv_psflist_Cal")


    print(directory)
    data_prep_loc_trans(directory)
    data_prep_loc_trans(directory, filters = 'F606W')
    run_xym2mat(directory)
    run_img2extract_wfc3uv_psflist(directory)
    run_img2extract_wfc3uv_psflist_Cal(directory)


def cmd_diagram(directory):

    copy_files(source=Path(directory).resolve() / "01.XYM" / "F814W", destination=Path(directory).resolve() / "02.CMD", extensions=[".02"])
    copy_files(source=Path(directory).resolve() / "01.XYM" / "F606W", destination=Path(directory).resolve() / "02.CMD", extensions=[".XYM"])
    subdir = Path(directory).resolve() / "02.CMD"
    #script_path = Path(directory).resolve() / "02.CMD" / "cmd_plot.py"
    #import numpy as np
    #import matplotlib.pyplot as plt
    
    #Load the MATCHUP files
    
    
    #xv, yv, mv are the x, y and magnitudes for the V band. Same logic for the I band
    
    xv, yv, mv = np.loadtxt(Path(directory).resolve() / "02.CMD" / "MATCHUP.F606W.XYM", unpack=True, usecols=(0,1,2))
    xi, yi, mi = np.loadtxt(Path(directory).resolve() / "02.CMD" / "MATCHUP.F814W.XYM.02", unpack=True, usecols=(0,1,2))
    
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
        log_file = Path(directory).resolve() / f"run_uvp2psf_simst_1.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "04.EXTRACT_PSF"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_uvp2psf_simst_1.src' if f == "F814W" else 'run_uvp2psf_simstV_1.src'
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_1")
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
                
    def run_uvp2psf_simst(directory, script='run_uvp2psf_simst_1.src'):
        """Finds a sky value for each star in each exposure using the pixels between 8.5 and 13.5 pixels of the center."""
        log_file = Path(directory).resolve() / f"run_uvp2psf_simst_1.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "04.EXTRACT_PSF"
                filters = ['F814W', 'F606W']
                for f in filters:
                    subdir = base_dir / f
                    script = 'run_uvp2psf_simst_2.src' if f == "F814W" else 'run_uvp2psf_simstV_2.src'
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_2")
        
    prepare_data(good_psf, directory)
    prepare_data(good_psf, directory, f='F606W')
    run_uvp2psf_simst(directory)


def tri_fit_dataprep_twostar(directory, f = 'F814W'):

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
    output_file_dir = base_dir / '06.FIT' / f / '2star-fit'
    output_file = os.path.join(output_file_dir, filename)
    
    with open(output_file, "w") as f:
        for line in content:
            f.write(line.rstrip() + "\n")

    print(f"File '{filename}' successfully created.")



def tri_fit_dataprep_threestar(directory, f = 'F814W'):

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


def tri_fit_dataprep_onestar(directory, f = 'F814W'):

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


def tri_fit_final_F814W_opt(source, directory):
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
        log_file = Path(directory).resolve() / f"uvp2tri_scon_fs_asym_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script = script
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_1")

        
    def run_mcmc_expand_average_814W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / f"run_mcmc_expand_average.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_mcmc_expand_average")
        

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename = "mcmc_expand_average.xOg")

    copy_entire_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", filename = "uvp2tri_scon_fs_asym_mcmc.xOg")
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "03.LOC_TRANS" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit",  extensions=[".gz"])
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "02.CMD", extensions=[".XYIVB_targ"], destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit")
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "06.FIT" / "F814W" / "3star-fit", extensions=[".fits"])
    copy_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F606W", destination=Path(directory).resolve() / "06.FIT" / "F606W" / "3star-fit", extensions=[".fits"])

    
    run_uvp2psf_simst_1(directory)
    #run_uvp2psf_simst_2(directory)
    run_mcmc_expand_average_814W(directory)
    #run_mcmc_expand_average_606W(directory)



def tri_fit_final_F814W(directory):
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
        log_file = Path(directory).resolve() / f"uvp2tri_scon_fs_asym_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                folders = ['1star-fit', '2star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script = script
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_1")

        
    def run_mcmc_expand_average_814W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / f"run_mcmc_expand_average.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F814W"
                folders = ['1star-fit', '2star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_mcmc_expand_average")
        
    
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
    #run_uvp2psf_simst_2(directory)
    run_mcmc_expand_average_814W(directory)
    #run_mcmc_expand_average_606W(directory)




def tri_fit_final_F606W(directory):
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
        log_file = Path(directory).resolve() / f"uvp2tri_scon_fs_asym_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                folders = ['1star-fit', '2star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_1")
        
    def run_mcmc_expand_average_606W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / f"run_mcmc_expand_average.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                folders = ['1star-fit', '2star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_mcmc_expand_average")
        
    run_uvp2psf_simst_2(directory)
    #run_uvp2psf_simst_2(directory)
    run_mcmc_expand_average_606W(directory)
    #run_mcmc_expand_average_606W(directory)


def tri_fit_final_F606W_opt(directory):
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
        log_file = Path(directory).resolve() / f"uvp2tri_scon_fs_asym_mcmc.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                folders = ['1star-fit', '2star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_uvp2psf_simst_1")
        
    def run_mcmc_expand_average_606W(directory, script='run_mcmc_expand_average.src'):
        log_file = Path(directory).resolve() / f"run_mcmc_expand_average.log"
        with open(log_file, "w") as logf:
            try:
                base_dir = Path(directory).resolve() / "06.FIT" / "F606W"
                folders = ['3star-fit']
                for f in folders:
                    subdir = base_dir / f
                    script_path = base_dir / f / script
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
        print(f"Finished run_mcmc_expand_average")
        
    run_uvp2psf_simst_2(directory)
    #run_uvp2psf_simst_2(directory)
    run_mcmc_expand_average_606W(directory)
    #run_mcmc_expand_average_606W(directory)


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

    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "MATCHUP.F606W.XYM")
    copy_entire_files(source=Path(directory).resolve() / "04.EXTRACT_PSF" / "F814W", destination=Path(directory).resolve() / "07.CALIBRATION" , filename = "MATCHUP.F814W.XYM.02")

    def psf_star_mags_mcmc(directory, script='run_psf_star_Imags_mcmc.src'):
        
        log_file = Path(directory).resolve() / f"run_psf_star_Imags_mcmc.log"
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

        log_file = Path(directory).resolve() / f"run_psf_star_Vmags_mcmc.log"
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
        
        log_file = Path(directory).resolve() / f"run_cal_star_num_2_MATCHUP.log"
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

    #psf_star_mags_mcmc(directory)
    cal_star_num(directory)
    #VI_HST_ogle_man_match4(directory)
    #fit_HST_IV_ogle_col_1(directory)
    
    

def calibration_hst_ogle_match(directory):
    """
    The goal here is to calibrate the HST photometry to the OGLE-III database.
    """

    def VI_HST_ogle_man_match4(directory, script='run_VI_HST_ogle_man_match4.src'):
        
        log_file = Path(directory).resolve() / f"run_VI_HST_ogle_man_match4.log"
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
        
        log_file = Path(directory).resolve() / f"run_fit_VI_HST_ogle_man_match4.log"
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
     
    
