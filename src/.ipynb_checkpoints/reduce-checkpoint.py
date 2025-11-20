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


"""
The HST fly star reduction pipeline depends on a certain directory structure.

<root_dir>/
    00.DATA/
    01.XYM/
"""

fcode_dir = '/Users/tmbhadra/Documents/Work/NASA/moira'

def copy_files(source, destination, extensions=[".fits"]):
    """
    Copy files from one folder to another.
    """
    if not os.path.exists(source):
        raise FileNotFoundError()

    for f in os.listdir(source):
        path = os.path.join(source, f)
        if os.path.isfile(path) and (any(f.endswith(e) for e in extensions)):
            shutil.copy2(source / f, os.path.join(destination, f))

def data_prep_early(source, destination):
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "00.DATA" / "F814W", extensions=[".xgf"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "00.DATA" / "F606W", extensions=[".xgf"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "01.XYM" / "F814W", extensions=[".xgf"])
    copy_files(source=Path(source).resolve() / "src" / "fortran_compile", destination=Path(destination).resolve() / "01.XYM" / "F606W", extensions=[".xgf"])


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

    # Make sure script exists
    script_path = base_dir / script

    if not script_path.exists():
        raise FileNotFoundError(f'Conversion script not found: {script_path}')

    # Subdirectories to process
    filters = ['F814W', 'F606W']

    for f in filters:
        subdir = base_dir / f
        if not subdir.exists():
            print(f'Missing subdirectory.')
            continue

        # Get list of _flc files
        flc_files = list(subdir.glob('*_flc.fits'))

        if not flc_files:
            print(f'Missing flc files')
            continue

        # Run the script inside the subdir
        subprocess.run(
            ['csh', str(script_path)],  # assumes csh script
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
        log_file = Path(directory).resolve() /'"run_img2xym_wfc3uv.log'
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
                    print("Done", f)
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

        print(f"Finished run_img2xym_wfc3uv")
        

    def run_xym2mat(directory, filters = ['F814W'], script='run_xym2mat_1.src'):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        log_file = Path(directory).resolve() / f'run_xym2mat_{script.replace('.src','')}.log'

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

        print(f"Finished run_xym2mat")
        

    def run_xym2bar(directory, filters = ['F814W'], script='run_xym2bar_1.src'):
        """
        Get a final list of photometry that allowed small zeropoint shifts for each exposure
        """
        log_file = Path(directory).resolve() / f'run_xym2bar_{script.replace('.src','')}.log'

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
        print(f"Finished run_xym2bar")

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


    #copy files
    #Repeat for F606W


def run_output_stack(directory, script='run_img2sam_wfc3uv_379.src'):
    """
    Create a stack of the scene in the reference frame.
    """

    log_file = Path(directory).resolve() / f'output_stack.log'

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
                    text=True,
                    check=False
                )
        finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__        
    print(f"Created output stacks")

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
    script_path = Path(directory).resolve() / "02.CMD" / "cmd_plot.py"
    subprocess.run(
                    ["python", str(script_path)],
                    cwd=subdir,
                )


    

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
    run_uvp2psf_simst(directory)
        
def extract_psf_2(directory):
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
        
    run_uvp2psf_simst(directory)



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


    
    
