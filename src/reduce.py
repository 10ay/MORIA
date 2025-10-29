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
        print(subdir)
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

    base_dir = Path(directory).resolve()


    filters = ['F814W', 'F606W']

    for f in filters:
        subdir = base_dir / f
        in_img2sam_wfc3uv = 'IN.img2sam_wfc3uv'
        in_xym2bar_1 = 'IN.xym2bar.1'
        in_xym2bar_2 = 'IN.xym2bar.2'
        in_xym2mat_1 = 'IN.xym2mat.1'
        in_xym2mat_2 = 'IN.xym2mat.2'

        base_dir_one = base_dir / '01.XYM'/ f
        files = sorted([f for f in os.listdir(base_dir_one) if f.endswith('WJ2.xym')])

        output_file_dir = base_dir / '01.XYM' / f

        output_file_img2sam = os.path.join(output_file_dir, in_img2sam_wfc3uv)
        output_file_xym2mat1 = os.path.join(output_file_dir, in_xym2mat_1)
        output_file_xym2mat2 = os.path.join(output_file_dir, in_xym2mat_2)
        output_file_xym2bar1 = os.path.join(output_file_dir, in_xym2bar_1)
        output_file_xym2bar2 = os.path.join(output_file_dir, in_xym2bar_2)

        

        with open(output_file_img2sam, "w") as f:
            for i, filename in enumerate(files, start=1):
                if f == 'F606W':
                    f.write(f"{i:02d} \"{filename}\" 6 0\n")
                else:
                    f.write(f"{i:02d} \"{filename}\" 8 0\n")

        with open(output_file_xym2bar1, "w") as f:
            if f == 'F606W':
                f.write("00 MATCHUP.F814W.XYM.02\n")
                for i, filename in enumerate(files, start=1):
                    f.write(f"{i:02d} {filename} c8 f8 z0\n")
            else:
                f.write(f"{i:02d} \"{filename}\" c8 f8 z0\n")

        with open(output_file_xym2bar2, "w") as f:
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 z0\n")

        with open(output_file_xym2mat1, "w") as f:
            if f == 'F606W':
                f.write("#00 MATCHUP.XYM.02 c0\n")
            else:
                 f.write("#00 MATCHUP.XYM.01 c0\n")
            for i, filename in enumerate(files, start=1):
                if f == 'F606W':
                    f.write(f"{i:02d} {filename} c8 f6 \"m-14.75,-5.5\" \n")
                else:
                    #if i == 1:
                        #f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")
                    f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")

        with open(output_file_xym2mat2, "w") as f:
            f.write("00 MATCHUP.XYM.01 c0\n")
            for i, filename in enumerate(files, start=1):
                f.write(f"{i:02d} {filename} c8 f8 \"m-13.75,-8.5\" \n")




    return

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

    def run_xym2mat(directory, script='run_xym2mat_1.src'):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        log_file = Path(directory).resolve() / f'run_xym2mat_{script.replace('.src','')}.log'

        with open(log_file, "w") as logf:
            sys.stdout = sys.stderr = logf
            try:
                base_dir = Path(directory).resolve() / "01.XYM"
                filters = ['F814W', 'F606W']

                for f in filters:
                    subdir = base_dir / f
                    script_to_use = script if f == 'F814W' else 'run_xym2mat_VI.src'
                    script_path = subdir / script_to_use

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
                    if f == "F814W":
                        copy_files(base_dir / "F814W", base_dir / "F606W", extensions=[".XYMEEE"])
                        copy_files(base_dir / "F814W", base_dir / "F606W", extensions=[".02"])
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

        print(f"Finished run_xym2mat")

    def run_xym2bar(directory, script='run_xym2bar_1.src'):
        """
        Get a final list of photometry that allowed small zeropoint shifts for each exposure
        """
        log_file = Path(directory).resolve() / f'run_xym2bar_{script.replace('.src','')}.log'

        with open(log_file, "w") as logf:
            sys.stdout = sys.stderr = logf
            try:
                base_dir = Path(directory).resolve() / "01.XYM"
                filters = ['F814W', 'F606W']

                for f in filters:
                    subdir = base_dir / f
                    script_to_use = script if f == 'F814W' else 'run_xym2bar.src'
                    script_path = subdir / script_to_use

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
        print(f"Finished run_xym2bar")

    run_img2xym(directory)
    run_xym2mat(directory)
    run_xym2bar(directory)
    run_xym2mat(directory, script='run_xym2mat_2.src')
    run_xym2bar(directory, script='run_xym2bar_2.src')


def run_output_stack(directory, script='run_img2sam_wfc3uv_379.src'):
    """
    Create a stack of the scene in the reference frame.
    """

    base_dir = Path(directory).resolve() / "01.XYM"
    filters = ['F814W', 'F606W']

    for f in filters:
        subdir = base_dir / f
        script_path = subdir / script if f == "F814W" else subdir / script

        print(base_dir)
        print(script_path)

        subprocess.run(
            ["csh", str(script_path)],  # assumes csh script
            cwd=subdir,
            text=True,
            check=False
        )
        
    print(f"Created output stacks")

def loc_trans(directory):
    """
    Run the local transformation scripts in F814W and F606W subdirectories to extract 
    the pixels from each exposure and accurately transform their locations into the 
    reference frame so that we can use them to solve for a PSF and then use this PSF to model the target star.

    Parameters
    ----------
    directory - The root directory to operate on. There is a specific directory
    structure that is expected within <dir>:
        00.DATA/
        01.XYM/

    Returns
    -------
    Pixels for PSF generation
    """
    copy_files(source, destination, extensions=[".xym"])

    def run_xym2mat(directory, script='run_xym2mat.src'):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        base_dir = Path(directory).resolve()/"03.LOC_TRANS_final"

        filters = ['F814W', 'F606W']

        for f in filters:
            subdir = base_dir / f
            if f == "F814W":
                script_path = base_dir/"F814W"/script
            else:
                script_path = base_dir/"F606W"/script

            print(base_dir)
            print(script_path)

            subprocess.run(
                ["csh", str(script_path)],  # assumes csh script
                cwd=subdir,
                check=False
            )

    def run_img2extract_wfc3uv_psflist(directory):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        base_dir = Path(directory).resolve()/"03.LOC_TRANS_final"

        filters = ['F814W', 'F606W']

        for f in filters:
            subdir = base_dir / f
            if f == "F814W":
                script = 'run_img2extract_wfc3uv_psflist_simst.src'
                script_path = base_dir/"F814W"/script
            else:
                script = 'run_img2extract_wfc3uv_psflist_simstV.src'
                script_path = base_dir/"F606W"/script

            print(base_dir)
            print(script_path)

            subprocess.run(
                ["csh", str(script_path)],  # assumes csh script
                cwd=subdir,
                check=False
            )

    def run_img2extract_wfc3uv_psflist_Cal(directory):
        """
        Produces TRANS.xym2mat, as well as 16 MAT.0 files.
        """
        base_dir = Path(directory).resolve()/"03.LOC_TRANS_final"

        filters = ['F814W', 'F606W']

        for f in filters:
            subdir = base_dir / f
            if f == "F814W":
                script = 'run_img2extract_wfc3uv_psflist_Cal.src'
                script_path = base_dir/"F814W"/script
            else:
                script = 'run_img2extract_wfc3uv_psflist_CalV.src'
                script_path = base_dir/"F606W"/script

            print(base_dir)
            print(script_path)

            subprocess.run(
                ["csh", str(script_path)],  # assumes csh script
                cwd=subdir,
                check=False
            )

    copy_files(source = Path(directory).resolve()/"02.XYM"/"F814W", destination = Path(directory).resolve()/"03.LOC_TRANS_final"/"F814W")
    copy_files(source = Path(directory).resolve()/"02.XYM"/"F606W", destination = Path(directory).resolve()/"03.LOC_TRANS_final"/"F606W")
    run_xym2mat(directory)
    run_img2extract_wfc3uv_psflist(directory)
    run_img2extract_wfc3uv_psflist_Cal(directory)
    

    
    
