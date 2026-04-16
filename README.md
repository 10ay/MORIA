# MORIA: Microlensing Object high-Resolution Imaging Analysis (MORIA)

To be able to run

## How to Install

The easiest way to install this code is to first clone it. You can do this by running the command below:

    git clone https://github.com/10ay/MORIA.git

Once installed, navigate to inside the MORIA folder and run:

    pip install -e .

Git clone the repository, and run 'pip install -e .' in the src folder. Don't forget the dot after the '-e'. 

## Tutorial

We have notebooks in the 'demo/notebook' folder to show how to run each step of the pipeline. The order of notebooks is as follows:

    output_stacks.ipynb
    cmd_diagram.ipynb
    creating_psf.ipynb
    fitting_psfs.ipynb
    calibrations.ipynb

The notebook was run on the target KMT-BLG-2019-0253. 

Lastly, there is an optional notebook to facilitate a coordinate transformation between the highest resolution HST and Keck images. This can be used if we want to use the Keck analysis to constrain the lens-source separation in the HST images. The notebook below gives instructions on how to run keck_calib.ipynb. 

    keck_calib.ipynb

Note that keck_calib.ipynb was not run with KMT-BLG-2019-0253 but a different target.

## Dependencies

To be able to run the notebooks and the pipeline, we recommend installing a fortran compiler along with the following packages:

    matplotlib
    pandas
    pathlib
    astropy
    numpy 

The entire pipeline is also summarized in detail [here](https://docs.google.com/document/d/1t8rLScKMqQ0oCqxvHKvAT6LXi1fT8aqSIPO8eTue8SQ/edit?tab=t.0)

![Logo](https://github.com/10ay/MORIA/blob/main/new_logo.png)
