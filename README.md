# MORIA: Microlensing Object high-Resolution Imaging Analysis (MORIA)

## How to Install

The easiest way to install this code is to first clone it. You can do this by running the command below:

    git clone https://github.com/10ay/MORIA.git

Once installed, navigate to inside the MORIA folder and run:

    pip install -e

Git clone the repository, and run 'pip install -e' in the src folder. 

## Tutorial

We have notebooks in the 'demo/notebook' folder to show how to run each step of the pipeline. The order of notebooks is as follows:

    output_stacks.ipynb
    cmd_diagram.ipynb
    creating_psf.ipynb
    fitting_psfs.ipynb
    calibrations.ipynb

Lastly, there is an optional notebook to facilitate a coordinate transformation between the highest resolution HST and Keck images. This can be used if we want to use the Keck analysis to constrain the lens-source separation in the HST images. The notebook below gives instructions on how to run keck_calib.ipynb. 

    keck_calib.ipynb

The entire pipeline is also summarized in detail [here](https://docs.google.com/document/d/1t8rLScKMqQ0oCqxvHKvAT6LXi1fT8aqSIPO8eTue8SQ/edit?tab=t.0)

