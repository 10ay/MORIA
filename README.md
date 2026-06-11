# MORIA: Microlensing Object high-Resolution Imaging Analysis (MORIA)

## Website

[Website Interface for MORIA](https://10ay.github.io/MORIA)

## How to Install

The easiest way to install MORIA is to first clone it. You can do this by running the command below:

    git clone https://github.com/10ay/MORIA.git

Once installed, open MORIA from the terminal using 'cd MORIA'. Then run:

    pip install -e .

Don't forget the dot after the '-e'. 

## Tutorial

We have notebooks in the 'demo/notebook' folder to show how to run each step of the pipeline. 

The order of notebooks is as follows:

    00.output_stacks.ipynb
    01.cmd_diagram.ipynb
    02.creating_psf.ipynb
    03.fitting_psfs.ipynb
    04.calibrations.ipynb

The recommended way to run MORIA is to copy everything in demo/notebook to the same repository level at which you want to perform the data analysis. 

The demo notebooks were run on the target KMT-BLG-2019-0253. In the 'demo' folder, we store both: Jupyter notebooks (in demo/notebook) and their respective PDFs (in demo/pdf_notebook)

Lastly, there is an optional notebook to facilitate a coordinate transformation between the highest resolution HST and Keck images. This can be used if you want to use the Keck analysis to constrain the lens-source separation in the HST images. This notebook is called:

    opt.keck_calib.ipynb

Note that keck_calib.ipynb was not run on KMT-BLG-2019-0253 but on a different target.

## Dependencies

To be able to run the notebooks and the pipeline, we recommend installing a fortran compiler along with the following packages:

    matplotlib
    pandas
    pathlib
    astropy
    numpy 

The entire pipeline is also summarized in detail [here](https://docs.google.com/document/d/1t8rLScKMqQ0oCqxvHKvAT6LXi1fT8aqSIPO8eTue8SQ/edit?tab=t.0)

## Acknowledgment 

This following scientific publication should be acknowledged by a citation to the works relevant to your study:

[Bhadra, T. D. et al. “You Shall Not Pass (Without Modeling): High-Resolution Analysis of KMT-2019-BLG-0253 using MORIA.”](https://arxiv.org/abs/2605.08340)

    @article{bhadra2026you,
            title={You Shall Not Pass (Without Modeling): High-Resolution Analysis of KMT-2019-BLG-0253 using MORIA},
            author={Bhadra, T and Terry, Sean K and Bennett, David P and Bhattacharya, Aparna and Bond, Ian A and Hulberg, Jon and Silva, Stela Ishitani and Mr{\u{A}}{\l}z, Przemek and Vandorou, Aikaterini},
            journal={arXiv preprint arXiv:2605.08340},
            year={2026}
            }

The Zenodo Citation for this package is:

[![Zenodo DOI: 10.5281/zenodo.20648200]](https://doi.org/10.5281/zenodo.20648200)


![Logo](https://github.com/10ay/MORIA/blob/main/new_logo.png)
