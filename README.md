[![GitHub Release](https://img.shields.io/github/v/release/eippertlab/spinalcord-fmri-7T)](https://github.com/eippertlab/spinalcord-fmri-7T/releases/tag/v1.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# spinalcord-fmri-7T

This repository is associated with the following manuscript (LINK). The data for this project will be made publicly available via OpenNeuro and OSF upon acceptance of the manuscript.

If you have any questions regarding this code, please feel free to reach out to uhorn@cbs.mpg.de.

### Content

The code is organized into the sections:
  dataset1 (containing processing scripts for the fMRI and physiological data in the screening dataset),
  dataset2_3 (containing processing scripts for the fMRI and physiological data in the discovery and validation datasets),
  dataset4 (containing processing scripts for the physiological data in the behavioral dataset), and
  figures (containing extra analyses and scripts to create the figures of the abovementioned manuscript)

In each of these folders there are bash and python scripts that can be run in the order in which they are named. Some of them need specific extra functions that are included in the subfolder helper_functions within each folder.

### Required software

* [Spinal Cord Toolbox (SCT)](https://spinalcordtoolbox.com/index.html)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)
* [MATLAB](https://de.mathworks.com/products/matlab.html)
* [ANTS](https://github.com/ANTsX/ANTs) (only for N4BiasFieldCorrection of anatomical images and shift in z direction of functional images)
* [AFNI](https://afni.nimh.nih.gov/) (only for 3dFWHMx command)
* [Python](https://www.python.org/) + packages listed in code

### Additional repos we rely on

* https://github.com/NYU-DiffusionMRI/mppca_denoise  (for denoising of functional images)
* https://sites.google.com/site/pierrickcoupe/softwares/super-resolution/monomodal  (for shift in z-direction of functional images)
* https://sites.google.com/site/pierrickcoupe/softwares/denoising/mri-denoising/mri-denoising-software  (for denoising of anatomical images)

### For preprocessing and calculation of results the following software was used:

* Bash: GNU bash, version 5.2.15(1)-release (x86_64-pc-linux-gnu)
* Spinal cord toolbox: 6.1
* FSL version 6.0.3
* ANTs version 2.3.1
* Python version 3.12
