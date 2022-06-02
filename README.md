# CensGMR
CensGMR verion: 1.0

Authors: Ganzhong Tian <gavin.tian@emory.edu>; Benjamin Risk <benjamin.risk@emory.edu>

This project is an extension of the Gaussian Mixture Regression (GMR) model to handel censored multivariate responses.

The **MixCenMVReg_EM.R** is an EM algorithm implementation of the Censored GMR model wuth multivariate responses, whereas **MixCenUVReg_EM.R** is 
a univariate response Version.

The dependent R packages are listed as below:

* **matrixStats**  # dependent for MixCenMVReg_EM.R & MixCenUVReg_EM.R
* **condMVNorm**   # dependent for MixCenMVReg_EM.R
* **MomTrunc**     # dependent for MixCenMVReg_EM.R
* **truncnorm**    # dependent for MixCenUVReg_EM.R

Please see the jupyter-notebook **Example.ipynb** for instructions on how to use this software.

