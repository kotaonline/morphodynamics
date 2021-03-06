Please read the following carefully to understand how to run this code successfully.

## Runtime Environment

Requires MatPIV and its subdirectories to be accessible during runtime. This is currently hard-coded in morphodynamics.m. If you change the directory organization of the src files, please ensure that paths are set for the MatPIV files.

## Quality of input data

The quality of output depends on the quality of input -- Garbage in Garbage out! More often than not, stage drift becomes a factor during image acquisition. This drift causes undesirable artifacts in PIV analysis. For instance, very large values for velocities are artifically assumed due to stage movement even when negligible cell movement is present between successive frames. To avoid problems due to stage drift, an optional drift correction function (imgReg_wrapper.m) is included in the src directory. This function uses dftregistration.m developed by Manuel Guizar for sub-pixel registration of two images.

Depending on the cell type being analyzed, it might be prudent to determine the optimal frame rate for image acquisition. We found that 3-5 min between frames works well for many epithelial cell types that we have tested.

## Command-line for analysis

> morphodynamics('/input/tif/file','/input/parameter/file');

where input tif file is a multi-page tiff containing one frame of the timelapse movie per page. Sample input tif and parameter files are provided in the example directory. The parameter file contains basic information about the experiment such as camera scaling and frame rate. Please see the params.txt in the example directory for keywords and input format.

## Output tree

```
 % parentDir - directory in which input tif file is present
 %  |
 %  ---> expDir - has the same name as input tif file without ext
 %  |
 %  ---> Preprocess
 %        |
 %        ---> iframes - individual frames in the multiframe tiff
 %        |
 %        ---> overlay - tifs overlaid with detected scratch edge
 %        |
 %        ---> maskmat - mat file masks removing areas with no cells
 %        |
 %        ---> masktif - same mask files as in mat, but in tif format
 %  |
 %  ---> PIV - contains all PIV analysis
 %        |
 %        ---> velFields - contains images with quiver plots
 %        |
 %        ---> OrderParam - contains images with quiver plots
 %  |
 %  ---> Postprocess
 %        |
 %        ---> Kymograph - contains velocity and orderparam kymographs
 %        |
 %        ---> DistCorr  - contains spatial correlation data and plots
 %        |
 %        ---> BOD - contains biorthogonal decomposition data and plots
 %  |
 %  ---> log - temporary directory to hold all log files in parallel
 %  loops. Note that this directory is deleted once log files are
 %  consolidated.
 ```
    
