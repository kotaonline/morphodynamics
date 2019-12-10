Please read the following carefully to understand how to run this code successfully.

## Setting the environment

Make sure the src directory and the MatPIV subdirectory are in your MatLab path. Alternatively, you could run the code from the src directory with relative links to input directories

## Quality of input data

The quality of output depends on the quality of input -- Garbage in Garbage out! More often than not, stage drift becomes a factor during image acquisition. This drift causes undesirable artifacts in PIV analysis. For instance, very large values for velocities are artifically assumed due to stage movement even when negligible cell movement is present between successive frames. To avoid problems due to stage drift, an optional drift correction function (imgReg_wrapper.m) is included in the src directory. This function uses dftregistration.m developed by Manuel Guizar for sub-pixel registration of two images.

## Command-line for analysis

