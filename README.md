# Morphodynamics

## Morphodynamics -- what is it for?

Morphodynamics was developed for analysis of collective cell migration. This code was developed and tested using Matlab R2014b.

### How to cite this work?

If you end up using all or part of this code, please cite:

Kota, P., Terrell, E.M., Ritt, D.A., Insinna, C., Westlake, C.J., and Morrison, D.K. M-Ras/Shoc2 signaling modulates E-cadherin turnover and cell–cell adhesion during collective cell migration. Proc. Natl. Acad. Sci. USA 116:(9)3536-3545 (2019)

### How to run this code?

Please see the example directory for sample results
Please note that this package comes with MatPIV 1.7 which was used for developing the code. A few changes were made to the MatPIV code to fit our requirements. We strongly recommend using the MatPIV version provided with this package. Please cite MatPIV along with the citation to our publication provided above.

### What do I need to run this code?

For starts, you need MatLab (Duh!). Here are additional dependencies

>> [fList, pList] = matlab.codetools.requiredFilesAndProducts('morphodynamics.m');
>> {pList.Name}'
ans = 
    'MATLAB'
    'Image Processing Toolbox'
    'Statistics Toolbox'
