# Purpose

This repository contains the data used to generate the results published in XXX.


# Code
The processed data output by the imaging device is contained in Imaging_Data.zip.
Being very heavy, the unprocessed data (raw images, contour files) cannot be uploded here: please contact the corresponding author of the article if interested.

The data have been analyzed individually (measurement by measurement) to identify outliers; spherical objects as well as aggregates are removed by default.
More entries are removed based on the criteria indicated in the Supporting Information.
The procedure is run by removeOutliers.m, which is the function that loads the raw data, removes the outliers, and automatically saves the polished data as .mat file.

The contour levels for plotting each population are calculated by running saveContourLevels.m, which appends this data to the .mat files previously generated.

The results are displayed by running plotPopulations.m

Dependencies:
To run "plotPopulations.m", the "Violinplot-Matlab" library needs to be installed. It can be found under the following link:
https://github.com/bastibe/Violinplot-Matlab
