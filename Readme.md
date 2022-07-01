For the measurements performed for each batch of platelets the uncompressed videos have been recorded to allow storing the images.
The function runQuickBurst.m has been used for this purpose.

Later, the sequence of actions to generate the data is:

1.  Extract the images from the videos. This has been done with the function avi2png which is a wrapper of the
    extractAvi function part of the KIDDO software.
    Only the extracted images have been stored to save space and the videos have been permanently deleted.

2.  The images are processed by running EXEprocessPng.m
    The function relies on processPng.m which runs the C++ code staticAnalysis, which is functionally equal to BURST_pump.
    Based on the known particle size, some outliers are removed already at this stage. In particular,
    particles with contour lower than 2*round(L2+L3) and larger than round(2*1.2*(sqrt(L3^2+L2^2)+L1)) are discarded
    (minimum possible observable contour and maximum observable possible contour x1.2 as safety margin.

3.  The data is processed by running EXEprocessData.m which is a wrapper for processMuDISCOData.

4.  The data has been analyzed individually (measurement by measurement) to identify outliers;
    Spherical objects as well as aggregates are removed by default.
    More entries are removed based on
        a. Visual hull
        b. Asbsolute L1, L2, and L3 sizes
    The thresholds are summarized in removeOutliers.m, which is the function that has to be run to load the raw data,
    remove the outliers, and automatically save the polished data as .mat file.

5.  The contour levels for each population are calculated by running saveContourLevels.m, which appends this data
    to the .mat files previously generated.

6.  The results are displayed by running plotPopulations.m

Dependancies:
To run "plotPopulations.m", the "Violinplot-Matlab" library needs to be installed. It can be found under the following link:
https://github.com/bastibe/Violinplot-Matlab