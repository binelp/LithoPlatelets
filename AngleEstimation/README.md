# Angle Estimation
This folder and it's subfolders contains all codes related to the angle distribution estimation. A short description of how the important, higher level scripts and functions interact is given here. To see the purposes of each individual file please check the file headers. The files which are not mentioned here just contain auxilliary functions which are called by the scripts and functions listed here. Some scripts require data produced and saved by other scripts to run, which is indicated in brackets behind the filename. In these cases, the scripts also require entering the path and filenames of the saved results to function.

### 1. Obtaining Orientation Distributions 
**From Experimental Images**
For estimating the angle distributions from real images pbtained from experiments, 
the main/wrapper script and their internal implementation functions are:

*wrapperReal.m  &rarr; angleEstimReal.m &rarr; getAngleCuboidQuat.m*

Executing *wrapperReal.m* requires a folder with subfolders containing the images used for the angle estimation.

**From the Simulation**
For estimating the accuracy of the angle estimation procedure by simulating 
the imaging process with cuboids, and comparing their real orientation with 
the orientation obtained by estimating it only from the cuboids projections, 
the main/wrapper script and their internal implementation functions are:

*simulateRotationVTB.m &rarr; VTB_main_rotation.m &rarr; technique_IA_rotation.m  &rarr; getAngleCuboidQuat.m*

## 2. Post-processsing and Plotting
- *plotting/physAnglesPlot.m* (requires data produced by running *wrapperReal.m* and *simulateRotationVTB.m*)
- *plotting/analyzeAngleDists.m* (requires data produced by running *wrapperReal.m*)
- *plotting/analyzeRotVTB.m* (requires data produced by running *simulateRotationVTB.m*)
- *combineDists.m* (requires data produced by running *analyzeAngleDists.m*)

