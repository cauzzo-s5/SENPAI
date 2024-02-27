# SENPAI
 From tissue to segmentation: a modular framework for multi-scale neuron isolation
 
%%%%%%%%%%%%
Manuscript: "From tissue to segmentation: a modular framework for multi-scale neuron isolation", Cauzzo et al.
%%%%%%%%%%%%

!!!! System Requirements and Dependencies:
This code is written in the Matlab Environment and intended for use in Matlab. The code will be made freely available to the scientific community upon publication.
The toolbox requires the installation of Matlab (MathWorks Inc.). Additional required Matlab packages are: 
- 'Statistics and Machine Learning Toolbox'
- 'Image Processing Toolbox'
The toolbox has been written and tested using Matlab version R2022b (version 9.13), 'Statistics and Machine Learning Toolbox' version 12.4, and 'Image Processing Toolbox' 11.6. 
The toolbox does not require Matlab to be installed on a specific OS, and has been tested on Microsoft Windows, Linux Ubuntu and Apple iOS.

!!!! Installation guide: 
once Matlab is installed, the toolbox only requires to download the present folder and add it to the Matlab path. No further installation is required.

!!!! Demo:
The script "test_script.m" can be run from within this folder to test the segmentation and parcellation step on a 40x confocal microscopy image stack.
NOTE WELL: CD IN THE DOWNLOADED FOLDER BEFORE RUNNING THE TEST SCRIPT!
"somas.mat" and "markers.mat" are files that are needed to run the test script.
"somas.mat" is a logical 3D matrix containing markers for somas, that can be produced with any automatic, semiautomatic or manual segmentation routine.
"markers.mat" is a logical 3D matrix containing additional markers produced with senpai_prune, a posteriori parcellation correction routine.
The test script prompts a 3D rendering of 10 selected close-by segmented neurons.
The expected output can be seen in Figure 3, upper-right panel, in the main manuscript.
Expected run time on a "normal" desktop computer is around 10 minutes.


STRUCTURE OF THE TOOLBOX:

MAIN FUNCTIONS:

%%%%
Function "senpai_seg_core_v4" is the core of the segmentation algorithm, which takes in input a .tif image stack and provides a binary segmentation.

%%%%
Function "senpai_separator" is intended for parcellation of segmentations that failed in separating neurons from one another.

%%%%
Function "senpai_spinecatch" is intended for parcellation of segmentations that failed in connecting structures belonging to the same neuron.

%%%%
Function "senpai_prune" runs a GUI for a posteriori correction of parcellations.

OTHER FUNCTIONS:

%%%%
Function "senpai_estimateK" is intended for estimation of the optimal K parameter of the K-means embedded in senpai_seg_core.

%%%%
Function "senpai_skeletonize" is intended for skeletonization of segmentations, and provides a tree structure and a matrix in SWC notation.

%%%%
Function "senpai_strahlerord" extracts from the skeleton statistics based on Strahler Ordering.


!!!!!!
Detailed instructions for each function are reported in the header of the .m file.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
For any question please contact authors at simone.cauzzo@ing.unipi.it / cauzzo.simone@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
