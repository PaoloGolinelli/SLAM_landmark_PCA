# SLAM_landmark_PCA
This repository contains the Matlab code to perform a simplified version of the SLAM (Simultaneous Localization And Mapping) algorithm. 

The algorithm works by matching consecutive Lidar scans using natural landmarks in the environment (edges and discontinuities) obtained through the use of PCA (Principal Component Analysis).
To match the landmarks from different scans the Ransac algorithm is applyed. 
Finally the map and the motion of the robot is reconstracted by integrating the roto-translation used to match the consecutive scans. 


# required toolbox
- Navigation Tolbox 


# Description of files
- lidarScans: list of consecutive Lidar scans provided by 
- keyP2: contains the function extract the landmarks in a scan and returns the coordinades and a description of the landmarks
- rototranslation_PCA: is a function that uses the PCA to find the rototranslation that matches the principal directions of the two sets of points it takes as input
- main_scan_randsac: perform the matching between two lidar scans
- main_SLAM: perform the SLAM given a series of consecutive Lidar scans
