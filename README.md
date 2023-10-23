
Geodesic Models with convexity shape prior - A scheme that can compute curvature-penalized minimal paths (or shortest paths) in an orientation-lifted space, whose physical projection curves are simple closed and convex. Furthermore, this scheme also provides a minimal paths solution to the region-based image segmentation problem, featuring both convexity shape prior and curvature penalization, in a curve evolution manner. As a consequence, this scheme gives an efficient and reliable solution to the region-based curvature-penalized active contour problem.

Copyright (C) 2023 Da Chen, Shandong Artificial Intelligence Institute, Jinan, China. dachen.cn 'at' hotmail 'dot' com.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Geodesic Models with convexity shape prior

## Program description
This code implements the curvature-penalized minimal paths framework for solving the interactive image segmentation problem and the region-based active contour problem. Three types of curvature-penalized terms are considered, which are derived from the forward variant of the Reeds-Shepp car model, the Euler-Mumford elastica model and the Dubins car model. It computes simple closed and convex curves as the boundaries of target regions. The minimal paths are tracked by estimating the viscosity solution to the Hamiltonian-Jacobi-Bellman equation, using state-of-art Hamiltonian Fast-Marching method. The code is implemented in MATLAB in conjunction with C++.

## Compilation instructions
- Most of the C++ parts of the code are adapted from the Hamiltonian Fast-Marching which is written in C++17. You will need a recent compiler, such as gcc 8.0, VS 15, or a recent clang. (Compilation is known to fail on gcc 7.4 and VS 14.)

## Steps for the installation of the MATLAB usage
- Step 1: Go to the folder “MexFiles_SrcFiles”.
- Step 2: Run the “mexAllHFMFilters.m” files. This will generate four MEX files, and will store them in the folder “MexFiles_Binary”.
- Step 3: Go to the folder “ComputeCircularGeodesicPaths/CostConstruction/BhattacharyyaCoefficient/”, and run the file “mexBhattacharyyaCoefficient.m”. This is used to compile the region-based appearance model based on the Bhattacharyya Coefficient of two histograms. This is optional.
- Step 4: Go to the main directory and run the files “example_LandMarkPts_ConvexityDB” and “example_Scribbles_ConvexityDB” for examples.

## Fair use
If you use this program for an academic or commercial project, then please cite at least two of the following papers.

- Da Chen, Jean-Marie Mirebeau, Minglei Shu, Xuecheng Tai and Laurent D. Cohen, Geodesic models with convexity shape prior, Vol. 45(7):8433-8452, IEEE Trans. on PAMI, 2023.

- Da Chen, Jean-Marie Mirebeau, Xuecheng Tai and Laurent D. Cohen, An elastica geodesic approach with convexity shape prior, in Proc. ICCV 2021.

- Jean-Marie Mirebeau, Fast Marching methods for Curvature Penalized Shortest Paths, Journal of Mathematical Imaging and Vision, Vol.60:784–815, 2018.







