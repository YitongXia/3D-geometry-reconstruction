# 3D geometry reconstruction

 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/

This program is to reconstruct 3D geometry from corresponding image points of two images, it is implemented in C++ with Easy3D library.

The algorithm consists of the following step:
* Estimate fundamental matrix F using normalized 8-point algorithm;
* Recover relative pose from the fundamental matrix;
* Determine the 3D coordinates for all corresponding image points.

the initial two images are shown below:

<div align=center>  <img src="https://user-images.githubusercontent.com/75926656/169708252-38e2ea27-4da8-4c6b-a6ee-b03986d97a92.png" width = "500" alt="initial geometries" align=center/> </div>

the reconstructed geometries are shown below:

<div align=center>  <img src="https://user-images.githubusercontent.com/75926656/169708254-3554622f-dc41-42a6-974b-4f1930819d41.png" width = "500" alt="initial geometries" align=center/> </div>
