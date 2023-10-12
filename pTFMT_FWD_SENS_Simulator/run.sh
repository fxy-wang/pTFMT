#!/bin/bash
g++ -std=c++0x -O3 -fopenmp main.cpp -o pTFMT

#./pTFMT ../setup/3D_mouse_x10mm_y10mm_z5mm_4s_256d_homo_sensitivity.ini
#./pTFMT ../setup/3D_mouse_x10mm_y10mm_z5mm_4s_1024d_homo_sensitivity.ini
./pTFMT ../setup/3D_mouse_x10mm_y10mm_z5mm_4s_1024d_1obj_fwd.ini

