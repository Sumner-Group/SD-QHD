# SD-QHD
This is C++ code that will compute electron trajectories in molecular orbitals using spin-dependent quantum hydrodynamics (SD-QHD). The code assumes that the molecular orbitals are linear combinations of atom-centered, Cartesian Gaussian-type orbitals (GTOs). It requires that the wavefunction information is written in an AIMPAC WFX format.

Included is a shell script that will assist in writing the input file. It requires the following information:

electron mass (1 in atomic units); number of x grid points; number of y grid points; number of z grid points; range in x-drection; range in y-direction; range in z-direction (note these ranges are in atomic units and are the distance beyond the x,y, and z position of tha atoms); molecular orbital number; input file name for SD-QHD code; WFX file

To create the input file:

./sd-qhd-gridinput.script //then follow prompts

To compile the code:

g++ -O2 -o sd-qhd-grid.exe sd-qhd-grid.cpp

To run the code:

./sd-qhd-grid.exe < inputfile >& outputfile

Sample inputs and outputs are given
