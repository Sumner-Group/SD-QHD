Included is a shell script that will assist in writing the input file, the c++ code, a sample input/output file, and a wfx file. 

The shell script that writes the input file requires the following information:

electron mass (1 in atomic units); number of x grid points; number of y grid points; number of z grid points; range in x-drection; range in y-direction; range in z-direction (note these ranges are in atomic units and are the distance beyond the x,y, and z position of the atoms); molecular orbital number; input file name for SD-QHD code; WFX file

To create the input file:

./sd-qhd-gridinput.script //then follow prompts

To compile the code:

g++ -O2 -o sd-qhd-grid.exe sd-qhd-grid.cpp

To run the code:

./sd-qhd-grid.exe < inputfile >& outputfile

Sample inputs and outputs are given for the fluorine molecule molecula orbitals, calculated a the HF/aug-cc-pvqz level of theory using the Gaussian09 suite of programs.  

Sample input is created from f2aug-cc-pvqz.wfx and sd-qhd-gridinput.script:

F2-MO8.inp

Sample outputs are:

F2-MO.out

==============Anticipated/Needed updates===========
1. Currently the code does not have a single option to switch between computing full trajectories, or points on a grid. This should change.
