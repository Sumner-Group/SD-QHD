Included is a shell script that will assist in writing the input file, the c++ code, a sample input/output file, and a wfx file.

The shell script that writes the input file requires the following information:
electron mass (1 in atomic units); time step (in atomic units); number of trajectories; range in x-drection; range in y-direction; range in z-direction (note these ranges are in atomic units and are the distance beyond the x,y, and z position of the atoms);  wrange - this is a constant function that should be grater than the maximum value of the MO density. This ensures that the trajectories are sampled according the the MO probability density; number of steps in the trajectory (the trajecotyr stops after one orbit); molecular orbital number; input file name for SD-QHD code; WFX file

To create the input file:
./sd-qhd-montecarloinput.script //then follow prompts

To compile the code:
g++ -O2 -o sd-qhd-montecarlo.exe sd-qhd-montecarlo.cpp

To run the code:
./sd-qhd-montecarlo.exe < inputfile >& outputfile

Sample inputs and outputs are given for the fluorine molecule molecula orbitals, calculated a the HF/aug-cc-pvqz level of theory using the Gaussian09 suite of programs.

Sample input is created from f2aug-cc-pvqz.wfx and sd-qhd-montecarloinput.script:
F2-MO8-monte.inp

Sample outputs are:
F2-MO8-monte.out

=========Important Note==========
To replicate the output, change the seed for the random number generator to 1. See lines 224 and 225 in the code
