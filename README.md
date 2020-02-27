# SD-QHD
This is C++ code that will compute electron trajectories in molecular orbitals using spin-dependent quantum hydrodynamics (SD-QHD). The code assumes that the molecular orbitals are linear combinations of atom-centered, Cartesian Gaussian-type orbitals (GTOs). It requires that the wavefunction information is written in an AIMPAC WFX format.

There are two folders. One folder (Grid) holds the code to compute trajectory information on a 3-dimensional Cartesian grid. The other folder (MonteCarlo) will compute full trajectories; the starting points are randomly assigned using acception/rejection method to ensure they are sampled according to the MO probability density. 

