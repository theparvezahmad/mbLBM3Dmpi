Development of an MPI-parallel 3D fluid flow solver based on multiblock lattice Boltzmann method(LBM).
It has been used to study flow around a cylinder, sphere and a NACA airfoil.
Developed a 3D fluid flow solver in Fortran.
Used D3Q19 velocity set, Multiple-relaxation-time.
Achieved 3D domain decomposition using MPI.
Standard LBM uses uniform grid.
Multiblock technique allows local refinement.
It creates new unknown points at the interface.
Bicubic spline interpolation employed.
