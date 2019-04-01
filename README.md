# ECLIPS3D
Public Fortran 90 code for linear wave and circulation calculations, developed originally for planetary atmospheres, with python scripts provided for data analysis. 

Version : March 28, 2019. Developed by Florian Debras, article submitted to A&A. 

This README aims at describing the code and how to use it. Simple examples are provided. 

4 setups of the code are provided: 2D_axi, 2D_shallow, 3D, 3D_steady. We first describe these setups before explaining how to run the code. 

-------------------------------------------------------
SETUPS OF THE CODE
-------------------------------------------------------
  2D_axi: eigenvector setup in spherical coordinates assuming axisymmetry around the axis of rotation. A longitudinal wavenumber, m, must therefore be provided. 

  2D_shallow: eigenvector setup for shallow water beta-plane. The latitude of the beta plane and characteristic height can be changed.

  3D: eigenvector setup in full 3D, spherical coordinates.

  3D_steady: linear circulation setup, hence matrix inversion. A forcing and a dissipation have to be implemented for a linear steady state to exist. In the current version of ECLIPS3D, the implemented forcing follows Komacek and Showman 2016 and the dissipation is either following Komacek and Showman 2016 or Iro et al. 2005 or is constant through the atmosphere.

-------------------------------------------------------
CONTENT OF THE SETUPS
-------------------------------------------------------

No matter the setup, the directories are then organised as follows:

1) src - contains the source file. The name of the master file depends on the setup (normal_modes_2D_para.F90, normal_modes_shallow.F90, normal_modes_3D.F90 and standing_solver_3D.F90 respectively) but the other files are the same, except for the steady state solution.
Namely: - mod_data.F90: contains the data public to all modules of the code.
        - mod_init_para.F90: initialises MPI
        - mod_init_matrix.F90: distributes the matrix on the processors
        - mod_fill_matrix.F90: fills the matrix with the equations
        - mod_eigenvalues.F90: calculates eigenvalues and eigenvectors of the matrix. Does not exist for 3D_steady, instead mod_solver.F90 inverts the matric to obtain the steady state linear circulation.

Two other files allow to select and write the eigenvectors: study_eigenvectors.F90 read the eigenvector file and select some according to the selection procedures described in Debras et al. 2019. write_eigenvectors.F90 write the selected eigenvectors in a new file to feed them into python. NOTE: in the 3D_steady version, both these files are contained in read_solution.F90, that just allows to write the steady solution in a format suited for the given python file.

2) bin - contains the Makefile, and the intermediate compilation files. 

3) run - contains the input files and the exe file created by bin. 

4) data - contains the input data to launch the code as well as the output.

5) python - contains typical python files to generate an initial set of data and study the output. We provide simple example in these python files, namely axisymmetric state at rest and an initial, steady baroclinically unstable jet. 

-------------------------------------------------------
OVERVIEW OF THE CODE
-------------------------------------------------------

ECLIPS3D first initialises the parallel computing, performed by MPI, and allocates the memory on each processor to store the matrix to invert/for eigenvector calculation. The size on each processor is set by the parameter nb, and the performances of the program are very dependent to nb. It seems that nb ~ ntot/sqrt(n_proc)/10 is a reasonable number for good performances. 

Once initialised, the program calls mod_fill_matrix in order to calculate the coefficients of the matrix. In this version, we have decided to generate input data file from python (which deals with interpolations very easily) and read them afterwards in  fortran in ECLIPS3D. Changing that is no problem.

We then calculate the coordinates of the grid, the drag and radiative timescales if needed (a tentative at diffusion is also proposed, although not tested yet), the sound speed and brunt vaisala frequency. We then write an output file with the initial state so that it can be used to reconstruct dimensional variables once the program has run. 

Afterwards, we fill the matrix. There is a loop on ntot which set to one only the ith point and zero all the other one, followed by an inner loop on the coefficient of the matrix so that the coefficient i,j is the impact on the jth point from the ith point only. The beginning of the loop ensures that some points are always zero, according to boundary condition. The default setup assumes north south symmetry, hence v must be zero at the equator. The usual condition on v is that vcosphi must be zero at the pole, hence we set to zero all the terms that do not cancel in the equation of v when multiplied by cos phi and setting cos phi = 0. 

Finally, the last routine either inverts the matrix, calculates the full spectrum of eigenvectors or simply calculates the selected eigenvectors. This involves combining numerous sclapack routine, as described in the submitted paper. 

-
IMPORTANT NOTE: In this version, all the matrices are filled with complex quantities. This is inspired by the 2D, axisymmetric case with an initial state at rest where there is a diffenrence of phase of Pi between u and $v$, hence one is pure imaginary when the other is real. THIS HAS TO BE CHANGED, as keeping a complex matrix increases the size of the memory to store the matrix, whereas the real part of each number is always zero in the 3D version. The interest lies in the fact that reading the complex eigenvectors of a real matrix is different than just reading the eigenvector of a complex matrix in SCALAPACK, and we have not implemented the former case in python.
-
                                        
                                        
-------------------------------------------------------
SIMPLE CASES
-------------------------------------------------------
                                        
                                        
