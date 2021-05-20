PROGRAM normal_modes_3D

  USE MPI
  USE mod_init_para
  USE mod_data
  USE mod_init_matrix
  USE mod_3D_fill_matrix
  USE mod_eigenvalues




  IMPLICIT NONE

  INTEGER ::  ntab, ndutab, ndvtab, ndwtab, ndptab, nqtab
  INTEGER :: num_args, n_input_file,n_timescale_file
  CHARACTER(LEN=300) :: c_namelist_file, c_timescale_file


  CALL init_para

  IF (myrow==0 .AND. mycol==0) THEN
    num_args=command_argument_count()
    IF (num_args .lt. 1) THEN
      WRITE(6,*) "Please supply a namelist file as a command line argument."
      STOP
    END IF
    CALL get_command_argument(1,c_namelist_file)
    CALL get_command_argument(2,c_timescale_file)
    WRITE(6,*)"parameter file name=",trim(c_namelist_file)
  END IF

  CALL mpi_bcast(c_namelist_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)
  CALL mpi_bcast(c_timescale_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)
  n_input_file = 65
  n_timescale_file = 75

  OPEN(UNIT=n_input_file,FILE=trim(c_namelist_file),STATUS="OLD")

  READ(n_input_file, resolution)
  READ(n_input_file,folder)
  READ(n_input_file,planet)

  CLOSE(n_input_file)

  OPEN(UNIT=n_timescale_file,FILE=trim(c_timescale_file),STATUS="OLD")
  READ(n_timescale_file,timescales)
  CLOSE(n_timescale_file)

  ntot=nlong*(2*nlat*nz+(nlat+1)*nz+2*nlat*(nz))
  kappa = gascons/cp

  ntab=6*nlong*(2*nlat*nz+(nlat+1)*nz+nlat*(nz+1))
  !6 variables (u,v,w,theta,p,rho), each of which have 4 different arrays in the data.dat file. Two of those arrays are of size nlong*nlat*nz (hence the nlong*2*nlat*nz term. Then one of them has size nlong*(nlat+1)*nz and the other has nlong*nlat*(nz+1).

  ndutab=nz*nlat*nlong*7
  ! 7 lots of derivatives on u

  ndvtab=nz*(nlat+1)*nlong*7
  ! 7 lots of derivatives on v, which have nlat+1 latitude points

  ndptab=nz*nlat*nlong*8
  ! 8 lots of derivatives on p

  ndwtab=(nz+1)*nlat*nlong*10
  ! 10 lots of derivatives on w, which each have nz+1 points in the vertical

  ! nqtab=nz*nlat*nlong+(nz+1)*nlat*nlong
  nqtab=2*nz*nlat*nlong
  ! 2 lots of heating arrays, each of size nz*nlat*nlong


  IF (myrow==0 .AND. mycol==0) THEN
    print *, 'n', ntab, ndutab, ndvtab, ndptab, ndwtab, nqtab
    print *, 'Job starting'
  END IF

  CALL init_matrix
  CALL fill_matrix(ntab, ndutab, ndvtab, ndptab, ndwtab, nqtab)
  CALL evalues

  IF (myrow==0 .AND. mycol==0) THEN
    print *, 'Job Finished'
  END IF


  CALL BLACS_GRIDEXIT( info_txt )
  CALL BLACS_EXIT(1)


  CALL MPI_FINALIZE(ierr)

END PROGRAM
