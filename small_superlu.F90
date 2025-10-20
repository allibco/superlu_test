program small_superlu

  #include "superlu_dist_config.fh"
  use superlu_mod
  use superlu_bindings
  use dist_spmv_mod
  use iso_c_binding
  use mpi
  

  ! Local matrix storage
  integer :: iam, nprow, npcol, nprocs, info, i, ierr
  integer :: m_loc, fst_row
  integer :: n = 4   ! Global size
  integer :: nrhs = 1
  integer :: nnz_loc
  !integer(kind=c_int), allocatable :: rowptr(:), colind(:)
  real(kind=c_double), allocatable :: nzval(:), b(:), berr(:), y(:), x(:)
  integer, allocatable :: rowptr(:), colind(:), task_row_starts(:)
  integer, allocatable, target:: row_to_proc(:)

  
  !superlu structures
 
  integer(c_int64_t) :: A, grid
  integer(c_int64_t) :: options
  integer(c_int64_t) :: ScalePermstruct
  integer(c_int64_t) :: LUstruct
  integer(c_int64_t) :: SOLVEstruct
  integer(c_int64_t) :: stat
  !integer(c_int64_t) :: row_to_proc_handle
  !integer(c_int64_t) :: gsmv_comm_handle
  !integer(c_int64_t) :: nnn = 4 !global size

 type(halo_t) :: halo

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)


! Create Fortran handles for the C structures used in SuperLU_DIST
  call f_create_gridinfo_handle(grid)
  call f_create_options_handle(options)
  call f_dcreate_ScalePerm_handle(ScalePermstruct)
  call f_dcreate_LUstruct_handle(LUstruct)
  call f_dcreate_SOLVEstruct_handle(SOLVEstruct)
  call f_create_SuperMatrix_handle(A)
  call f_create_SuperLUStat_handle(stat)

  !gsmv_comm_handle = f_pdgsmv_comm_create()

  !print *, 'rank', iam, 'handles: A=', A, 'grid=', grid, 'options=', options

  
! Check we have exactly 2 processes
  if (nprocs /= 2) then
     if (iam == 0) print *, "This example requires exactly 2 MPI processes"
     call MPI_Finalize(ierr)
     stop
  end if


  ! Create the process grid
  nprow = nprocs
  npcol = 1
  call f_superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid)
  !print *, "GRIDINIT done"
  !print *, 'rank', iam, 'after gridinit: grid=', grid

  
  ! Local matrix partition (2 rows per proc for 4×4)
  if (iam == 0) then
     m_loc = 2
     fst_row = 0
     nnz_loc = 4
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc), berr(nrhs),  x(m_loc),  y(m_loc))
     rowptr = [0, 2, 4] !0-based
     colind = [0, 1, 1, 2]
     nzval  = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
     b = [5.0d0, 6.0d0]
  else
     m_loc = 2
     fst_row = 2
     nnz_loc = 4
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc), berr(nrhs),  x(m_loc),  y(m_loc))
     rowptr = [0, 2, 4] ! 0-based
     colind = [2, 3, 3, 0]
     nzval  = [5.0d0, 6.0d0, 7.0d0, 8.0d0]
     b = [7.0d0, 8.0d0]
  end if

  allocate(task_row_starts(3))
  task_row_starts=[0, 2, 4];
  
  !allocate(row_to_proc(n))
  !row_to_proc=[ 0, 0, 1, 1]
  !row_to_proc_handle = transfer(c_loc(row_to_proc), row_to_proc_handle)

  call dist_spmv_init(m_loc, fst_row, n, nprocs, iam, rowptr, colind, task_row_starts, comm, halo, ierr)
  
  print *, "Rank", iam, "m_loc=", m_loc, "nnz_loc=", nnz_loc, "fst_row=", fst_row

  
  ! Create the distributed compressed row matrix pointed to by the F90 handle A
  call f_dCreate_CompRowLoc_Mat_dist(A, n, n, nnz_loc, m_loc, fst_row, &
       nzval, colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE)

  ! Set the default input options
  call f_set_default_options(options)

  ! Change one or more options
  !could also try ColPerm=NATURAL, RowPerm=NOROWPERM)
  call set_superlu_options(options,ColPerm=COLAMD)
  call set_superlu_options(options,RowPerm=LargeDiag_MC64)
  
  ! Initialize ScalePermstruct and LUstruct
  call get_SuperMatrix(A, nrow=n, ncol=n)
  call f_dScalePermstructInit(n, n, ScalePermstruct)
  call f_dLUstructInit(n, n, LUstruct)

  ! Initialize the statistics variables
  call f_PStatInit(stat)

  !copy b into x
  x = b

  ! ---- test matvec BEFORE solve
  ! Allocate the communication structure

  ! Initialize pdgsmv communication
  !call f_pdgsmv_init(A, row_to_proc_handle, grid, gsmv_comm_handle)

  ! Now you can use pdgsmv multiple times
  !call f_pdgsmv(0_c_int64_t, A, grid, gsmv_comm_handle, x, y)

  !print *, 'rank', iam,  'pdgsmv BEFORE solver: x = ', x
  !print *, 'rank', iam,  'pdgsmv BEFORE solver: b = ', b
  !print *, 'rank', iam, 'pdgsmv BEFORE solver: Ab = ', y

  !call f_pdgsmv_comm_destroy(gsmv_comm_handle)
  
  write(*,*) "calling pdgssvx"

  ! Call the linear equation solver (soln overwrites x)
  call f_pdgssvx(options, A, ScalePermstruct, x, m_loc, nrhs, &
       grid, LUstruct, SOLVEstruct, berr, stat, info)

  print *, 'rank', iam, 'AFTER solver x=', x

  
  if (info == 0 .and. iam == 1) then
     write (*,*) 'Backward error: ', (berr(i), i = 1, nrhs)
  else
     write(*,*) 'INFO from f_pdgssvx = ', info
  endif

  ! do y = A*x (now x is soln , so y = b)
  !call f_pdgsmv_init(A, row_to_proc_handle, grid, gsmv_comm_handle)

  !call f_pdgsmv(0_c_int64_t, A, grid, gsmv_comm_handle, x, y)

  !print *, 'rank', iam,  'pdgsmv AFTER solver, x = ', x
  !print *, 'rank', iam,  'pdgsmv AFTER solver, b = ', b
  !print *, 'rank', iam, 'pdgsmv AFTER solver, Ax = ', y
  !call f_pdgsmv_comm_destroy(gsmv_comm_handle)


  
  !this (y) should = b
  !print *, 'b = ', b
  !print *, 'Ax = ', y
  !print *, 'rank', iam, 'pdgsmv AFTER solver, y=', y
  
  !deallocate the storage allocated by SuperLU_DIST
  call f_PStatFree(stat)
  !do not call - tries to free the fortran-allocated rowptr,colind and nzval array
  !call f_Destroy_CompRowLoc_Mat_dist(A)
  call f_dScalePermstructFree(ScalePermstruct)
  call f_dDestroy_LU_SOLVE_struct(options, n, grid, LUstruct, SOLVEstruct)

! Release the SuperLU process grid
  call f_superlu_gridexit(grid)

! Deallocate the C structures pointed to by the Fortran handles
  call f_destroy_gridinfo_handle(grid)
  call f_destroy_options_handle(options)
  call f_destroy_ScalePerm_handle(ScalePermstruct)
  call f_destroy_LUstruct_handle(LUstruct)
  call f_destroy_SOLVEstruct_handle(SOLVEstruct)
  call f_destroy_SuperMatrix_handle(A)
  call f_destroy_SuperLUStat_handle(stat)


  
  deallocate(rowptr)
  deallocate(colind, nzval)
  deallocate(b, x, y, berr, row_to_proc)

  call MPI_Finalize(ierr)

end program
