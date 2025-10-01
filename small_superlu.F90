program small_superlu
  use iso_c_binding
  use superlu_mod   ! Assume you have interfaces for C bindings here
  implicit none

  include 'mpif.h'

  ! SuperLU types
  type(c_ptr) :: A, grid
  type(dScalePermstruct_t) :: ScalePermstruct
  type(dLUstruct_t) :: LUstruct
  type(superlu_dist_options_t) :: options
  type(SuperLUStat_t) :: stat


  ! Local matrix storage
  integer(C_INT) :: iam, nprow, npcol, nprocs, info, nrhs
  integer(C_INT) :: m_loc, fst_row
  integer(C_INT) :: n = 4   ! Global size
  integer(C_INT) :: nnz_loc
  integer(C_INT), allocatable, target :: rowptr(:), colind(:)
  real(C_DOUBLE), allocatable, target :: nzval(:), b(:), berr(:)

  call MPI_Init(info)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, info)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, info)

  ! Create 1D process grid
  nprow = nprocs
  npcol = 1
  call superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid)

  ! Local matrix partition (2 rows per proc for 4×4)
  if (iam == 0) then
     m_loc = 2
     fst_row = 0
     nnz_loc = 4
     allocate(rowptr(0:2), colind(4), nzval(4), b(2), berr(1))
     rowptr = [0, 2, 4]
     colind = [0, 1, 1, 2]
     nzval  = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
     b = [5.0d0, 6.0d0]
  else
     m_loc = 2
     fst_row = 2
     nnz_loc = 4
     allocate(rowptr(0:2), colind(4), nzval(4), b(2), berr(1))
     rowptr = [0, 2, 4]
     colind = [2, 3, 3, 0]
     nzval  = [5.0d0, 6.0d0, 7.0d0, 8.0d0]
     b = [7.0d0, 8.0d0]
  end if

  ! Create distributed matrix A
  call dCreate_CompRowLoc_Matrix_dist(A, n, n, m_loc, nnz_loc, &
       fst_row, c_loc(nzval), c_loc(colind), c_loc(rowptr), 0,1,0)

  ! Initialize SuperLU_DIST structures
  call dScalePermstructInit(n, n, ScalePermstruct)
  call dLUstructInit(n, LUstruct)
  call PStatInit(stat)

  ! Set options
  call set_default_options_dist(options)
  options%ColPerm = NATURAL
  options%RowPerm = NOROWPERM
  options%IterRefine = SLU_DOUBLE

  ! Solve Ax = b
  call pdgssvx(options, A, ScalePermstruct, c_loc(b), m_loc, 1, grid, LUstruct, c_loc(berr), stat, info)

  ! Print solution (just proc 0)
  if (iam == 0) then
     print *, "INFO =", info
     print *, "Solution x:"
     print *, b(1:2)
  else
     print *, "Rank", iam, "solution x:"
     print *, b(1:2)
  end if

  ! Cleanup
  
  call Destroy_SuperMatrix_Store_dist(A)
  call dScalePermstructFree(ScalePermstruct)
  call dLUstructFree(LUstruct)
  call dDestroy_LU(n,grid,LUstruct)
  call PStatFree(stat)
  call superlu_gridexit(grid)

  call MPI_Finalize()

end program
