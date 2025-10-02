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
  integer(C_INT) :: iam, nprow, npcol, nprocs, info
  integer(C_INT) :: m_loc, fst_row
  integer(C_INT) :: n = 4   ! Global size
  integer(C_INT) :: nrhs = 1
  integer(C_INT) :: nnz_loc
  integer(C_INT), allocatable, target :: rowptr(:), colind(:)
  real(C_DOUBLE), allocatable, target :: nzval(:), b(:), berr(:)



! SuperLU_DIST enum constants (from supermatrix.h)
  integer(C_INT), parameter :: SLU_NR_loc = 7
  integer(C_INT), parameter :: SLU_D = 1
  integer(C_INT), parameter :: SLU_GE = 0


  A = c_null_ptr
  
  call MPI_Init(info)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, info)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, info)

! Check we have exactly 2 processes
  if (nprocs /= 2) then
     if (iam == 0) print *, "This example requires exactly 2 MPI processes"
     call MPI_Finalize(info)
     stop
  end if
  
  ! Create 1D process grid
  nprow = nprocs
  npcol = 1
  call superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid)
  print *, "GRIDINIT"

  ! Local matrix partition (2 rows per proc for 4×4)
  if (iam == 0) then
     m_loc = 2
     fst_row = 0
     nnz_loc = 4
     allocate(rowptr(0:m_loc), colind(nnz_loc), nzval(nnz_loc), b(m_loc), berr(nrhs))
     rowptr = [0, 2, 4] !0-based
     colind = [0, 1, 1, 2]
     nzval  = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
     b = [5.0d0, 6.0d0]
  else
     m_loc = 2
     fst_row = 2
     nnz_loc = 4
     allocate(rowptr(0:m_loc), colind(nnz_loc), nzval(nnz_loc), b(m_loc), berr(nrhs))
     rowptr = [0, 2, 4] ! 0-based
     colind = [2, 3, 3, 0]
     nzval  = [5.0d0, 6.0d0, 7.0d0, 8.0d0]
     b = [7.0d0, 8.0d0]
  end if


    write(*,*) "rowptr=", rowptr(0:m_loc)
    write(*,*) "colind=", colind(1:nnz_loc)
    write(*,*) "values=", nzval(1:nnz_loc)

    print *, "Rank", iam, "m_loc=", m_loc, "nnz_loc=", nnz_loc, "fst_row=", fst_row

    
  ! Create distributed matrix A
    call dCreate_CompRowLoc_Matrix_dist(A, n, n, nnz_loc, m_loc, &
       fst_row, c_loc(nzval(1)), c_loc(colind(1)), c_loc(rowptr(0)), SLU_NR_loc, SLU_D, SLU_GE)
    write(*,*) "created A"
    if (.not. c_associated(A)) then
       write(*,*) "ERROR rank ", iam, ": A is NULL after creation"
    endif
    

  ! Initialize SuperLU_DIST structures
  call dScalePermstructInit(n, n, ScalePermstruct)
  call dLUstructInit(n, LUstruct)
  call PStatInit(stat)

  ! Set options
  call set_default_options_dist(options)
  options%ColPerm = 3
  options%RowPerm = 1


  write(*,*) "calling pdgssvx"

  ! Solve Ax = b
  call pdgssvx(options, A, ScalePermstruct, c_loc(b(1)), m_loc, 1, grid, LUstruct, c_null_ptr, c_loc(berr(1)), stat, info)

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

  deallocate(rowptr, colind, nzval, b, berr)

  
  call MPI_Finalize()

end program
