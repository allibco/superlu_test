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
  integer :: n = 14   ! Global size
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
  !integer(c_int64_t) :: nnn = 14 !global size

 type(halo_t) :: halo

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  
! Check we have exactly 2 processes
  if (nprocs /= 4) then
     if (iam == 0) print *, "This example requires exactly 2 MPI processes"
     call MPI_Finalize(ierr)
     stop
  end if
  
  ! Local matrix partition (2 rows per proc for 4×4)
  if (iam == 0) then
     m_loc = 4
     fst_row = 0
     nnz_loc = 10
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc),  x(m_loc),  y(m_loc))
     rowptr = [0, 2, 4, 7, 10] !0-based
     colind = [0, 1, 1, 2, 2, 3, 8, 3, 4, 13]
     nzval  = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 1.0d0, 8.0d0, 3.0d0, 4.0d0]
     b = [1.0d0, 2.0d0, 5.0d0, 6.0d0]
  else if (iam ==1) then
     m_loc = 4
     fst_row = 4
     nnz_loc = 9
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc),  x(m_loc),  y(m_loc))
     rowptr = [0, 2, 4, 6, 9] ! 0-based
     colind = [1,4, 5, 8, 6,7, 6, 7,9 ]
     nzval  = [3.0d0, 4.0d0, 4.0d0, 6.0d0, 4.0d0, 5.0d0, 3.0d0 ,4.0d0, 6.0d0]
     b = [7.0d0, 8.0d0, 5.0d0, 5.0d0]
  else if (iam ==2) then
     m_loc = 4
     fst_row = 8
     nnz_loc = 4
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc),  x(m_loc),  y(m_loc))
     rowptr = [0, 1, 2, 3, 4] ! 0-based
     colind = [8, 9, 10, 11 ]
     nzval  = [8.0d0, 7.0d0, 6.0d0, 3.0d0]
     b = [10.0d0, 10.0d0, 10.0d0, 10.0d0]
  else !3
     m_loc = 2
     fst_row = 10
     nnz_loc = 7
     allocate(rowptr(m_loc+1), colind(nnz_loc), nzval(nnz_loc), b(m_loc),  x(m_loc),  y(m_loc))
     rowptr = [0, 3, 7] ! 0-based
     colind = [11, 12, 13, 0,6, 12, 13 ]
     nzval  = [1.0d0, 2.0d0, 1.0d0, 1.0d0, 2.0d0 , 1.0d0, 2.0d0]
     b = [15.0d0, 1.0d0]

     
  end if

  allocate(task_row_starts(5))
  task_row_starts=[0, 4, 8, 12, 14];
  

  call dist_spmv_init(m_loc, fst_row, n, nprocs, iam, rowptr, colind, task_row_starts, MPI_COMM_WORLD, halo, ierr)

  print*,'DID spmv init: iam = ', iam

  call dist_spmv(rowptr, colind, nzval, b, y, halo, MPI_COMM_WORLD, ierr)

  print*,'AFTERs spmv: iam = ', iam, 'A*b = ', y

  
  
  
  deallocate(rowptr)
  deallocate(colind, nzval)
  deallocate(b, x, y, task_row_starts)

  call MPI_Finalize(ierr)

end program
