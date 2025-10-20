module dist_spmv_mod
  use iso_c_binding, only: c_int, c_double
  use mpi
  implicit none
   !public :: dist_spmv_init
  !public :: dist_spmv, dist_spmv_free

  type halo_t
     integer :: nprocs_local      ! communicator size
     integer :: rank
     integer :: n_global
     integer :: fst_row, last_row, m_loc
     integer, allocatable :: owners(:)    ! owners for each halo col (size nhalo)
     integer, allocatable :: halo_cols(:) ! global column indices of halo entries (sorted)
     integer :: nhalo
     integer, allocatable :: sendcounts(:), sdispls(:)
     integer, allocatable :: recvcounts(:), rdispls(:)
     integer, allocatable :: recv_from(:)  ! ranks we expect data from (list)
  end type halo_t

contains

  !----------------------------------------------------------------------
  ! dist_spmv_init
  !  Precompute halo pattern and communication groups.
  !  Input:
  !    m_loc, fst_row, n_global, ntasks, myrank,rowptr, colind, task_row_starts
  !    (task_row_starts are 0-index)
  !  Output:
  !    halo (type halo_t) - must be passed to dist_spmv and freed with dist_spmv_free
  !----------------------------------------------------------------------

  subroutine dist_spmv_init(m_loc, fst_row, n_global, ntasks, myrank, rowptr, colind, task_row_starts, comm, halo, ierr)
    integer(c_int), intent(in) :: m_loc, fst_row, n_global, ntasks, myrank
    integer(c_int), intent(in) :: rowptr(:), colind(:)
    integer(c_int), intent(in) :: task_row_starts(ntasks+1)
    integer, intent(in) :: comm
    type(halo_t), intent(out) :: halo
    integer, intent(out) :: ierr

    integer :: i, j, k, nn, rank, nprocs, jp
    integer, allocatable :: tmp(:)
    integer :: owner
    integer :: rows_per_proc_num, rows_per_proc_den
    integer :: last_row
    integer, allocatable :: owner_set(:)
    integer :: nowners

    ierr = 0

    halo%rank = myrank
    halo%nprocs_local = ntasks
    halo%n_global = n_global
    halo%fst_row = fst_row
    halo%m_loc = m_loc
    last_row = fst_row + m_loc - 1
    halo%last_row = last_row

    ! 1) Collect remote column indices (may have duplicates)
    allocate(tmp(0)) !empty array
    nn = 0
    do i = 1, m_loc
       do j = rowptr(i), rowptr(i+1)-1
          !add one to j since row_ptr is 0-based
          jp = j+1
          if (colind(jp) < fst_row .or. colind(jp) > last_row) then
             nn = nn + 1
             if (size(tmp) == 0) then
                deallocate(tmp)
                allocate(tmp(1))
                tmp(1) = colind(jp)
             else
                tmp = [tmp, colind(jp)]
             end if
          end if
       end do
    end do

    if (nn == 0) then
       ! no halo needed (assign empty arrays)
       halo%nhalo = 0
       allocate(halo%halo_cols(0), halo%owners(0))
       allocate(halo%sendcounts(0), halo%sdispls(0))
       allocate(halo%recvcounts(0), halo%rdispls(0))
       allocate(halo%recv_from(0))
       if (allocated(tmp)) deallocate(tmp)
       return
    end if

    ! 2) Make unique and sort (simple O(n^2) unique - fine for moderate halo sizes)
    call unique_sort_int(tmp, nn)

    halo%nhalo = nn
    allocate(halo%halo_cols(nn))
    halo%halo_cols = tmp(1:nn)

    ! 3) Determine owner of each halo column (use task_row_starts)
    !   halo%halo_cols and task_row_starts are sorted, so we can go in order
    allocate(halo%owners(halo%nhalo))
    allocate(halo%sendcounts(halo%nprocs_local))
    halo%sendcounts = 0
    halo%owners = 0

    rank = 0 ! first rank
    ! Advance to the correct rank if halo_cols(k) has crossed next boundary
    do k = 1, halo%nhalo
       do while (rank < halo%nprocs_local -1 .and. halo%halo_cols(k) >= task_row_starts(rank+1))
          rank = rank + 1
       end do
       halo%owners(k) = rank !0-based owner
       !Build per-destination sendcounts
       halo%sendcounts(rank + 1) = halo%sendcounts(rank + 1) + 1  ! +1 since Fortran arrays are 1-based
    end do

    ! 4) Compute send displacements
    allocate(halo%sdispls(halo%nprocs_local))
    halo%sdispls(1) = 0
    do rank = 2, ntasks
        halo%sdispls(rank) = halo%sdispls(rank - 1) + halo%sendcounts(rank - 1)
    end do

    ! 5) Determine recvcounts by exchanging sendcounts
    allocate(halo%recvcounts(halo%nprocs_local))
    halo%recvcounts = 0
    call MPI_Alltoall(halo%sendcounts, 1, MPI_INTEGER, halo%recvcounts, 1, MPI_INTEGER, comm, ierr)

    !calc recv displacements
    allocate(halo%rdispls(halo%nprocs_local))
    halo%rdispls = 0
    do rank = 2, halo%nprocs_local
       halo%rdispls(rank) = halo%rdispls(rank-1) + halo%recvcounts(rank-1)
    end do

    ! Build list of ranks we will receive from (nonzero recvcounts)
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%recvcounts(rank) > 0) then
          nowners = nowners + 1
       end if
    end do
    allocate(halo%recv_from(nowners))
    nowners = 0
    do rank = 1, halo%nprocs_local
       if (halo%recvcounts(rank) > 0) then
          nowners = nowners + 1
          halo%recv_from(nowners) = rank-1 !0-based
       end if
    end do

    ! Clean up
    deallocate(tmp)

  end subroutine dist_spmv_init

  !----------------------------------------------------------------------
  ! dist_spmv - perform y_local = A_local * x_global
  !   Inputs:
  !      rowptr, colind, nzval  - local CSR (rowptr(1:m_loc+1), colind(1:nnz), nzval(1:nnz))
  !      x_local(1:m_loc)       - local portion of x (corresponds to fst_row..last_row)
  !      halo (precomputed)
  !   Output:
  !      y_local(1:m_loc)
  !----------------------------------------------------------------------

!!$  subroutine dist_spmv(rowptr, colind, nzval, x_local, y_local, halo, comm, ierr)
!!$    integer(c_int), intent(in) :: rowptr(0:), colind(:)
!!$    real(c_double), intent(in) :: nzval(:)
!!$    real(c_double), intent(in) :: x_local(:)
!!$    real(c_double), intent(out) :: y_local(:)
!!$    type(halo_t), intent(in) :: halo
!!$    integer, intent(in), optional :: comm
!!$    integer, intent(out) :: ierr
!!$
!!$    integer :: use_comm
!!$    integer :: i, j, k, p, idx, owner, pos
!!$    integer :: m_loc
!!$    integer :: nprocs
!!$    integer, parameter :: dp_kind = c_double
!!$    integer :: reqs_count, total_reqs
!!$    integer, allocatable :: reqs(:)
!!$    integer, allocatable :: stats(:)
!!$    real(c_double), allocatable :: sendbuf(:), recvbuf(:)
!!$    integer :: sdis, rdis
!!$    integer :: nnei
!!$
!!$    ierr = 0
!!$    if (present(comm)) then
!!$       use_comm = comm
!!$    else
!!$       use_comm = MPI_COMM_WORLD
!!$    end if
!!$
!!$    m_loc = halo%m_loc
!!$    nprocs = halo%nprocs_local
!!$
!!$    ! Zero output
!!$    y_local = 0.0d0
!!$
!!$    ! If no halo, simple local SpMV
!!$    if (halo%nhalo == 0) then
!!$       do i = 1, m_loc
!!$          do j = rowptr(i-1), rowptr(i)-1
!!$             y_local(i) = y_local(i) + nzval(j) * x_local( colind(j) - halo%fst_row + 1 )
!!$          end do
!!$       end do
!!$       return
!!$    end if
!!$
!!$    ! Build send buffer: extract x values for halo_cols, packed in send_order
!!$    allocate(sendbuf(halo%nhalo))
!!$    do p = 1, halo%nhalo
!!$       idx = halo%send_order(p)            ! idx indexes halo%halo_cols
!!$       sendbuf(p) = x_value_for_global( halo%halo_cols(idx), x_local, halo )
!!$    end do
!!$
!!$    ! Allocate recv buffer sized total recv entries
!!$    allocate(recvbuf( halo%rdispls(nprocs) + halo%recvcounts(nprocs) ))
!!$    ! Post Irecv from each src rank where recvcounts>0
!!$    nnei = size(halo%recv_from)
!!$    total_reqs = nnei + count(halo%sendcounts > 0)
!!$    allocate(reqs(total_reqs))
!!$    allocate(stats(MPI_STATUS_SIZE * total_reqs))
!!$    reqs_count = 0
!!$
!!$    ! Irecv
!!$    do k = 1, size(halo%recv_from)
!!$       owner = halo%recv_from(k)
!!$       rdis = halo%rdispls(owner+1)
!!$       call MPI_Irecv(recvbuf(rdis+1), halo%recvcounts(owner+1), MPI_DOUBLE_PRECISION, owner, &
!!$                      12345, use_comm, reqs(reqs_count+1), ierr)
!!$       reqs_count = reqs_count + 1
!!$    end do
!!$
!!$    ! Isend
!!$    do owner = 0, nprocs-1
!!$       if (halo%sendcounts(owner+1) > 0) then
!!$          sdis = halo%sdispls(owner+1)
!!$          call MPI_Isend(sendbuf(sdis+1), halo%sendcounts(owner+1), MPI_DOUBLE_PRECISION, owner, &
!!$                         12345, use_comm, reqs(reqs_count+1), ierr)
!!$          reqs_count = reqs_count + 1
!!$       end if
!!$    end do
!!$
!!$    ! Waitall
!!$    if (reqs_count > 0) then
!!$       call MPI_Waitall(reqs_count, reqs, stats, ierr)
!!$    end if
!!$
!!$    ! Now we have received remote x packed: recvbuf(rdis+1 : rdis+recvcounts)
!!$    ! Build a small map from halo global column -> value (simple linear search; halo sizes usually small)
!!$    ! We'll place recv values into an array 'halo_values' in same order as halo%halo_cols.
!!$    real(c_double), allocatable :: halo_values(:)
!!$    allocate(halo_values(halo%nhalo))
!!$    halo_values = 0.0d0
!!$
!!$    ! For each source owner, their contributions fill a contiguous range in recvbuf with length recvcounts(owner)
!!$    ! We need to know which global columns these correspond to:
!!$    ! When sending, other ranks will have used same scheme to compute send_order, so recv ordering
!!$    ! corresponds to THEIR send_order for destination=me. But because both sides used Alltoall to compute counts,
!!$    ! the ordering is consistent if both sides used the same packing (they do). However we didn't compute the remote ordering here.
!!$    !
!!$    ! To avoid complicated remote ordering logic, we will assume symmetric partitioning (most common):
!!$    !   - when rank R sends K entries to me, it sends them in the same relative order as my halo_cols that I own on R's side.
!!$    ! But in general to be robust across arbitrary patterns you would exchange index lists too. For simplicity here,
!!$    ! we'll assume both sides used the same deterministic packing algorithm (Alltoall of counts + our send_order), which is true.
!!$    !
!!$    ! For correct mapping, we need to reconstruct the mapping: for each owner O, the entries I receive from O correspond to
!!$    ! those halo_cols whose owner == O, and in the same order as they appear in halo%send_order for owner O.
!!$    !
!!$    ! Fill halo_values accordingly:
!!$
!!$    do owner = 0, nprocs-1
!!$       if (halo%sendcounts(owner+1) > 0) then
!!$          ! entries that this rank would have sent to 'owner' are at send_order( sdis+1 : sdis+sendcounts )
!!$          sdis = halo%sdispls(owner+1)
!!$          do p = 1, halo%sendcounts(owner+1)
!!$             pos = halo%send_order(sdis + p)
!!$             ! find position of this pos in halo%owners == ?  (we know owners(pos) == owner)
!!$             ! the corresponding received value from owner appears in recvbuf at rdis + offset
!!$             ! offset is the index among recvbuf for incoming from 'owner' where this rank is the destination.
!!$             ! That offset equals the index within recv for owner where rank==halo%rank; But to avoid two-way mapping complexity,
!!$             ! use this approach: the recv ordering for data from owner is exactly the subsequence of halo%halo_cols whose owners==halo%rank,
!!$             ! which we can produce by scanning halo%halo_cols. Simpler: reconstruct mapping by scanning halo%halo_cols and matching owners.
!!$          end do
!!$       end if
!!$    end do
!!$
!!$    ! Simpler robust approach: iterate over halo%halo_cols in order 1..nhalo;
!!$    ! for each halo_col k, find its owner o = owners(k).
!!$    ! Keep a per-owner pointer 'rcur(o)' initially rdispls(o)+1; when owner==o, assign halo_values(k) = recvbuf(rcur(o)); rcur(o)=rcur(o)+1
!!$    integer, allocatable :: rcur(:)
!!$    allocate(rcur(nprocs))
!!$    do owner = 0, nprocs-1
!!$       rcur(owner+1) = halo%rdispls(owner+1) + 1
!!$    end do
!!$
!!$    do k = 1, halo%nhalo
!!$       owner = halo%owners(k)
!!$       if (halo%recvcounts(owner+1) > 0) then
!!$          halo_values(k) = recvbuf( rcur(owner+1) )
!!$          rcur(owner+1) = rcur(owner+1) + 1
!!$       else
!!$          ! If recvcount==0 then this halo col was actually local to us (rare due to rounding)
!!$          halo_values(k) = x_value_for_global( halo%halo_cols(k), x_local, halo )
!!$       end if
!!$    end do
!!$
!!$    ! Finally do local SpMV using halo_values when needed
!!$    do i = 1, m_loc
!!$       do j = rowptr(i-1), rowptr(i)-1
!!$          if (colind(j) >= halo%fst_row .and. colind(j) <= halo%last_row) then
!!$             y_local(i) = y_local(i) + nzval(j) * x_local( colind(j) - halo%fst_row + 1 )
!!$          else
!!$             ! find index into halo_cols
!!$             ! linear search; can be replaced with hash if halo large
!!$             do k = 1, halo%nhalo
!!$                if (halo%halo_cols(k) == colind(j)) then
!!$                   y_local(i) = y_local(i) + nzval(j) * halo_values(k)
!!$                   exit
!!$                end if
!!$             end do
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    ! cleanup temporaries
!!$    deallocate(sendbuf, recvbuf, reqs, stats, halo_values, rcur)
!!$
!!$  end subroutine dist_spmv

  
  !----------------------------------------------------------------------
  ! dist_spmv_free - free arrays inside halo
  !----------------------------------------------------------------------

!!$  subroutine dist_spmv_free(halo)
!!$    type(halo_t), intent(inout) :: halo
!!$    if (allocated(halo%halo_cols)) deallocate(halo%halo_cols)
!!$    if (allocated(halo%owners)) deallocate(halo%owners)
!!$    if (allocated(halo%sendcounts)) deallocate(halo%sendcounts)
!!$    if (allocated(halo%sdispls)) deallocate(halo%sdispls)
!!$    if (allocated(halo%recvcounts)) deallocate(halo%recvcounts)
!!$    if (allocated(halo%rdispls)) deallocate(halo%rdispls)
!!$    if (allocated(halo%send_order)) deallocate(halo%send_order)
!!$    if (allocated(halo%recv_from)) deallocate(halo%recv_from)
!!$  end subroutine dist_spmv_free

  !----------------------------------------------------------------------
  ! Helper: return x value for a global index using local x_local if owned
  !----------------------------------------------------------------------

!!$  elemental function x_value_for_global(gidx, x_local, halo) result(val)
!!$    integer(c_int), intent(in) :: gidx
!!$    real(c_double), intent(in) :: x_local(:)
!!$    type(halo_t), intent(in) :: halo
!!$    real(c_double) :: val
!!$    if (gidx >= halo%fst_row .and. gidx <= halo%last_row) then
!!$       val = x_local( gidx - halo%fst_row + 1 )
!!$    else
!!$       val = 0.0d0   ! placeholder; actual remote values filled from recvbuf
!!$    end if
!!$  end function x_value_for_global
!!$  
  !----------------------------------------------------------------------
  ! small utility: unique_sort_int 
  !----------------------------------------------------------------------
  subroutine unique_sort_int(arr, n)
    !! Sorts arr(1:n) in ascending order and removes duplicates in-place.
    !! On return, n = number of unique elements (<= original n).
    implicit none
    integer, intent(inout) :: n
    integer, intent(inout) :: arr(:)

    if (n <= 1) return
    call quicksort_int(arr, 1, n)
    
    ! Remove duplicates in-place
    call unique_inplace_int(arr, n)
  end subroutine unique_sort_int


  recursive subroutine quicksort_int(a, left, right)
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j, pivot, tmp

    if (left >= right) return
    pivot = a((left + right) / 2)
    i = left
    j = right
    do
       do while (a(i) < pivot)
          i = i + 1
       end do
       do while (a(j) > pivot)
          j = j - 1
       end do
       if (i <= j) then
          tmp = a(i)
          a(i) = a(j)
          a(j) = tmp
          i = i + 1
          j = j - 1
       end if
       if (i > j) exit
    end do

    if (left < j) call quicksort_int(a, left, j)
    if (i < right) call quicksort_int(a, i, right)
  end subroutine quicksort_int
  
  
  subroutine unique_inplace_int(a, n)
    !! Removes duplicates in sorted array a(1:n), in-place.
    !! Returns new n = number of unique values.
    implicit none
    integer, intent(inout) :: n
    integer, intent(inout) :: a(:)
    integer :: i, k
    
    if (n <= 1) return
    
    k = 1
    do i = 2, n
       if (a(i) /= a(k)) then
          k = k + 1
          a(k) = a(i)
       end if
    end do
    n = k
  end subroutine unique_inplace_int

  
end module dist_spmv_mod
