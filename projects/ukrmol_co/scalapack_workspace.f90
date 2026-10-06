! Query the installed ILP64 ScaLAPACK library without allocating a Hamiltonian.
! Run with mpirun -np (rows*columns) workspace dimension rows columns.
program workspace
  implicit none
  integer :: n, rows, columns, context, rank, processes, row, column
  integer :: nr, nc, descriptor(9), info, lwork, liwork
  integer, allocatable :: iwork(:)
  integer :: local_bytes
  integer :: numroc
  real(kind=8) :: a(1), z(1), eig_work(1), orth_work(1)
  real(kind=8) :: total_bytes, maximum_bytes, invalid_queries
  real(kind=8), allocatable :: eigenvalues(:)
  character(len=32) :: argument
  call get_command_argument(1, argument)
  read(argument, *) n
  call get_command_argument(2, argument)
  read(argument, *) rows
  call get_command_argument(3, argument)
  read(argument, *) columns
  if (min(n, rows, columns) < 1) stop 1
  call blacs_pinfo(rank, processes)
  if (processes /= rows*columns) stop 1
  call blacs_get(-1, 0, context)
  call blacs_gridinit(context, 'R', rows, columns)
  call blacs_gridinfo(context, rows, columns, row, column)
  nr = numroc(n, 64, row, 0, rows)
  nc = numroc(n, 64, column, 0, columns)
  call descinit(descriptor, n, n, 64, 64, 0, 0, context, max(1,nr), info)
  if (info /= 0) stop 2
  liwork = 7*n + 8*columns + 2
  allocate(eigenvalues(n), iwork(liwork))
  lwork = -1
  call pdsyevd('V', 'L', n, a, 1, 1, descriptor, eigenvalues, &
       z, 1, 1, descriptor, eig_work, lwork, iwork, liwork, info)
  if (info /= 0) stop 3
  call pdormtr('L', 'L', 'N', n, n, a, 1, 1, descriptor, eigenvalues, &
       z, 1, 1, descriptor, orth_work, -1, info)
  if (info /= 0) stop 4
  write(*,*) 'raw workspace queries (rank,eigensolver,orthogonal):', rank, eig_work(1), orth_work(1)
  invalid_queries = 0
  if (eig_work(1) < real(2*nr*nc, kind=8) .or. orth_work(1) <= 0) invalid_queries = 1
  call dgsum2d(context, 'All', ' ', 1, 1, invalid_queries, 1, -1, -1)
  if (invalid_queries > 0) then
    if (rank == 0) write(*,*) 'Invalid workspace query; no usable array-memory estimate.'
    call blacs_gridexit(context)
    call blacs_exit(0)
    stop 5
  end if
  ! Match the pinned UKRmol diagonalizer's extra PDORMTR workspace check.
  lwork = ceiling(max(eig_work(1), orth_work(1) + 2*n))
  local_bytes = 8*(2*nr*nc + lwork + liwork + n)
  ! Double represents these byte counts exactly and avoids the installed
  ! BLACS integer-reduction interface's 32-bit overflow for large allocations.
  total_bytes = real(local_bytes, kind=8)
  maximum_bytes = real(local_bytes, kind=8)
  call dgsum2d(context, 'All', ' ', 1, 1, total_bytes, 1, -1, -1)
  call dgamx2d(context, 'All', ' ', 1, 1, maximum_bytes, 1, row, column, -1, -1, -1)
  write(*,*) 'rank, local rows/cols, lwork/liwork, bytes:', rank, nr, nc, lwork, liwork, local_bytes
  if (rank == 0) write(*,*) 'dimension, grid, total bytes, max rank bytes:', &
       n, rows, columns, int(total_bytes,kind=8), int(maximum_bytes,kind=8)
  call blacs_gridexit(context)
  call blacs_exit(0)
end program workspace
