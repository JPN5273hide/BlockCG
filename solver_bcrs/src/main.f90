program main
    use openacc
    use solver_acc
    implicit none

    integer :: nnode, nnz, nblock, inode, iblock
    double precision, allocatable :: &
        kglobal_val(:, :, :), kglobal_diag_inv(:, :, :), &
        uglobal(:, :, :), fglobal(:, :, :)
    integer, allocatable :: kglobal_col(:), kglobal_ind(:)
    character(len=256) :: filepath
    character(len=*), parameter :: data_dir = "../bcrs_data/"

    print *, "Reading data..."
    call get_shape(data_dir, nnode, nnz, nblock)
    print *, nnode, nnz, nblock
    allocate (kglobal_val(3, 3, nnz), kglobal_col(nnz), kglobal_ind(nnode + 1), kglobal_diag_inv(3, 3, nnode))
    allocate (uglobal(nblock, 3, nnode), fglobal(nblock, 3, nnode))
    call read_data(data_dir, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

    print *, "Solving..."
    call block_conjugate_gradient(nnode, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

    ! print *, "Output data..."
    ! !  open (10, file='../data__/sol.dat', status='replace')
    ! !  do idof = 1, ndof
    ! !      do iblock = 1, nblock
    ! !          write (10, "(e20.10)", advance="no") uglobal(iblock, idof)
    ! !      end do
    ! !      write (10, *)
    ! !  end do

    ! filepath = trim(data_dir)//'sol.bin'
    ! open (10, file=filepath, status='replace', form='unformatted', access='stream')
    ! write (10) uglobal
    ! close (10)

    ! print *, "Validation..."
    ! call validate(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

contains

    subroutine get_shape(data_dir, nnode, nnz, nblock)
        implicit none

        integer, intent(out) :: nnode, nnz, nblock
        character(len=*), intent(in) :: data_dir

        character(len=256) :: filepath
        filepath = trim(data_dir)//'shape.dat'
        open (10, file=filepath, status='old')
        read (10, *) nnode, nnz, nblock
        close (10)
    end subroutine

    subroutine read_data(data_dir, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        implicit none
        character(len=*), intent(in) :: data_dir
        integer, intent(inout) :: kglobal_col(:), kglobal_ind(:)
        double precision, intent(inout) :: &
            kglobal_val(:, :, :), kglobal_diag_inv(:, :, :), uglobal(:, :, :), fglobal(:, :, :)

        character(len=256) :: filepath

        filepath = trim(data_dir)//'kglobal_val.bin'
        open (10, file=filepath, status='old', form='unformatted', access='stream')
        read (10) kglobal_val
        close (10)

        filepath = trim(data_dir)//'kglobal_col.bin'
        open (10, file=filepath, status='old', form='unformatted', access='stream')
        read (10) kglobal_col
        close (10)

        filepath = trim(data_dir)//'kglobal_ind.bin'
        open (10, file=filepath, status='old', form='unformatted', access='stream')
        read (10) kglobal_ind
        close (10)

        filepath = trim(data_dir)//'kglobal_diag_inv.bin'
        open (10, file=filepath, status='old', form='unformatted', access='stream')
        read (10) kglobal_diag_inv
        close (10)

        filepath = trim(data_dir)//'fglobal.bin'
        open (10, file=filepath, status='old', form='unformatted', access='stream')
        read (10) fglobal
        close (10)
        print *, sum(fglobal)
    end subroutine
end program
