program main
    use solver
    implicit none

    integer :: ndof, nnz, nblock, idof, iblock
    double precision, allocatable :: kglobal_val(:), kglobal_diag_inv(:), &
    uglobal(:, :), fglobal(:, :)
    integer, allocatable :: kglobal_col(:), kglobal_ind(:)

    print *, "Reading data..."
    call get_ndof_nnz_nblock(ndof, nnz, nblock)
    print *, ndof, nnz, nblock
    allocate(kglobal_val(nnz), kglobal_col(nnz), kglobal_ind(ndof + 1), kglobal_diag_inv(ndof))
    allocate(uglobal(nblock, ndof), fglobal(nblock, ndof))
    call read_data(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

    print *, "Solving..."
    call block_conjugate_gradient(ndof, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

    print *, "Output data..."
    open(10, file='../data/sol.dat', status='replace')
    do idof = 1, ndof
        do iblock = 1, nblock
            write (10, "(e20.10)", advance="no") uglobal(iblock, idof)
        end do
        write (10, *)
    end do
    close(10)

    print *, "Validation..."
    call validate(ndof, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

contains

    subroutine get_ndof_nnz_nblock(ndof, nnz, nblock)
        implicit none

        integer, intent(out) :: ndof, nnz, nblock

        open(10, file='../data/shape.dat', status='old')
        read (10, *) ndof, nnz, nblock
        close(10)
    end subroutine

    subroutine read_data(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        implicit none

        integer, intent(in) :: ndof, nnz, nblock
        integer, intent(inout) :: kglobal_col(:), kglobal_ind(:)
        double precision, intent(inout) :: kglobal_val(:), kglobal_diag_inv(:), uglobal(:, :), fglobal(:, :)

        integer :: inz, idof, iblock

        open(10, file='../data/kglobal_val.dat', status='old')
        read (10, *) kglobal_val
        close(10)

        open(10, file='../data/kglobal_col.dat', status='old')
        read (10, *) kglobal_col
        close(10)

        open(10, file='../data/kglobal_ind.dat', status='old')
        read (10, *) kglobal_ind
        close(10)

        open(10, file='../data/kglobal_diag_inv.dat', status='old')
        read (10, *) kglobal_diag_inv
        close(10)

        ! open(10, file='../data/uinit.dat', status='old')
        ! read (10, *) uglobal
        ! close(10)

        open(10, file='../data/fglobal.dat', status='old')
        read (10, *) fglobal
        close(10)
        print *, sum(fglobal)
    end subroutine
end program