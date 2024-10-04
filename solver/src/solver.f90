module solver
    implicit none
    private :: spmatvec_block, norm_block, axpby_block, preconditioner_block, dot_block
    public :: block_conjugate_gradient, validate
contains
    subroutine spmatvec_block(ndof, nblock, val, col, ind, vec, res)
        implicit none
        integer, intent(in) :: ndof, nblock, col(:), ind(:)
        double precision, intent(in) :: val(:), vec(:, :)
        double precision, intent(inout) :: res(:, :)

        double precision :: tmp(nblock)
        integer :: idof, iblock, j

        do idof = 1, ndof
            do iblock = 1, nblock
                tmp(iblock) = 0d0
            end do
            do j = ind(idof), ind(idof + 1) - 1
                do iblock = 1, nblock
                    tmp(iblock) = tmp(iblock) + val(j)*vec(iblock, col(j))
                end do
            end do
            do iblock = 1, nblock
                res(iblock, idof) = tmp(iblock)
            end do
        end do
    end subroutine

    subroutine norm_block(ndof, nblock, vec, norm)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: vec(:, :)
        double precision, intent(inout) :: norm(:)

        integer :: idof, iblock

        do iblock = 1, nblock
            norm(iblock) = 0d0
        end do
        do idof = 1, ndof
            do iblock = 1, nblock
                norm(iblock) = norm(iblock) + vec(iblock, idof)*vec(iblock, idof)
            end do
        end do
        do iblock = 1, nblock
            norm(iblock) = sqrt(norm(iblock))
        end do
    end

    subroutine axpby_block(ndof, nblock, xvec, yvec, zvec, a, b)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: xvec(:, :), yvec(:, :)
        double precision, intent(inout) :: zvec(:, :)
        double precision, intent(in) :: a(:), b(:)

        integer :: idof, iblock

        do idof = 1, ndof
            do iblock = 1, nblock
                zvec(iblock, idof) = a(iblock)*xvec(iblock, idof) + b(iblock)*yvec(iblock, idof)
            end do
        end do
    end

    subroutine preconditioner_block(ndof, nblock, rvec, zvec, diag_inv)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: rvec(:, :), diag_inv(:)
        double precision, intent(inout) :: zvec(:, :)

        integer :: idof, iblock

        do idof = 1, ndof
            do iblock = 1, nblock
                zvec(iblock, idof) = rvec(iblock, idof)*diag_inv(idof)
            end do
        end do
    end

    subroutine dot_block(ndof, nblock, xvec, yvec, dot)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: xvec(:, :), yvec(:, :)
        double precision, intent(inout) :: dot(:)

        integer :: idof, iblock

        do iblock = 1, nblock
            dot(iblock) = 0d0
        end do
        do idof = 1, ndof
            do iblock = 1, nblock
                dot(iblock) = dot(iblock) + xvec(iblock, idof)*yvec(iblock, idof)
            end do
        end do
    end

    subroutine block_conjugate_gradient &
        (ndof, nblock, amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
        implicit none
        double precision, intent(in) :: amat_val(:), amat_diag_inv(:), &
            bvec(:, :)
        double precision, intent(inout) :: uvec(:, :)
        integer, intent(in) :: ndof, nblock, amat_col(:), amat_ind(:)

        double precision :: rvec(nblock, ndof), &
            pvec(nblock, ndof), &
            qvec(nblock, ndof), &
            zvec(nblock, ndof), &
            bnorm(nblock), rnorm(nblock), &
            alpha(nblock), beta(nblock), rho(nblock), &
            ones(nblock)
        integer :: iblock, iter

        ! initialize
        do iblock = 1, nblock
            alpha(iblock) = 0d0
            beta(iblock) = 0d0
            rho(iblock) = 0d0
            ones(iblock) = 1d0
        end do
        iter = 1
        uvec = 0d0

        print *, "start solving..."

        ! calculate bnorm
        call norm_block(ndof, nblock, bvec, bnorm)

        ! q <- A u
        call spmatvec_block(ndof, nblock, amat_val, amat_col, amat_ind, uvec, qvec)

        ! r <- b - q
        call axpby_block(ndof, nblock, bvec, qvec, rvec, ones, -ones)

        ! calculate rnorm
        call norm_block(ndof, nblock, rvec, rnorm)

        ! print *, "iter = ", iter, "err = ", maxval(rnorm / bnorm)
        do while (maxval(rnorm/bnorm) > 10d0**(-8))
            ! z <- M^-1 r
            call preconditioner_block(ndof, nblock, rvec, zvec, amat_diag_inv)

            if (iter > 1) then
                ! beta <- (r, z) / rho
                call dot_block(ndof, nblock, rvec, zvec, beta)
                do iblock = 1, nblock
                    beta(iblock) = beta(iblock)/rho(iblock)
                end do
            end if

            ! p <- z + beta p
            call axpby_block(ndof, nblock, zvec, pvec, pvec, ones, beta)

            ! q <- A p
            call spmatvec_block(ndof, nblock, amat_val, amat_col, amat_ind, pvec, qvec)

            ! rho <- (r, z)
            call dot_block(ndof, nblock, rvec, zvec, rho)

            ! alpha <- rho / (p, q)
            call dot_block(ndof, nblock, pvec, qvec, alpha)
            do iblock = 1, nblock
                alpha(iblock) = rho(iblock)/alpha(iblock)
            end do

            ! r <- r - alpha q
            call axpby_block(ndof, nblock, rvec, qvec, rvec, ones, -alpha)

            ! u <- u + alpha p
            call axpby_block(ndof, nblock, uvec, pvec, uvec, ones, alpha)

            ! calculate rnorm
            call norm_block(ndof, nblock, rvec, rnorm)

            ! print *, "iter = ", iter, "err = ", maxval(rnorm / bnorm)
            iter = iter + 1
        end do
        print *, "maxiter = ", iter
    end subroutine

    subroutine validate(ndof, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        integer, intent(in) :: ndof, nblock, kglobal_col(:), kglobal_ind(:)
        double precision, intent(in) :: kglobal_val(:), kglobal_diag_inv(:), uglobal(:, :), fglobal(:, :)

        integer :: iblock, idof
        double precision :: res(nblock, ndof)

        call spmatvec_block(ndof, nblock, kglobal_val, kglobal_col, kglobal_ind, uglobal, res)

        do idof = 1, ndof
            do iblock = 1, nblock
                res(iblock, idof) = (res(iblock, idof) - fglobal(iblock, idof))**2
            end do
        end do

        print *, "L2 norm of residual: ", sqrt(sum(res))

    end
end module
