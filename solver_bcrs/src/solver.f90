module solver
    implicit none
    private :: spmatvec_block, norm_block, axpby_block, preconditioner_block, dot_block
    public :: block_conjugate_gradient, validate
contains
    subroutine spmatvec_block(nnode, nblock, val, col, ind, vec, res)
        implicit none
        integer, intent(in) :: nnode, nblock, col(:), ind(:)
        double precision, intent(in) :: val(:, :, :), vec(:, :, :)
        double precision, intent(inout) :: res(:, :, :)

        ! double precision :: tmp(nblock, 3)
        double precision :: tmp
        integer :: inode, idim, jdim, iblock, j

        ! do inode = 1, nnode
        !     do idim = 1, 3
        !         do iblock = 1, nblock
        !             tmp(iblock, idim) = 0d0
        !         end do
        !     end do
        !     do j = ind(inode), ind(inode + 1) - 1
        !         do iblock = 1, nblock
        !             do idim = 1, 3
        !                 do jdim = 1, 3
        !                     tmp(iblock, idim) = tmp(iblock, idim) + val(idim, jdim, j)*vec(iblock, jdim, col(j))
        !                 end do
        !             end do
        !         end do
        !     end do
        !     do iblock = 1, nblock
        !         do idim = 1, 3
        !             res(iblock, idim, inode) = tmp(iblock, idim)
        !         end do
        !     end do
        ! end do
        do inode = 1, nnode
            do idim = 1, 3
                do iblock = 1, nblock
                    tmp = 0d0
                    do j = ind(inode), ind(inode + 1) - 1
                        tmp = tmp &
                              + val(idim, 1, j)*vec(iblock, 1, col(j)) &
                              + val(idim, 2, j)*vec(iblock, 2, col(j)) &
                              + val(idim, 3, j)*vec(iblock, 3, col(j))
                    end do
                    res(iblock, idim, inode) = tmp
                end do
            end do
        end do
    end subroutine

    subroutine norm_block(nnode, nblock, vec, norm)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: vec(:, :, :)
        double precision, intent(inout) :: norm(:)

        integer :: iblock, idim, inode

        do iblock = 1, nblock
            norm(iblock) = 0d0
        end do
        do inode = 1, nnode
            do idim = 1, 3
                do iblock = 1, nblock
                    norm(iblock) = norm(iblock) + vec(iblock, idim, inode)*vec(iblock, idim, inode)
                end do
            end do
        end do
        do iblock = 1, nblock
            norm(iblock) = sqrt(norm(iblock))
        end do
    end

    subroutine axpby_block(nnode, nblock, xvec, yvec, zvec, a, b)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: xvec(:, :, :), yvec(:, :, :)
        double precision, intent(inout) :: zvec(:, :, :)
        double precision, intent(in) :: a(:), b(:)

        integer :: inode, idim, iblock

        do inode = 1, nnode
            do idim = 1, 3
                do iblock = 1, nblock
                    zvec(iblock, idim, inode) = a(iblock)*xvec(iblock, idim, inode) + b(iblock)*yvec(iblock, idim, inode)
                end do
            end do
        end do
    end

    subroutine preconditioner_block(nnode, nblock, rvec, zvec, diag_inv)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: rvec(:, :, :), diag_inv(:, :, :)
        double precision, intent(inout) :: zvec(:, :, :)

        integer :: inode, iblock, idim, jdim

        do inode = 1, nnode
            do iblock = 1, nblock
                do idim = 1, 3
                    zvec(iblock, idim, inode) = 0d0
                    do jdim = 1, 3
                        zvec(iblock, idim, inode) = zvec(iblock, idim, inode) + &
                                                    diag_inv(idim, jdim, inode)*rvec(iblock, jdim, inode)
                    end do
                end do
            end do
        end do
    end

    subroutine dot_block(nnode, nblock, xvec, yvec, dot)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: xvec(:, :, :), yvec(:, :, :)
        double precision, intent(inout) :: dot(:)

        integer :: inode, idim, iblock

        do iblock = 1, nblock
            dot(iblock) = 0d0
        end do
        do inode = 1, nnode
            do idim = 1, 3
                do iblock = 1, nblock
                    dot(iblock) = dot(iblock) + xvec(iblock, idim, inode)*yvec(iblock, idim, inode)
                end do
            end do
        end do
    end

    subroutine block_conjugate_gradient &
        (nnode, nnz, nblock, amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
        implicit none
        double precision, intent(in) :: amat_val(:, :, :), amat_diag_inv(:, :, :), &
            bvec(:, :, :)
        double precision, intent(inout) :: uvec(:, :, :)
        integer, intent(in) :: nnode, nnz, nblock, amat_col(:), amat_ind(:)

        double precision :: rvec(nblock, 3, nnode), &
            pvec(nblock, 3, nnode), &
            qvec(nblock, 3, nnode), &
            zvec(nblock, 3, nnode), &
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
        call norm_block(nnode, nblock, bvec, bnorm)

        ! q <- A u
        call spmatvec_block(nnode, nblock, amat_val, amat_col, amat_ind, uvec, qvec)

        ! r <- b - q
        call axpby_block(nnode, nblock, bvec, qvec, rvec, ones, -ones)

        ! calculate rnorm
        call norm_block(nnode, nblock, rvec, rnorm)

        print *, "iter = ", iter, "err = ", maxval(rnorm/bnorm)
        do while (maxval(rnorm/bnorm) > 10d0**(-8))
            ! z <- M^-1 r
            call preconditioner_block(nnode, nblock, rvec, zvec, amat_diag_inv)

            if (iter > 1) then
                ! beta <- (r, z) / rho
                call dot_block(nnode, nblock, rvec, zvec, beta)
                do iblock = 1, nblock
                    beta(iblock) = beta(iblock)/rho(iblock)
                end do
            end if

            ! p <- z + beta p
            call axpby_block(nnode, nblock, zvec, pvec, pvec, ones, beta)

            ! q <- A p
            call spmatvec_block(nnode, nblock, amat_val, amat_col, amat_ind, pvec, qvec)

            ! rho <- (r, z)
            call dot_block(nnode, nblock, rvec, zvec, rho)

            ! alpha <- rho / (p, q)
            call dot_block(nnode, nblock, pvec, qvec, alpha)
            do iblock = 1, nblock
                alpha(iblock) = rho(iblock)/alpha(iblock)
            end do

            ! r <- r - alpha q
            call axpby_block(nnode, nblock, rvec, qvec, rvec, ones, -alpha)

            ! u <- u + alpha p
            call axpby_block(nnode, nblock, uvec, pvec, uvec, ones, alpha)

            ! calculate rnorm
            call norm_block(nnode, nblock, rvec, rnorm)

            print *, "iter = ", iter, "err = ", maxval(rnorm/bnorm)
            iter = iter + 1
        end do
        print *, "maxiter = ", iter
    end subroutine

    subroutine validate(nnode, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        integer, intent(in) :: nnode, nblock, kglobal_col(:), kglobal_ind(:)
        double precision, intent(in) :: kglobal_val(:, :, :), kglobal_diag_inv(:, :, :), uglobal(:, :, :), fglobal(:, :, :)

        integer :: iblock, idim, inode
        double precision :: res(nblock, 3, nnode)

        call spmatvec_block(nnode, nblock, kglobal_val, kglobal_col, kglobal_ind, uglobal, res)

        do inode = 1, nnode
            do idim = 1, 3
                do iblock = 1, nblock
                    res(iblock, idim, inode) = (res(iblock, idim, inode) - fglobal(iblock, idim, inode))**2
                end do
            end do
        end do

        print *, "L2 norm of residual: ", sqrt(sum(res))

    end
end module
