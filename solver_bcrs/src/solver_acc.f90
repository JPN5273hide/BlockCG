module solver_acc
    use omp_lib
    implicit none
    private :: spmatvec_block, norm_block, axpby_block, preconditioner_block, dot_block
    public :: block_conjugate_gradient, validate
contains
    subroutine spmatvec_block(nnode, nblock, val, col, ind, vec, res)
        implicit none
        integer, intent(in) :: nnode, nblock, col(:), ind(:)
        double precision, intent(in) :: val(:, :, :), vec(:, :, :)
        double precision, intent(inout) :: res(:, :, :)

        double precision :: tmp1, tmp2, tmp3, vec1, vec2, vec3
        integer :: inode, idim, jdim, iblock, j

        !$acc parallel present(val, col, ind, vec, res)
        !$acc loop private(tmp1, tmp2, tmp3, vec1, vec2, vec3) &
        !$acc collapse(2)
        do inode = 1, nnode
            do iblock = 1, nblock
                tmp1 = 0d0
                tmp2 = 0d0
                tmp3 = 0d0
                !$acc loop seq
                do j = ind(inode), ind(inode + 1) - 1
                    vec1 = vec(iblock, 1, col(j))
                    vec2 = vec(iblock, 2, col(j))
                    vec3 = vec(iblock, 3, col(j))
                    tmp1 = tmp1 &
                           + val(1, 1, j)*vec1 + val(2, 1, j)*vec2 + val(3, 1, j)*vec3
                    tmp2 = tmp2 &
                           + val(1, 2, j)*vec1 + val(2, 2, j)*vec2 + val(3, 2, j)*vec3
                    tmp3 = tmp3 &
                           + val(1, 3, j)*vec1 + val(2, 3, j)*vec2 + val(3, 3, j)*vec3
                end do
                res(iblock, 1, inode) = tmp1
                res(iblock, 2, inode) = tmp2
                res(iblock, 3, inode) = tmp3
            end do
        end do
        !$acc end parallel
    end subroutine

    subroutine norm_block(nnode, nblock, vec, norm)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: vec(:, :, :)
        double precision, intent(inout) :: norm(:)

        integer :: iblock, idim, inode
        double precision :: &
            norm01, norm02, norm03, norm04, norm05, norm06, norm07, norm08, &
            norm09, norm10, norm11, norm12, norm13, norm14, norm15, norm16

        norm01 = 0d0; norm02 = 0d0; norm03 = 0d0; norm04 = 0d0
        norm05 = 0d0; norm06 = 0d0; norm07 = 0d0; norm08 = 0d0
        norm09 = 0d0; norm10 = 0d0; norm11 = 0d0; norm12 = 0d0
        norm13 = 0d0; norm14 = 0d0; norm15 = 0d0; norm16 = 0d0
        !$acc kernels present(vec)
        !$acc loop independent &
        !$acc reduction(+:norm01, norm02, norm03, norm04, norm05, norm06, norm07, norm08, &
        !$acc             norm09, norm10, norm11, norm12, norm13, norm14, norm15, norm16)
        do inode = 1, nnode
            norm01 = norm01 + vec(1, 1, inode)*vec(1, 1, inode) &
                     + vec(1, 2, inode)*vec(1, 2, inode) &
                     + vec(1, 3, inode)*vec(1, 3, inode)
            norm02 = norm02 + vec(2, 1, inode)*vec(2, 1, inode) &
                     + vec(2, 2, inode)*vec(2, 2, inode) &
                     + vec(2, 3, inode)*vec(2, 3, inode)
            norm03 = norm03 + vec(3, 1, inode)*vec(3, 1, inode) &
                     + vec(3, 2, inode)*vec(3, 2, inode) &
                     + vec(3, 3, inode)*vec(3, 3, inode)
            norm04 = norm04 + vec(4, 1, inode)*vec(4, 1, inode) &
                     + vec(4, 2, inode)*vec(4, 2, inode) &
                     + vec(4, 3, inode)*vec(4, 3, inode)
            norm05 = norm05 + vec(5, 1, inode)*vec(5, 1, inode) &
                     + vec(5, 2, inode)*vec(5, 2, inode) &
                     + vec(5, 3, inode)*vec(5, 3, inode)
            norm06 = norm06 + vec(6, 1, inode)*vec(6, 1, inode) &
                     + vec(6, 2, inode)*vec(6, 2, inode) &
                     + vec(6, 3, inode)*vec(6, 3, inode)
            norm07 = norm07 + vec(7, 1, inode)*vec(7, 1, inode) &
                     + vec(7, 2, inode)*vec(7, 2, inode) &
                     + vec(7, 3, inode)*vec(7, 3, inode)
            norm08 = norm08 + vec(8, 1, inode)*vec(8, 1, inode) &
                     + vec(8, 2, inode)*vec(8, 2, inode) &
                     + vec(8, 3, inode)*vec(8, 3, inode)
            norm09 = norm09 + vec(9, 1, inode)*vec(9, 1, inode) &
                     + vec(9, 2, inode)*vec(9, 2, inode) &
                     + vec(9, 3, inode)*vec(9, 3, inode)
            norm10 = norm10 + vec(10, 1, inode)*vec(10, 1, inode) &
                     + vec(10, 2, inode)*vec(10, 2, inode) &
                     + vec(10, 3, inode)*vec(10, 3, inode)
            norm11 = norm11 + vec(11, 1, inode)*vec(11, 1, inode) &
                     + vec(11, 2, inode)*vec(11, 2, inode) &
                     + vec(11, 3, inode)*vec(11, 3, inode)
            norm12 = norm12 + vec(12, 1, inode)*vec(12, 1, inode) &
                     + vec(12, 2, inode)*vec(12, 2, inode) &
                     + vec(12, 3, inode)*vec(12, 3, inode)
            norm13 = norm13 + vec(13, 1, inode)*vec(13, 1, inode) &
                     + vec(13, 2, inode)*vec(13, 2, inode) &
                     + vec(13, 3, inode)*vec(13, 3, inode)
            norm14 = norm14 + vec(14, 1, inode)*vec(14, 1, inode) &
                     + vec(14, 2, inode)*vec(14, 2, inode) &
                     + vec(14, 3, inode)*vec(14, 3, inode)
            norm15 = norm15 + vec(15, 1, inode)*vec(15, 1, inode) &
                     + vec(15, 2, inode)*vec(15, 2, inode) &
                     + vec(15, 3, inode)*vec(15, 3, inode)
            norm16 = norm16 + vec(16, 1, inode)*vec(16, 1, inode) &
                     + vec(16, 2, inode)*vec(16, 2, inode) &
                     + vec(16, 3, inode)*vec(16, 3, inode)
        end do
        !$acc end kernels
        norm(1) = sqrt(norm01); norm(2) = sqrt(norm02); norm(3) = sqrt(norm03); norm(4) = sqrt(norm04)
        norm(5) = sqrt(norm05); norm(6) = sqrt(norm06); norm(7) = sqrt(norm07); norm(8) = sqrt(norm08)
        norm(9) = sqrt(norm09); norm(10) = sqrt(norm10); norm(11) = sqrt(norm11); norm(12) = sqrt(norm12)
        norm(13) = sqrt(norm13); norm(14) = sqrt(norm14); norm(15) = sqrt(norm15); norm(16) = sqrt(norm16)
    end

    subroutine axpby_block(nnode, nblock, xvec, yvec, zvec, a, b)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: xvec(:, :, :), yvec(:, :, :)
        double precision, intent(inout) :: zvec(:, :, :)
        double precision, intent(in) :: a(:), b(:)

        integer :: inode, idim, iblock

        !$acc kernels present(xvec, yvec, zvec)
        !$acc loop independent collapse(2)
        do inode = 1, nnode
            do iblock = 1, nblock
                zvec(iblock, 1, inode) = a(iblock)*xvec(iblock, 1, inode) + b(iblock)*yvec(iblock, 1, inode)
                zvec(iblock, 2, inode) = a(iblock)*xvec(iblock, 2, inode) + b(iblock)*yvec(iblock, 2, inode)
                zvec(iblock, 3, inode) = a(iblock)*xvec(iblock, 3, inode) + b(iblock)*yvec(iblock, 3, inode)
            end do
        end do
        !$acc end kernels
    end

    subroutine preconditioner_block(nnode, nblock, rvec, zvec, diag_inv)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: rvec(:, :, :), diag_inv(:, :, :)
        double precision, intent(inout) :: zvec(:, :, :)

        integer :: inode, iblock

        !$acc kernels present(rvec, diag_inv, zvec)
        !$acc loop independent collapse(2)
        do inode = 1, nnode
            do iblock = 1, nblock
                zvec(iblock, 1, inode) = diag_inv(1, 1, inode)*rvec(iblock, 1, inode) + &
                                         diag_inv(2, 1, inode)*rvec(iblock, 2, inode) + &
                                         diag_inv(3, 1, inode)*rvec(iblock, 3, inode)
                zvec(iblock, 2, inode) = diag_inv(1, 2, inode)*rvec(iblock, 1, inode) + &
                                         diag_inv(2, 2, inode)*rvec(iblock, 2, inode) + &
                                         diag_inv(3, 2, inode)*rvec(iblock, 3, inode)
                zvec(iblock, 3, inode) = diag_inv(1, 3, inode)*rvec(iblock, 1, inode) + &
                                         diag_inv(2, 3, inode)*rvec(iblock, 2, inode) + &
                                         diag_inv(3, 3, inode)*rvec(iblock, 3, inode)
            end do
        end do
        !$acc end kernels
    end

    subroutine dot_block(nnode, nblock, xvec, yvec, dot)
        implicit none
        integer, intent(in) :: nnode, nblock
        double precision, intent(in) :: xvec(:, :, :), yvec(:, :, :)
        double precision, intent(inout) :: dot(:)

        integer :: inode, iblock
        double precision :: &
            dot01, dot02, dot03, dot04, dot05, dot06, dot07, dot08, &
            dot09, dot10, dot11, dot12, dot13, dot14, dot15, dot16

        dot01 = 0d0; dot02 = 0d0; dot03 = 0d0; dot04 = 0d0
        dot05 = 0d0; dot06 = 0d0; dot07 = 0d0; dot08 = 0d0
        dot09 = 0d0; dot10 = 0d0; dot11 = 0d0; dot12 = 0d0
        dot13 = 0d0; dot14 = 0d0; dot15 = 0d0; dot16 = 0d0

        !$acc kernels present(xvec, yvec)
        !$acc loop independent &
        !$acc reduction(+:dot01, dot02, dot03, dot04, &
        !$acc             dot05, dot06, dot07, dot08, &
        !$acc             dot09, dot10, dot11, dot12, &
        !$acc             dot13, dot14, dot15, dot16)
        do inode = 1, nnode
            dot01 = dot01 + xvec(1, 1, inode)*yvec(1, 1, inode) + &
                    xvec(1, 2, inode)*yvec(1, 2, inode) + &
                    xvec(1, 3, inode)*yvec(1, 3, inode)
            dot02 = dot02 + xvec(2, 1, inode)*yvec(2, 1, inode) + &
                    xvec(2, 2, inode)*yvec(2, 2, inode) + &
                    xvec(2, 3, inode)*yvec(2, 3, inode)
            dot03 = dot03 + xvec(3, 1, inode)*yvec(3, 1, inode) + &
                    xvec(3, 2, inode)*yvec(3, 2, inode) + &
                    xvec(3, 3, inode)*yvec(3, 3, inode)
            dot04 = dot04 + xvec(4, 1, inode)*yvec(4, 1, inode) + &
                    xvec(4, 2, inode)*yvec(4, 2, inode) + &
                    xvec(4, 3, inode)*yvec(4, 3, inode)
            dot05 = dot05 + xvec(5, 1, inode)*yvec(5, 1, inode) + &
                    xvec(5, 2, inode)*yvec(5, 2, inode) + &
                    xvec(5, 3, inode)*yvec(5, 3, inode)
            dot06 = dot06 + xvec(6, 1, inode)*yvec(6, 1, inode) + &
                    xvec(6, 2, inode)*yvec(6, 2, inode) + &
                    xvec(6, 3, inode)*yvec(6, 3, inode)
            dot07 = dot07 + xvec(7, 1, inode)*yvec(7, 1, inode) + &
                    xvec(7, 2, inode)*yvec(7, 2, inode) + &
                    xvec(7, 3, inode)*yvec(7, 3, inode)
            dot08 = dot08 + xvec(8, 1, inode)*yvec(8, 1, inode) + &
                    xvec(8, 2, inode)*yvec(8, 2, inode) + &
                    xvec(8, 3, inode)*yvec(8, 3, inode)
            dot09 = dot09 + xvec(9, 1, inode)*yvec(9, 1, inode) + &
                    xvec(9, 2, inode)*yvec(9, 2, inode) + &
                    xvec(9, 3, inode)*yvec(9, 3, inode)
            dot10 = dot10 + xvec(10, 1, inode)*yvec(10, 1, inode) + &
                    xvec(10, 2, inode)*yvec(10, 2, inode) + &
                    xvec(10, 3, inode)*yvec(10, 3, inode)
            dot11 = dot11 + xvec(11, 1, inode)*yvec(11, 1, inode) + &
                    xvec(11, 2, inode)*yvec(11, 2, inode) + &
                    xvec(11, 3, inode)*yvec(11, 3, inode)
            dot12 = dot12 + xvec(12, 1, inode)*yvec(12, 1, inode) + &
                    xvec(12, 2, inode)*yvec(12, 2, inode) + &
                    xvec(12, 3, inode)*yvec(12, 3, inode)
            dot13 = dot13 + xvec(13, 1, inode)*yvec(13, 1, inode) + &
                    xvec(13, 2, inode)*yvec(13, 2, inode) + &
                    xvec(13, 3, inode)*yvec(13, 3, inode)
            dot14 = dot14 + xvec(14, 1, inode)*yvec(14, 1, inode) + &
                    xvec(14, 2, inode)*yvec(14, 2, inode) + &
                    xvec(14, 3, inode)*yvec(14, 3, inode)
            dot15 = dot15 + xvec(15, 1, inode)*yvec(15, 1, inode) + &
                    xvec(15, 2, inode)*yvec(15, 2, inode) + &
                    xvec(15, 3, inode)*yvec(15, 3, inode)
            dot16 = dot16 + xvec(16, 1, inode)*yvec(16, 1, inode) + &
                    xvec(16, 2, inode)*yvec(16, 2, inode) + &
                    xvec(16, 3, inode)*yvec(16, 3, inode)
        end do
        !$acc end kernels
        dot(1) = dot01; dot(2) = dot02; dot(3) = dot03; dot(4) = dot04
        dot(5) = dot05; dot(6) = dot06; dot(7) = dot07; dot(8) = dot08
        dot(9) = dot09; dot(10) = dot10; dot(11) = dot11; dot(12) = dot12
        dot(13) = dot13; dot(14) = dot14; dot(15) = dot15; dot(16) = dot16
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
        integer :: iblock, iter, maxiter
        double precision :: st_time, en_time

        maxiter = 50

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
        st_time = omp_get_wtime()

        !$acc data copyin(amat_val, amat_col, amat_ind, amat_diag_inv, bvec, uvec) &
        !$acc create(rvec, pvec, qvec, zvec) &
        !$acc copyout(uvec)
        ! calculate bnorm
        call norm_block(nnode, nblock, bvec, bnorm)

        ! q <- A u
        call spmatvec_block(nnode, nblock, amat_val, amat_col, amat_ind, uvec, qvec)

        ! r <- b - q
        call axpby_block(nnode, nblock, bvec, qvec, rvec, ones, -ones)

        ! calculate rnorm
        call norm_block(nnode, nblock, rvec, rnorm)

        print *, "iter = ", iter, "err = ", maxval(rnorm/bnorm)
        do while (maxval(rnorm/bnorm) > 10d0**(-8) .and. iter <= maxiter)
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
        !$acc end data
        en_time = omp_get_wtime()
        print *, "maxiter = ", iter
        print *, "elapsed time = ", en_time - st_time
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
