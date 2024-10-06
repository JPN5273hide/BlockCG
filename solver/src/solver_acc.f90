module solver_acc
    implicit none
    private :: spmatvec_block, norm_block, xpby_block, preconditioner_block, dot_block
    public :: block_conjugate_gradient, validate
contains
    subroutine spmatvec_block(ndof, nnz, nblock, val, col, ind, vec, res, tmp)
        implicit none
        integer, intent(in) :: ndof, nnz, nblock, col(:), ind(:)
        double precision, intent(in) :: val(:), vec(:, :)
        double precision, intent(inout) :: res(:, :), tmp

        integer :: idof, iblock, j

        !$acc kernels present(res, vec, val, col, ind)
        !$acc loop independent private(tmp) collapse(2)
        do idof = 1, ndof
            do iblock = 1, nblock
                tmp = 0.0
                do j = ind(idof), ind(idof + 1) - 1
                    tmp = tmp + val(j)*vec(iblock, col(j))
                end do
                res(iblock, idof) = tmp
            end do
        end do
        !$acc end kernels
    end subroutine

    subroutine norm_block(ndof, nblock, vec, &
                          norm01, norm02, norm03, norm04, &
                          norm05, norm06, norm07, norm08, &
                          norm09, norm10, norm11, norm12, &
                          norm13, norm14, norm15, norm16)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: vec(:, :)
        double precision, intent(inout) :: &
            norm01, norm02, norm03, norm04, &
            norm05, norm06, norm07, norm08, &
            norm09, norm10, norm11, norm12, &
            norm13, norm14, norm15, norm16

        integer :: idof, iblock

        norm01 = 0d0
        norm02 = 0d0
        norm03 = 0d0
        norm04 = 0d0
        norm05 = 0d0
        norm06 = 0d0
        norm07 = 0d0
        norm08 = 0d0
        norm09 = 0d0
        norm10 = 0d0
        norm11 = 0d0
        norm12 = 0d0
        norm13 = 0d0
        norm14 = 0d0
        norm15 = 0d0
        norm16 = 0d0
        !$acc kernels present(vec)
        !$acc loop independent &
        !$acc reduction(+:norm01, norm02, norm03, norm04, &
        !$acc             norm05, norm06, norm07, norm08, &
        !$acc             norm09, norm10, norm11, norm12, &
        !$acc             norm13, norm14, norm15, norm16)
        do idof = 1, ndof
            norm01 = norm01 + vec(1, idof)*vec(1, idof)
            norm02 = norm02 + vec(2, idof)*vec(2, idof)
            norm03 = norm03 + vec(3, idof)*vec(3, idof)
            norm04 = norm04 + vec(4, idof)*vec(4, idof)
            norm05 = norm05 + vec(5, idof)*vec(5, idof)
            norm06 = norm06 + vec(6, idof)*vec(6, idof)
            norm07 = norm07 + vec(7, idof)*vec(7, idof)
            norm08 = norm08 + vec(8, idof)*vec(8, idof)
            norm09 = norm09 + vec(9, idof)*vec(9, idof)
            norm10 = norm10 + vec(10, idof)*vec(10, idof)
            norm11 = norm11 + vec(11, idof)*vec(11, idof)
            norm12 = norm12 + vec(12, idof)*vec(12, idof)
            norm13 = norm13 + vec(13, idof)*vec(13, idof)
            norm14 = norm14 + vec(14, idof)*vec(14, idof)
            norm15 = norm15 + vec(15, idof)*vec(15, idof)
            norm16 = norm16 + vec(16, idof)*vec(16, idof)
        end do
        !$acc end kernels
        norm01 = sqrt(norm01)
        norm02 = sqrt(norm02)
        norm03 = sqrt(norm03)
        norm04 = sqrt(norm04)
        norm05 = sqrt(norm05)
        norm06 = sqrt(norm06)
        norm07 = sqrt(norm07)
        norm08 = sqrt(norm08)
        norm09 = sqrt(norm09)
        norm10 = sqrt(norm10)
        norm11 = sqrt(norm11)
        norm12 = sqrt(norm12)
        norm13 = sqrt(norm13)
        norm14 = sqrt(norm14)
        norm15 = sqrt(norm15)
        norm16 = sqrt(norm16)
    end subroutine

    subroutine xpby_block(ndof, nblock, xvec, yvec, zvec, &
                          b01, b02, b03, b04, b05, b06, b07, b08, &
                          b09, b10, b11, b12, b13, b14, b15, b16)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: xvec(:, :), yvec(:, :)
        double precision, intent(inout) :: zvec(:, :)
        double precision, intent(in) :: &
            b01, b02, b03, b04, b05, b06, b07, b08, &
            b09, b10, b11, b12, b13, b14, b15, b16

        integer :: idof, iblock

        !$acc kernels present(xvec, yvec, zvec)
        !$acc loop independent
        do idof = 1, ndof
            zvec(1, idof) = xvec(1, idof) + b01*yvec(1, idof)
            zvec(2, idof) = xvec(2, idof) + b02*yvec(2, idof)
            zvec(3, idof) = xvec(3, idof) + b03*yvec(3, idof)
            zvec(4, idof) = xvec(4, idof) + b04*yvec(4, idof)
            zvec(5, idof) = xvec(5, idof) + b05*yvec(5, idof)
            zvec(6, idof) = xvec(6, idof) + b06*yvec(6, idof)
            zvec(7, idof) = xvec(7, idof) + b07*yvec(7, idof)
            zvec(8, idof) = xvec(8, idof) + b08*yvec(8, idof)
            zvec(9, idof) = xvec(9, idof) + b09*yvec(9, idof)
            zvec(10, idof) = xvec(10, idof) + b10*yvec(10, idof)
            zvec(11, idof) = xvec(11, idof) + b11*yvec(11, idof)
            zvec(12, idof) = xvec(12, idof) + b12*yvec(12, idof)
            zvec(13, idof) = xvec(13, idof) + b13*yvec(13, idof)
            zvec(14, idof) = xvec(14, idof) + b14*yvec(14, idof)
            zvec(15, idof) = xvec(15, idof) + b15*yvec(15, idof)
            zvec(16, idof) = xvec(16, idof) + b16*yvec(16, idof)
        end do
        !$acc end kernels
    end

    subroutine preconditioner_block(ndof, nblock, rvec, zvec, diag_inv)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: rvec(:, :), diag_inv(:)
        double precision, intent(inout) :: zvec(:, :)

        integer :: idof, iblock

        !$acc kernels present(rvec, zvec, diag_inv)
        !$acc loop independent
        do idof = 1, ndof
            zvec(1, idof) = rvec(1, idof)*diag_inv(idof)
            zvec(2, idof) = rvec(2, idof)*diag_inv(idof)
            zvec(3, idof) = rvec(3, idof)*diag_inv(idof)
            zvec(4, idof) = rvec(4, idof)*diag_inv(idof)
            zvec(5, idof) = rvec(5, idof)*diag_inv(idof)
            zvec(6, idof) = rvec(6, idof)*diag_inv(idof)
            zvec(7, idof) = rvec(7, idof)*diag_inv(idof)
            zvec(8, idof) = rvec(8, idof)*diag_inv(idof)
            zvec(9, idof) = rvec(9, idof)*diag_inv(idof)
            zvec(10, idof) = rvec(10, idof)*diag_inv(idof)
            zvec(11, idof) = rvec(11, idof)*diag_inv(idof)
            zvec(12, idof) = rvec(12, idof)*diag_inv(idof)
            zvec(13, idof) = rvec(13, idof)*diag_inv(idof)
            zvec(14, idof) = rvec(14, idof)*diag_inv(idof)
            zvec(15, idof) = rvec(15, idof)*diag_inv(idof)
            zvec(16, idof) = rvec(16, idof)*diag_inv(idof)
        end do
        !$acc end kernels
    end

    subroutine dot_block(ndof, nblock, xvec, yvec, &
                         dot01, dot02, dot03, dot04, &
                         dot05, dot06, dot07, dot08, &
                         dot09, dot10, dot11, dot12, &
                         dot13, dot14, dot15, dot16)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: xvec(:, :), yvec(:, :)
        double precision, intent(inout) :: &
            dot01, dot02, dot03, dot04, &
            dot05, dot06, dot07, dot08, &
            dot09, dot10, dot11, dot12, &
            dot13, dot14, dot15, dot16

        integer :: idof, iblock

        dot01 = 0d0
        dot02 = 0d0
        dot03 = 0d0
        dot04 = 0d0
        dot05 = 0d0
        dot06 = 0d0
        dot07 = 0d0
        dot08 = 0d0
        dot09 = 0d0
        dot10 = 0d0
        dot11 = 0d0
        dot12 = 0d0
        dot13 = 0d0
        dot14 = 0d0
        dot15 = 0d0
        dot16 = 0d0
        !$acc kernels present(xvec, yvec)
        !$acc loop independent &
        !$acc reduction(+:dot01, dot02, dot03, dot04, &
        !$acc             dot05, dot06, dot07, dot08, &
        !$acc             dot09, dot10, dot11, dot12, &
        !$acc             dot13, dot14, dot15, dot16)
        do idof = 1, ndof
            dot01 = dot01 + xvec(1, idof)*yvec(1, idof)
            dot02 = dot02 + xvec(2, idof)*yvec(2, idof)
            dot03 = dot03 + xvec(3, idof)*yvec(3, idof)
            dot04 = dot04 + xvec(4, idof)*yvec(4, idof)
            dot05 = dot05 + xvec(5, idof)*yvec(5, idof)
            dot06 = dot06 + xvec(6, idof)*yvec(6, idof)
            dot07 = dot07 + xvec(7, idof)*yvec(7, idof)
            dot08 = dot08 + xvec(8, idof)*yvec(8, idof)
            dot09 = dot09 + xvec(9, idof)*yvec(9, idof)
            dot10 = dot10 + xvec(10, idof)*yvec(10, idof)
            dot11 = dot11 + xvec(11, idof)*yvec(11, idof)
            dot12 = dot12 + xvec(12, idof)*yvec(12, idof)
            dot13 = dot13 + xvec(13, idof)*yvec(13, idof)
            dot14 = dot14 + xvec(14, idof)*yvec(14, idof)
            dot15 = dot15 + xvec(15, idof)*yvec(15, idof)
            dot16 = dot16 + xvec(16, idof)*yvec(16, idof)
        end do
        !$acc end kernels
    end subroutine

    subroutine div1_block(a01, a02, a03, a04, &
                          a05, a06, a07, a08, &
                          a09, a10, a11, a12, &
                          a13, a14, a15, a16, &
                          b01, b02, b03, b04, &
                          b05, b06, b07, b08, &
                          b09, b10, b11, b12, &
                          b13, b14, b15, b16)
        implicit none
        double precision, intent(in) :: &
            b01, b02, b03, b04, &
            b05, b06, b07, b08, &
            b09, b10, b11, b12, &
            b13, b14, b15, b16
        double precision, intent(inout) :: &
            a01, a02, a03, a04, &
            a05, a06, a07, a08, &
            a09, a10, a11, a12, &
            a13, a14, a15, a16

        a01 = a01/b01
        a02 = a02/b02
        a03 = a03/b03
        a04 = a04/b04
        a05 = a05/b05
        a06 = a06/b06
        a07 = a07/b07
        a08 = a08/b08
        a09 = a09/b09
        a10 = a10/b10
        a11 = a11/b11
        a12 = a12/b12
        a13 = a13/b13
        a14 = a14/b14
        a15 = a15/b15
        a16 = a16/b16
    end subroutine

    subroutine div2_block(a01, a02, a03, a04, &
                          a05, a06, a07, a08, &
                          a09, a10, a11, a12, &
                          a13, a14, a15, a16, &
                          b01, b02, b03, b04, &
                          b05, b06, b07, b08, &
                          b09, b10, b11, b12, &
                          b13, b14, b15, b16)
        implicit none
        double precision, intent(in) :: &
            b01, b02, b03, b04, &
            b05, b06, b07, b08, &
            b09, b10, b11, b12, &
            b13, b14, b15, b16
        double precision, intent(inout) :: &
            a01, a02, a03, a04, &
            a05, a06, a07, a08, &
            a09, a10, a11, a12, &
            a13, a14, a15, a16

        a01 = b01/a01
        a02 = b02/a02
        a03 = b03/a03
        a04 = b04/a04
        a05 = b05/a05
        a06 = b06/a06
        a07 = b07/a07
        a08 = b08/a08
        a09 = b09/a09
        a10 = b10/a10
        a11 = b11/a11
        a12 = b12/a12
        a13 = b13/a13
        a14 = b14/a14
        a15 = b15/a15
        a16 = b16/a16
    end subroutine
    subroutine err_block(err, &
                         rnorm01, rnorm02, rnorm03, rnorm04, &
                         rnorm05, rnorm06, rnorm07, rnorm08, &
                         rnorm09, rnorm10, rnorm11, rnorm12, &
                         rnorm13, rnorm14, rnorm15, rnorm16, &
                         bnorm01, bnorm02, bnorm03, bnorm04, &
                         bnorm05, bnorm06, bnorm07, bnorm08, &
                         bnorm09, bnorm10, bnorm11, bnorm12, &
                         bnorm13, bnorm14, bnorm15, bnorm16)
        implicit none
        double precision, intent(out) :: err
        double precision, intent(in) :: &
            rnorm01, rnorm02, rnorm03, rnorm04, &
            rnorm05, rnorm06, rnorm07, rnorm08, &
            rnorm09, rnorm10, rnorm11, rnorm12, &
            rnorm13, rnorm14, rnorm15, rnorm16, &
            bnorm01, bnorm02, bnorm03, bnorm04, &
            bnorm05, bnorm06, bnorm07, bnorm08, &
            bnorm09, bnorm10, bnorm11, bnorm12, &
            bnorm13, bnorm14, bnorm15, bnorm16

        err = max(rnorm01/bnorm01, rnorm02/bnorm02, rnorm03/bnorm03, rnorm04/bnorm04, &
                  rnorm05/bnorm05, rnorm06/bnorm06, rnorm07/bnorm07, rnorm08/bnorm08, &
                  rnorm09/bnorm09, rnorm10/bnorm10, rnorm11/bnorm11, rnorm12/bnorm12, &
                  rnorm13/bnorm13, rnorm14/bnorm14, rnorm15/bnorm15, rnorm16/bnorm16)
    end subroutine

    subroutine block_conjugate_gradient &
        (ndof, nnz, nblock, amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
        implicit none
        double precision, intent(in) :: amat_val(:), amat_diag_inv(:), &
            bvec(:, :)
        double precision, intent(inout) :: uvec(:, :)
        integer, intent(in) :: ndof, nnz, nblock, amat_col(:), amat_ind(:)

        double precision :: rvec(nblock, ndof), &
            pvec(nblock, ndof), &
            qvec(nblock, ndof), &
            zvec(nblock, ndof), &
            bnorm01, bnorm02, bnorm03, bnorm04, &
            bnorm05, bnorm06, bnorm07, bnorm08, &
            bnorm09, bnorm10, bnorm11, bnorm12, &
            bnorm13, bnorm14, bnorm15, bnorm16, &
            rnorm01, rnorm02, rnorm03, rnorm04, &
            rnorm05, rnorm06, rnorm07, rnorm08, &
            rnorm09, rnorm10, rnorm11, rnorm12, &
            rnorm13, rnorm14, rnorm15, rnorm16, &
            alpha01, alpha02, alpha03, alpha04, &
            alpha05, alpha06, alpha07, alpha08, &
            alpha09, alpha10, alpha11, alpha12, &
            alpha13, alpha14, alpha15, alpha16, &
            beta01, beta02, beta03, beta04, &
            beta05, beta06, beta07, beta08, &
            beta09, beta10, beta11, beta12, &
            beta13, beta14, beta15, beta16, &
            rho01, rho02, rho03, rho04, &
            rho05, rho06, rho07, rho08, &
            rho09, rho10, rho11, rho12, &
            rho13, rho14, rho15, rho16, &
            tmp, err
        integer :: iblock, iter

        ! initialize
        beta01 = 0d0; beta02 = 0d0; beta03 = 0d0; beta04 = 0d0
        beta05 = 0d0; beta06 = 0d0; beta07 = 0d0; beta08 = 0d0
        beta09 = 0d0; beta10 = 0d0; beta11 = 0d0; beta12 = 0d0
        beta13 = 0d0; beta14 = 0d0; beta15 = 0d0; beta16 = 0d0
        iter = 1

        print *, "start solving..."

        !$acc data copyin(amat_val, amat_col, amat_ind, amat_diag_inv, bvec, uvec) &
        !$acc create(rvec, pvec, qvec, zvec) &
        !$acc copyout(uvec)
        ! calculate bnorm
        call norm_block(ndof, nblock, bvec, &
                        bnorm01, bnorm02, bnorm03, bnorm04, &
                        bnorm05, bnorm06, bnorm07, bnorm08, &
                        bnorm09, bnorm10, bnorm11, bnorm12, &
                        bnorm13, bnorm14, bnorm15, bnorm16)

        ! q <- A u
        call spmatvec_block(ndof, nnz, nblock, amat_val, amat_col, amat_ind, uvec, qvec, tmp)

        ! r <- b - q
        call xpby_block(ndof, nblock, bvec, qvec, rvec, &
                        -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, &
                        -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, -1d0)

        ! calculate rnorm
        call norm_block(ndof, nblock, rvec, &
                        rnorm01, rnorm02, rnorm03, rnorm04, &
                        rnorm05, rnorm06, rnorm07, rnorm08, &
                        rnorm09, rnorm10, rnorm11, rnorm12, &
                        rnorm13, rnorm14, rnorm15, rnorm16)

        ! calculate err
        call err_block(err, &
                       rnorm01, rnorm02, rnorm03, rnorm04, &
                       rnorm05, rnorm06, rnorm07, rnorm08, &
                       rnorm09, rnorm10, rnorm11, rnorm12, &
                       rnorm13, rnorm14, rnorm15, rnorm16, &
                       bnorm01, bnorm02, bnorm03, bnorm04, &
                       bnorm05, bnorm06, bnorm07, bnorm08, &
                       bnorm09, bnorm10, bnorm11, bnorm12, &
                       bnorm13, bnorm14, bnorm15, bnorm16)

        do while (err > 10d0**(-8))
            ! z <- M^-1 r
            call preconditioner_block(ndof, nblock, rvec, zvec, amat_diag_inv)

            if (iter > 1) then
                ! beta <- (r, z) / rho
                call dot_block(ndof, nblock, rvec, zvec, &
                               beta01, beta02, beta03, beta04, &
                               beta05, beta06, beta07, beta08, &
                               beta09, beta10, beta11, beta12, &
                               beta13, beta14, beta15, beta16)
                call div1_block(beta01, beta02, beta03, beta04, &
                                beta05, beta06, beta07, beta08, &
                                beta09, beta10, beta11, beta12, &
                                beta13, beta14, beta15, beta16, &
                                rho01, rho02, rho03, rho04, &
                                rho05, rho06, rho07, rho08, &
                                rho09, rho10, rho11, rho12, &
                                rho13, rho14, rho15, rho16)
            end if

            ! p <- z + beta p
            call xpby_block(ndof, nblock, zvec, pvec, pvec, &
                            beta01, beta02, beta03, beta04, beta05, beta06, beta07, beta08, &
                            beta09, beta10, beta11, beta12, beta13, beta14, beta15, beta16)

            ! q <- A p
            call spmatvec_block(ndof, nnz, nblock, amat_val, amat_col, amat_ind, pvec, qvec, tmp)

            ! rho <- (r, z)
            call dot_block(ndof, nblock, rvec, zvec, &
                           rho01, rho02, rho03, rho04, &
                           rho05, rho06, rho07, rho08, &
                           rho09, rho10, rho11, rho12, &
                           rho13, rho14, rho15, rho16)

            ! alpha <- rho / (p, q)
            call dot_block(ndof, nblock, pvec, qvec, &
                           alpha01, alpha02, alpha03, alpha04, &
                           alpha05, alpha06, alpha07, alpha08, &
                           alpha09, alpha10, alpha11, alpha12, &
                           alpha13, alpha14, alpha15, alpha16)
            call div2_block(alpha01, alpha02, alpha03, alpha04, &
                            alpha05, alpha06, alpha07, alpha08, &
                            alpha09, alpha10, alpha11, alpha12, &
                            alpha13, alpha14, alpha15, alpha16, &
                            rho01, rho02, rho03, rho04, &
                            rho05, rho06, rho07, rho08, &
                            rho09, rho10, rho11, rho12, &
                            rho13, rho14, rho15, rho16)

            ! r <- r - alpha q
            call xpby_block(ndof, nblock, rvec, qvec, rvec, &
                            -alpha01, -alpha02, -alpha03, -alpha04, -alpha05, -alpha06, -alpha07, -alpha08, &
                            -alpha09, -alpha10, -alpha11, -alpha12, -alpha13, -alpha14, -alpha15, -alpha16)

            ! u <- u + alpha p
            call xpby_block(ndof, nblock, uvec, pvec, uvec, &
                            alpha01, alpha02, alpha03, alpha04, alpha05, alpha06, alpha07, alpha08, &
                            alpha09, alpha10, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16)

            ! calculate rnorm
            call norm_block(ndof, nblock, rvec, &
                            rnorm01, rnorm02, rnorm03, rnorm04, &
                            rnorm05, rnorm06, rnorm07, rnorm08, &
                            rnorm09, rnorm10, rnorm11, rnorm12, &
                            rnorm13, rnorm14, rnorm15, rnorm16)

            call err_block(err, &
                           rnorm01, rnorm02, rnorm03, rnorm04, &
                           rnorm05, rnorm06, rnorm07, rnorm08, &
                           rnorm09, rnorm10, rnorm11, rnorm12, &
                           rnorm13, rnorm14, rnorm15, rnorm16, &
                           bnorm01, bnorm02, bnorm03, bnorm04, &
                           bnorm05, bnorm06, bnorm07, bnorm08, &
                           bnorm09, bnorm10, bnorm11, bnorm12, &
                           bnorm13, bnorm14, bnorm15, bnorm16)

            iter = iter + 1

            ! if (iter == 20) then
            !     print *, "maxiter reached"
            !     exit
            ! end if
            print *, "iter = ", iter, "err = ", err
        end do
        !$acc end data
        print *, "maxiter = ", iter
    end subroutine

    subroutine validate(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        integer, intent(in) :: ndof, nnz, nblock, kglobal_col(:), kglobal_ind(:)
        double precision, intent(in) :: kglobal_val(:), kglobal_diag_inv(:), uglobal(:, :), fglobal(:, :)

        integer :: iblock, idof
        double precision :: res(nblock, ndof), tmp, err

        !$acc data copyin(uglobal, res, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, fglobal) &
        !$acc create(tmp, res, err) &
        !$acc copyout(err)
        call spmatvec_block(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, uglobal, res, tmp)

        !$acc kernels present(res, fglobal)
        !$acc loop independent collapse(2) reduction(+:err)
        do idof = 1, ndof
            do iblock = 1, nblock
                err = err + (res(iblock, idof) - fglobal(iblock, idof))**2
            end do
        end do
        !$acc end kernels
        !$acc end data

        print *, "L2 norm of residual: ", sqrt(err)

    end
end module
