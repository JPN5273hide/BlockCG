module solver_acc_breakfree
    implicit none
    private :: spmatvec, calc_norm, xpby, dot
    public :: block_conjugate_gradient, validate
contains
    subroutine calc_norm(ndof, nblock, vec, norm)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: vec(:, :)
        double precision, intent(inout) :: norm(:)
        double precision :: &
            norm01, norm02, norm03, norm04, norm05, norm06, norm07, norm08, &
            norm09, norm10, norm11, norm12, norm13, norm14, norm15, norm16
        integer :: idof

        norm01 = 0d0; norm02 = 0d0; norm03 = 0d0; norm04 = 0d0
        norm05 = 0d0; norm06 = 0d0; norm07 = 0d0; norm08 = 0d0
        norm09 = 0d0; norm10 = 0d0; norm11 = 0d0; norm12 = 0d0
        norm13 = 0d0; norm14 = 0d0; norm15 = 0d0; norm16 = 0d0

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

        norm(1) = sqrt(norm01)
        norm(2) = sqrt(norm02)
        norm(3) = sqrt(norm03)
        norm(4) = sqrt(norm04)
        norm(5) = sqrt(norm05)
        norm(6) = sqrt(norm06)
        norm(7) = sqrt(norm07)
        norm(8) = sqrt(norm08)
        norm(9) = sqrt(norm09)
        norm(10) = sqrt(norm10)
        norm(11) = sqrt(norm11)
        norm(12) = sqrt(norm12)
        norm(13) = sqrt(norm13)
        norm(14) = sqrt(norm14)
        norm(15) = sqrt(norm15)
        norm(16) = sqrt(norm16)
    end subroutine

    subroutine spmatvec(ndof, nnz, nblock, val, col, ind, vec, res)
        implicit none
        integer, intent(in) :: ndof, nnz, nblock, col(:), ind(:)
        double precision, intent(in) :: val(:), vec(:, :)
        double precision, intent(inout) :: res(:, :)

        double precision :: tmp
        integer :: idof, iblock, j

        !$acc kernels present(val, col, ind, vec, res)
        !$acc loop independent collapse(2) private(tmp)
        do idof = 1, ndof
            do iblock = 1, nblock
                tmp = 0d0
                do j = ind(idof), ind(idof + 1) - 1
                    tmp = tmp + val(j)*vec(iblock, col(j))
                end do
                res(iblock, idof) = tmp
            end do
        end do
        !$acc end kernels
    end subroutine

    ! subroutine dot(ndof, nblock, vec1, vec2, res)
    !     implicit none
    !     integer, intent(in) :: ndof, nblock
    !     double precision, intent(in) :: vec1(:, :), vec2(:, :)
    !     double precision, intent(inout) :: res(:, :)
    !     integer :: idof, iblock, jblock

    !     do jblock = 1, nblock
    !         do iblock = 1, nblock
    !             res(iblock, jblock) = 0d0
    !         end do
    !     end do
    !     !$acc kernels present(vec1, vec2)
    !     !$acc loop independent collapse(3)
    !     do idof = 1, ndof
    !         do jblock = 1, nblock
    !             do iblock = 1, nblock
    !                 !$acc atomic
    !                 res(iblock, jblock) = res(iblock, jblock) + vec1(iblock, idof)*vec2(jblock, idof)
    !             end do
    !         end do
    !     end do
    !     !$acc end kernels
    ! end subroutine

    subroutine dot(ndof, nblock, vec1, vec2, res)
        implicit none
        integer, intent(in) :: ndof, nblock
        double precision, intent(in) :: vec1(:, :), vec2(:, :)
        double precision, intent(inout) :: res(:, :)
        integer :: idof, iblock, jblock
        double precision :: &
            res0101, res0102, res0103, res0104, res0105, res0106, res0107, res0108, &
            res0109, res0110, res0111, res0112, res0113, res0114, res0115, res0116, &
            res0201, res0202, res0203, res0204, res0205, res0206, res0207, res0208, &
            res0209, res0210, res0211, res0212, res0213, res0214, res0215, res0216, &
            res0301, res0302, res0303, res0304, res0305, res0306, res0307, res0308, &
            res0309, res0310, res0311, res0312, res0313, res0314, res0315, res0316, &
            res0401, res0402, res0403, res0404, res0405, res0406, res0407, res0408, &
            res0409, res0410, res0411, res0412, res0413, res0414, res0415, res0416, &
            res0501, res0502, res0503, res0504, res0505, res0506, res0507, res0508, &
            res0509, res0510, res0511, res0512, res0513, res0514, res0515, res0516, &
            res0601, res0602, res0603, res0604, res0605, res0606, res0607, res0608, &
            res0609, res0610, res0611, res0612, res0613, res0614, res0615, res0616, &
            res0701, res0702, res0703, res0704, res0705, res0706, res0707, res0708, &
            res0709, res0710, res0711, res0712, res0713, res0714, res0715, res0716, &
            res0801, res0802, res0803, res0804, res0805, res0806, res0807, res0808, &
            res0809, res0810, res0811, res0812, res0813, res0814, res0815, res0816, &
            res0901, res0902, res0903, res0904, res0905, res0906, res0907, res0908, &
            res0909, res0910, res0911, res0912, res0913, res0914, res0915, res0916, &
            res1001, res1002, res1003, res1004, res1005, res1006, res1007, res1008, &
            res1009, res1010, res1011, res1012, res1013, res1014, res1015, res1016, &
            res1101, res1102, res1103, res1104, res1105, res1106, res1107, res1108, &
            res1109, res1110, res1111, res1112, res1113, res1114, res1115, res1116, &
            res1201, res1202, res1203, res1204, res1205, res1206, res1207, res1208, &
            res1209, res1210, res1211, res1212, res1213, res1214, res1215, res1216, &
            res1301, res1302, res1303, res1304, res1305, res1306, res1307, res1308, &
            res1309, res1310, res1311, res1312, res1313, res1314, res1315, res1316, &
            res1401, res1402, res1403, res1404, res1405, res1406, res1407, res1408, &
            res1409, res1410, res1411, res1412, res1413, res1414, res1415, res1416, &
            res1501, res1502, res1503, res1504, res1505, res1506, res1507, res1508, &
            res1509, res1510, res1511, res1512, res1513, res1514, res1515, res1516, &
            res1601, res1602, res1603, res1604, res1605, res1606, res1607, res1608, &
            res1609, res1610, res1611, res1612, res1613, res1614, res1615, res1616

        res0101 = 0d0; res0102 = 0d0; res0103 = 0d0; res0104 = 0d0
        res0105 = 0d0; res0106 = 0d0; res0107 = 0d0; res0108 = 0d0
        res0109 = 0d0; res0110 = 0d0; res0111 = 0d0; res0112 = 0d0
        res0113 = 0d0; res0114 = 0d0; res0115 = 0d0; res0116 = 0d0
        res0201 = 0d0; res0202 = 0d0; res0203 = 0d0; res0204 = 0d0
        res0205 = 0d0; res0206 = 0d0; res0207 = 0d0; res0208 = 0d0
        res0209 = 0d0; res0210 = 0d0; res0211 = 0d0; res0212 = 0d0
        res0213 = 0d0; res0214 = 0d0; res0215 = 0d0; res0216 = 0d0
        res0301 = 0d0; res0302 = 0d0; res0303 = 0d0; res0304 = 0d0
        res0305 = 0d0; res0306 = 0d0; res0307 = 0d0; res0308 = 0d0
        res0309 = 0d0; res0310 = 0d0; res0311 = 0d0; res0312 = 0d0
        res0313 = 0d0; res0314 = 0d0; res0315 = 0d0; res0316 = 0d0
        res0401 = 0d0; res0402 = 0d0; res0403 = 0d0; res0404 = 0d0
        res0405 = 0d0; res0406 = 0d0; res0407 = 0d0; res0408 = 0d0
        res0409 = 0d0; res0410 = 0d0; res0411 = 0d0; res0412 = 0d0
        res0413 = 0d0; res0414 = 0d0; res0415 = 0d0; res0416 = 0d0
        res0501 = 0d0; res0502 = 0d0; res0503 = 0d0; res0504 = 0d0
        res0505 = 0d0; res0506 = 0d0; res0507 = 0d0; res0508 = 0d0
        res0509 = 0d0; res0510 = 0d0; res0511 = 0d0; res0512 = 0d0
        res0513 = 0d0; res0514 = 0d0; res0515 = 0d0; res0516 = 0d0
        res0601 = 0d0; res0602 = 0d0; res0603 = 0d0; res0604 = 0d0
        res0605 = 0d0; res0606 = 0d0; res0607 = 0d0; res0608 = 0d0
        res0609 = 0d0; res0610 = 0d0; res0611 = 0d0; res0612 = 0d0
        res0613 = 0d0; res0614 = 0d0; res0615 = 0d0; res0616 = 0d0
        res0701 = 0d0; res0702 = 0d0; res0703 = 0d0; res0704 = 0d0
        res0705 = 0d0; res0706 = 0d0; res0707 = 0d0; res0708 = 0d0
        res0709 = 0d0; res0710 = 0d0; res0711 = 0d0; res0712 = 0d0
        res0713 = 0d0; res0714 = 0d0; res0715 = 0d0; res0716 = 0d0
        res0801 = 0d0; res0802 = 0d0; res0803 = 0d0; res0804 = 0d0
        res0805 = 0d0; res0806 = 0d0; res0807 = 0d0; res0808 = 0d0
        res0809 = 0d0; res0810 = 0d0; res0811 = 0d0; res0812 = 0d0
        res0813 = 0d0; res0814 = 0d0; res0815 = 0d0; res0816 = 0d0
        res0901 = 0d0; res0902 = 0d0; res0903 = 0d0; res0904 = 0d0
        res0905 = 0d0; res0906 = 0d0; res0907 = 0d0; res0908 = 0d0
        res0909 = 0d0; res0910 = 0d0; res0911 = 0d0; res0912 = 0d0
        res0913 = 0d0; res0914 = 0d0; res0915 = 0d0; res0916 = 0d0
        res1001 = 0d0; res1002 = 0d0; res1003 = 0d0; res1004 = 0d0
        res1005 = 0d0; res1006 = 0d0; res1007 = 0d0; res1008 = 0d0
        res1009 = 0d0; res1010 = 0d0; res1011 = 0d0; res1012 = 0d0
        res1013 = 0d0; res1014 = 0d0; res1015 = 0d0; res1016 = 0d0
        res1101 = 0d0; res1102 = 0d0; res1103 = 0d0; res1104 = 0d0
        res1105 = 0d0; res1106 = 0d0; res1107 = 0d0; res1108 = 0d0
        res1109 = 0d0; res1110 = 0d0; res1111 = 0d0; res1112 = 0d0
        res1113 = 0d0; res1114 = 0d0; res1115 = 0d0; res1116 = 0d0
        res1201 = 0d0; res1202 = 0d0; res1203 = 0d0; res1204 = 0d0
        res1205 = 0d0; res1206 = 0d0; res1207 = 0d0; res1208 = 0d0
        res1209 = 0d0; res1210 = 0d0; res1211 = 0d0; res1212 = 0d0
        res1213 = 0d0; res1214 = 0d0; res1215 = 0d0; res1216 = 0d0
        res1301 = 0d0; res1302 = 0d0; res1303 = 0d0; res1304 = 0d0
        res1305 = 0d0; res1306 = 0d0; res1307 = 0d0; res1308 = 0d0
        res1309 = 0d0; res1310 = 0d0; res1311 = 0d0; res1312 = 0d0
        res1313 = 0d0; res1314 = 0d0; res1315 = 0d0; res1316 = 0d0
        res1401 = 0d0; res1402 = 0d0; res1403 = 0d0; res1404 = 0d0
        res1405 = 0d0; res1406 = 0d0; res1407 = 0d0; res1408 = 0d0
        res1409 = 0d0; res1410 = 0d0; res1411 = 0d0; res1412 = 0d0
        res1413 = 0d0; res1414 = 0d0; res1415 = 0d0; res1416 = 0d0
        res1501 = 0d0; res1502 = 0d0; res1503 = 0d0; res1504 = 0d0
        res1505 = 0d0; res1506 = 0d0; res1507 = 0d0; res1508 = 0d0
        res1509 = 0d0; res1510 = 0d0; res1511 = 0d0; res1512 = 0d0
        res1513 = 0d0; res1514 = 0d0; res1515 = 0d0; res1516 = 0d0
        res1601 = 0d0; res1602 = 0d0; res1603 = 0d0; res1604 = 0d0
        res1605 = 0d0; res1606 = 0d0; res1607 = 0d0; res1608 = 0d0
        res1609 = 0d0; res1610 = 0d0; res1611 = 0d0; res1612 = 0d0
        res1613 = 0d0; res1614 = 0d0; res1615 = 0d0; res1616 = 0d0

        !$acc kernels present(vec1, vec2)
        !$acc loop independent &
        !$acc reduction(+:res0101, res0102, res0103, res0104, &
        !$acc             res0105, res0106, res0107, res0108, &
        !$acc             res0109, res0110, res0111, res0112, &
        !$acc             res0113, res0114, res0115, res0116, &
        !$acc             res0201, res0202, res0203, res0204, &
        !$acc             res0205, res0206, res0207, res0208, &
        !$acc             res0209, res0210, res0211, res0212, &
        !$acc             res0213, res0214, res0215, res0216, &
        !$acc             res0301, res0302, res0303, res0304, &
        !$acc             res0305, res0306, res0307, res0308, &
        !$acc             res0309, res0310, res0311, res0312, &
        !$acc             res0313, res0314, res0315, res0316, &
        !$acc             res0401, res0402, res0403, res0404, &
        !$acc             res0405, res0406, res0407, res0408, &
        !$acc             res0409, res0410, res0411, res0412, &
        !$acc             res0413, res0414, res0415, res0416, &
        !$acc             res0501, res0502, res0503, res0504, &
        !$acc             res0505, res0506, res0507, res0508, &
        !$acc             res0509, res0510, res0511, res0512, &
        !$acc             res0513, res0514, res0515, res0516, &
        !$acc             res0601, res0602, res0603, res0604, &
        !$acc             res0605, res0606, res0607, res0608, &
        !$acc             res0609, res0610, res0611, res0612, &
        !$acc             res0613, res0614, res0615, res0616, &
        !$acc             res0701, res0702, res0703, res0704, &
        !$acc             res0705, res0706, res0707, res0708, &
        !$acc             res0709, res0710, res0711, res0712, &
        !$acc             res0713, res0714, res0715, res0716, &
        !$acc             res0801, res0802, res0803, res0804, &
        !$acc             res0805, res0806, res0807, res0808, &
        !$acc             res0809, res0810, res0811, res0812, &
        !$acc             res0813, res0814, res0815, res0816, &
        !$acc             res0901, res0902, res0903, res0904, &
        !$acc             res0905, res0906, res0907, res0908, &
        !$acc             res0909, res0910, res0911, res0912, &
        !$acc             res0913, res0914, res0915, res0916, &
        !$acc             res1001, res1002, res1003, res1004, &
        !$acc             res1005, res1006, res1007, res1008, &
        !$acc             res1009, res1010, res1011, res1012, &
        !$acc             res1013, res1014, res1015, res1016, &
        !$acc             res1101, res1102, res1103, res1104, &
        !$acc             res1105, res1106, res1107, res1108, &
        !$acc             res1109, res1110, res1111, res1112, &
        !$acc             res1113, res1114, res1115, res1116, &
        !$acc             res1201, res1202, res1203, res1204, &
        !$acc             res1205, res1206, res1207, res1208, &
        !$acc             res1209, res1210, res1211, res1212, &
        !$acc             res1213, res1214, res1215, res1216, &
        !$acc             res1301, res1302, res1303, res1304, &
        !$acc             res1305, res1306, res1307, res1308, &
        !$acc             res1309, res1310, res1311, res1312, &
        !$acc             res1313, res1314, res1315, res1316, &
        !$acc             res1401, res1402, res1403, res1404, &
        !$acc             res1405, res1406, res1407, res1408, &
        !$acc             res1409, res1410, res1411, res1412, &
        !$acc             res1413, res1414, res1415, res1416, &
        !$acc             res1501, res1502, res1503, res1504, &
        !$acc             res1505, res1506, res1507, res1508, &
        !$acc             res1509, res1510, res1511, res1512, &
        !$acc             res1513, res1514, res1515, res1516, &
        !$acc             res1601, res1602, res1603, res1604, &
        !$acc             res1605, res1606, res1607, res1608, &
        !$acc             res1609, res1610, res1611, res1612, &
        !$acc             res1613, res1614, res1615, res1616)
        do idof = 1, ndof
            res0101 = res0101 + vec1(1, idof)*vec2(1, idof)
            res0102 = res0102 + vec1(1, idof)*vec2(2, idof)
            res0103 = res0103 + vec1(1, idof)*vec2(3, idof)
            res0104 = res0104 + vec1(1, idof)*vec2(4, idof)
            res0105 = res0105 + vec1(1, idof)*vec2(5, idof)
            res0106 = res0106 + vec1(1, idof)*vec2(6, idof)
            res0107 = res0107 + vec1(1, idof)*vec2(7, idof)
            res0108 = res0108 + vec1(1, idof)*vec2(8, idof)
            res0109 = res0109 + vec1(1, idof)*vec2(9, idof)
            res0110 = res0110 + vec1(1, idof)*vec2(10, idof)
            res0111 = res0111 + vec1(1, idof)*vec2(11, idof)
            res0112 = res0112 + vec1(1, idof)*vec2(12, idof)
            res0113 = res0113 + vec1(1, idof)*vec2(13, idof)
            res0114 = res0114 + vec1(1, idof)*vec2(14, idof)
            res0115 = res0115 + vec1(1, idof)*vec2(15, idof)
            res0116 = res0116 + vec1(1, idof)*vec2(16, idof)
            res0201 = res0201 + vec1(2, idof)*vec2(1, idof)
            res0202 = res0202 + vec1(2, idof)*vec2(2, idof)
            res0203 = res0203 + vec1(2, idof)*vec2(3, idof)
            res0204 = res0204 + vec1(2, idof)*vec2(4, idof)
            res0205 = res0205 + vec1(2, idof)*vec2(5, idof)
            res0206 = res0206 + vec1(2, idof)*vec2(6, idof)
            res0207 = res0207 + vec1(2, idof)*vec2(7, idof)
            res0208 = res0208 + vec1(2, idof)*vec2(8, idof)
            res0209 = res0209 + vec1(2, idof)*vec2(9, idof)
            res0210 = res0210 + vec1(2, idof)*vec2(10, idof)
            res0211 = res0211 + vec1(2, idof)*vec2(11, idof)
            res0212 = res0212 + vec1(2, idof)*vec2(12, idof)
            res0213 = res0213 + vec1(2, idof)*vec2(13, idof)
            res0214 = res0214 + vec1(2, idof)*vec2(14, idof)
            res0215 = res0215 + vec1(2, idof)*vec2(15, idof)
            res0216 = res0216 + vec1(2, idof)*vec2(16, idof)
            res0301 = res0301 + vec1(3, idof)*vec2(1, idof)
            res0302 = res0302 + vec1(3, idof)*vec2(2, idof)
            res0303 = res0303 + vec1(3, idof)*vec2(3, idof)
            res0304 = res0304 + vec1(3, idof)*vec2(4, idof)
            res0305 = res0305 + vec1(3, idof)*vec2(5, idof)
            res0306 = res0306 + vec1(3, idof)*vec2(6, idof)
            res0307 = res0307 + vec1(3, idof)*vec2(7, idof)
            res0308 = res0308 + vec1(3, idof)*vec2(8, idof)
            res0309 = res0309 + vec1(3, idof)*vec2(9, idof)
            res0310 = res0310 + vec1(3, idof)*vec2(10, idof)
            res0311 = res0311 + vec1(3, idof)*vec2(11, idof)
            res0312 = res0312 + vec1(3, idof)*vec2(12, idof)
            res0313 = res0313 + vec1(3, idof)*vec2(13, idof)
            res0314 = res0314 + vec1(3, idof)*vec2(14, idof)
            res0315 = res0315 + vec1(3, idof)*vec2(15, idof)
            res0316 = res0316 + vec1(3, idof)*vec2(16, idof)
            res0401 = res0401 + vec1(4, idof)*vec2(1, idof)
            res0402 = res0402 + vec1(4, idof)*vec2(2, idof)
            res0403 = res0403 + vec1(4, idof)*vec2(3, idof)
            res0404 = res0404 + vec1(4, idof)*vec2(4, idof)
            res0405 = res0405 + vec1(4, idof)*vec2(5, idof)
            res0406 = res0406 + vec1(4, idof)*vec2(6, idof)
            res0407 = res0407 + vec1(4, idof)*vec2(7, idof)
            res0408 = res0408 + vec1(4, idof)*vec2(8, idof)
            res0409 = res0409 + vec1(4, idof)*vec2(9, idof)
            res0410 = res0410 + vec1(4, idof)*vec2(10, idof)
            res0411 = res0411 + vec1(4, idof)*vec2(11, idof)
            res0412 = res0412 + vec1(4, idof)*vec2(12, idof)
            res0413 = res0413 + vec1(4, idof)*vec2(13, idof)
            res0414 = res0414 + vec1(4, idof)*vec2(14, idof)
            res0415 = res0415 + vec1(4, idof)*vec2(15, idof)
            res0416 = res0416 + vec1(4, idof)*vec2(16, idof)
            res0501 = res0501 + vec1(5, idof)*vec2(1, idof)
            res0502 = res0502 + vec1(5, idof)*vec2(2, idof)
            res0503 = res0503 + vec1(5, idof)*vec2(3, idof)
            res0504 = res0504 + vec1(5, idof)*vec2(4, idof)
            res0505 = res0505 + vec1(5, idof)*vec2(5, idof)
            res0506 = res0506 + vec1(5, idof)*vec2(6, idof)
            res0507 = res0507 + vec1(5, idof)*vec2(7, idof)
            res0508 = res0508 + vec1(5, idof)*vec2(8, idof)
            res0509 = res0509 + vec1(5, idof)*vec2(9, idof)
            res0510 = res0510 + vec1(5, idof)*vec2(10, idof)
            res0511 = res0511 + vec1(5, idof)*vec2(11, idof)
            res0512 = res0512 + vec1(5, idof)*vec2(12, idof)
            res0513 = res0513 + vec1(5, idof)*vec2(13, idof)
            res0514 = res0514 + vec1(5, idof)*vec2(14, idof)
            res0515 = res0515 + vec1(5, idof)*vec2(15, idof)
            res0516 = res0516 + vec1(5, idof)*vec2(16, idof)
            res0601 = res0601 + vec1(6, idof)*vec2(1, idof)
            res0602 = res0602 + vec1(6, idof)*vec2(2, idof)
            res0603 = res0603 + vec1(6, idof)*vec2(3, idof)
            res0604 = res0604 + vec1(6, idof)*vec2(4, idof)
            res0605 = res0605 + vec1(6, idof)*vec2(5, idof)
            res0606 = res0606 + vec1(6, idof)*vec2(6, idof)
            res0607 = res0607 + vec1(6, idof)*vec2(7, idof)
            res0608 = res0608 + vec1(6, idof)*vec2(8, idof)
            res0609 = res0609 + vec1(6, idof)*vec2(9, idof)
            res0610 = res0610 + vec1(6, idof)*vec2(10, idof)
            res0611 = res0611 + vec1(6, idof)*vec2(11, idof)
            res0612 = res0612 + vec1(6, idof)*vec2(12, idof)
            res0613 = res0613 + vec1(6, idof)*vec2(13, idof)
            res0614 = res0614 + vec1(6, idof)*vec2(14, idof)
            res0615 = res0615 + vec1(6, idof)*vec2(15, idof)
            res0616 = res0616 + vec1(6, idof)*vec2(16, idof)
            res0701 = res0701 + vec1(7, idof)*vec2(1, idof)
            res0702 = res0702 + vec1(7, idof)*vec2(2, idof)
            res0703 = res0703 + vec1(7, idof)*vec2(3, idof)
            res0704 = res0704 + vec1(7, idof)*vec2(4, idof)
            res0705 = res0705 + vec1(7, idof)*vec2(5, idof)
            res0706 = res0706 + vec1(7, idof)*vec2(6, idof)
            res0707 = res0707 + vec1(7, idof)*vec2(7, idof)
            res0708 = res0708 + vec1(7, idof)*vec2(8, idof)
            res0709 = res0709 + vec1(7, idof)*vec2(9, idof)
            res0710 = res0710 + vec1(7, idof)*vec2(10, idof)
            res0711 = res0711 + vec1(7, idof)*vec2(11, idof)
            res0712 = res0712 + vec1(7, idof)*vec2(12, idof)
            res0713 = res0713 + vec1(7, idof)*vec2(13, idof)
            res0714 = res0714 + vec1(7, idof)*vec2(14, idof)
            res0715 = res0715 + vec1(7, idof)*vec2(15, idof)
            res0716 = res0716 + vec1(7, idof)*vec2(16, idof)
            res0801 = res0801 + vec1(8, idof)*vec2(1, idof)
            res0802 = res0802 + vec1(8, idof)*vec2(2, idof)
            res0803 = res0803 + vec1(8, idof)*vec2(3, idof)
            res0804 = res0804 + vec1(8, idof)*vec2(4, idof)
            res0805 = res0805 + vec1(8, idof)*vec2(5, idof)
            res0806 = res0806 + vec1(8, idof)*vec2(6, idof)
            res0807 = res0807 + vec1(8, idof)*vec2(7, idof)
            res0808 = res0808 + vec1(8, idof)*vec2(8, idof)
            res0809 = res0809 + vec1(8, idof)*vec2(9, idof)
            res0810 = res0810 + vec1(8, idof)*vec2(10, idof)
            res0811 = res0811 + vec1(8, idof)*vec2(11, idof)
            res0812 = res0812 + vec1(8, idof)*vec2(12, idof)
            res0813 = res0813 + vec1(8, idof)*vec2(13, idof)
            res0814 = res0814 + vec1(8, idof)*vec2(14, idof)
            res0815 = res0815 + vec1(8, idof)*vec2(15, idof)
            res0816 = res0816 + vec1(8, idof)*vec2(16, idof)
            res0901 = res0901 + vec1(9, idof)*vec2(1, idof)
            res0902 = res0902 + vec1(9, idof)*vec2(2, idof)
            res0903 = res0903 + vec1(9, idof)*vec2(3, idof)
            res0904 = res0904 + vec1(9, idof)*vec2(4, idof)
            res0905 = res0905 + vec1(9, idof)*vec2(5, idof)
            res0906 = res0906 + vec1(9, idof)*vec2(6, idof)
            res0907 = res0907 + vec1(9, idof)*vec2(7, idof)
            res0908 = res0908 + vec1(9, idof)*vec2(8, idof)
            res0909 = res0909 + vec1(9, idof)*vec2(9, idof)
            res0910 = res0910 + vec1(9, idof)*vec2(10, idof)
            res0911 = res0911 + vec1(9, idof)*vec2(11, idof)
            res0912 = res0912 + vec1(9, idof)*vec2(12, idof)
            res0913 = res0913 + vec1(9, idof)*vec2(13, idof)
            res0914 = res0914 + vec1(9, idof)*vec2(14, idof)
            res0915 = res0915 + vec1(9, idof)*vec2(15, idof)
            res0916 = res0916 + vec1(9, idof)*vec2(16, idof)
            res1001 = res1001 + vec1(10, idof)*vec2(1, idof)
            res1002 = res1002 + vec1(10, idof)*vec2(2, idof)
            res1003 = res1003 + vec1(10, idof)*vec2(3, idof)
            res1004 = res1004 + vec1(10, idof)*vec2(4, idof)
            res1005 = res1005 + vec1(10, idof)*vec2(5, idof)
            res1006 = res1006 + vec1(10, idof)*vec2(6, idof)
            res1007 = res1007 + vec1(10, idof)*vec2(7, idof)
            res1008 = res1008 + vec1(10, idof)*vec2(8, idof)
            res1009 = res1009 + vec1(10, idof)*vec2(9, idof)
            res1010 = res1010 + vec1(10, idof)*vec2(10, idof)
            res1011 = res1011 + vec1(10, idof)*vec2(11, idof)
            res1012 = res1012 + vec1(10, idof)*vec2(12, idof)
            res1013 = res1013 + vec1(10, idof)*vec2(13, idof)
            res1014 = res1014 + vec1(10, idof)*vec2(14, idof)
            res1015 = res1015 + vec1(10, idof)*vec2(15, idof)
            res1016 = res1016 + vec1(10, idof)*vec2(16, idof)
            res1101 = res1101 + vec1(11, idof)*vec2(1, idof)
            res1102 = res1102 + vec1(11, idof)*vec2(2, idof)
            res1103 = res1103 + vec1(11, idof)*vec2(3, idof)
            res1104 = res1104 + vec1(11, idof)*vec2(4, idof)
            res1105 = res1105 + vec1(11, idof)*vec2(5, idof)
            res1106 = res1106 + vec1(11, idof)*vec2(6, idof)
            res1107 = res1107 + vec1(11, idof)*vec2(7, idof)
            res1108 = res1108 + vec1(11, idof)*vec2(8, idof)
            res1109 = res1109 + vec1(11, idof)*vec2(9, idof)
            res1110 = res1110 + vec1(11, idof)*vec2(10, idof)
            res1111 = res1111 + vec1(11, idof)*vec2(11, idof)
            res1112 = res1112 + vec1(11, idof)*vec2(12, idof)
            res1113 = res1113 + vec1(11, idof)*vec2(13, idof)
            res1114 = res1114 + vec1(11, idof)*vec2(14, idof)
            res1115 = res1115 + vec1(11, idof)*vec2(15, idof)
            res1116 = res1116 + vec1(11, idof)*vec2(16, idof)
            res1201 = res1201 + vec1(12, idof)*vec2(1, idof)
            res1202 = res1202 + vec1(12, idof)*vec2(2, idof)
            res1203 = res1203 + vec1(12, idof)*vec2(3, idof)
            res1204 = res1204 + vec1(12, idof)*vec2(4, idof)
            res1205 = res1205 + vec1(12, idof)*vec2(5, idof)
            res1206 = res1206 + vec1(12, idof)*vec2(6, idof)
            res1207 = res1207 + vec1(12, idof)*vec2(7, idof)
            res1208 = res1208 + vec1(12, idof)*vec2(8, idof)
            res1209 = res1209 + vec1(12, idof)*vec2(9, idof)
            res1210 = res1210 + vec1(12, idof)*vec2(10, idof)
            res1211 = res1211 + vec1(12, idof)*vec2(11, idof)
            res1212 = res1212 + vec1(12, idof)*vec2(12, idof)
            res1213 = res1213 + vec1(12, idof)*vec2(13, idof)
            res1214 = res1214 + vec1(12, idof)*vec2(14, idof)
            res1215 = res1215 + vec1(12, idof)*vec2(15, idof)
            res1216 = res1216 + vec1(12, idof)*vec2(16, idof)
            res1301 = res1301 + vec1(13, idof)*vec2(1, idof)
            res1302 = res1302 + vec1(13, idof)*vec2(2, idof)
            res1303 = res1303 + vec1(13, idof)*vec2(3, idof)
            res1304 = res1304 + vec1(13, idof)*vec2(4, idof)
            res1305 = res1305 + vec1(13, idof)*vec2(5, idof)
            res1306 = res1306 + vec1(13, idof)*vec2(6, idof)
            res1307 = res1307 + vec1(13, idof)*vec2(7, idof)
            res1308 = res1308 + vec1(13, idof)*vec2(8, idof)
            res1309 = res1309 + vec1(13, idof)*vec2(9, idof)
            res1310 = res1310 + vec1(13, idof)*vec2(10, idof)
            res1311 = res1311 + vec1(13, idof)*vec2(11, idof)
            res1312 = res1312 + vec1(13, idof)*vec2(12, idof)
            res1313 = res1313 + vec1(13, idof)*vec2(13, idof)
            res1314 = res1314 + vec1(13, idof)*vec2(14, idof)
            res1315 = res1315 + vec1(13, idof)*vec2(15, idof)
            res1316 = res1316 + vec1(13, idof)*vec2(16, idof)
            res1401 = res1401 + vec1(14, idof)*vec2(1, idof)
            res1402 = res1402 + vec1(14, idof)*vec2(2, idof)
            res1403 = res1403 + vec1(14, idof)*vec2(3, idof)
            res1404 = res1404 + vec1(14, idof)*vec2(4, idof)
            res1405 = res1405 + vec1(14, idof)*vec2(5, idof)
            res1406 = res1406 + vec1(14, idof)*vec2(6, idof)
            res1407 = res1407 + vec1(14, idof)*vec2(7, idof)
            res1408 = res1408 + vec1(14, idof)*vec2(8, idof)
            res1409 = res1409 + vec1(14, idof)*vec2(9, idof)
            res1410 = res1410 + vec1(14, idof)*vec2(10, idof)
            res1411 = res1411 + vec1(14, idof)*vec2(11, idof)
            res1412 = res1412 + vec1(14, idof)*vec2(12, idof)
            res1413 = res1413 + vec1(14, idof)*vec2(13, idof)
            res1414 = res1414 + vec1(14, idof)*vec2(14, idof)
            res1415 = res1415 + vec1(14, idof)*vec2(15, idof)
            res1416 = res1416 + vec1(14, idof)*vec2(16, idof)
            res1501 = res1501 + vec1(15, idof)*vec2(1, idof)
            res1502 = res1502 + vec1(15, idof)*vec2(2, idof)
            res1503 = res1503 + vec1(15, idof)*vec2(3, idof)
            res1504 = res1504 + vec1(15, idof)*vec2(4, idof)
            res1505 = res1505 + vec1(15, idof)*vec2(5, idof)
            res1506 = res1506 + vec1(15, idof)*vec2(6, idof)
            res1507 = res1507 + vec1(15, idof)*vec2(7, idof)
            res1508 = res1508 + vec1(15, idof)*vec2(8, idof)
            res1509 = res1509 + vec1(15, idof)*vec2(9, idof)
            res1510 = res1510 + vec1(15, idof)*vec2(10, idof)
            res1511 = res1511 + vec1(15, idof)*vec2(11, idof)
            res1512 = res1512 + vec1(15, idof)*vec2(12, idof)
            res1513 = res1513 + vec1(15, idof)*vec2(13, idof)
            res1514 = res1514 + vec1(15, idof)*vec2(14, idof)
            res1515 = res1515 + vec1(15, idof)*vec2(15, idof)
            res1516 = res1516 + vec1(15, idof)*vec2(16, idof)
            res1601 = res1601 + vec1(16, idof)*vec2(1, idof)
            res1602 = res1602 + vec1(16, idof)*vec2(2, idof)
            res1603 = res1603 + vec1(16, idof)*vec2(3, idof)
            res1604 = res1604 + vec1(16, idof)*vec2(4, idof)
            res1605 = res1605 + vec1(16, idof)*vec2(5, idof)
            res1606 = res1606 + vec1(16, idof)*vec2(6, idof)
            res1607 = res1607 + vec1(16, idof)*vec2(7, idof)
            res1608 = res1608 + vec1(16, idof)*vec2(8, idof)
            res1609 = res1609 + vec1(16, idof)*vec2(9, idof)
            res1610 = res1610 + vec1(16, idof)*vec2(10, idof)
            res1611 = res1611 + vec1(16, idof)*vec2(11, idof)
            res1612 = res1612 + vec1(16, idof)*vec2(12, idof)
            res1613 = res1613 + vec1(16, idof)*vec2(13, idof)
            res1614 = res1614 + vec1(16, idof)*vec2(14, idof)
            res1615 = res1615 + vec1(16, idof)*vec2(15, idof)
            res1616 = res1616 + vec1(16, idof)*vec2(16, idof)
        end do
        !$acc end kernels

        res(1, 1) = res0101; res(1, 2) = res0102; res(1, 3) = res0103; res(1, 4) = res0104
        res(1, 5) = res0105; res(1, 6) = res0106; res(1, 7) = res0107; res(1, 8) = res0108
        res(1, 9) = res0109; res(1, 10) = res0110; res(1, 11) = res0111; res(1, 12) = res0112
        res(1, 13) = res0113; res(1, 14) = res0114; res(1, 15) = res0115; res(1, 16) = res0116

        res(2, 1) = res0201; res(2, 2) = res0202; res(2, 3) = res0203; res(2, 4) = res0204
        res(2, 5) = res0205; res(2, 6) = res0206; res(2, 7) = res0207; res(2, 8) = res0208
        res(2, 9) = res0209; res(2, 10) = res0210; res(2, 11) = res0211; res(2, 12) = res0212
        res(2, 13) = res0213; res(2, 14) = res0214; res(2, 15) = res0215; res(2, 16) = res0216

        res(3, 1) = res0301; res(3, 2) = res0302; res(3, 3) = res0303; res(3, 4) = res0304
        res(3, 5) = res0305; res(3, 6) = res0306; res(3, 7) = res0307; res(3, 8) = res0308
        res(3, 9) = res0309; res(3, 10) = res0310; res(3, 11) = res0311; res(3, 12) = res0312
        res(3, 13) = res0313; res(3, 14) = res0314; res(3, 15) = res0315; res(3, 16) = res0316

        res(4, 1) = res0401; res(4, 2) = res0402; res(4, 3) = res0403; res(4, 4) = res0404
        res(4, 5) = res0405; res(4, 6) = res0406; res(4, 7) = res0407; res(4, 8) = res0408
        res(4, 9) = res0409; res(4, 10) = res0410; res(4, 11) = res0411; res(4, 12) = res0412
        res(4, 13) = res0413; res(4, 14) = res0414; res(4, 15) = res0415; res(4, 16) = res0416

        res(5, 1) = res0501; res(5, 2) = res0502; res(5, 3) = res0503; res(5, 4) = res0504
        res(5, 5) = res0505; res(5, 6) = res0506; res(5, 7) = res0507; res(5, 8) = res0508
        res(5, 9) = res0509; res(5, 10) = res0510; res(5, 11) = res0511; res(5, 12) = res0512
        res(5, 13) = res0513; res(5, 14) = res0514; res(5, 15) = res0515; res(5, 16) = res0516

        res(6, 1) = res0601; res(6, 2) = res0602; res(6, 3) = res0603; res(6, 4) = res0604
        res(6, 5) = res0605; res(6, 6) = res0606; res(6, 7) = res0607; res(6, 8) = res0608
        res(6, 9) = res0609; res(6, 10) = res0610; res(6, 11) = res0611; res(6, 12) = res0612
        res(6, 13) = res0613; res(6, 14) = res0614; res(6, 15) = res0615; res(6, 16) = res0616

        res(7, 1) = res0701; res(7, 2) = res0702; res(7, 3) = res0703; res(7, 4) = res0704
        res(7, 5) = res0705; res(7, 6) = res0706; res(7, 7) = res0707; res(7, 8) = res0708
        res(7, 9) = res0709; res(7, 10) = res0710; res(7, 11) = res0711; res(7, 12) = res0712
        res(7, 13) = res0713; res(7, 14) = res0714; res(7, 15) = res0715; res(7, 16) = res0716

        res(8, 1) = res0801; res(8, 2) = res0802; res(8, 3) = res0803; res(8, 4) = res0804
        res(8, 5) = res0805; res(8, 6) = res0806; res(8, 7) = res0807; res(8, 8) = res0808
        res(8, 9) = res0809; res(8, 10) = res0810; res(8, 11) = res0811; res(8, 12) = res0812
        res(8, 13) = res0813; res(8, 14) = res0814; res(8, 15) = res0815; res(8, 16) = res0816

        res(9, 1) = res0901; res(9, 2) = res0902; res(9, 3) = res0903; res(9, 4) = res0904
        res(9, 5) = res0905; res(9, 6) = res0906; res(9, 7) = res0907; res(9, 8) = res0908
        res(9, 9) = res0909; res(9, 10) = res0910; res(9, 11) = res0911; res(9, 12) = res0912
        res(9, 13) = res0913; res(9, 14) = res0914; res(9, 15) = res0915; res(9, 16) = res0916

        res(10, 1) = res1001; res(10, 2) = res1002; res(10, 3) = res1003; res(10, 4) = res1004
        res(10, 5) = res1005; res(10, 6) = res1006; res(10, 7) = res1007; res(10, 8) = res1008
        res(10, 9) = res1009; res(10, 10) = res1010; res(10, 11) = res1011; res(10, 12) = res1012
        res(10, 13) = res1013; res(10, 14) = res1014; res(10, 15) = res1015; res(10, 16) = res1016

        res(11, 1) = res1101; res(11, 2) = res1102; res(11, 3) = res1103; res(11, 4) = res1104
        res(11, 5) = res1105; res(11, 6) = res1106; res(11, 7) = res1107; res(11, 8) = res1108
        res(11, 9) = res1109; res(11, 10) = res1110; res(11, 11) = res1111; res(11, 12) = res1112
        res(11, 13) = res1113; res(11, 14) = res1114; res(11, 15) = res1115; res(11, 16) = res1116

        res(12, 1) = res1201; res(12, 2) = res1202; res(12, 3) = res1203; res(12, 4) = res1204
        res(12, 5) = res1205; res(12, 6) = res1206; res(12, 7) = res1207; res(12, 8) = res1208
        res(12, 9) = res1209; res(12, 10) = res1210; res(12, 11) = res1211; res(12, 12) = res1212
        res(12, 13) = res1213; res(12, 14) = res1214; res(12, 15) = res1215; res(12, 16) = res1216

        res(13, 1) = res1301; res(13, 2) = res1302; res(13, 3) = res1303; res(13, 4) = res1304
        res(13, 5) = res1305; res(13, 6) = res1306; res(13, 7) = res1307; res(13, 8) = res1308
        res(13, 9) = res1309; res(13, 10) = res1310; res(13, 11) = res1311; res(13, 12) = res1312
        res(13, 13) = res1313; res(13, 14) = res1314; res(13, 15) = res1315; res(13, 16) = res1316

        res(14, 1) = res1401; res(14, 2) = res1402; res(14, 3) = res1403; res(14, 4) = res1404
        res(14, 5) = res1405; res(14, 6) = res1406; res(14, 7) = res1407; res(14, 8) = res1408
        res(14, 9) = res1409; res(14, 10) = res1410; res(14, 11) = res1411; res(14, 12) = res1412
        res(14, 13) = res1413; res(14, 14) = res1414; res(14, 15) = res1415; res(14, 16) = res1416

        res(15, 1) = res1501; res(15, 2) = res1502; res(15, 3) = res1503; res(15, 4) = res1504
        res(15, 5) = res1505; res(15, 6) = res1506; res(15, 7) = res1507; res(15, 8) = res1508
        res(15, 9) = res1509; res(15, 10) = res1510; res(15, 11) = res1511; res(15, 12) = res1512
        res(15, 13) = res1513; res(15, 14) = res1514; res(15, 15) = res1515; res(15, 16) = res1516

        res(16, 1) = res1601; res(16, 2) = res1602; res(16, 3) = res1603; res(16, 4) = res1604
        res(16, 5) = res1605; res(16, 6) = res1606; res(16, 7) = res1607; res(16, 8) = res1608
        res(16, 9) = res1609; res(16, 10) = res1610; res(16, 11) = res1611; res(16, 12) = res1612
        res(16, 13) = res1613; res(16, 14) = res1614; res(16, 15) = res1615; res(16, 16) = res1616
    end subroutine

    subroutine xpby(nblock, ndof, xvec, yvec, b, zvec)
        implicit none
        integer, intent(in) :: nblock, ndof
        double precision, intent(in) :: xvec(:, :), yvec(:, :), b(:, :)
        double precision, intent(inout) :: zvec(:, :)

        integer :: idof, iblock, jblock
        double precision :: tmp_b(nblock), tmp

        !$acc kernels present(xvec, yvec, zvec)
        !$acc loop independent private(tmp, tmp_b, &
        !$acc     iblock, jblock)
        do idof = 1, ndof
            !$acc loop seq
            do iblock = 1, nblock
                tmp = 0d0
                !$acc loop seq
                do jblock = 1, nblock
                    tmp = tmp + b(jblock, iblock)*yvec(jblock, idof)
                end do
                tmp_b(iblock) = tmp
            end do
            !$acc loop seq
            do iblock = 1, nblock
                zvec(iblock, idof) = xvec(iblock, idof) + tmp_b(iblock)
            end do
        end do
        !$acc end kernels
    end subroutine

    subroutine gram_schmidt(nblock, ndof, vec, tol)
        implicit none
        integer, intent(in) :: nblock, ndof
        double precision, intent(inout) :: vec(:, :), tol
        double precision :: v1(ndof), v2(ndof), tmp, tmp_b(nblock), &
            tmp_b01, tmp_b02, tmp_b03, tmp_b04, tmp_b05, tmp_b06, tmp_b07, tmp_b08, &
            tmp_b09, tmp_b10, tmp_b11, tmp_b12, tmp_b13, tmp_b14, tmp_b15, tmp_b16

        integer :: idof, iblock, jblock

        !$acc data create(v1, v2)
        do iblock = 1, nblock
            ! v1 <- vec(iblock, :)
            ! v2 <- vec(iblock, :)
            !$acc kernels present(vec, v1, v2)
            !$acc loop independent
            do idof = 1, ndof
                v1(idof) = vec(iblock, idof)
                v2(idof) = vec(iblock, idof)
            end do
            !$acc end kernels

            tmp_b01 = 0d0; tmp_b02 = 0d0; tmp_b03 = 0d0; tmp_b04 = 0d0
            tmp_b05 = 0d0; tmp_b06 = 0d0; tmp_b07 = 0d0; tmp_b08 = 0d0
            tmp_b09 = 0d0; tmp_b10 = 0d0; tmp_b11 = 0d0; tmp_b12 = 0d0
            tmp_b13 = 0d0; tmp_b14 = 0d0; tmp_b15 = 0d0; tmp_b16 = 0d0
            !$acc kernels present(vec, v1)
            !$acc loop independent &
            !$acc reduction(+:tmp_b01, tmp_b02, tmp_b03, tmp_b04, &
            !$acc             tmp_b05, tmp_b06, tmp_b07, tmp_b08, &
            !$acc             tmp_b09, tmp_b10, tmp_b11, tmp_b12, &
            !$acc             tmp_b13, tmp_b14, tmp_b15, tmp_b16)
            do idof = 1, ndof
                tmp_b01 = tmp_b01 + v1(idof)*vec(1, idof)
                tmp_b02 = tmp_b02 + v1(idof)*vec(2, idof)
                tmp_b03 = tmp_b03 + v1(idof)*vec(3, idof)
                tmp_b04 = tmp_b04 + v1(idof)*vec(4, idof)
                tmp_b05 = tmp_b05 + v1(idof)*vec(5, idof)
                tmp_b06 = tmp_b06 + v1(idof)*vec(6, idof)
                tmp_b07 = tmp_b07 + v1(idof)*vec(7, idof)
                tmp_b08 = tmp_b08 + v1(idof)*vec(8, idof)
                tmp_b09 = tmp_b09 + v1(idof)*vec(9, idof)
                tmp_b10 = tmp_b10 + v1(idof)*vec(10, idof)
                tmp_b11 = tmp_b11 + v1(idof)*vec(11, idof)
                tmp_b12 = tmp_b12 + v1(idof)*vec(12, idof)
                tmp_b13 = tmp_b13 + v1(idof)*vec(13, idof)
                tmp_b14 = tmp_b14 + v1(idof)*vec(14, idof)
                tmp_b15 = tmp_b15 + v1(idof)*vec(15, idof)
                tmp_b16 = tmp_b16 + v1(idof)*vec(16, idof)
            end do
            !$acc end kernels

            tmp_b(1) = tmp_b01; tmp_b(2) = tmp_b02; tmp_b(3) = tmp_b03; tmp_b(4) = tmp_b04
            tmp_b(5) = tmp_b05; tmp_b(6) = tmp_b06; tmp_b(7) = tmp_b07; tmp_b(8) = tmp_b08
            tmp_b(9) = tmp_b09; tmp_b(10) = tmp_b10; tmp_b(11) = tmp_b11; tmp_b(12) = tmp_b12
            tmp_b(13) = tmp_b13; tmp_b(14) = tmp_b14; tmp_b(15) = tmp_b15; tmp_b(16) = tmp_b16
            !$acc kernels present(vec, v2)
            !$acc loop independent private(jblock)
            do idof = 1, ndof
                do jblock = 1, iblock - 1
                    v2(idof) = v2(idof) - tmp_b(jblock)*vec(jblock, idof)
                end do
            end do
            !$acc end kernels

            ! tmp <- norm(v2)
            tmp = 0d0
            !$acc kernels present(v2)
            !$acc loop independent reduction(+:tmp)
            do idof = 1, ndof
                tmp = tmp + v2(idof)**2
            end do
            !$acc end kernels
            tmp = sqrt(tmp)

            if (tmp > tol) then
                ! vec(iblock, :) <- v2 / tmp
                !$acc kernels present(vec, v2)
                !$acc loop independent
                do idof = 1, ndof
                    vec(iblock, idof) = v2(idof)/tmp
                end do
                !$acc end kernels
            else
                ! vec(iblock, :) <- 0
                !$acc kernels present(vec, v2)
                !$acc loop independent
                do idof = 1, ndof
                    vec(iblock, idof) = 0d0
                end do
                !$acc end kernels
            end if
        end do
        !$acc end data
    end subroutine

    ! Calculate pseudo-inverse of vec
    subroutine pinv(nblock, vec, tol)
        implicit none
        integer, intent(in) :: nblock
        double precision, intent(in) :: tol
        double precision, intent(inout) :: vec(:, :)
        integer :: info, i, j, k
        double precision :: U(nblock, nblock), S(nblock), &
            VT(nblock, nblock), &
            tmp, work(5*nblock)

        ! Call LAPACK subroutine for SVD
        call dgesvd('A', 'A', nblock, nblock, &
                    vec, nblock, S, U, nblock, VT, nblock, &
                    work, 5*nblock, info)

        ! Compute S^(-1)
        do i = 1, nblock
            if (S(i) > tol) then
                S(i) = 1.0d0/S(i)
            end if
        end do

        ! Calculate A_pseudo_inv = VT^T * S_inv * U^T
        do j = 1, nblock
            do i = 1, nblock
                vec(i, j) = 0.0d0
                do k = 1, nblock
                    vec(i, j) = vec(i, j) + VT(k, i)*S(k)*U(j, k)
                end do
            end do
        end do
    end subroutine

    subroutine block_conjugate_gradient &
        (ndof, nnz, nblock, amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
        implicit none
        double precision, intent(in) :: amat_val(:), amat_diag_inv(:), &
            bvec(:, :)
        double precision, intent(inout) :: uvec(:, :)
        integer, intent(in) :: ndof, nnz, nblock, amat_col(:), amat_ind(:)

        double precision :: &
            rvec(nblock, ndof), pvec(nblock, ndof), &
            qvec(nblock, ndof), zvec(nblock, ndof), &
            alpha(nblock, nblock), beta(nblock, nblock), &
            rho(nblock, nblock), bnorm(nblock), rnorm(nblock), &
            tmp, tmp_bb(nblock, nblock), err, tol
        integer :: idof, iblock, jblock, kblock, iter, max_iter

        tol = 1d-8
        max_iter = 20

        !$acc data copyin(amat_val, amat_col, amat_ind, amat_diag_inv, bvec) &
        !$acc create(uvec, rvec, pvec, qvec, zvec) &
        !$acc copyout(uvec)
        ! initializations
        do iblock = 1, nblock
            do jblock = 1, nblock
                beta(iblock, jblock) = 0d0
            end do
        end do
        iter = 1
        call calc_norm(ndof, nblock, bvec, bnorm)

        ! u <- u0
        !$acc kernels present(uvec)
        !$acc loop independent private(iblock)
        do idof = 1, ndof
            do iblock = 1, nblock
                uvec(iblock, idof) = 0d0
            end do
        end do
        !$acc end kernels

        ! q <- A u
        call spmatvec(ndof, nnz, nblock, amat_val, amat_col, amat_ind, uvec, qvec)

        ! r <- b - q
        !$acc kernels present(bvec, qvec, rvec)
        !$acc loop independent private(iblock)
        do idof = 1, ndof
            do iblock = 1, nblock
                rvec(iblock, idof) = bvec(iblock, idof) - qvec(iblock, idof)
            end do
        end do
        !$acc end kernels
        call calc_norm(ndof, nblock, rvec, rnorm)
        err = maxval(rnorm/bnorm)

        do while (err > tol)
            ! z <- M^-1 r
            !$acc kernels present(rvec, zvec, amat_diag_inv)
            !$acc loop independent private(iblock)
            do idof = 1, ndof
                do iblock = 1, nblock
                    zvec(iblock, idof) = amat_diag_inv(idof)*rvec(iblock, idof)
                end do
            end do
            !$acc end kernels

            if (iter > 1) then
                ! TODO: beta <- rho^-1 * (R, Z)
                call pinv(nblock, rho, tol)
                call dot(ndof, nblock, rvec, zvec, tmp_bb)
                do iblock = 1, nblock
                    do jblock = 1, nblock
                        tmp = 0d0
                        do kblock = 1, nblock
                            tmp = tmp - rho(iblock, kblock)*tmp_bb(kblock, jblock)
                        end do
                    end do
                end do
                
                !deleted: beta <- -rho * (q, z)
                !deleted: call dot(ndof, nblock, qvec, zvec, tmp_bb)
                !deleted: do iblock = 1, nblock
                !deleted:     do jblock = 1, nblock
                !deleted:         tmp = 0d0
                !deleted:         do kblock = 1, nblock
                !deleted:             tmp = tmp - rho(iblock, kblock)*tmp_bb(kblock, jblock)
                !deleted:         end do
                !deleted:         beta(iblock, jblock) = tmp
                !deleted:     end do
                !deleted: end do
            end if

            ! p <- z + beta p
            call xpby(nblock, ndof, zvec, pvec, beta, pvec)

            ! q <- A p
            call spmatvec(ndof, nnz, nblock, amat_val, amat_col, amat_ind, pvec, qvec)

            ! rho <- (R, Z)
            call dot(ndof, nblock, rvec, zvec, rho)
            ! deleted: call dot(ndof, nblock, pvec, qvec, rho)
            ! deleted: call pinv(nblock, rho, tol)

            ! alpha <- (P, Q)^{-1} * rho
            call dot(ndof, nblock, pvec, qvec, tmp_bb)
            call pinv(nblock, tmp_bb, tol)
            do iblock=1, nblock
                do jblock = 1, nblock
                    tmp = 0d0
                    do kblock = 1, nblock
                        tmp = tmp + tmp_bb(iblock, kblock)*rho(iblock, jblock)
                    end do
                    alpha(iblock, jblock) = tmp
                end do
            end do
            ! alpha <- rho * (p, r)
            ! call dot(ndof, nblock, pvec, rvec, tmp_bb)
            ! do iblock = 1, nblock
            !     do jblock = 1, nblock
            !         tmp = 0d0
            !         do kblock = 1, nblock
            !             tmp = tmp + rho(iblock, kblock)*tmp_bb(kblock, jblock)
            !         end do
            !         alpha(iblock, jblock) = tmp
            !     end do
            ! end do

            ! r <- r - alpha q
            call xpby(nblock, ndof, rvec, qvec, -alpha, rvec)

            ! u <- u + alpha p
            call xpby(nblock, ndof, uvec, pvec, alpha, uvec)

            call calc_norm(ndof, nblock, rvec, rnorm)
            do iblock = 1, nblock
                if (rnorm(iblock)/bnorm(iblock) < tol) then
                    rnorm(iblock) = 0d0
                    !$acc kernels present(rvec, qvec, pvec)
                    !$acc loop independent
                    do idof = 1, ndof
                        rvec(iblock, idof) = 0d0
                        qvec(iblock, idof) = 0d0
                        pvec(iblock, idof) = 0d0
                    end do
                    !$acc end kernels
                end if
            end do
            err = maxval(rnorm/bnorm)

            iter = iter + 1
            print *, "iter = ", iter, "err = ", err

            if (iter > max_iter) then
                print *, "maximum iteration reached"
                exit
            end if
        end do
        !$acc end data
    end subroutine

    subroutine validate(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
        integer, intent(in) :: ndof, nnz, nblock, kglobal_col(:), kglobal_ind(:)
        double precision, intent(in) :: kglobal_val(:), kglobal_diag_inv(:), uglobal(:, :), fglobal(:, :)

        integer :: iblock, idof
        double precision :: res(nblock, ndof), tmp, err

        res = 0d0
        call spmatvec(ndof, nnz, nblock, kglobal_val, kglobal_col, kglobal_ind, uglobal, res)
        err = 0d0
        do idof = 1, ndof
            do iblock = 1, nblock
                err = err + (res(iblock, idof) - fglobal(iblock, idof))**2
            end do
        end do
        print *, "L2 norm of residual: ", sqrt(err)

    end
end module
