using DelimitedFiles
using WriteVTK
using LinearAlgebra

function gen_model()
    ds = 0.2
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 2.0
    zmin = 0.0
    zmax = 4.0

    # number of elements in each direction
    mx = Int((xmax - xmin) / ds)
    my = Int((ymax - ymin) / ds)
    mz = Int((zmax - zmin) / ds)

    # coordinates
    coor = zeros(Float64, 3, (mx + 1) * (my + 1) * (mz + 1))
    for iz = 1:mz+1
        for iy = 1:my+1
            for ix = 1:mx+1
                i = (iz - 1) * (mx + 1) * (my + 1) + (iy - 1) * (mx + 1) + ix
                coor[1, i] = xmin + (ix - 1) * ds
                coor[2, i] = ymin + (iy - 1) * ds
                coor[3, i] = zmin + (iz - 1) * ds
            end
        end
    end

    # connectivity
    cny = zeros(Int, 8, mx * my * mz)
    for iz = 1:mz
        for iy = 1:my
            for ix = 1:mx
                i = (iz - 1) * mx * my + (iy - 1) * mx + ix
                cny[1, i] = (iz - 1) * (mx + 1) * (my + 1) + (iy - 1) * (mx + 1) + ix
                cny[2, i] = (iz - 1) * (mx + 1) * (my + 1) + (iy - 1) * (mx + 1) + ix + 1
                cny[3, i] = (iz - 1) * (mx + 1) * (my + 1) + iy * (mx + 1) + ix + 1
                cny[4, i] = (iz - 1) * (mx + 1) * (my + 1) + iy * (mx + 1) + ix
                cny[5, i] = iz * (mx + 1) * (my + 1) + (iy - 1) * (mx + 1) + ix
                cny[6, i] = iz * (mx + 1) * (my + 1) + (iy - 1) * (mx + 1) + ix + 1
                cny[7, i] = iz * (mx + 1) * (my + 1) + iy * (mx + 1) + ix + 1
                cny[8, i] = iz * (mx + 1) * (my + 1) + iy * (mx + 1) + ix
            end
        end
    end

    @show "nelem: ", size(cny, 2)
    @show "nnode: ", size(coor, 2)
    return cny, coor
end

function write(cny, coor, kglobal_val, kglobal_col, kglobal_ind, uglobal, fglobal)
    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, cny[:, i]) for i = 1:size(cny, 2)]

    ux = uglobal[1:3:end]
    uy = uglobal[2:3:end]
    uz = uglobal[3:3:end]
    vtk_grid("../data/mesh", coor, cells) do vtk
        vtk["ux"] = ux
        vtk["uy"] = uy
        vtk["uz"] = uz
    end

    writedlm("../data/cny.dat", cny)
    writedlm("../data/coor.dat", coor)
    writedlm("../data/kglobal_val.dat", kglobal_val)
    writedlm("../data/kglobal_col.dat", kglobal_col)
    writedlm("../data/kglobal_ind.dat", kglobal_ind)
    writedlm("../data/uinit.dat", uglobal)
    writedlm("../data/fglobal.dat", fglobal)
end

function spmatvec!(val, col, ind, vec, res)
    for i in 1:length(ind)-1
        tmp = 0.0
        for j in ind[i]:ind[i+1]-1
            tmp += val[j] * vec[col[j]]
        end
        res[i] = tmp
    end
end

# solve A u = b
function ConjugateGradient!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
    rvec = zeros(size(bvec))
    pvec = zeros(size(bvec))
    qvec = zeros(size(bvec))
    zvec = zeros(size(bvec))

    bnorm = sqrt(dot(bvec, bvec))

    # initial guess
    uvec .= bvec .* amat_diag_inv

    # r <- b - A u
    # qvec .= amat * uvec
    spmatvec!(amat_val, amat_col, amat_ind, uvec, qvec)
    rvec .= bvec .- qvec
    rnorm = sqrt(dot(rvec, rvec))

    beta = 0.0
    rho = 0.0
    iter = 1
    @show rnorm / bnorm
    while rnorm / bnorm > 1.0e-8
        # z <- M^-1 r
        zvec .= rvec .* amat_diag_inv

        if iter > 1
            # beta <- (r, z) / rho
            beta = dot(rvec, zvec) / rho
        end

        # p <- z + beta p
        pvec .= zvec .+ beta .* pvec

        # q <- A p
        # qvec .= amat * pvec
        spmatvec!(amat_val, amat_col, amat_ind, pvec, qvec)

        # rho <- (r, z)
        rho = dot(rvec, zvec)

        # alpha <- rho / (p, q)
        alpha = rho / dot(pvec, qvec)

        # q <- -alpha q
        qvec .*= -alpha

        # r <- r + q
        rvec .+= qvec

        # u <- u + alpha p
        uvec .+= alpha .* pvec

        rnorm = sqrt(dot(rvec, rvec))
        @show rnorm / bnorm
        iter += 1
    end
end

function gen_matvec(cny, coor)
    # material properties
    young = 10.0
    nyu = 0.3
    lam = young * nyu / (1.0 + nyu) / (1.0 - 2.0 * nyu)
    mu = young / 2.0 / (1.0 + nyu)

    nnode = size(coor, 2)
    nelement = size(cny, 2)
    kglobal_tmp_val = [[] for _ in 1:3*nnode]
    kglobal_tmp_col = [[] for _ in 1:3*nnode]
    fglobal = zeros(Float64, 3 * nnode)
    uglobal = zeros(Float64, 3 * nnode)

    for ie = 1:nelement
        node_id = zeros(Int, 8)
        for inode = 1:8
            node_id[inode] = cny[inode, ie]
        end
        xnode = zeros(Float64, 3, 8)
        for inode = 1:8
            xnode[1, inode] = coor[1, node_id[inode]]
            xnode[2, inode] = coor[2, node_id[inode]]
            xnode[3, inode] = coor[3, node_id[inode]]
        end

        kelement = zeros(Float64, 24, 24)
        # gauss quadrature
        for r1 in [-1.0 / sqrt(3), 1.0 / sqrt(3)]
            for r2 in [-1.0 / sqrt(3), 1.0 / sqrt(3)]
                for r3 in [-1.0 / sqrt(3), 1.0 / sqrt(3)]
                    nvec = zeros(Float64, 8)
                    nvec[1] = 0.125 * (1.0 - r1) * (1.0 - r2) * (1.0 - r3)
                    nvec[2] = 0.125 * (1.0 + r1) * (1.0 - r2) * (1.0 - r3)
                    nvec[3] = 0.125 * (1.0 + r1) * (1.0 + r2) * (1.0 - r3)
                    nvec[4] = 0.125 * (1.0 - r1) * (1.0 + r2) * (1.0 - r3)
                    nvec[5] = 0.125 * (1.0 - r1) * (1.0 - r2) * (1.0 + r3)
                    nvec[6] = 0.125 * (1.0 + r1) * (1.0 - r2) * (1.0 + r3)
                    nvec[7] = 0.125 * (1.0 + r1) * (1.0 + r2) * (1.0 + r3)
                    nvec[8] = 0.125 * (1.0 - r1) * (1.0 + r2) * (1.0 + r3)

                    nmat = zeros(Float64, 3, 24)
                    for inode = 1:8
                        nmat[1, 3*inode-2] = nvec[inode]
                        nmat[2, 3*inode-1] = nvec[inode]
                        nmat[3, 3*inode] = nvec[inode]
                    end

                    dndr = zeros(Float64, 3, 8)
                    dndr[1, 1] = -0.125 * (1.0 - r2) * (1.0 - r3)
                    dndr[1, 2] = 0.125 * (1.0 - r2) * (1.0 - r3)
                    dndr[1, 3] = 0.125 * (1.0 + r2) * (1.0 - r3)
                    dndr[1, 4] = -0.125 * (1.0 + r2) * (1.0 - r3)
                    dndr[1, 5] = -0.125 * (1.0 - r2) * (1.0 + r3)
                    dndr[1, 6] = 0.125 * (1.0 - r2) * (1.0 + r3)
                    dndr[1, 7] = 0.125 * (1.0 + r2) * (1.0 + r3)
                    dndr[1, 8] = -0.125 * (1.0 + r2) * (1.0 + r3)

                    dndr[2, 1] = -0.125 * (1.0 - r1) * (1.0 - r3)
                    dndr[2, 2] = -0.125 * (1.0 + r1) * (1.0 - r3)
                    dndr[2, 3] = 0.125 * (1.0 + r1) * (1.0 - r3)
                    dndr[2, 4] = 0.125 * (1.0 - r1) * (1.0 - r3)
                    dndr[2, 5] = -0.125 * (1.0 - r1) * (1.0 + r3)
                    dndr[2, 6] = -0.125 * (1.0 + r1) * (1.0 + r3)
                    dndr[2, 7] = 0.125 * (1.0 + r1) * (1.0 + r3)
                    dndr[2, 8] = 0.125 * (1.0 - r1) * (1.0 + r3)

                    dndr[3, 1] = -0.125 * (1.0 - r1) * (1.0 - r2)
                    dndr[3, 2] = -0.125 * (1.0 + r1) * (1.0 - r2)
                    dndr[3, 3] = -0.125 * (1.0 + r1) * (1.0 + r2)
                    dndr[3, 4] = -0.125 * (1.0 - r1) * (1.0 + r2)
                    dndr[3, 5] = 0.125 * (1.0 - r1) * (1.0 - r2)
                    dndr[3, 6] = 0.125 * (1.0 + r1) * (1.0 - r2)
                    dndr[3, 7] = 0.125 * (1.0 + r1) * (1.0 + r2)
                    dndr[3, 8] = 0.125 * (1.0 - r1) * (1.0 + r2)

                    dxdr = zeros(Float64, 3, 3)
                    for inode in 1:8
                        for i in 1:3
                            for j in 1:3
                                dxdr[i, j] += dndr[i, inode] * xnode[j, inode]
                            end
                        end
                    end
                    jac = det(dxdr)
                    drdx = inv(dxdr)

                    dndx = zeros(Float64, 3, 8)
                    for inode in 1:8
                        for i in 1:3
                            for j in 1:3
                                dndx[i, inode] += drdx[j, i] * dndr[j, inode]
                            end
                        end
                    end

                    bmat = zeros(Float64, 6, 24)
                    for inode in 1:8
                        bmat[1, 3*(inode-1)+1] = dndx[1, inode]
                        bmat[2, 3*(inode-1)+2] = dndx[2, inode]
                        bmat[3, 3*(inode-1)+3] = dndx[3, inode]
                        bmat[4, 3*(inode-1)+1] = dndx[2, inode]
                        bmat[4, 3*(inode-1)+2] = dndx[1, inode]
                        bmat[5, 3*(inode-1)+2] = dndx[3, inode]
                        bmat[5, 3*(inode-1)+3] = dndx[2, inode]
                        bmat[6, 3*(inode-1)+1] = dndx[3, inode]
                        bmat[6, 3*(inode-1)+3] = dndx[1, inode]
                    end

                    dmat = zeros(Float64, 6, 6)
                    dmat[1, 1] = lam + 2.0 * mu
                    dmat[2, 2] = lam + 2.0 * mu
                    dmat[3, 3] = lam + 2.0 * mu
                    dmat[4, 4] = mu
                    dmat[5, 5] = mu
                    dmat[6, 6] = mu
                    dmat[1, 2] = lam
                    dmat[1, 3] = lam
                    dmat[2, 1] = lam
                    dmat[2, 3] = lam
                    dmat[3, 1] = lam
                    dmat[3, 2] = lam

                    kelement .+= transpose(bmat) * dmat * bmat * jac
                end
            end
        end
        for inode in 1:8
            for jnode in 1:8
                for i in 1:3
                    idof = 3 * (node_id[inode] - 1) + i
                    for j in 1:3
                        jdof = 3 * (node_id[jnode] - 1) + j
                        knode = findfirst(x -> x == jdof, kglobal_tmp_col[idof])
                        if isnothing(knode)
                            push!(kglobal_tmp_col[idof], jdof)
                            push!(kglobal_tmp_val[idof], kelement[3*(inode-1)+i, 3*(jnode-1)+j])
                        else
                            kglobal_tmp_val[idof][knode] += kelement[3*(inode-1)+i, 3*(jnode-1)+j]
                        end
                    end
                end
            end
        end
    end

    nnz = sum(length, kglobal_tmp_col)
    kglobal_col = zeros(Int, nnz)
    kglobal_val = zeros(Float64, nnz)
    kglobal_ind = zeros(Int, 3 * nnode + 1)
    kglobal_ind[1] = 1
    for idof in 1:3*nnode
        kglobal_ind[idof+1] = kglobal_ind[idof] + length(kglobal_tmp_col[idof])
        kglobal_col[kglobal_ind[idof]:kglobal_ind[idof+1]-1] = kglobal_tmp_col[idof]
        kglobal_val[kglobal_ind[idof]:kglobal_ind[idof+1]-1] = kglobal_tmp_val[idof]
    end

    kglobal_diag = zeros(Float64, 3 * nnode)
    kglobal_diag_inv = zeros(Float64, 3 * nnode)

    # inpulse force
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        if norm([x, y, z] .- [1.0, 1.0, 4.0]) < 1.0e-6
            fglobal[3*(inode-1)+3] = -1.0
        end
    end

    # boundary conditions
    bcflag = zeros(Bool, 3 * nnode)
    uglobal = zeros(Float64, 3 * nnode)
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        if abs(z) < 1.0e-6
            bcflag[3*(inode-1)+1] = true
            bcflag[3*(inode-1)+2] = true
            bcflag[3*(inode-1)+3] = true
            uglobal[3*(inode-1)+1] = 0.0
            uglobal[3*(inode-1)+2] = 0.0
            uglobal[3*(inode-1)+3] = 0.0
        end
    end

    for idof in 1:3*nnode
        if bcflag[idof]
            for jdof in 1:kglobal_ind[idof+1]-kglobal_ind[idof]
                if (kglobal_col[kglobal_ind[idof]+jdof-1] == idof)
                    kglobal_val[kglobal_ind[idof]+jdof-1] = 1.0
                    kglobal_diag[idof] = 1.0
                else
                    kglobal_val[kglobal_ind[idof]+jdof-1] = 0.0
                end
            end
            fglobal[idof] = uglobal[idof]
        else
            for jdof in 1:kglobal_ind[idof+1]-kglobal_ind[idof]
                if bcflag[kglobal_col[kglobal_ind[idof]+jdof-1]]
                    fglobal[idof] -= kglobal_val[kglobal_ind[idof]+jdof-1] * uglobal[kglobal_col[kglobal_ind[idof]+jdof-1]]
                    kglobal_val[kglobal_ind[idof]+jdof-1] = 0.0
                end
                if (kglobal_col[kglobal_ind[idof]+jdof-1] == idof)
                    kglobal_diag[idof] = kglobal_val[kglobal_ind[idof]+jdof-1]
                end
            end
        end
    end

    for idof in 1:3*nnode
        kglobal_diag_inv[idof] = 1.0 / kglobal_diag[idof]
    end

    # solve linear system
    ConjugateGradient!(kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

    return kglobal_val, kglobal_col, kglobal_ind, uglobal, fglobal
end

cny, coor = gen_model()

kglobal_val, kglobal_col, kglobal_ind, uglobal, fglobal = gen_matvec(cny, coor)

write(cny, coor, kglobal_val, kglobal_col, kglobal_ind, uglobal, fglobal)
