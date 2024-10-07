using Random
using DelimitedFiles
using WriteVTK
using LinearAlgebra

function gen_model(ds, xmin, xmax, ymin, ymax, zmin, zmax)

    # number of elements in each direction
    mx = Int32((xmax - xmin) / ds)
    my = Int32((ymax - ymin) / ds)
    mz = Int32((zmax - zmin) / ds)

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
    cny = zeros(Int32, 8, mx * my * mz)
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

function write_data(data_dir, ndof, nnz, nblock, cny, coor, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
    open(joinpath(data_dir, "cny.bin"), "w") do io
        write(io, cny)
    end
    open(joinpath(data_dir, "coor.bin"), "w") do io
        write(io, coor)
    end
    open(joinpath(data_dir, "kglobal_val.bin"), "w") do io
        write(io, kglobal_val)
    end
    open(joinpath(data_dir, "kglobal_col.bin"), "w") do io
        write(io, kglobal_col)
    end
    open(joinpath(data_dir, "kglobal_ind.bin"), "w") do io
        write(io, kglobal_ind)
    end
    open(joinpath(data_dir, "kglobal_diag_inv.bin"), "w") do io
        write(io, kglobal_diag_inv)
    end
    open(joinpath(data_dir, "uinit.bin"), "w") do io
        write(io, uglobal)
    end
    open(joinpath(data_dir, "fglobal.bin"), "w") do io
        write(io, fglobal)
    end
    writedlm(joinpath(data_dir, "shape.dat"), [ndof, nnz, nblock])
end

function write_data_bcrs(data_dir, nnode, nnz, nblock, cny, coor, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
    open(joinpath(data_dir, "cny.bin"), "w") do io
        write(io, cny)
    end
    open(joinpath(data_dir, "coor.bin"), "w") do io
        write(io, coor)
    end
    open(joinpath(data_dir, "kglobal_val.bin"), "w") do io
        write(io, kglobal_val)
    end
    open(joinpath(data_dir, "kglobal_col.bin"), "w") do io
        write(io, kglobal_col)
    end
    open(joinpath(data_dir, "kglobal_ind.bin"), "w") do io
        write(io, kglobal_ind)
    end
    open(joinpath(data_dir, "kglobal_diag_inv.bin"), "w") do io
        write(io, kglobal_diag_inv)
    end
    open(joinpath(data_dir, "uinit.bin"), "w") do io
        write(io, uglobal)
    end
    open(joinpath(data_dir, "fglobal.bin"), "w") do io
        write(io, fglobal)
    end
    writedlm(joinpath(data_dir, "shape.dat"), [nnode, nnz, nblock])
end

function write_mesh(data_dir, nblock)
    cny = reshape(reinterpret(Int32, read(joinpath(data_dir, "cny.bin"))), (8, :))
    coor = reshape(reinterpret(Float64, read(joinpath(data_dir, "coor.bin"))), (3, :))
    uglobal = reshape(reinterpret(Float64, read(joinpath(data_dir, "sol.bin"))), (nblock, :))
    @show size(uglobal)

    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, cny[:, i]) for i = 1:size(cny, 2)]

    ux = uglobal[1, 1:3:end]
    uy = uglobal[1, 2:3:end]
    uz = uglobal[1, 3:3:end]
    @show size(ux), size(uy), size(uz)
    vtk_grid(joinpath(data_dir, "mesh"), coor, cells) do vtk
        vtk["ux"] = ux
        vtk["uy"] = uy
        vtk["uz"] = uz
    end
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

function spmatvec_block!(val, col, ind, vec, res)
    nblock, ndof = size(vec)
    tmp = zeros(Float64, nblock)
    for i in 1:ndof
        for iblock in 1:nblock
            tmp[iblock] = 0.0
        end
        for j in ind[i]:ind[i+1]-1
            for iblock in 1:nblock
                tmp[iblock] += val[j] * vec[iblock, col[j]]
            end
        end
        for iblock in 1:nblock
            res[iblock, i] = tmp[iblock]
        end
    end
end

function spmatvec_block_bcrs!(val, col, ind, vec, res)
    nblock, _, nnode = size(vec)
    tmp = zeros(Float64, 3, nblock)
    for inode in 1:nnode
        for iblock in 1:nblock
            tmp[:, iblock] .= 0.0
        end
        for j in ind[inode]:ind[inode+1]-1
            for iblock in 1:nblock
                tmp[:, iblock] .+= val[:, :, j] * vec[iblock, :, col[j]]
            end
        end
        for iblock in 1:nblock
            res[iblock, :, inode] .= tmp[:, iblock]
        end
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

function pinv_(A)
    F = svd(A)
    S_inv = Diagonal(1.0 ./ F.S)
    return (F.Vt)' * S_inv * (F.U)'
end

# solve A u = b
function BlockConjugateGradient!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
    nblock, ndof = size(bvec)
    rvec = zeros(Float64, size(bvec))
    pvec = zeros(Float64, size(bvec))
    qvec = zeros(Float64, size(bvec))
    zvec = zeros(Float64, size(bvec))

    bnorm = zeros(Float64, nblock)
    rnorm = zeros(Float64, nblock)

    bnorm .= 0.0
    for iblock in 1:nblock
        bnorm[iblock] = norm(bvec[iblock, :])
    end

    # initial guess
    for idof in 1:ndof
        for iblock in 1:nblock
            uvec[iblock, idof] = bvec[iblock, idof] * amat_diag_inv[idof]
        end
    end

    # q <- A u
    spmatvec_block!(amat_val, amat_col, amat_ind, uvec, qvec)

    # r <- b - q
    rvec .= bvec .- qvec
    for iblock in 1:nblock
        rnorm[iblock] = norm(rvec[iblock, :])
    end

    alpha = zeros(Float64, nblock, nblock)
    beta = zeros(Float64, nblock, nblock)
    rho = zeros(Float64, nblock, nblock)
    iter = 1
    @show rnorm[1] ./ bnorm[1]
    while maximum(rnorm / bnorm) > 1.0e-8
        # z <- M^-1 r
        for idof in 1:ndof
            for iblock in 1:nblock
                zvec[iblock, idof] = rvec[iblock, idof] * amat_diag_inv[idof]
            end
        end

        if iter > 1
            # beta <- (r, z) / rho
            beta .= pinv(rho) * (rvec * transpose(zvec))
        end

        # p <- z + beta p
        pvec .= zvec .+ transpose(beta) * pvec


        # q <- A p
        spmatvec_block!(amat_val, amat_col, amat_ind, pvec, qvec)

        # rho <- (r, z)
        rho .= rvec * transpose(zvec)

        # alpha <- rho / (p, q)
        alpha .= pinv(pvec * transpose(qvec)) * rho

        # q <- -alpha q
        qvec .= -transpose(alpha) * qvec

        # r <- r + q
        rvec .= rvec .+ qvec

        # u <- u + alpha p
        uvec .= uvec .+ transpose(alpha) * pvec

        rnorm .= 0.0
        for iblock in 1:nblock
            rnorm[iblock] = norm(rvec[iblock, :])
        end
        @show rnorm[1] ./ bnorm[1]

        iter += 1
        if (iter > 500)
            break
        end
    end
    @show "max iter: ", iter
end


# solve A u = b
function BlockConjugateGradient_diag!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
    nblock, ndof = size(bvec)
    rvec = zeros(Float64, size(bvec))
    pvec = zeros(Float64, size(bvec))
    qvec = zeros(Float64, size(bvec))
    zvec = zeros(Float64, size(bvec))

    bnorm = zeros(Float64, nblock)
    rnorm = zeros(Float64, nblock)

    bnorm .= 0.0
    for iblock in 1:nblock
        bnorm[iblock] = norm(bvec[iblock, :])
    end

    # initial guess
    for idof in 1:ndof
        for iblock in 1:nblock
            uvec[iblock, idof] = bvec[iblock, idof] * amat_diag_inv[idof]
        end
    end

    # q <- A u
    spmatvec_block!(amat_val, amat_col, amat_ind, uvec, qvec)

    # r <- b - q
    rvec .= bvec .- qvec
    for iblock in 1:nblock
        rnorm[iblock] = norm(rvec[iblock, :])
    end

    alpha = zeros(Float64, nblock)
    beta = zeros(Float64, nblock)
    rho = zeros(Float64, nblock)
    iter = 1
    @show (rnorm./bnorm)[1]
    while maximum(rnorm / bnorm) > 1.0e-8
        # z <- M^-1 r
        for idof in 1:ndof
            for iblock in 1:nblock
                zvec[iblock, idof] = rvec[iblock, idof] * amat_diag_inv[idof]
            end
        end

        if iter > 1
            # beta <- (r, z) / rho
            # beta .= inv(rho) * (rvec * transpose(zvec))
            for iblock in 1:nblock
                beta[iblock] = dot(rvec[iblock, :], zvec[iblock, :]) / rho[iblock]
            end
        end

        # p <- z + beta p
        for iblock in 1:nblock
            pvec[iblock, :] = zvec[iblock, :] .+ beta[iblock] * pvec[iblock, :]
        end

        # q <- A p
        spmatvec_block!(amat_val, amat_col, amat_ind, pvec, qvec)

        # rho <- (r, z)
        for iblock in 1:nblock
            rho[iblock] = dot(rvec[iblock, :], zvec[iblock, :])
        end

        # alpha <- rho / (p, q)
        for iblock in 1:nblock
            alpha[iblock] = rho[iblock] / dot(pvec[iblock, :], qvec[iblock, :])
        end

        # q <- -alpha q
        for iblock in 1:nblock
            qvec[iblock, :] .= -alpha[iblock] .* qvec[iblock, :]
        end

        # r <- r + q
        rvec .= rvec .+ qvec

        # u <- u + alpha p
        for iblock in 1:nblock
            uvec[iblock, :] .= uvec[iblock, :] .+ alpha[iblock] .* pvec[iblock, :]
        end

        rnorm .= 0.0
        for iblock in 1:nblock
            rnorm[iblock] = norm(rvec[iblock, :])
        end
        @show (rnorm./bnorm)[1]

        iter += 1
    end
    @show "max iter: ", iter
end

# solve A u = b
function BlockConjugateGradient_diag_bcrs!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
    nblock, _, nnode = size(bvec)
    @show nblock, nnode
    rvec = zeros(Float64, size(bvec))
    pvec = zeros(Float64, size(bvec))
    qvec = zeros(Float64, size(bvec))
    zvec = zeros(Float64, size(bvec))

    bnorm = zeros(Float64, nblock)
    rnorm = zeros(Float64, nblock)

    bnorm .= 0.0
    for iblock in 1:nblock
        bnorm[iblock] = norm(bvec[iblock, :, :])
    end

    # initial guess
    for inode in 1:nnode
        for iblock in 1:nblock
            uvec[iblock, :, inode] = amat_diag_inv[:, :, inode] * bvec[iblock, :, inode]
        end
    end

    # q <- A u
    spmatvec_block_bcrs!(amat_val, amat_col, amat_ind, uvec, qvec)

    # r <- b - q
    rvec .= bvec .- qvec
    for iblock in 1:nblock
        rnorm[iblock] = norm(rvec[iblock, :, :])
    end

    alpha = zeros(Float64, nblock)
    beta = zeros(Float64, nblock)
    rho = zeros(Float64, nblock)
    iter = 1
    @show (rnorm./bnorm)[1]
    while maximum(rnorm / bnorm) > 1.0e-8
        # z <- M^-1 r
        # for idof in 1:ndof
        #     for iblock in 1:nblock
        #         zvec[iblock, idof] = rvec[iblock, idof] * amat_diag_inv[idof]
        #     end
        # end
        for inode in 1:nnode
            for iblock in 1:nblock
                zvec[iblock, :, inode] .= amat_diag_inv[:, :, inode] * rvec[iblock, :, inode] 
            end
        end

        if iter > 1
            # beta <- (r, z) / rho
            for iblock in 1:nblock
                beta[iblock] = 0.0
                for inode in 1:nnode
                    beta[iblock] += dot(rvec[iblock, :, inode], zvec[iblock, :, inode])
                end
                beta[iblock] = beta[iblock] / rho[iblock]
            end
        end

        # p <- z + beta p
        for iblock in 1:nblock
            pvec[iblock, :, :] = zvec[iblock, :, :] .+ beta[iblock] * pvec[iblock, :, :]
        end

        # q <- A p
        spmatvec_block_bcrs!(amat_val, amat_col, amat_ind, pvec, qvec)

        # rho <- (r, z)
        for iblock in 1:nblock
            rho[iblock] = 0.0
            for inode in 1:nnode
                rho[iblock] += dot(rvec[iblock, :, inode], zvec[iblock, :, inode])
            end
        end

        # alpha <- rho / (p, q)
        for iblock in 1:nblock
            alpha[iblock] = 0.0
            for inode in 1:nnode
                alpha[iblock] += dot(pvec[iblock, :, inode], qvec[iblock, :, inode])
            end
            alpha[iblock] = rho[iblock] / alpha[iblock]
        end

        # q <- -alpha q
        for iblock in 1:nblock
            qvec[iblock, :, :] .= -alpha[iblock] .* qvec[iblock, :, :]
        end

        # r <- r + q
        rvec .= rvec .+ qvec

        # u <- u + alpha p
        for iblock in 1:nblock
            uvec[iblock, :, :] .= uvec[iblock, :, :] .+ alpha[iblock] .* pvec[iblock, :, :]
        end

        rnorm .= 0.0
        for iblock in 1:nblock
            rnorm[iblock] = norm(rvec[iblock, :, :])
        end
        @show (rnorm./bnorm)[1]

        iter += 1
    end
    @show "max iter: ", iter
end
function gen_matvec_bcrs(cny, coor, ds, xmin, xmax, ymin, ymax, zmin, zmax, nblock, mblock)
    # material properties
    young = 10.0
    nyu = 0.3
    lam = young * nyu / (1.0 + nyu) / (1.0 - 2.0 * nyu)
    mu = young / 2.0 / (1.0 + nyu)

    nnode = size(coor, 2)
    nelement = size(cny, 2)
    kglobal_tmp_val = [[] for _ in 1:nnode]
    kglobal_tmp_col = [[] for _ in 1:nnode]

    for ie = 1:nelement
        if ie % 1000 == 0
            @show "element ", ie
        end
        node_id = zeros(Int32, 8)
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

                    kelement .+= transpose(bmat) * dmat * bmat .* jac
                end
            end
        end
        for inode in 1:8
            inode_g = node_id[inode]
            for jnode in 1:8
                jnode_g = node_id[jnode]
                knode = findfirst(x -> x == jnode_g, kglobal_tmp_col[inode_g])
                if isnothing(knode)
                    push!(kglobal_tmp_col[inode_g], jnode_g)
                    push!(kglobal_tmp_val[inode_g], kelement[3*(inode-1)+1:3*inode, 3*(jnode-1)+1:3*jnode])
                else
                    kglobal_tmp_val[inode_g][knode] .+= kelement[3*(inode-1)+1:3*inode, 3*(jnode-1)+1:3*jnode]
                end
            end
        end
    end

    nnz = sum(length, kglobal_tmp_col)
    kglobal_col = zeros(Int32, nnz)
    kglobal_val = zeros(Float64, 3, 3, nnz)
    kglobal_ind = zeros(Int32, nnode + 1)

    kglobal_ind[1] = 1
    for inode in 1:nnode
        kglobal_ind[inode+1] = kglobal_ind[inode] + length(kglobal_tmp_col[inode])
        kglobal_col[kglobal_ind[inode]:kglobal_ind[inode+1]-1] = kglobal_tmp_col[inode]
        for jnode in 1:length(kglobal_tmp_col[inode])
            kglobal_val[:, :, kglobal_ind[inode]+jnode-1] = kglobal_tmp_val[inode][jnode]
        end
    end

    kglobal_diag = zeros(Float64, 3, 3, nnode)
    kglobal_diag_inv = zeros(Float64, 3, 3, nnode)

    # inpulse force
    fglobal = zeros(Float64, nblock, 3, nnode)
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        for iblock in 1:nblock
            iy, ix = divrem(iblock - 1, mblock) .+ 1
            xt = (xmax - xmin) / mblock * (ix - 0.5)
            yt = (ymax - ymin) / mblock * (iy - 0.5)
            zt = zmax
            if norm([x, y, z] .- [xt, yt, zt]) < 1.0e-6
                @show "force at node ", inode
                fglobal[iblock, 3, inode] = 1.0
            end
        end
    end

    # boundary conditions
    bcflag = zeros(Bool, nnode)
    uglobal = zeros(Float64, nblock, 3, nnode)
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        if abs(z - zmin) < 1.0e-6
            bcflag[inode] = true
            for iblock in 1:nblock
                uglobal[iblock, 1, inode] = 0.0
                uglobal[iblock, 2, inode] = 0.0
                uglobal[iblock, 3, inode] = 0.0
            end
        end
    end

    for inode in 1:nnode
        if bcflag[inode]
            for jnode in 1:length(kglobal_tmp_col[inode])
                if kglobal_col[kglobal_ind[inode]+jnode-1] == inode
                    kglobal_val[:, :, kglobal_ind[inode]+jnode-1] .= 0.0
                    kglobal_val[1, 1, kglobal_ind[inode]+jnode-1] = 1.0
                    kglobal_val[2, 2, kglobal_ind[inode]+jnode-1] = 1.0
                    kglobal_val[3, 3, kglobal_ind[inode]+jnode-1] = 1.0
                    kglobal_diag[:, :, inode] .= 0.0
                    kglobal_diag[1, 1, inode] = 1.0
                    kglobal_diag[2, 2, inode] = 1.0
                    kglobal_diag[3, 3, inode] = 1.0
                else
                    kglobal_val[:, :, kglobal_ind[inode]+jnode-1] .= 0.0
                end
            end
            for iblock in 1:nblock
                fglobal[iblock, :, inode] .= uglobal[iblock, :, inode]
            end
        else
            for jnode in 1:length(kglobal_tmp_col[inode])
                if bcflag[kglobal_col[kglobal_ind[inode]+jnode-1]]
                    for iblock in 1:nblock
                        fglobal[iblock, :, inode] .-= kglobal_val[:, :, kglobal_ind[inode]+jnode-1] * uglobal[iblock, :, kglobal_col[kglobal_ind[inode]+jnode-1]]
                    end
                    kglobal_val[:, :, kglobal_ind[inode]+jnode-1] .= 0.0
                end
                if (kglobal_col[kglobal_ind[inode]+jnode-1] == inode)
                    kglobal_diag[:, :, inode] .= kglobal_val[:, :, kglobal_ind[inode]+jnode-1]
                end
            end
        end
    end

    # for idof in 1:3*nnode
    #     kglobal_diag_inv[idof] = 1.0 / kglobal_diag[idof]
    # end
    for inode in 1:nnode
        kglobal_diag_inv[:, :, inode] .= inv(kglobal_diag[:, :, inode])
    end

    return kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal
end

function gen_matvec(cny, coor, ds, xmin, xmax, ymin, ymax, zmin, zmax, nblock, mblock)
    # material properties
    young = 10.0
    nyu = 0.3
    lam = young * nyu / (1.0 + nyu) / (1.0 - 2.0 * nyu)
    mu = young / 2.0 / (1.0 + nyu)

    nnode = size(coor, 2)
    nelement = size(cny, 2)
    kglobal_tmp_val = [[] for _ in 1:3*nnode]
    kglobal_tmp_col = [[] for _ in 1:3*nnode]

    for ie = 1:nelement
        if ie % 1000 == 0
            @show "element ", ie
        end
        node_id = zeros(Int32, 8)
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

                    kelement .+= transpose(bmat) * dmat * bmat .* jac
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
    kglobal_col = zeros(Int32, nnz)
    kglobal_val = zeros(Float64, nnz)
    kglobal_ind = zeros(Int32, 3 * nnode + 1)

    kglobal_ind[1] = 1
    for idof in 1:3*nnode
        kglobal_ind[idof+1] = kglobal_ind[idof] + length(kglobal_tmp_col[idof])
        kglobal_col[kglobal_ind[idof]:kglobal_ind[idof+1]-1] = kglobal_tmp_col[idof]
        kglobal_val[kglobal_ind[idof]:kglobal_ind[idof+1]-1] = kglobal_tmp_val[idof]
    end

    kglobal_diag = zeros(Float64, 3 * nnode)
    kglobal_diag_inv = zeros(Float64, 3 * nnode)

    # inpulse force
    fglobal = zeros(Float64, nblock, 3 * nnode)
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        for iblock in 1:nblock
            iy, ix = divrem(iblock - 1, mblock) .+ 1
            # xt = 2.0 * (x - 1.0) / 2.0
            xt = (xmax - xmin) / mblock * (ix - 0.5)
            yt = (ymax - ymin) / mblock * (iy - 0.5)
            zt = zmax
            if norm([x, y, z] .- [xt, yt, zt]) < 1.0e-6
                @show "force at node ", inode
                @show xt, yt, zt
                fglobal[iblock, 3*(inode-1)+3] = 1.0
            end
        end
    end

    # boundary conditions
    bcflag = zeros(Bool, 3 * nnode)
    uglobal = zeros(Float64, nblock, 3 * nnode)
    for inode in 1:nnode
        x, y, z = coor[:, inode]
        if abs(z - zmin) < 1.0e-6
            bcflag[3*(inode-1)+1] = true
            bcflag[3*(inode-1)+2] = true
            bcflag[3*(inode-1)+3] = true
            for iblock in 1:nblock
                uglobal[iblock, 3*(inode-1)+1] = 0.0
                uglobal[iblock, 3*(inode-1)+2] = 0.0
                uglobal[iblock, 3*(inode-1)+3] = 0.0
            end
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
            for iblock in 1:nblock
                fglobal[iblock, idof] = uglobal[iblock, idof]
            end
        else
            for jdof in 1:kglobal_ind[idof+1]-kglobal_ind[idof]
                if bcflag[kglobal_col[kglobal_ind[idof]+jdof-1]]
                    for iblock in 1:nblock
                        fglobal[iblock, idof] -= kglobal_val[kglobal_ind[idof]+jdof-1] * uglobal[iblock, kglobal_col[kglobal_ind[idof]+jdof-1]]
                    end
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

    return kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal
end


# nblock for block conjugate gradient
# nblock needs to be square of an integer
mblock = 4
nblock = mblock * mblock

data_dir = "../bcrs_data/"
if !isdir(data_dir)
    mkpath(data_dir)
end

ds = 0.03125
xmin = 0.0
xmax = 4.0
ymin = 0.0
ymax = 4.0
zmin = 0.0
zmax = 4.0
cny, coor = gen_model(ds, xmin, xmax, ymin, ymax, zmin, zmax)
kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal = gen_matvec_bcrs(cny, coor, ds, xmin, xmax, ymin, ymax, zmin, zmax, nblock, mblock)
# BlockConjugateGradient_diag_bcrs!(kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
nblock, _, nnode = size(uglobal)
nnz = length(kglobal_col)
write_data_bcrs(data_dir, nnode, nnz, nblock, cny, coor, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)


# kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal = gen_matvec(cny, coor, ds, xmin, xmax, ymin, ymax, zmin, zmax, nblock, mblock)
# BlockConjugateGradient_diag!(kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)
# nblock, ndof = size(uglobal)
# nnz = length(kglobal_col)
# write_data(data_dir, ndof, nnz, nblock, cny, coor, kglobal_val, kglobal_col, kglobal_ind, kglobal_diag_inv, uglobal, fglobal)

# write_mesh(data_dir, nblock)


# # 正定値行列を作成するための関数
# function generate_positive_definite_matrix(ndof)
#     A = randn(ndof, ndof)  # ランダムな行列を生成
#     return A' * A           # A^T * A で正定値行列を作成
# end

# # 正定値行列をCRS形式に変換する関数
# function convert_to_crs(matrix)
#     val = []
#     col = []
#     ind = [1]  # CRSの行インデックス（1-indexed）

#     for i in 1:size(matrix, 1)
#         for j in 1:size(matrix, 2)
#             if matrix[i, j] != 0
#                 push!(val, matrix[i, j])
#                 push!(col, j)
#             end
#         end
#         push!(ind, length(val) + 1)
#     end

#     return val, col, ind
# end

# # テスト行列の作成
# Random.seed!(0)
# ndof = 100   # 自由度（行列のサイズ）
# nblock = 16 # ブロック数を16に設定

# # 正定値行列の生成
# A = generate_positive_definite_matrix(ndof)
# amat_val, amat_col, amat_ind = convert_to_crs(A)

# # 対角成分の逆数を計算（前処理用）
# amat_diag_inv = 1.0 ./ diag(A)  # 対角成分を抽出して逆数を取る

# # 右辺ベクトルを作成
# # ここでブロック数に合わせてベクトルを初期化
# bvec = randn(nblock, ndof)  # 16x4 の右辺ベクトルをランダムに生成

# # 解ベクトルの初期化
# uvec = zeros(Float64, nblock, ndof)

# uvec = uvec[1:16, :]
# bvec = bvec[1:16, :]
# # BlockConjugateGradient! 関数をテスト
# BlockConjugateGradient!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)
# uvec .= 0.0
# BlockConjugateGradient_diag!(amat_val, amat_col, amat_ind, amat_diag_inv, uvec, bvec)


# # --- uvec が正しいか確認するコードを追加 ---
# # A * uvec を計算し、bvec に近いか確認
# function verify_solution(A, uvec, bvec)
#     nblock = size(bvec, 1)
#     ndof = size(bvec, 2)

#     # A * uvec を計算
#     b_approx = uvec * A

#     # 各ブロックの誤差を計算
#     for iblock in 1:nblock
#         error = norm(b_approx[iblock, :] - bvec[iblock, :])
#         println("ブロック $iblock の誤差: ", error)
#     end
# end

# # 解の検証
# verify_solution(A, uvec, bvec)
