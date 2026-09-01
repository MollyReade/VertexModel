#
#  Model.jl
#  VertexModel
#
#  Function to calculate force vector on vertex k from cell i (Fᵢₖ) for all vertices and update vertex positions.

module Model

# Julia packages
using LinearAlgebra
using StaticArrays
using UnPack
using SparseArrays
using .Threads
using FromFile 
using DrWatson
using FastBroadcast

# Local modules
@from "SpatialData.jl" using SpatialData

function model!(du, u, p, t)

    params, matrices = p
    @unpack A,
        B,
        Ā,
        B̄,
        areaJacobian,
        cellαᵢs,
        cellTensions,
        cellPressures,
        edgeLengths,
        edgeTangents,
        edgeTensions,
        F,
        FEdges,
        externalF,
        ϵ,
        boundaryVertices,
        boundaryEdges,
        edgeCellNormals,
        normEdges,
        vertexAreas = matrices
    @unpack nVerts,
        nCells,
        nEdges,
        Pₘ,
        κ,
        P₀,
        Λ,
        Amp,
        ipModel,
        boundaryToggle,
        peripheralTension,
        energyModel = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    dR = reinterpret(SVector{2,Float64}, du)

    spatialData!(R, params, matrices)

    fill!(F, @SVector zeros(2))
    dropzeros!(F)
    fill!(FEdges, @SVector zeros(2))
    dropzeros!(FEdges)
    fill!(externalF, @SVector zeros(2))

    peripheryLength = sum(boundaryEdges .* edgeLengths)

    if ipModel == "sinusoidal"
        Pᵢₚ = P₀ + Amp * sin(t)
    else
        Pᵢₚ = P₀
    end

    for k = 1:nVerts
        for j = 1:nEdges
            dR[k] += - A[j, k] * edgeTensions[j] * (normEdges[j])
            for i = 1:nCells
                dR[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])
            end
            
        end
        for j in nzrange(A,k)
            externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j],k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
        end
        dR[k] += externalF[k]
    end
    for i = 1:nCells
        for k = 1:nVerts
            dR[k] += areaJacobian[i,k] * Pₘ
        end
    end
    #preMult = κ * Diagonal(vertexAreas) + areaJacobian' * Λ * Diagonal(cellαᵢs) * areaJacobian
    preMult = zeros(Float64, nVerts, nVerts)
    for k = 1:nVerts
        preMult[k,k] += κ * vertexAreas[k] 
        for l = 1:nVerts
            for i = 1:nCells
                preMult[k,l] += Λ * cellαᵢs[i] * (areaJacobian[i,k]' * areaJacobian[i,l])
            end
        end
    end
    
    dR = inv(preMult) * dR
    cellPressures .= Pₘ * ones(nCells) - Λ * Diagonal(cellαᵢs) * [sum(dot(areaJacobian[i,k], dR[k]) for k in axes(areaJacobian, 2)) for i in axes(areaJacobian, 1)]

    if boundaryToggle == 1
        for k = 1:nVerts
            if boundaryVertices[k] != 1
                dR[k] = @SVector zeros(2)
            end
        end
    end

    
    # dR accesses the same underlying data as du, so by altering dR we have already updated du appropriately
    return du
end

function calculateDrag(params, matrices)

    @unpack A,
        edgeTangents,
        edgeLengths,
        vertexAreas = matrices
    @unpack nVerts,
        nEdges,
        edgeDissToggle,
        vertexDissToggle = params

    μ = 1.0
    I2 = Matrix{Float64}(I, 2, 2)  # 2×2 identity
    W_blocks = [ sparse(Matrix(t*t')) / (norm(t)^2) for t in edgeTangents ]
    W = blockdiag(W_blocks...)
    R_edges = μ * kron(A', I2) * W * kron(A, I2)  # Drag matrix in vertex space
    D = Diagonal(vertexAreas)
    R_verts = kron(D, I2)  # Drag matrix in vertex space
    R = zeros(size(R_edges))
    if edgeDissToggle == 1
        R += R_edges
    end
    if vertexDissToggle == 1
        R += R_verts
    end

    if any(isnan, R) || any(isinf, R)
        @show "NaN or Inf in R"
        @show edgeTangents
        @show norm.(edgeTangents)
        @show vertexAreas
        error("R contains NaNs or Infs")
    end
    return R
end


export model!

end
