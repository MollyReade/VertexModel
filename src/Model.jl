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

# Local modules
@from "SpatialData.jl" using SpatialData

function model!(du, u, p, t)

    params, matrices = p
    @unpack A,
        B,
        Ā,
        B̄,
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
        vertexAreas = matrices
    @unpack nVerts,
        nCells,
        nEdges,
        pressureExternal,
        boundaryToggle,
        peripheralTension,
        dissipationToggle,
        edgeDissToggle,
        vertexDissToggle,
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

    

    if startswith(energyModel, "ventilation")
        for k = 1:nVerts
            for j in nzrange(A, k)
                for i in nzrange(B, rowvals(A)[j])
                    # Force components from cell pressure perpendicular to edge tangents 
                    F[k, rowvals(B)[i]] -= 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                    externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                end
                #Force component from edge tension on vertex k
                FEdges[k] -= edgeTensions[rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from peripheral tension
                externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            
            end
            dR[k] = (sum(@view F[k, :]) .+ externalF[k] .+ FEdges[k])
        end
    else
        for k = 1:nVerts
            for j in nzrange(A, k)
                for i in nzrange(B, rowvals(A)[j])
                    # Force components from cell pressure perpendicular to edge tangents 
                    F[k, rowvals(B)[i]] -= 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                    # Force components from cell membrane tension parallel to edge tangents 
                    F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                    #F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                    # Force on vertex from external pressure 
                    externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                    #externalF[k] -= peripheralTension
                end
                # Force on vertex from peripheral tension
                externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                #externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            end
            
            dR[k] = (sum(@view F[k, :]) .+ externalF[k])
        end
    end

    if boundaryToggle == 1
        for k = 1:nVerts
            if boundaryVertices[k] != 1
                dR[k] = @SVector zeros(2)
            end
        end
    end

    if startswith(energyModel,"ventilation_rational") && dissipationToggle == 1
        # Rational vertex weighting
        R = calculateDrag( params, matrices)
        Fvec = reinterpret(Float64, dR)
        v = R \ Fvec
        dR .= reinterpret(SVector{2,Float64}, v)
    else
        dissipationToggle == 1 ? dR ./= vertexAreas : nothing 
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
