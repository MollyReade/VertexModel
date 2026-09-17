#
#  SpatialData.jl
#  VertexModel
#
#  Function to calculate spatial data including tangents, lengths, midpoints, tensions etc from incidence and vertex position matrices.

module SpatialData

# Julia packages
using FromFile
using LinearAlgebra
using StaticArrays
using UnPack
using FastBroadcast
using SparseArrays
using GeometryBasics

@from "OrderAroundCell.jl" using OrderAroundCell
# @from "AnalysisFunctions.jl" using AnalysisFunctions

function spatialData!(R,params,matrices)

    @unpack A,
        B,
        Ā,
        Āᵀ,
        B̄,
        Bᵀ,
        C,
        ϵ,
        areaJacobian,
        boundaryEdges,
        cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        cellPerimeters,
        cellOrientedAreas,
        cellShapeTensor,
        cellAreas,
        cellL₀s,
        cellA₀s,
        cellTensions,
        cellPressures,
        edgeCellNormals,
        avgEdgeCellNormals,
        edgeLengths,
        edgeTangents,
        edgeTensions,
        edgeMidpoints,
        edgeMidpointLinks,
        normEdges,
        vertexAreas,
        ∂𝒜∂r,
        T = matrices
    @unpack nCells,
        nEdges,
        nVerts,
        energyModel,
        l₀,
        currentTime,
        peripheralTension = params

    cellPositions  .= C*R./cellEdgeCount
    
    edgeTangents   .= A*R
    
    @.. thread=false edgeLengths .= norm.(edgeTangents)

    normEdges .= edgeTangents ./ edgeLengths

    fill!(edgeCellNormals,  SVector{2,Float64}(zeros(2)))
    dropzeros!(edgeCellNormals)
    # ----------- NOT EQUIVALENT
    for j in 1:nEdges
        for p in B.colptr[j]:(B.colptr[j+1]-1)
            i=B.rowval[p]
            edgeCellNormals[i,j] = -ϵ * B.nzval[p] * edgeTangents[j]
        end
    end

    # for i = 1:nCells
    #     for j = 1:nEdges
    #         edgeCellNormals[i, j] = - ϵ * B[i,j] * edgeTangents[j]
    #     end
    # end
    # ---------------

    fill!(avgEdgeCellNormals, SVector{2,Float64}(zeros(2)))

    for k in 1:nVerts
        n = @SVector zeros(2)

        for i in findnz(@view C[:,k])[1]
            for j in cellEdgeOrders[i]
                if A[j,k] != 0
                    n +=  edgeCellNormals[i,j]
                end
            end
        end

        avgEdgeCellNormals[k] = n/ norm(n)
    end
    # for k = 1:nVerts
    #     for i = 1:nCells
    #         for j = 1:nEdges
    #             #FIXME Heavy line
    #             avgEdgeCellNormals[k] += abs(A[j,k]) * edgeCellNormals[i,j] * C[i,k]
    #         end
    #     end
    #     avgEdgeCellNormals[k] = avgEdgeCellNormals[k] ./ norm(avgEdgeCellNormals[k])
    # end


    fill!(areaJacobian, SVector{2,Float64}(zeros(2)))
    dropzeros!(areaJacobian)
    for i in 1:nCells
        for j in cellEdgeOrders[i]
            normal = 0.5 * edgeCellNormals[i,j]
            for k in findnz(@view Āᵀ[:,j])[1]
                areaJacobian[i,k] += normal * Ā[j,k]
            end
        end
    end

    # for  i = 1:nCells
    #     for k = 1:nVerts
    #         for j = 1:nEdges
    #             #FIXME Heaviest line, make faster
    #             areaJacobian[i, k] += 0.5 * edgeCellNormals[i,j] * Ā[j,k]
    #         end
    #     end
    # end

    fill!(∂𝒜∂r,  SVector{2,Float64}(zeros(2)))

    for i in 1:nCells
        for k in 1:nVerts
            ∂𝒜∂r[k] += areaJacobian[i,k]
        end
    end
    # for k = 1:nVerts
    #     ∂𝒜∂r[k] = @SVector zeros(2)
    #     for i = 1:nCells
    #         for j = 1:nEdges
    #             #FIXME Heavy line
    #             ∂𝒜∂r[k] += 0.5* edgeCellNormals[i,j] * Ā[j,k]
    #         end
    #     end
    # end

    edgeMidpoints  .= 0.5.*Ā*R
    
    fill!(edgeMidpointLinks, SVector{2,Float64}(zeros(2)))
    dropzeros!(edgeMidpointLinks)
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1], nzC[2])
    for (i, k) in ikPairs
        for j in cellEdgeOrders[i]
            edgeMidpointLinks[i, k] = edgeMidpointLinks[i, k] .+ 0.5 .* B[i, j] .* edgeTangents[j] .* Ā[j, k]
        end
    end

    # Find vertex areas, with special consideration of peripheral vertices with 1 or 2 adjacent cells
    for k = 1:nVerts
        k_is = findall(x -> x != 0, @view C[:, k])
        if length(k_is) == 1
            k_js = findall(x -> x != 0, A[:, k])
            vertexAreas[k] = 0.5^3 * norm([edgeTangents[k_js[1]]..., 0.0] × [edgeTangents[k_js[2]]..., 0.0])
        elseif length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x -> x != 0, B[k_is[1], :] .* A[:, k])
            vertexAreas[k] = 0.5^3 * norm([edgeTangents[edgesSharedBy_i1_And_k[1]]..., 0.0] × [edgeTangents[edgesSharedBy_i1_And_k[2]]..., 0.0])
            edgesSharedBy_i2_And_k = findall(x -> x != 0, B[k_is[2], :] .* A[:, k])
            vertexAreas[k] += 0.5^3 * norm([edgeTangents[edgesSharedBy_i2_And_k[1]]..., 0.0] × [edgeTangents[edgesSharedBy_i2_And_k[2]]..., 0.0])
        else
            vertexAreas[k] = 0.5 * norm([edgeMidpointLinks[k_is[1], k]..., 0.0] × [edgeMidpointLinks[k_is[2], k]..., 0.0])
        end
    end

    cellPerimeters .= B̄ * edgeLengths

    # Find cell areas and shape tensors 
    for i = 1:nCells
        cellAreas[i] = abs(area(Point{2,Float64}.(R[cellVertexOrders[i]])))

        Rα = [R[kk].-matrices.cellPositions[i] for kk in cellVertexOrders[i]]
        cellShapeTensor[i] = sum(Rα.*transpose.(Rα))./cellEdgeCount[i]
    end

    # Calculate cell pressures and tensions according to energy model choice 
    if energyModel == "ventilation"
        # Ventilation energy model
        # Calculate cell boundary tensions
        @.. thread = false edgeTensions .=  (edgeLengths .- 1)   
        # Calculate cell internal pressures
        #@.. thread = false cellPressures .= 0.0
    elseif energyModel == "ventilation_rational"
        # Ventilation energy model with rational tension law
        #Aα^n + Bα^-m + C
        #A = 1/n(n+m), B = 1/m(n+m), C = -1/mn
        # Calculate cell boundary tensions
        @.. thread = false edgeTensions .=  ((edgeLengths/l₀).^2 .- (edgeLengths/l₀).^(-3)) ./ 5  
        # Calculate cell internal pressures
        #@.. thread = false cellPressures .= cellPᵢs
    else
        #something here
    end

    edgeTensions[boundaryEdges .==1 ] .+= peripheralTension

    fill!(T,@SVector zeros(2))
    for j = 1:nEdges
        T[j] = edgeTensions[j] * normEdges[j]
    end

    return nothing

end


export spatialData!

end


# Calculate oriented cell areas
# fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
# for i=1:nCells
#     for j in nzrange(Bᵀ,i)
#         cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
#     end
#     cellAreas[i] = cellOrientedAreas[i][1,2]
# end
