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

function modelT!(du, u, p, t)
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
        dR[k] = @SVector zeros(2)
        for j = 1:nEdges
            dR[k] += - A[j, k] * edgeTensions[j] * (normEdges[j])
            for i = 1:nCells
                #dR[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])
            end
            
        end
        
    end
    
    #preMult = κ * Diagonal(vertexAreas) + areaJacobian' * Λ * Diagonal(cellαᵢs) * areaJacobian
    preMult = zeros(Float64, nVerts, nVerts)
    for k = 1:nVerts
        preMult[k,k] += κ * vertexAreas[k] 
        preMult[k,k] = κ
    end
    
    #println(
    #    "t = ", t,
    #    ", min vertex area = ", minimum(vertexAreas),
    #    ", max vertex area = ", maximum(vertexAreas)
    #)
    dR .= preMult \ dR

    if boundaryToggle == 1
        for k = 1:nVerts
            if boundaryVertices[k] != 1
                dR[k] = @SVector zeros(2)
            end
        end
    end

    for i = 1:nCells
        cellPressures[i] = 0
        for k = 1:nVerts
            cellPressures[i] += dot(areaJacobian[i,k], (vertexAreas[k]^-1 * sum(A[j,k]* edgeTensions[j] * normEdges[j] for j = 1:nEdges)))  
            cellPressures[i] -=  dot(areaJacobian[i,k], (vertexAreas[k]^-1 * (Pᵢₚ - Pₘ) * sum(sum(0.5*edgeCellNormals[i,j]*abs(A[j,k]) for j = 1:nEdges) for i = 1:nCells)))
            cellPressures[i] +=  κ * cellαᵢs[i]^-1 * Λ^-1 * Pₘ
        end
        cellPressures[i] = cellPressures[i]/(κ*cellαᵢs[i]^-1 * Λ^-1 + sum(dot(areaJacobian[i,k], (vertexAreas[k]^-1 * sum(areaJacobian[i,k] for i = 1:nCells))) for k = 1:nVerts))
    end
    

    # dR accesses the same underlying data as du, so by altering dR we have already updated du appropriately
    return du
end

function modelP!(du,u,p,t)
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
    if ipModel == "sinusoidal"
        Pᵢₚ = P₀ + Amp * sin(t)
    else
        Pᵢₚ = P₀
    end
    R = reinterpret(SVector{2,Float64}, u)
    dR = reinterpret(SVector{2,Float64}, du)

    spatialData!(R, params, matrices)
    #cellPressures .= 0
    for k = 1:nVerts
        dR[k] = @SVector zeros(2)
        for i = 1:nCells
            dR[k] += areaJacobian[i,k] * cellPressures[i]
        end
    end
    for k = 1:nVerts
        for j = 1:nEdges
            for i = 1:nCells
                dR[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])
            end
        end
    end
    
    #preMult = κ * Diagonal(vertexAreas) + areaJacobian' * Λ * Diagonal(cellαᵢs) * areaJacobian
    preMult = zeros(Float64, nVerts, nVerts)
    for k = 1:nVerts
        preMult[k,k] += κ * vertexAreas[k] 
        preMult[k,k] = κ
    end

    # println(
    #     "t = ", t,
    #     ", min vertex area = ", minimum(vertexAreas),
    #     ", max vertex area = ", maximum(vertexAreas)
    # )
    
    dR .= preMult \ dR #problematic line

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
        dR[k] = @SVector zeros(2)
        for j = 1:nEdges
            dR[k] += - A[j, k] * edgeTensions[j] * (normEdges[j])
            for i = 1:nCells
                dR[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])
            end
            
        end
        for j in nzrange(A,k)
            externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j],k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
        end
        #dR[k] += externalF[k]
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
        
        dR .= preMult \ dR
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

function gradₐ(areaJacobian, vertexAreas, P)
    _, nVerts = size(areaJacobian)
    return [sum(areaJacobian[i,k] *P[i] for i in axes(areaJacobian, 1)) / vertexAreas[k] for k in 1:nVerts]
end

function divₐ(areaJacobian, cellAreas, V)
    nCells, _ = size(areaJacobian)
    return [-sum(dot(areaJacobian[i,k], V[k]) for k in axes(areaJacobian, 2)) / cellAreas[i] for i in 1:nCells]
end


function tensionComponent!(u,p,t)
     params, matrices = p
    @unpack A,
         B,
         edgeTensions,
         boundaryVertices,
         edgeCellNormals,
         normEdges,
         vertexAreas,
         dR⁰,
         areaJacobian,
         cellPressures = matrices
    @unpack nVerts,
         nCells,
         nEdges,
         Pₘ,
         κ,
         P₀,
         Λ,
         Amp,
         ipModel,
         boundaryCondition = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    fill!(dR⁰, @SVector zeros(2))

    spatialData!(R, params, matrices) 

    if boundaryCondition == "force"
        if ipModel == "sinusoidal"
            Pᵢₚ = P₀ + Amp * sin(t)
        else
            Pᵢₚ = P₀
        end
    else
        Pᵢₚ = 0 # Set pleural pressure to zero if using displacement boundary condition
    end

    #----------------- THESE TWO LOOPS PRODUCE VERY DIFFERENT RESULTS
    # for k = 1:nVerts
    #     for j = 1:nEdges
    #         dR⁰[k] += - A[j, k] * edgeTensions[j] * (normEdges[j]) # Tension contribution
    #         for i = 1:nCells
    #             #TODO Can i replace this triple loop?
    #             #TODO Replace abs A with already existing abs a
    #             #Bottleneck line
    #             dR⁰[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])  # Contribution from prescibed pressure at mouth and pleural pressure
    #         end
            
    #     end
        
    # end

    for k = 1:nVerts
        for p in A.colptr[k]:(A.colptr[k+1]-1)
            j = A.rowval[p]
            Aⱼₖ = A.nzval[p]
            dR⁰[k]+= - Aⱼₖ * edgeTensions[j] * normEdges[j]
            for q in B.colptr[j]:(B.colptr[j+1]-1)
                i = B.rowval[q]
                dR⁰[k]+= (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j]
            end
        end
    end

    # for k in 1:nVerts
    #     edges = findnz(@view A[:,k])[1]
    #     for a in 1:length(edges)
    #         j = edges[a]
    #         dR⁰[k] += - A[j,k] * edgeTensions[j] * normEdges[j]
    #         cells = findnz(@view B[:,j])[1]
    #         for b in 1:length(cells)
    #             i = cells[b]
    #             dR⁰[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j]
    #         end
    #     end
    # end
    # ------------------------------
    for k = 1:nVerts
        for i = 1:nCells
            dR⁰[k] += areaJacobian[i,k] * cellPressures[i] # Force due to internal presuures
        end
    end
    
    if boundaryCondition == "displacement"
        for k in 1:nVerts
            if boundaryVertices[k] == 1
                dR⁰[k] = @SVector zeros(2)
            end
        end
    end
    
    dR⁰ ./= κ .* vertexAreas

    
end

function inversePressure!(p,t)
    params, matrices = p
    @unpack A,
        areaJacobian,
        vertexAreas,
        cellAreas,
        edgeCellNormals,
        edgeTensions,
        normEdges,
        cellαᵢs,
        cellPressures,
        R⁺,
        ∂𝒜∂r,
        T,
        laplacian = matrices
    @unpack κ,
        nCells,
        nEdges,
        nVerts,
        P₀,
        Amp,
        ipModel,
        Pₘ,
        boundaryCondition = params

    spatialData!(R⁺, params, matrices)

    gradₗ = A'
    
    if boundaryCondition == "force"
        if ipModel == "sinusoidal" 
            Pᵢₚ = P₀ + Amp * sin(t)
        else
            Pᵢₚ = P₀
        end
    else
        Pᵢₚ = 0
    end

    fill!(laplacian,0)
    dropzeros!(laplacian)

    # ------ THESE TWO LOOPS DONT MATCH EITHER

    # for i in 1:nCells
    #     for m in 1:nCells
    #         for k in 1:nVerts
    #             laplacian[i,m] += dot(areaJacobian[i,k],
    #                         areaJacobian[m,k]) / vertexAreas[k]    # Define the laplacian
    #                         #TODO Main bottleneck now
    #         end
    #     end
    # end

    for k = 1:nVerts
        invDₖ = 1/vertexAreas[k]
        rng = areaJacobian.colptr[k]:(areaJacobian.colptr[k+1]-1)
        n = length(rng)
        for a in 1:n
            pa = rng[a]
            i = areaJacobian.rowval[pa]
            vi = areaJacobian.nzval[pa]
            for b in a:n
                pb = rng[b]
                m=areaJacobian.rowval[pb]
                vm = areaJacobian.nzval[pb]
                Lim = dot(vi, vm) * invDₖ
                laplacian[i,m] += Lim
                if i != m 
                    laplacian[m,i] += Lim
                end
            end
        end
    end

    # for k in 1:nVerts
    #     cells = findnz(@view areaJacobian[:,k])[1]
    #     for a in 1:length(cells)
    #         i = cells[a]
    #         for b in a:length(cells)
    #             m = cells[b]
    #             Lᵢₘ = dot(areaJacobian[i,k],areaJacobian[m,k]) / vertexAreas[k]
    #             laplacian[i,m] += Lᵢₘ
    #             if i != m
    #                 laplacian[m,i] += Lᵢₘ
    #             end
    #         end
    #     end
    # end

    # --------------------------------------------

    # TODO swap out Diagonals?
    RHS = -spdiagm(cellAreas)*divₐ(areaJacobian, cellAreas,(spdiagm(1 ./ vertexAreas)*gradₗ*T))
    RHS += κ*spdiagm(1 ./ cellαᵢs)*Pₘ*ones(nCells)
    RHS += spdiagm(cellAreas)*divₐ(areaJacobian,cellAreas,spdiagm(1 ./ vertexAreas)*(Pᵢₚ - Pₘ)*∂𝒜∂r)
    cellPressures .= (κ*spdiagm(1 ./ cellαᵢs) + laplacian) \ (RHS)   # Invert for cell pressures

end

function pressureComponent!(p,t,oldJacobian,oldPressures)
     params, matrices = p
    @unpack areaJacobian,
        cellPressures,
        boundaryVertices,
        vertexAreas,
        R⁺,
        dR⁺ = matrices
    @unpack nVerts,
        nCells,
        boundaryCondition,
        κ = params

    # Reinterpret state vector as a vector of SVectors 
    fill!(dR⁺, @SVector zeros(2))

    for k = 1:nVerts
        for i = 1:nCells
            dR⁺[k] += areaJacobian[i,k] * cellPressures[i] -oldJacobian[i,k] * oldPressures[i] # Force due to internal presuures
        end
    end

    if boundaryCondition == "displacement"
        for k in 1:nVerts
            if boundaryVertices[k] == 1
                dR⁺[k] = @SVector zeros(2) # If the boundary is prescibed, no force balance on it
            end
        end
    end
    
    dR⁺ ./= κ .* vertexAreas # Dissipation

end

function splitStep(r⁰, u0, p,t, Δt)
    params, matrices = p
    @unpack A, boundaryVertices, avgEdgeCellNormals, cellPressures, areaJacobian, cellαᵢs, dR⁰, R⁺, dR⁺, R¹, cellAreas = matrices
    @unpack boundaryCondition, Amp, ω, Pₘ, nCells, nVerts = params

    initialR = reinterpret(SVector{2,Float64}, u0)
    R⁰ = reinterpret(SVector{2,Float64}, r⁰)

    tensionComponent!(r⁰, p, t) # Updates dR0 using tension terms, Pm and Pip
    #println("a₀: ", cellAreas)
    oldJacobian = copy(areaJacobian)
    oldPressures = copy(cellPressures)

    # Tension Step
    R⁺ .= R⁰ + Δt * dR⁰ # Calculate intermediate position
    if boundaryCondition == "displacement"
        for k in 1:length(R⁺)
            if boundaryVertices[k] == 1
                R⁺[k] = initialR[k] .+ (Amp/ω)*(1-cos(ω*(t+Δt)))*avgEdgeCellNormals[k] # Prescibe boundary
            end
        end
    end

    # Intermediate Pressure

    inversePressure!(p, t) # Updates cell pressures using R+

    #println("a⁺: ", cellAreas)

    # Pressure step

    pressureComponent!(p,t,oldJacobian,oldPressures) # Updates dR+ using internal pressure term

    R¹ .= R⁺ .+ Δt * dR⁺ # Calculate final position

    if boundaryCondition == "displacement"
        for k in 1:length(R⁺)
            if boundaryVertices[k] == 1
                R¹[k] = initialR[k] .+ (Amp/ω)*(1-cos(ω*(t+2*Δt)))*avgEdgeCellNormals[k] # Prescibe boundary
            end
        end
    end

    
    daᵢdt = zeros(nCells)
    for i in 1:nCells
        daᵢdt[i] = sum(dot(areaJacobian[i,k],dR⁺[k]+dR⁰[k]) for k in 1:nVerts)
    end

    #println("Pressure check: ",cellPressures .- (Pₘ*ones(nCells) .- Diagonal(cellαᵢs)*daᵢdt) )

    return reinterpret(Float64, R¹)
end



export modelP!, modelT!, model!, splitStep

end
