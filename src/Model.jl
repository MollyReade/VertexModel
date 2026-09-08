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


function tensionComponent(u,p,t)
     params, matrices = p
    @unpack A,
         edgeTensions,
         boundaryVertices,
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
         boundaryCondition = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    dR = spzeros(SVector{2,Float64}, nVerts)

    spatialData!(R, params, matrices)

    if boundaryCondition == "force"
        if ipModel == "sinusoidal"
            Pᵢₚ = P₀ + Amp * sin(t)
        else
            Pᵢₚ = P₀
        end
    else
        Pᵢₚ =0
    end

    for k = 1:nVerts
        dR[k] = @SVector zeros(2)
        for j = 1:nEdges
            dR[k] += - A[j, k] * edgeTensions[j] * (normEdges[j])
            for i = 1:nCells
                dR[k] += (Pᵢₚ - Pₘ) * 0.5 * edgeCellNormals[i,j] * abs(A[j,k])
            end
            
        end
        
    end
    
    preMult = zeros(Float64, nVerts, nVerts)
    for k = 1:nVerts
        preMult[k,k] += κ * vertexAreas[k] 
    end
    if boundaryCondition == "displacement"
        for k in 1:nVerts
            if boundaryVertices[k] == 1
                dR[k] = @SVector zeros(2)
            end
        end
    end
    
    dR .= preMult \ dR

    return dR
end

function inversePressure!(R,p,t)
    params, matrices = p
    @unpack A,
        areaJacobian,
        vertexAreas,
        cellAreas,
        edgeCellNormals,
        edgeTensions,
        normEdges,
        cellαᵢs,
        cellPressures = matrices
    @unpack κ,
        nCells,
        nEdges,
        nVerts,
        P₀,
        Amp,
        ipModel,
        Pₘ,
        boundaryCondition = params

    spatialData!(R, params, matrices)

    #gradₐ = Diagonal(vertexAreas)^(-1) * areaJacobian'
    #divₐ = - Diagonal(cellAreas)^(-1) * areaJacobian
    gradₗ = A'
    ∂𝒜∂r = spzeros(SVector{2,Float64}, nVerts, 1)
    T = spzeros(SVector{2,Float64}, nEdges, 1)

    if boundaryCondition == "force"
        if ipModel == "sinusoidal" 
            Pᵢₚ = P₀ + Amp * sin(t)
        else
            Pᵢₚ = P₀
        end
    else
        Pᵢₚ = 0
    end

    for k = 1:nVerts
        ∂𝒜∂r[k] = @SVector zeros(2)
        for i = 1:nCells
            for j = 1:nEdges
                ∂𝒜∂r[k] += 0.5* edgeCellNormals[i,j] * abs(A[j,k])
            end
        end
    end

    for j = 1:nEdges
        T[j] = edgeTensions[j] * normEdges[j]
    end

    L = zeros(Float64, nCells, nCells)

    for i in 1:nCells
        for m in 1:nCells
            for k in 1:nVerts
                L[i,m] += dot(areaJacobian[i,k],
                            areaJacobian[m,k]) / vertexAreas[k]
            end
        end
    end

    RHS = -Diagonal(cellAreas)*divₐ(areaJacobian, cellAreas,(Diagonal(vertexAreas)^(-1)*gradₗ*T))
    RHS += κ*Diagonal(cellαᵢs)^(-1)*Pₘ*ones(nCells)
    RHS += Diagonal(cellAreas)*divₐ(areaJacobian,cellAreas,Diagonal(vertexAreas)^(-1)*(Pᵢₚ - Pₘ)*∂𝒜∂r)
    cellPressures .= (κ*Diagonal(cellαᵢs)^(-1) + L) \ (RHS)

    # println("L range           : ", extrema(L))
    # println("cond(M)            : ", cond(Matrix(
    #     κ*Diagonal(1.0 ./ cellαᵢs) + L
    # )))
    # println("RHS norm           : ", norm(RHS))
    # println("pressure range     : ", extrema(cellPressures))

end

function pressureComponent(u,p,t)
     params, matrices = p
    @unpack areaJacobian,
        cellPressures,
        boundaryVertices,
        vertexAreas = matrices
    @unpack nVerts,
        nCells,
        boundaryCondition,
        κ = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    dR = spzeros(SVector{2,Float64}, nVerts)

    spatialData!(R, params, matrices)

    for k = 1:nVerts
        dR[k] = @SVector zeros(2)
        for i = 1:nCells
            dR[k] += areaJacobian[i,k] * cellPressures[i]
        end
    end
    
    preMult = zeros(Float64, nVerts, nVerts)
    for k = 1:nVerts
        preMult[k,k] += κ * vertexAreas[k] 
    end

    if boundaryCondition == "displacement"
        for k in 1:nVerts
            if boundaryVertices[k] == 1
                dR[k] = @SVector zeros(2)
            end
        end
    end
    
    dR .= preMult \ dR 

    return dR
end

function splitStep(r⁰, u0, p,t, Δt)
    params, matrices = p
    @unpack A, boundaryVertices, avgEdgeCellNormals, cellPressures, areaJacobian, cellαᵢs = matrices
    @unpack boundaryCondition, Amp, ω, Pₘ, nCells, nVerts = params

    initialR = reinterpret(SVector{2,Float64}, u0)
    R⁰ = reinterpret(SVector{2,Float64}, r⁰)

    dR⁰ = tensionComponent(r⁰, p, t)

    # Tension Step
    R⁺ = R⁰ + Δt *  dR⁰
    if boundaryCondition == "displacement"
        for k in 1:length(R⁺)
            if boundaryVertices[k] == 1
                R⁺[k] = initialR[k] .+ (Amp/ω)*(1-cos(ω*(t+Δt)))*avgEdgeCellNormals[k]
            end
        end
    end

    # Intermediate Pressure

    inversePressure!(R⁺, p, t)

    # Pressure step

    dR⁺ = pressureComponent(R⁺,p,t)

    R¹ = R⁺ .+ Δt * dR⁺

    if boundaryCondition == "displacement"
        for k in 1:length(R⁺)
            if boundaryVertices[k] == 1
                R¹[k] = initialR[k] .+ (Amp/ω)*(1-cos(ω*(t+2*Δt)))*avgEdgeCellNormals[k]
            end
        end
    end

    
    # daᵢdt = zeros(nCells)
    # for i in 1:nCells
    #     daᵢdt[i] = sum(dot(areaJacobian[i,k],dR⁺[k]+dR⁰[k]) for k in 1:nVerts)
    # end

    # println("Pressure check: ",cellPressures .- (Pₘ*ones(nCells) .- Diagonal(cellαᵢs)*daᵢdt) )

    return reinterpret(Float64, R¹)
end



export modelP!, modelT!, model!, splitStep

end
