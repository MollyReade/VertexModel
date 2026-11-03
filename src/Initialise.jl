#
#  Initialise.jl
#  VertexModel
#
#  Function to initialise vertex model system matrices and derived parameters

module Initialise

# Julia packages
using SparseArrays
using StaticArrays
using JLD2
using UnPack
using FromFile
using DrWatson
using Random
using Distributions
using Dates
using CircularArrays
using LinearAlgebra

# Local modules
@from "initialSystemLayout.jl" using InitialSystemLayout
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

function initialise(; initialSystem = "new",
        generationMode = "squares", # "squares" or "Matern" or "Gibbs"
        nCycles = 1,
        realCycleTime = 86400.0,
        realTimetMax = nCycles*realCycleTime,
        κ = 0.2,
        Λ = 0.5,
        γ = 0.2,
        L₀ = 0.75,
        l₀ = 0.15,
        A₀ = 1.0,
        Pᵢ = 0.2,
        Pₘ = 0.0,
        P₀ = 0.0,
        ω = π/(2*realCycleTime),
        Amp = 3,
        ipModel = "sinusoidal",
        viscousTimeScale = 1000.0,
        boundaryCondition = "displacement",
        boundaryToggle = 0,
        edgeCellsToggle = 0,
        disorderedInitial= false,
        outputTotal = 100,
        t1Threshold = 0.05,
        peripheralTension = 0.0,
        randomSeed = 0,
        nRows = 9,
        energyModel = "ventilation",
        resistanceModel = "Linear",
        R_in= spzeros(2),
        A_in= spzeros(2),
        B_in= spzeros(2),
        initialEdgeLength = 9/6,
        unstressedSetup = "initial" # "one" or "initial" - whether to set all edge unstressed lengths to the same value (one) or to the initial edge lengths (initial)
    )

    # Calculate derived parameters
    tMax = realTimetMax / viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval = tMax / (outputTotal-1)     # Time interval for storing system data (non dimensionalised)
    λ = -2.0 * L₀ * γ
    nonDimCycleTime = realCycleTime / viscousTimeScale # Non dimensionalised cell cycle time

    # Set random seed value and allocate random number generator
    # Random seed set from current unix time, 
    # unless non zero value of randomSeed is passed, in which case random seed is passed value of randomSeed
    seed = (randomSeed == 0 ? floor(Int64, datetime2unix(now())) : randomSeed)
    Random.seed!(seed)
    rng = MersenneTwister(seed)

    # Initialise system matrices from function or file
    if initialSystem in ["one"]
        A,B,R = initialSmallConfig("one")
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "new"
        isodd(nRows) && (nRows>1)  ? nothing : throw("nRows must be an odd number greater than 1.")
        A, B, R = initialSystemLayout(nRows=nRows,generationMode=generationMode,boundaryToggle=boundaryToggle, initialEdgeLength=initialEdgeLength, edgeCellsToggle=edgeCellsToggle, disorderedInitial = disorderedInitial)
        println("INITIAL R:")
        println("x range = ", extrema(r[1] for r in R))
        println("y range = ", extrema(r[2] for r in R))

        rawLengths = norm.(A * R)

        println("INITIAL EDGE LENGTHS:")
        println(extrema(rawLengths))

    elseif initialSystem == "argument"
        R = R_in
        A = A_in
        B = B_in
    else
        # Import system matrices from final state of previous run
        importedData = load("$initialSystem"; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
        @unpack A,B = importedData["matrices"]
        R = importedData["R"]
    end

    nCells = size(B, 1)
    nEdges = size(A, 1)
    nVerts = size(A, 2)
    

    # Fill preallocated matrices into struct for convenience
    matrices = MatricesContainer(
        A                 = A,
        B                 = B,
        Aᵀ                = spzeros(Int64, nVerts, nEdges),
        Ā                 = spzeros(Int64, nEdges, nVerts),
        Āᵀ                = spzeros(Int64, nVerts, nEdges),
        Bᵀ                = spzeros(Int64, nEdges, nCells),
        B̄                 = spzeros(Int64, nCells, nEdges),
        B̄ᵀ                = spzeros(Int64, nEdges, nCells),
        C                 = spzeros(Int64, nCells, nVerts),
        areaJacobian      = spzeros(SVector{2,Float64}, nCells, nVerts),
        cellEdgeCount     = zeros(Int64, nCells),
        cellVertexOrders  = fill(CircularVector(Int64[]), nCells),
        cellEdgeOrders    = fill(CircularVector(Int64[]), nCells),
        boundaryVertices  = zeros(Int64, nVerts),
        boundaryEdges     = zeros(Int64, nEdges),
        cellPositions     = fill(SVector{2,Float64}(zeros(2)), nCells),
        cellPerimeters    = zeros(nCells),
        cellOrientedAreas = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        cellAreas         = zeros(nCells),
        cellA₀s           = fill(A₀, nCells),
        cellαᵢs           = ones(nCells),
        cellL₀s           = fill(L₀, nCells),
        cellTensions      = zeros(nCells),
        cellPressures     = zeros(nCells),
        μ                 = ones(nCells),
        Γ                 = γ.*ones(nCells), 
        avgEdgeCellNormals = fill(SVector{2,Float64}(zeros(2)), nVerts),
        edgeCellNormals     = spzeros(SVector{2,Float64}, nCells, nEdges),
        edgeLengths       = zeros(nEdges),
        edgeTangents      = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgeTensions      = zeros(nEdges),
        edgeMidpoints     = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgel₀s          = fill(l₀, nEdges),
        edgeMidpointLinks = spzeros(SVector{2,Float64}, nCells, nVerts),
        timeSinceT1       = zeros(nEdges),
        normEdges         = fill(SVector{2,Float64}(zeros(2)), nEdges),
        vertexAreas       = ones(nVerts),
        F                 = spzeros(SVector{2,Float64}, nVerts, nCells),
        FEdges            = spzeros(SVector{2,Float64}, nVerts),
        externalF         = fill(SVector{2,Float64}(zeros(2)), nVerts),
        totalF            = fill(SVector{2,Float64}(zeros(2)), nVerts),
        ϵ                 = SMatrix{2, 2, Float64}([
                                0.0 1.0
                                -1.0 0.0
                            ]),
        cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        dR⁰ = fill(SVector{2,Float64}(zeros(2)), nVerts),
        R⁺ = fill(SVector{2,Float64}(zeros(2)), nVerts), 
        dR⁺ = fill(SVector{2,Float64}(zeros(2)), nVerts), 
        R¹ = fill(SVector{2,Float64}(zeros(2)), nVerts),
        ∂𝒜∂r = fill(SVector{2,Float64}(zeros(2)),nVerts),
        T = fill(SVector{2,Float64}(zeros(2)),nEdges),
        laplacian = spzeros(Float64,nCells,nCells)
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(
        initialSystem     = initialSystem,
        nCells            = nCells,
        nEdges            = nEdges,
        nVerts            = nVerts,
        γ                 = γ,
        λ                 = λ,
        κ                = κ,
        Λ                = Λ,
        L₀                = L₀,
        A₀                = A₀,
        l₀                = l₀,
        Pᵢ                = Pᵢ,
        Pₘ               = Pₘ,
        P₀                = P₀,
        Amp                 = Amp,
        ipModel           = ipModel,
        ω = ω,
        boundaryCondition   = boundaryCondition,
        boundaryToggle    = boundaryToggle,
        edgeCellsToggle  = edgeCellsToggle,
        outputTotal       = outputTotal,
        outputInterval    = outputInterval,
        viscousTimeScale  = viscousTimeScale,
        realTimetMax      = realTimetMax,
        tMax              = tMax,
        realCycleTime     = realCycleTime,
        nCycles           = nCycles,
        nRows            = nRows,
        nonDimCycleTime   = nonDimCycleTime,
        t1Threshold       = t1Threshold,
        peripheralTension = peripheralTension,
        seed              = seed,
        rng               = rng,
        distLogNormal     = LogNormal(0.0, 0.2),
        energyModel       = energyModel,
        resistanceModel   = resistanceModel,
        folderName        = "",
        currentTime        = 0.0,
    )

    # Initial evaluation of matrices based on system topology

    topologyChange!(matrices, params)
    spatialData!(R, params, matrices)

    if resistanceModel == "Linear"
        matrices.cellαᵢs = ones(size(matrices.B, 1))
    elseif resistanceModel == "One Blocked"
        matrices.cellαᵢs = ones(size(matrices.B, 1))
        i=1
        # choose a 'double landlocked' cell to block, i.e. a cell which has no boundary edges and is not adjacent to any cell with boundary edges. For each edge on cell i, check if boundaryEdges is 1 or 0, if 1, try next cell, if zero check that other cell adjacent to that edge also has no boundary edges, if so, block that cell, if not try next cell. If no double landlocked cells exist, throw error.
        while true
            # cellEdgeOrders[i] gives edges of cell i in a CircularArray, so for each edge, find the other cell adjacent to that edge and check if that cell has boundaryEdges[cellEdgeOrders] one or zero. If all edges have other cells with boundaryEdges zero, then cell i is double landlocked and can be blocked.
            if iszero(matrices.boundaryEdges[matrices.cellEdgeOrders[i]]) 
                break
            else
                i += 1
                if i > size(matrices.B, 1)
                    throw("No double landlocked cells found to block.")
                end
            end
        end

        println("Blocking cell ", i)
        matrices.cellαᵢs[i] = 1000000000
        println("Cell αᵢs = ", matrices.cellαᵢs)
    elseif resistanceModel == "Lobes"
        # Currently splits the monolayer into four quadrants with different airway resistances
        # TODO Want to have more random lobes and to identify edges between different lobes
        matrices.cellαᵢs = ones(size(matrices.B, 1))
        for i in 1:size(matrices.B, 1)
            if matrices.cellPositions[i][2] > 0.0 && matrices.cellPositions[i][1] > 0.0
                matrices.cellαᵢs[i] = 1.5
            elseif matrices.cellPositions[i][2] > 0.0 && matrices.cellPositions[i][1] <= 0.0
                matrices.cellαᵢs[i] = 1
            elseif matrices.cellPositions[i][2] <= 0.0 && matrices.cellPositions[i][1] > 0.0
                matrices.cellαᵢs[i] = 0.5
            elseif matrices.cellPositions[i][2] <= 0.0 && matrices.cellPositions[i][1] <= 0.0
                matrices.cellαᵢs[i] = 100
            end
        end

        # Identify edges that are between cells with different αᵢs and store them in a vector
        edgeLobeLinks = []
        for j in 1:size(matrices.B, 2)
            cells = findall(!iszero, matrices.B[:,j])
            if length(cells) == 2
                if matrices.cellαᵢs[cells[1]] != matrices.cellαᵢs[cells[2]]
                    push!(edgeLobeLinks, j)
                end
            end
        end

        # Stiffen edges that are lobe walls
        matrices.edgel₀s[edgeLobeLinks] .= matrices.edgeLengths[edgeLobeLinks] * 4

    else
        throw("Invalid resistance model specified. Must be either 'Linear' or 'One Blocked'.")
    end
    println("Cell αᵢs = ", matrices.cellαᵢs)

    if unstressedSetup == "initial"
        matrices.edgel₀s .= matrices.edgeLengths/2
    end

    # Convert vector of SVectors to flat vector of Float64
    u0 = Float64[]
    for r in R
        append!(u0, r)
    end

    return u0, params, matrices

end

export initialise

end
