#
#  VertexModel.jl
#  VertexModel
#
#

module VertexModel

# Julia packages
using PrecompileTools
using DrWatson
using FromFile
using OrdinaryDiffEq
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf
using OrdinaryDiffEqSDIRK
using DifferentialEquations

# Local modules
@from "CreateRunDirectory.jl" using CreateRunDirectory
@from "Visualise.jl" using Visualise
@from "Initialise.jl" using Initialise
@from "SpatialData.jl" using SpatialData
@from "PlotSetup.jl" using PlotSetup
@from "Model.jl" using Model
@from "T1Transitions.jl" using T1Transitions
@from "TopologyChange.jl" using TopologyChange
@from "SenseCheck.jl" using SenseCheck

function vertexModel(;
    initialSystem = "new",
    nRows = 7,
    nCycles = 1,
    realCycleTime = 16.640,
    realTimetMax = nCycles*realCycleTime,
    Λ = 0.1,
    κ = 0.05,
    γ = 0.2,
    L₀ = 0.75,
    l₀ = 0.45,
    A₀ = 1.0,
    Pᵢ = 0.8,
    Pₘ = 0.0,
    P₀ = -0.0,
    Amp = 0.4,
    ipModel = "sinusoidal",
    viscousTimeScale = 1.0,
    peripheralTension = 1.0,
    t1Threshold = 0.00,
    boundaryToggle = 1,
    edgeCellsToggle = 1,
    solver = Tsit5(),
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 100,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 1,
    printToggle = 1,
    videoToggle = 1,
    plotCells = 1,
    scatterEdges = 1,
    scatterVertices = 1,
    scatterCells = 1,
    plotForces = 1,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
    energyModel = "ventilation_rational",
    resistanceModel = "Linear",
    dissipationToggle = 0,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2),
    existingMov = nothing,
    existingFig = nothing,
    existingAx = nothing,
    returnPlots = false
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem = initialSystem,
        realTimetMax = realTimetMax,
        γ = γ,
        Λ = Λ,
        κ = κ,
        L₀ = L₀,
        A₀ = A₀,
        l₀ = l₀,
        Pᵢ = Pᵢ,
        Pₘ = Pₘ,
        P₀ = P₀,
        Amp = Amp,
        ipModel = ipModel,
        viscousTimeScale = viscousTimeScale,
        boundaryToggle = boundaryToggle,
        edgeCellsToggle = edgeCellsToggle,
        outputTotal = outputTotal,
        t1Threshold = t1Threshold,
        realCycleTime = realCycleTime,
        peripheralTension = peripheralTension,
        randomSeed = randomSeed,
        nRows = nRows,
        energyModel = energyModel,
        resistanceModel = resistanceModel,
        R_in = R_in,
        A_in = A_in,
        B_in = B_in,
    )
    fig = existingFig !== nothing ? existingFig : nothing
    ax  = existingAx  !== nothing ? existingAx  : nothing
    mov = existingMov !== nothing ? existingMov : nothing
    # Create directory in which to store date. Save parameters and store directory name for later use.
    if outputToggle == 1
        subFolder=energyModel
        folderName = createRunDirectory(params,subFolder)
        params.folderName = folderName
        
        # Create plot object for later use 
        if (frameImageToggle==1 || videoToggle==1) && isnothing(mov)
            fig, ax, mov = plotSetup()
        elseif (frameImageToggle==1 || videoToggle==1)
            fig= existingFig
            ax= existingAx
            mov= existingMov
        end
    end

    # Set up ODE integrator 
    prob = ODEProblem(model!, u0, (0.0, Inf), (params, matrices))
    alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
    if energyModel == "log"
        integrator = init(prob, Tsit5(), tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    else
        integrator = init(prob, solver, tstops=alltStops, abstol=abstol, reltol=reltol)
    end
    outputCounter = [1]

    # Iterate until integrator time reaches max system time 
    while integrator.t <= params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
        
        # Reinterpret state vector as a vector of SVectors 
        R = reinterpret(SVector{2,Float64}, integrator.u)
        # Note that reinterpreting accesses the same underlying data, so changes to R will update integrator.u and vice versa 

        # Output data to file 
        if integrator.t == alltStops[outputCounter[1]]
            # Update progress on command line 
            printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
            if frameDataToggle == 1
                # Save system data to file 
                
                jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
            end
            if frameImageToggle == 1 || videoToggle == 1
                # Render visualisation of system and add frame to movie
                visualise(R, integrator.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
            end
            # Save still image of this time step 
            frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
            outputCounter[1] += 1
        end

        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        params.currentTime = integrator.t
        spatialData!(R, params, matrices)

        # Check system for T1 transitions 
        if t1Transitions!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices, params) # Update system matrices after T1 transition
            spatialData!(R, params, matrices) # Update spatial data after T1 transition  
        end
        # Update cell ages with (variable) timestep used in integration step
        matrices.timeSinceT1 .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    (outputToggle == 1 && videoToggle == 1 && isnothing(existingMov)) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

    if returnPlots
        return integrator, fig, ax, mov
    else
        return integrator
    end
end

# Function to load previously saved simulation data 
function loadData(relativePath; outputNumber=100)
    data = load(projectdir(relativePath, "frameData", "systemData$(@sprintf("%03d", outputNumber)).jld2"))
    return data["R"], data["matrices"], data["params"]
end

# Ensure code is precompiled
@compile_workload begin
    vertexModel(nCycles=0.01, outputToggle=0, frameDataToggle=0, frameImageToggle=0, printToggle=0, videoToggle=0)
end

export vertexModel
export loadData 

end
