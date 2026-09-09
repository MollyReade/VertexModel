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
    nRows = 5,
    nCycles = 1,
    realCycleTime = 8.640,
    realTimetMax = nCycles*realCycleTime,
    Λ = 1,
    κ = 0.05,
    γ = 0.2,
    L₀ = 1,
    l₀ = 1,
    A₀ = 1.0,
    Pᵢ = 0.0,
    Pₘ = 0.0,
    P₀ = 0.0,
    ω = π/(2*realCycleTime),
    Amp = 0.0,
    ipModel = "sinusoidal",
    viscousTimeScale = 1.0,
    peripheralTension = 0.0,
    t1Threshold = 0.00,
    boundaryCondition = "displacement", # "displacement" or "force"
    boundaryToggle = 0,
    edgeCellsToggle = 0,
    solveMethod = "Manual",  #Manual or ODEProblem
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
        ω = ω,
        Amp = Amp,
        ipModel = ipModel,
        viscousTimeScale = viscousTimeScale,
        boundaryToggle = boundaryToggle,
        boundaryCondition = boundaryCondition,
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

    if solveMethod == "Manual"
        t = 0.0
        alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
        outputCounter = [1]
        u = u0
        while t < params.tMax
            dt = params.outputInterval/10000.0
            R = reinterpret(SVector{2,Float64}, u)
            if abs((t - alltStops[outputCounter[1]])) < 1e-8
                # Update progress on command line 
                printToggle == 1 ? println("$(@sprintf("%.2f", t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
                if frameDataToggle == 1
                    # Save system data to file 
                    jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
                end
                if frameImageToggle == 1 || videoToggle == 1
                    # Render visualisation of system and add frame to movie
                    visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
                end
                # Save still image of this time step 
                frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
                outputCounter[1] += 1
            end
            u = splitStep(u, u0,(params,matrices), t, dt/2)
            t += dt

            params.currentTime = t
            spatialData!(R, params, matrices)
            # println("\n========== DEBUG ==========")

            # println("edge length       : ", extrema(matrices.edgeLengths))
            # println("edge length/l₀    : ", extrema(matrices.edgeLengths ./ l₀))
            # println("edge tension      : ", extrema(matrices.edgeTensions))
            # println("cell area         : ", extrema(matrices.cellAreas))
            # println("vertex area D     : ", extrema(matrices.vertexAreas))
            # println("cell alpha        : ", extrema(matrices.cellαᵢs))

            # println("============================\n")
        end
        (outputToggle == 1 && videoToggle == 1 && isnothing(existingMov)) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing
        if returnPlots
            return u, fig, ax, mov
        else
            return u
        end

    elseif solveMethod == "ODEProblem"
        # Set up ODE integrator 
        #prob = ODEProblem(model!, u0, (0.0, Inf), (params, matrices))
        probT = ODEProblem(modelT!, u0, (0.0, Inf), (params, matrices))
        probP = ODEProblem(modelP!, u0, (0.0, Inf), (params, matrices))
        alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
        dt = params.outputInterval
        #if energyModel == "log"
        #    integrator = init(prob, Tsit5(), tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
        #else
        #    integrator = init(prob, solver, tstops=alltStops, abstol=abstol, reltol=reltol)
        #end
        #integrator = LieIntegrator(probT, probP, copy(u0),0.0,dt,nothing, solver)
        integratorT = init(probT, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
        integratorP = init(probP, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
        outputCounter = [1]

        u = copy(u0)
        t=0.0

        while t <= params.tMax
            dtStep = min(params.outputInterval, params.tMax - t)

            probT = remake(probT; u0=u, tspan=(t, t+dtStep))
            solT = solve(probT, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
            u = solT.u[end]

            R = reinterpret(SVector{2,Float64}, u)
            params.currentTime = t+dtStep
            spatialData!(R, params, matrices)

            probP = remake(probP; u0=u, tspan=(t, t+dtStep))
            solP = solve(probP, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
            u = solP.u[end]
            
        
            R = reinterpret(SVector{2,Float64}, u)
            spatialData!(R, params, matrices)
            println("t = ", t, ", alltStops[outputCounter[1]] = ", alltStops[outputCounter[1]], ", outputCounter[1] = ", outputCounter[1])
            if (t - alltStops[outputCounter[1]]) < 1e-8
                # Update progress on command line 
                printToggle == 1 ? println("$(@sprintf("%.2f", t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
                if frameDataToggle == 1
                    # Save system data to file 
                    
                    jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
                end
                if frameImageToggle == 1 || videoToggle == 1
                    # Render visualisation of system and add frame to movie
                    visualise(R, integratorP.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
                end
                # Save still image of this time step 
                frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
                outputCounter[1] += 1
            end
            t += dtStep
            params.currentTime = t
        end

        #Iterate until integrator time reaches max system time 
        # while integratorP.t <= params.tMax 
            
        #     # Reinterpret state vector as a vector of SVectors 
        #     R = reinterpret(SVector{2,Float64}, integratorP.u)
        #     # Note that reinterpreting accesses the same underlying data, so changes to R will update integrator.u and vice versa 

        #     # Output data to file 
        #     if integratorP.t == alltStops[outputCounter[1]] 
        #         # Update progress on command line 
        #         printToggle == 1 ? println("$(@sprintf("%.2f", integratorP.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
        #         if frameDataToggle == 1
        #             # Save system data to file 
                    
        #             jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
        #         end
        #         if frameImageToggle == 1 || videoToggle == 1
        #             # Render visualisation of system and add frame to movie
        #             visualise(R, integratorP.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
        #         end
        #         # Save still image of this time step 
        #         frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
        #         outputCounter[1] += 1
        #     end

        #     # Step integrator forwards in time to update vertex positions 
        #     step!(integratorP)
        #     spatialData!(R, params, matrices)
        #     #step!(integratorP)

        #     # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        #     params.currentTime = integratorP.t
        #     #spatialData!(R, params, matrices)

        #     # Check system for T1 transitions 
        #     #if t1Transitions!(integrator, params, matrices) > 0
        #     #    u_modified!(integrator, true)
        #     #    # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
        #     #    topologyChange!(matrices, params) # Update system matrices after T1 transition
        #     #    spatialData!(R, params, matrices) # Update spatial data after T1 transition  
        #     #end
        #     # Update cell ages with (variable) timestep used in integration step
        #     matrices.timeSinceT1 .+= integratorP.dt
        # end

        # If outputToggle==1, save animation object and save final system matrices
        (outputToggle == 1 && videoToggle == 1 && isnothing(existingMov)) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

        if returnPlots
            return integratorP, fig, ax, mov
        else
            return integratorP
        end
    end
end

# Function to load previously saved simulation data 
function loadData(relativePath; outputNumber=100)
    data = load(projectdir(relativePath, "frameData", "systemData$(@sprintf("%03d", outputNumber)).jld2"))
    return data["R"], data["matrices"], data["params"]
end

function lieStep(u,t,dt,params,matrices,solver,abstol,reltol)
    # Tension step

    probT = ODEProblem(modelT!, u, (t, t+dt), (params, matrices))

    solT = solve(probT, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    u = solT.u[end]

    R  = reinterpret(SVector{2,Float64}, u)
    params.currentTime = t+dt
    spatialData!(R, params, matrices)

    # Pressure step
    probP = ODEProblem(modelP!, u, (t, t+dt), (params, matrices))
    solP = solve(probP, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    u = solP.u[end]
    return u
end

mutable struct LieIntegrator{IT, IP, S}
    probT::IT
    probP::IP
    u
    t
    dt
    sol
    solver::S
end

function initLieIntegrator(probT,probP,solver;dt,abstol=1e-7,reltol=1e-4)
    intT = init(probT, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    intP = init(probP, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)

    return LieIntegrator(intT, intP,copy(probT.u0),first(probT.tspan),dt,nothing,solver)
end
function lieStep!(u,t,dt,params,matrices,solver,abstol,reltol)
    probT = ODEProblem(modelT!, u, (t, t+dt), (params, matrices))
    solT = solve(probT, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    u = solT.u[end]

    R = reinterpret(SVector{2,Float64}, u)
    params.currentTime = t+dt
    spatialData!(R, params, matrices)

    probP = ODEProblem(modelP!, u, (t, t+dt), (params, matrices))
    solP = solve(probP, solver, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    u = solP.u[end]
    return u
end
# function step!(I::LieIntegrator)

#     # T step
#     probT = remake(I.probT, u=I.u, tspan=(I.t,I.t+I.dt))
#     solT = solve(probT, I.solver; save_on=false, save_start=false, save_end=true)
#     I.u = solT.u[end]

#     # P step
#     probP = remake(I.probP, u=I.u, tspan=(I.t,I.t+I.dt))
#     solP = solve(probP, I.solver; save_on=false, save_start=false, save_end=true)
#     I.u = solP.u[end]

#     # Advance time
#     I.t += I.dt
#     I.sol = solP
#     return nothing
# end
# Ensure code is precompiled
# @compile_workload begin
#     vertexModel(nCycles=0.01, outputToggle=0, frameDataToggle=0, frameImageToggle=0, printToggle=0, videoToggle=0)
# end


export vertexModel
export loadData 

end
