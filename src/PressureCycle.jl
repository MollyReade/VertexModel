#
# PressureCycle.jl
# VertexModel
#
# Function to run a sinusoidal pressure cycle simulation, with discrete steps saving system state at each pressure value
#

module PressureCycle

using VertexModel
using UnPack
using Statistics
using FromFile
using DrWatson
using CairoMakie
using Printf

@from "PlotSetup.jl" using PlotSetup
@from "VertexModelContainers.jl" using VertexModelContainers



function pressureCycle(; nSamples=9, nR=5, initialSystem = "new", boundaryToggle = 0,
        edgeCellsToggle = 1,)
    #Set up plot for movie of all runs
    fig, ax, mov = PlotSetup.plotSetup()

    # Create sinusoidal pressure cycle
    points = LinRange(0,12*π,6*nSamples)
    pressures = 0.5*sin.(points) .+ 0.8

    edgeLengths_all = Float64[]
    cellAreas_all = Float64[]

    # Run first simulation from new initial conditions
    integ, _,_,_ = vertexModel(initialSystem = initialSystem,
        nRows = nR,
        realCycleTime = 86400.0,
        l₀ = 0.15,
        A₀ = 1.0,
        Pᵢ = pressures[1],
        boundaryToggle = boundaryToggle,
        edgeCellsToggle = edgeCellsToggle,
        outputToggle = 1,
        outputTotal = 20,
        frameDataToggle = 1,
        frameImageToggle = 1,
        printToggle = 1,
        videoToggle = 1,
        energyModel = "ventilation_rational",
        existingMov = mov,
        existingFig = fig,
        existingAx = ax,
        returnPlots = true)

    (params, matrices) = integ.p
    @unpack edgeLengths, cellAreas = matrices
    @unpack folderName, outputTotal = params

    push!(edgeLengths_all, mean(edgeLengths))
    push!(cellAreas_all, mean(cellAreas))

    # Set lastfile to load for next simulation
    lastfile = "data\\"* folderName * "\\frameData\\systemData$(@sprintf("%03d", outputTotal)).jld2"


    for i in 2:6*nSamples
        # Run simulation from last saved file
        integ, _,_,_ = vertexModel(initialSystem = lastfile,
            nRows = nR,
            realCycleTime = 26400.0,
            l₀ = 0.15,
            A₀ = 1.0,
            Pᵢ = pressures[i],
            boundaryToggle = boundaryToggle,
            edgeCellsToggle = edgeCellsToggle,
            outputToggle = 1,
            outputTotal = 5,
            frameDataToggle = 1,
            frameImageToggle = 1,
            printToggle = 1,
            videoToggle = 1,
            energyModel = "ventilation_rational",
            existingMov = mov,
            existingFig = fig,
            existingAx = ax,
            returnPlots = true)

        (params, matrices) = integ.p
        @unpack edgeLengths, cellAreas = matrices
        @unpack folderName, outputTotal = params

        push!(edgeLengths_all, mean(edgeLengths))
        push!(cellAreas_all, mean(cellAreas)) 

        lastfile = "data\\"* folderName * "\\frameData\\systemData$(@sprintf("%03d", outputTotal)).jld2"
    end
    # Save combined movie
    CairoMakie.save("data\\" * folderName * "\\" * "$(splitpath(folderName)[end]).mp4", mov)

    println("Saved movie to: data\\" * folderName * "\\" * "$(splitpath(folderName)[end]).mp4")

    return pressures, edgeLengths_all, cellAreas_all
end

export pressureCycle

end 
