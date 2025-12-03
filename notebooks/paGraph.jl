using VertexModel
using CairoMakie
using UnPack
using JLD2
using FromFile
using DrWatson
using Printf

nRows = 7
initialP = 0.5*sin(1.25*0)+0.6

integInitial,_,_,_ = vertexModel(
    initialSystem = "new",
    nRows = nRows,
    nCycles = 1,
    realCycleTime = 864.0,
    γ = 0.2,
    L₀ = 0.6*6,
    l₀ = 0.15,
    A₀ = 1.0,
    Pᵢ = initialP,
    viscousTimeScale = 1.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    t1Threshold = 0.00,
    divisionToggle = 0,
    boundaryToggle = 0,
    edgeCellsToggle = 1,
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 2,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 0,
    printToggle = 1,
    videoToggle = 0,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
    energyModel = "ventilation_rational",
    dissipationToggle = 1,
    edgeDissToggle = 0,
    vertexDissToggle = 1
)

(paramsInitial, matricesInitial) = integInitial.p
@unpack edgeLengths= matricesInitial
@unpack folderName, outputTotal = paramsInitial

lastfile = "data\\"* folderName * "\\frameData\\systemData$(@sprintf("%03d", outputTotal)).jld2"


integ,_,_,_ = vertexModel(
    initialSystem = lastfile,
    nRows = nRows,
    nCycles = 1,
    realCycleTime = 106.40,
    γ = 0.2,
    L₀ = 0.6*6,
    l₀ = 0.15,
    A₀ = 1.0,
    Pᵢ = initialP,
    viscousTimeScale = 1.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    t1Threshold = 0.00,
    divisionToggle = 0,
    boundaryToggle = 0,
    edgeCellsToggle = 1,
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 500,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 0,
    printToggle = 1,
    videoToggle = 1,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
    energyModel = "ventilation_rational_cycle",
    dissipationToggle = 1,
    edgeDissToggle = 0,
    vertexDissToggle = 1
)

(params, matrices) = integ.p

@unpack folderName = params

println("Loading from folder: ", folderName)




#@from "$(srcdir("PressureCycle.jl"))" using PressureCycle



#pressures, edgeLengths_all, cellAreas_all = PressureCycle.pressureCycle(nSamples=19, nR=7)

#@save "cellAreas_allDrag.jld2" cellAreas_all
#@save "edgeLengths_allDrag.jld2" edgeLengths_all
#@save "pressures.jld2" pressures
