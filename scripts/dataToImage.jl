using DrWatson
using FromFile
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf
using Statistics

@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("PlotSetup.jl"))" using PlotSetup

# Select the folder name where the data is stored, how many frames there are and desired features
folderName = "sims\\Full-2D-Experiment\\ventilation_rational\\26-10-21-11-13-29_Pᵢ=0.0_l₀=1.0_nCells=61_realTimetMax=62.8"
outputTotal = 100
plotCells = 1
scatterEdges = 1
scatterVertices = 1
scatterCells = 1
plotForces = 1
plotEdgeMidpointLinks = 1

fig, ax, mov = plotSetup()

xlimits = Float64[0, 0]
ylimits = Float64[0, 0]
pressureLimits = Float64[0, 0]
strainLimits = Float64[0, 2]
for dataFile in 1:outputTotal
    if isfile(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", dataFile)).jld2")) == false
        continue
    end
    data = jldopen(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", dataFile)).jld2"), "r") do file
        return read(file, "matrices"), read(file, "params"), read(file, "R")
    end
    matrices, params, R = data
    if pressureLimits[1] > extrema(matrices.cellPressures)[1]
        pressureLimits[1] = extrema(matrices.cellPressures)[1]
    end
    if pressureLimits[2] < extrema(matrices.cellPressures)[2]
        pressureLimits[2] = extrema(matrices.cellPressures)[2]
    end
    if strainLimits[2] < extrema(matrices.edgeLengths ./ matrices.edgel₀s)[2]
        strainLimits[2] = extrema(matrices.edgeLengths ./ matrices.edgel₀s)[2]
    end
    if ylimits[1] > extrema(r[2] for r in R)[1]
        extreme = extrema(r[2] for r in R)
        ylimits[1] = extreme[1]
    end
    if ylimits[2] < extrema(r[2] for r in R)[2]
        ylimits[2] = extrema(r[2] for r in R)[2]
    end
    if xlimits[1] > extrema(r[1] for r in R)[1]
        xlimits[1] = extrema(r[1] for r in R)[1]
    end
    if xlimits[2] < extrema(r[1] for r in R)[2]
        xlimits[2] = extrema(r[1] for r in R)[2]
    end
end

for dataFile in 1:outputTotal
    if isfile(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", dataFile)).jld2")) == false
        println("File $(dataFile) does not exist, skipping...")
        continue
    end
    data = jldopen(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", dataFile)).jld2"), "r") do file
        return read(file, "matrices"), read(file, "params"), read(file, "R"), read(file, "t")
    end
    matrices, params, R, t = data
    # Render visualisation of system and add frame to movie
    visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks,strainLimits,pressureLimits, limits = (xlimits[1], xlimits[2], ylimits[1], ylimits[2]))
    save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", dataFile)).png"), fig)
    
end