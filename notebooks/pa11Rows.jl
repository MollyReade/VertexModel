# Using set up with 11 rows, κ = 0.6, P\_0 = 0.7, sinusoidal pleural pressure with amplitude 0.3 and with edge cells
# Shows that motion propogates from the outer alveoli to the inner alveoli with decreasing amplitude of motion - inside, the attained pressures are of smaller size and there are less extreme changes in area


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
# jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)

dataFolder = "data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-09-28-13-10-59_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=100.0\\frameData"

innerPressures = zeros(99)
innerAreas = zeros(99)
outerPressures = zeros(99)
outerAreas = zeros(99)

for i in 1:99
    matrices, params = load("$(dataFolder)\\systemData$(lpad(string(i),3,'0')).jld2", "matrices", "params")
    cellPressures = matrices.cellPressures
    cellAreas = matrices.cellAreas
    innerPressures[i] = cellPressures[convert(Int,(size(cellPressures,1)-1)/2)] # should pick out the central alveoli
    innerAreas[i] = cellAreas[convert(Int,(size(cellAreas,1)-1)/2)]
    outerPressures[i] = cellPressures[params.nRows-1] # an alveoli one layer in, to avoid the non-hexagonal edge cells
    outerAreas[i] = cellAreas[params.nRows-1]
end

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Cell Pressure Pᵢ",
    ylabel = "Mean Cell Area Aᵢ",
    title = "Inner vs Outer Alveoli"
)

lines!(ax, innerPressures, innerAreas; color = :red)
lines!(ax, outerPressures, outerAreas; color = :blue)

fig