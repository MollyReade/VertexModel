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


# Files
#data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-08-12-05-47_Pᵢ=0.0_l₀=1.0_nCells=37_realTimetMax=50.0\\frameData
#data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-08-11-34-03_Pᵢ=0.0_l₀=1.0_nCells=37_realTimetMax=50.0\\frameData
#data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-08-11-26-17_Pᵢ=0.0_l₀=1.0_nCells=37_realTimetMax=50.0\\frameData

dataFolder = "data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-08-13-21-21_Pᵢ=0.0_l₀=1.0_nCells=37_realTimetMax=8.0\\frameData"
innerPressures = zeros(99)
innerAreas = zeros(99)
outerPressures = zeros(99)
outerAreas = zeros(99)
appliedPressures = zeros(99)

matrices, params = load("$(dataFolder)\\systemData$(lpad(string(1),3,'0')).jld2", "matrices", "params")
alltStops = collect(0.0:params.outputInterval:params.tMax)

for i in 1:99
    matrices, params = load("$(dataFolder)\\systemData$(lpad(string(i),3,'0')).jld2", "matrices", "params")
    if params.ipModel == "sinusoidal" 
        appliedPressures[i] = params.P₀ + params.Amp * sin(alltStops[i]*params.ω)
    else
        appliedPressures[i] = params.P₀
    end

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
    ylabel = "Cell Area Aᵢ",
    title = "Inner vs Outer Alveoli"
)

fig2 = Figure()
ax2 = Axis(fig2[1, 1];
    xlabel = "Applied Pleural Pressure Pₚ",
    ylabel = "Cell Area Aᵢ",
    title = "Applied Pleural Pressure vs Inner and Outer Alveoli Area"
)

fig3 = Figure()
ax3 = Axis(fig3[1, 1];
    xlabel = "Time t",
    ylabel = "Cell Pressure Pᵢ",
    title = "Inner and Outer Alveoli Pressure vs Time"
)

lines!(ax3, innerPressures[3:end]; color = :red)
lines!(ax3, appliedPressures[3:end]; color = :blue)

lines!(ax2, appliedPressures[3:end], innerAreas[3:end]; color = :red)
lines!(ax2, appliedPressures[3:end], outerAreas[3:end]; color = :blue)

lines!(ax, innerPressures[3:end], innerAreas[3:end]; color = :red)
lines!(ax, outerPressures[3:end], outerAreas[3:end]; color = :blue)


println("Loaded parameters: ", params, "Alpha: ", matrices.cellαᵢs)
fig3