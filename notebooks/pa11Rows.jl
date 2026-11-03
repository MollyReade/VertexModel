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


# Files for small alpha, small kappa
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-05-11-46-57_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-05-12-11-03_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-05-12-45-05_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-06-07-36-48_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-06-09-15-48_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-06-09-29-20_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData
# C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-06-10-05-52_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData low P_p and Amp

dataFolder = "C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-10-06-09-29-20_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=50.0\\frameData" #"data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-09-28-13-10-59_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=100.0\\frameData"
innerPressures = zeros(99)
innerAreas = zeros(99)
outerPressures = zeros(99)
outerAreas = zeros(99)
appliedPressures = zeros(99)
totalAreas = zeros(99)

matrices, params = load("$(dataFolder)\\systemData$(lpad(string(1),3,'0')).jld2", "matrices", "params")
alltStops = collect(0.0:params.outputInterval:params.tMax)

for i in 1:99
    matrices, params = load("$(dataFolder)\\systemData$(lpad(string(i),3,'0')).jld2", "matrices", "params")
    if params.ipModel == "sinusoidal" 
        appliedPressures[i] = params.P₀ + params.Amp * sin(alltStops[i])
    else
        appliedPressures[i] = params.P₀
    end

    cellPressures = matrices.cellPressures
    cellAreas = matrices.cellAreas
    innerPressures[i] = cellPressures[7] # should pick out the central alveoli
    innerAreas[i] = cellAreas[7]
    outerPressures[i] = cellPressures[66] # an alveoli one layer in, to avoid the non-hexagonal edge cells
    outerAreas[i] = cellAreas[66]
    totalAreas[i] = sum(cellAreas)
end

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Inner Alveoli Pressure Pᵢ",
    ylabel = "Outer Alveoli Pressure Pₒ",
    title = "Inner vs Outer Alveoli",
    xgridvisible = false,
    ygridvisible = false
)

fig2 = Figure()
ax2 = Axis(fig2[1, 1];
    xlabel = "Applied Pleural Pressure Pₚ",
    ylabel = "Cell Area Aᵢ",
    title = "Applied Pleural Pressure vs Inner and Outer Alveoli Area",
    xgridvisible = false,
    ygridvisible = false
)

fig3 = Figure()
ax3 = Axis(fig3[1, 1];
    xlabel = "Time t",
    ylabel = "Cell Pressure Pᵢ",
    title = "Inner and Outer Alveoli Pressure vs Time",
    xgridvisible = false,
    ygridvisible = false
)

lines!(ax3, innerPressures[3:end]; color = :red)
lines!(ax3, outerPressures[3:end]; color = :blue)
#axislegend(ax3, ["Inner Alveoli Pressure" "Outer Alveoli Pressure"]; position = :rt)

lines!(ax2, appliedPressures[3:end], innerAreas[3:end]; color = :red)
lines!(ax2, appliedPressures[3:end], outerAreas[3:end]; color = :blue)
lines!(ax2, appliedPressures[3:end], totalAreas[3:end]; color = :green)
#axislegend(ax2, ["Inner Alveoli Area" "Outer Alveoli Area" "Total Area"]; position = :rt)

lines!(ax, innerPressures[3:end], outerPressures[3:end]; color = :red)


println("Loaded parameters: ", params, "Alpha: ", matrices.cellαᵢs)
fig3