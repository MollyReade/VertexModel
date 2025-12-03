using VertexModel
using JLD2
using UnPack
using FromFile
using Statistics
using Printf
using CairoMakie


meanCellAreas = Float64[]
meanCellPressures = Float64[]

folderName = "sims\\molly-main\\ventilation_rational_cycle\\25-12-03-12-41-47_Pᵢ=0.6_l₀=0.15_nCells=37_realTimetMax=56.4"

for i in 1:500
    lastfile = "data\\"* folderName * "\\frameData\\systemData$(@sprintf("%03d", i)).jld2"
    @load lastfile matrices params
    @unpack cellAreas, cellPressures = matrices
    push!(meanCellAreas, mean(cellAreas))
    push!(meanCellPressures, mean(cellPressures))
end

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Cell Pressure Pᵢ",
    ylabel = "Mean Cell Area Aᵢ",
    title = "Mean Cell Area vs Cell Pressure"
)

lines!(ax, meanCellPressures, meanCellAreas; color = :red)
fig

EVmeanCellAreas = Float64[]
NomeanCellAreas = Float64[]
VmeanCellAreas = Float64[]
time = range(1,104.40,500)

folderEVDiss = "sims\\molly-main\\ventilation_rational_cycle\\25-12-03-11-27-40_Pᵢ=0.6_l₀=0.15_nCells=37_realTimetMax=106.0"
folderNoDiss = "sims\\molly-main\\ventilation_rational_cycle\\25-12-03-11-05-46_Pᵢ=0.6_l₀=0.15_nCells=37_realTimetMax=106.0"
folverVDiss = "sims\\molly-main\\ventilation_rational_cycle\\25-12-03-13-16-18_Pᵢ=0.6_l₀=0.15_nCells=37_realTimetMax=106.0"

for i in 1:500
    lastfileEVDiss = "data\\"* folderEVDiss * "\\frameData\\systemData$(@sprintf("%03d", i)).jld2"
    @load lastfileEVDiss matrices params
    @unpack cellAreas = matrices
    push!(EVmeanCellAreas, mean(cellAreas))

    lastfileNoDiss = "data\\"* folderNoDiss * "\\frameData\\systemData$(@sprintf("%03d", i)).jld2"
    @load lastfileNoDiss matrices params
    @unpack cellAreas = matrices
    push!(NomeanCellAreas, mean(cellAreas))

    lastfileVDiss = "data\\"* folverVDiss * "\\frameData\\systemData$(@sprintf("%03d", i)).jld2"
    @load lastfileVDiss matrices params
    @unpack cellAreas = matrices
    push!(VmeanCellAreas, mean(cellAreas))
end

 

fig2 = Figure()
ax = Axis(fig2[1, 1];
    xlabel = "Time",
    ylabel = "Mean Cell Area Aᵢ",
    title = "Mean Cell Area vs Time"
)

lines!(ax,time, EVmeanCellAreas; color = :blue, label = "With Energy Dissipation")
lines!(ax,time, NomeanCellAreas; color = :red, label = "Without Energy Dissipation")
lines!(ax,time, VmeanCellAreas; color = :green, label = "With Vertex Dissipation")
axislegend(ax; position = :rt)
display(fig2)