# Sweep values of omega (period) and plot different graphs for each value. Want ideally one well sampled oscillation cycle for each value, so need to adjust run time
using VertexModel
using CairoMakie
using JLD2
using DrWatson

P₀Values = [0.2, 0.4, 0.8, 1.6, 3.2]

nRows = 11
κ = 0.2
l₀ = 1
Amp = 0.2
ω = 0.05
realCycleTime = 2π/ω
boundaryModel = "sinusoidal"
relaxFirst = true

#Redo simulations
# for P₀ in P₀Values
#     nCycles = 3
#     vertexModel(nRows=nRows, κ=κ, l₀=l₀, P₀=P₀, Amp=Amp, ω=ω, boundaryModel=boundaryModel, realCycleTime=realCycleTime, nCycles = nCycles, relaxFirst=relaxFirst, frameImageToggle=false, videoToggle=false, subFolder="P₀Sweep",outputTotal=150)
# end 

P₀DataPaths = ["data\\sims\\Full-2D-Experiment\\P₀Sweep\\26-10-15-11-22-51_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=377.0\\frameData", "data\\sims\\Full-2D-Experiment\\P₀Sweep\\26-10-15-11-36-19_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=377.0\\frameData","data\\sims\\Full-2D-Experiment\\P₀Sweep\\26-10-15-11-54-51_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=377.0\\frameData", "data\\sims\\Full-2D-Experiment\\P₀Sweep\\26-10-15-12-18-31_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=377.0\\frameData", "data\\sims\\Full-2D-Experiment\\P₀Sweep\\26-10-15-12-51-29_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=377.0\\frameData"]


P₀SweepFig = Figure(resolution = (2000, 2000)) # Will have 3 rows and 5 columns, one for each value of ω, and each row will be a different plot type (applied pressure vs inner/outer areas, alveolar pressures over time, inner vs outer alveolar pressures)

for j in 1:length(P₀Values)
    folder = P₀DataPaths[j]
    P₀ = P₀Values[j]
    innerPressures = zeros(149)
    outerPressures = zeros(149)
    innerAreas = zeros(149)
    outerAreas = zeros(149)
    appliedPressures = zeros(149)
    totalAreas = zeros(149)
    alltStops = collect(0.0:0.1:9.8)

    # Load axis for applied pressure vs inner/outer areas
    ax1 = Axis(P₀SweepFig[1, j];
        xlabel = "Applied Pressure Pₚ",
        ylabel = "Alveolar Areas",
        title = "Applied Pressure vs Alveolar Areas for P₀ = $(P₀)",
        xgridvisible = false,
        ygridvisible = false
    )

    # Load axis for alveolar pressures over time
    ax2 = Axis(P₀SweepFig[2, j];
        xlabel = "Time t",
        ylabel = "Alveolar Pressures",
        title = "Alveolar Pressures vs Time for P₀ = $(P₀)",
        xgridvisible = false,
        ygridvisible = false
    )

    # Load axis for inner vs outer alveolar pressures   
    ax3 = Axis(P₀SweepFig[3, j];
        xlabel = "Inner Alveolar Pressure Pᵢ",
        ylabel = "Outer Alveolar Pressure Pₒ",
        title = "Inner vs Outer Alveolar Pressures for P₀ = $(P₀)",
        xgridvisible = false,
        ygridvisible = false
    )

    for i in 1:149
        matrices, params = load("$(folder)\\systemData$(lpad(string(i),3,'0')).jld2", "matrices", "params")
        if params.ipModel == "sinusoidal" 
            appliedPressures[i] = params.P₀ + params.Amp * sin(params.ω * (i-1) * params.outputInterval)
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
        alltStops = collect(0.0:params.outputInterval:params.tMax)
    end

    lines!(ax1, appliedPressures, innerAreas; color = :red, label = "Inner Alveoli Area")
    lines!(ax1, appliedPressures, outerAreas; color = :blue, label = "Outer Alveoli Area")
    #lines!(ax1, appliedPressures, totalAreas; color = :green, label = "Total Alveolar Area")
    axislegend(ax1; position = :rt)

    lines!(ax2, alltStops[1:149], innerPressures; color = :red, label = "Inner Alveoli Pressure")
    lines!(ax2, alltStops[1:149], outerPressures; color = :blue, label = "Outer Alveoli Pressure")
    axislegend(ax2; position = :rt)

    lines!(ax3, innerPressures, outerPressures; color = :purple, label = "Inner vs Outer Alveoli Pressures")
end

display(P₀SweepFig)


