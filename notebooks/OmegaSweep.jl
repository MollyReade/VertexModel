# Sweep values of omega (period) and plot different graphs for each value. Want ideally one well sampled oscillation cycle for each value, so need to adjust run time
using VertexModel
using CairoMakie
using JLD2
using DrWatson

#ωValues = [ 0.001, 0.01, 0.1, 1.0, 10.0]
ωValues = [ 0.01, 0.1, 1.0, 10.0]
# small omega case stops running with no obvious singularities
 
realCycleTimes = [2π/ω for ω in ωValues]

nRows = 11
κ = 0.2
l₀ = 1
P₀ = 0.2
Amp = 0.1
boundaryModel = "sinusoidal"
relaxFirst = true

#Redo simulations
for ω in ωValues
    realCycleTime = 2π/ω
    nCycles = 5
    vertexModel(nRows=nRows, κ=κ, l₀=l₀, P₀=P₀, Amp=Amp, ω=ω, boundaryModel=boundaryModel, realCycleTime=realCycleTime, nCycles = nCycles, relaxFirst=relaxFirst, frameImageToggle=false, videoToggle=false, subFolder="ωSweep",outputTotal=550)
end 
#data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-15-09-08-56_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=126.0\\frameData
#data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-52-04_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=12600.0\\frameData
ωDataPaths = ["data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-15-08-58-12_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=1260.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-15-10-02-26_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=251.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-15-09-20-39_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=12.6\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-15-09-31-59_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=1.26\\frameData"]
# ωDataPaths = [ "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-10-04-17_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=1260.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-10-08-13_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=126.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-10-12-06_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=12.6\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-10-15-47_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=1.26\\frameData"]
# #ωDataPaths = ["data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-25-54_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6910.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-33-18_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=691.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-36-16_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=69.1\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-39-14_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6.91\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-09-42-06_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=0.691\\frameData"]
# #ωDataPaths = ["data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-25-51_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6280.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-32-49_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=628.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-35-15_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=62.8\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-37-49_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6.28\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-40-20_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=0.628\\frameData"]
# #ωDataPaths = ["data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-10-38_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6280.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-12-55_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=628.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-15-44_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=62.8\\frameData","data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-18-37_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6.28\\frameData","data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-13-08-21-07_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=0.628\\frameData"]
# #ωDataPaths = ["data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-12-14-39-23_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6280.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-12-13-34-04_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=628.0\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-12-14-44-00_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=62.8\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-12-13-40-53_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=6.28\\frameData", "data\\sims\\Full-2D-Experiment\\ωSweep\\26-10-12-14-47-33_Pᵢ=0.0_l₀=1.0_nCells=91_realTimetMax=0.628\\frameData"]
ωSweepFig = Figure(resolution = (2000, 2000)) # Will have 3 rows and 5 columns, one for each value of ω, and each row will be a different plot type (applied pressure vs inner/outer areas, alveolar pressures over time, inner vs outer alveolar pressures)

for j in 1:length(ωValues)
    folder = ωDataPaths[j]
    ω = ωValues[j]
    innerPressures = zeros(149)
    outerPressures = zeros(149)
    innerAreas = zeros(149)
    outerAreas = zeros(149)
    appliedPressures = zeros(149)
    totalAreas = zeros(149)
    alltStops = collect(0.0:0.1:9.8)

   

    # Load axis for applied pressure vs inner/outer areas
    ax1 = Axis(ωSweepFig[1, j];
        xlabel = "Applied Pressure Pₚ",
        ylabel = "Alveolar Areas",
        title = "ω = $(ω)",
        xgridvisible = false,
        ygridvisible = false,
        titlesize = 20,
        xlabelsize = 18,
        ylabelsize = 18
    )

    # Load axis for alveolar pressures over time
    ax2 = Axis(ωSweepFig[3, j];
        xlabel = "Time t",
        ylabel = "Alveolar Pressures",
        title = "ω = $(ω)",
        xgridvisible = false,
        ygridvisible = false,
        titlesize = 20,
        xlabelsize = 18,
        ylabelsize = 18
    )

    # Load axis for inner vs outer alveolar pressures   
    ax3 = Axis(ωSweepFig[5, j];
        xlabel = "Inner Alveolar Pressure Pᵢ",
        ylabel = "Outer Alveolar Pressure Pₒ",
        title = "ω = $(ω)",
        xgridvisible = false,
        ygridvisible = false,
        titlesize = 20,
        xlabelsize = 18,
        ylabelsize = 18
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
    axislegend(ax1; position = :rt, fontsize = 17)

    lines!(ax2, alltStops[1:149], innerPressures; color = :red, label = "Inner Alveoli Pressure")
    lines!(ax2, alltStops[1:149], outerPressures; color = :blue, label = "Outer Alveoli Pressure")
    axislegend(ax2; position = :rt, fontsize = 17)

    lines!(ax3, innerPressures, outerPressures; color = :purple, label = "Inner vs Outer Alveoli Pressures")
end
supertitle1 = Label(ωSweepFig[0,:], "Applied Pressure vs Alveolar Areas", fontsize = 28, justification = :center)
supertitle2 = Label(ωSweepFig[2,:], "Alveolar Pressures vs Time", fontsize = 28, justification = :center)
supertitle3 = Label(ωSweepFig[4,:], "Inner vs Outer Alveolar Pressures", fontsize = 28, justification = :center)

display(ωSweepFig)


