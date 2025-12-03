using JLD2
using CairoMakie

@load "cellAreas_allDrag.jld2" cellAreas_all
@load "edgeLengths_allDrag.jld2" edgeLengths_all
@load "pressures.jld2" pressures

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Cell Pressure Pᵢ",
    ylabel = "Mean Alveolar Area Aᵢ",
    title = "Mean Alveolar Area vs Cell Pressure"
)
lines!(ax, pressures, cellAreas_all; color = :blue)
scatter!(ax, pressures[1:19], cellAreas_all[1:19]; color = :red)
scatter!(ax, pressures[20:38], cellAreas_all[20:38]; color = :purple)
scatter!(ax, pressures[39:57], cellAreas_all[39:57]; color = :yellow)

#axislegend(ax; position = :rt)
display(fig)

fig2 = Figure()
ax = Axis(fig2[1, 1];
    xlabel = "Cell Pressure Pᵢ",
    ylabel = "Mean Edge Length lⱼ",
    title = "Mean Edge Length vs Cell Pressure"
)
lines!(ax, pressures, edgeLengths_all; color = :blue)
scatter!(ax, pressures[1:19], edgeLengths_all[1:19]; color = :red)
scatter!(ax, pressures[20:38], edgeLengths_all[20:38]; color = :purple)
scatter!(ax, pressures[39:57], edgeLengths_all[39:57]; color = :yellow)

#axislegend(ax; position = :rt)
display(fig2)


fig3 = Figure()
ax = Axis(fig3[1, 1];
    xlabel = "Time",
    ylabel = "Mean Alveolar Area Aᵢ",
    title = "Mean Alveolar Area vs Time"
)
lines!(ax, cellAreas_all; color = :blue)
scatter!(ax, cellAreas_all[1:19]; color = :red)
scatter!(ax, cellAreas_all[20:38]; color = :purple)
scatter!(ax, cellAreas_all[39:57]; color = :yellow)

#axislegend(ax; position = :rt)
display(fig3)

UVᵢ(Aᵢ,lᵢ,l₀,Pᵢ,Tᵢ) = Pᵢ * Aᵢ + 1/(2) * Tᵢ

tens_sum = B̄*((edgeLengths_all.^2 .- edgeLengths_all.^(-3))./6)
for i = 1:nCells
    energyPerCell[i] = UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPᵢs[i],tens_sum[i])
end