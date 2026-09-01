#
#  Energy.jl
#  VertexModel
#
#  Function to calculate system energy

module Energy

# Julia packages
using LinearAlgebra
using UnPack

# Energy per Cowley et al. 2024 Section 2a
𝒰(θ) = θ*(log(θ)-1.0)
Uᵢ(Aᵢ, A₀, Lᵢ, L₀, μᵢ, Γᵢ) = μᵢ*(𝒰(Aᵢ/A₀) + Γᵢ*L₀^2*𝒰(Lᵢ/L₀))

#Energy for ventilation model with pressure, drag and tension

UVᵢ(Aᵢ,lᵢ,l₀,Pᵢ,Tᵢ) = Pᵢ * Aᵢ + 1/(2) * Tᵢ

function energy(params,matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        B̄,
        edgeLengths,
        cellPressures,
        edgel₀s = matrices
    @unpack energyModel,
    l₀ = params
    print(stacktrace())
    if energyModel == "ventilation"
        # Ventilation energy
        tens_sum = B̄*(edgeLengths.-1).^2
        for i = 1:nCells
            
            energyTotal += UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPressures[i],tens_sum[i])
        end
    elseif energyModel == "ventilation_rational"
        # Ventilation energy with rational tension law

        tens_sum = B̄*((edgeLengths.^2 .- edgeLengths.^(-3))./6)
        for i = 1:nCells
            energyTotal += UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPressures[i],tens_sum[i])
        end
    end

    return energyTotal
end

function energyCells(params, matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        B̄,
        edgeLengths,
        cellPressures,
        edgel₀s = matrices
    @unpack energyModel,
    l₀ = params

    nCells = length(cellAreas)
    energyPerCell = zeros(nCells)

    if energyModel == "ventilation"
        # Ventilation energy
        tens_sum = B̄*(edgeLengths.-1).^2
        for i = 1:nCells
            energyPerCell[i] = UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPressures[i],tens_sum[i])
        end
    elseif energyModel == "ventilation_rational"
        # Ventilation energy with rational tension law

        tens_sum = B̄*((edgeLengths.^2 .- edgeLengths.^(-3))./5)
        for i = 1:nCells
            energyPerCell[i] = UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPressures[i],tens_sum[i])
        end
    end

    return energyPerCell
end
function tempMeanEnergy(edgeLengths, cellAreas, pressures, nCells)
    
    tens_sum = B̄*((edgeLengths.^2 .- edgeLengths.^(-3))./6)
    
    energyPerCell = UVᵢ.(cellAreas, edgeLengths, 1, pressures,tens_sum)
return energyPerCell
export energy, tempMeanEnergy, energyCells
end

end
