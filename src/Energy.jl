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
        cellL₀s,
        μ,
        Γ ,
        B̄,
        edgeLengths,
        cellPressures,
        edgel₀s,
        cellPᵢs= matrices
    @unpack energyModel,
    l₀ = params
    print(stacktrace())
    if energyModel == "log"
        # Logarithmic energy
        energyTotal
        for i = 1:nCells
            energyTotal += Uᵢ.(cellAreas[i], cellA₀s[i], cellPerimeters[i], cellL₀s[i], μ[i], Γ[i])
        end
    elseif energyModel == "ventilation"
        # Ventilation energy
        tens_sum = B̄*(edgeLengths.-1).^2
        energyTotal 
        for i = 1:nCells
            
            energyTotal += UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPᵢs[i],tens_sum[i])
        end
    elseif energyModel == "ventilation_rational"
        # Ventilation energy with rational tension law

        tens_sum = B̄*((edgeLengths.^2 .- edgeLengths.^(-3))./6)
        for i = 1:nCells
            energyTotal += UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPᵢs[i],tens_sum[i])
        end
    else
        # Quadratic energy
        energyTotal = sum(μ.*(0.5 .* (cellAreas .- cellA₀s).^2 .+ 0.5 .* Γ .* (cellPerimeters .- cellL₀s).^2))
    end

    return energyTotal
end

function energyCells(params, matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        cellL₀s,
        μ,
        Γ ,
        B̄,
        edgeLengths,
        cellPressures,
        edgel₀s,
        cellPᵢs= matrices
    @unpack energyModel,
    l₀ = params

    nCells = length(cellAreas)
    energyPerCell = zeros(nCells)

    if energyModel == "log"
        # Logarithmic energy
        for i = 1:nCells
            energyPerCell[i] = Uᵢ.(cellAreas[i], cellA₀s[i], cellPerimeters[i], cellL₀s[i], μ[i], Γ[i])
        end
    elseif energyModel == "ventilation"
        # Ventilation energy
        tens_sum = B̄*(edgeLengths.-1).^2
        for i = 1:nCells
            energyPerCell[i] = UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPᵢs[i],tens_sum[i])
        end
    elseif energyModel == "ventilation_rational"
        # Ventilation energy with rational tension law

        tens_sum = B̄*((edgeLengths.^2 .- edgeLengths.^(-3))./6)
        for i = 1:nCells
            energyPerCell[i] = UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPᵢs[i],tens_sum[i])
        end
    else
        # Quadratic energy
        for i = 1:nCells
            energyPerCell[i] = μ[i]*(0.5 *(cellAreas[i] - cellA₀s[i])^2 + 0.5 * Γ[i] * (cellPerimeters[i] - cellL₀s[i])^2)
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
