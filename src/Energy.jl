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
W(x) = ((x^3)/3 - (x^(-2)/2))/5
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
    l₀,
    nEdges,
    P₀ = params
    if energyModel == "ventilation"
        # Ventilation energy
        tens_sum = B̄*(edgeLengths.-1).^2
        energyTotal = 0
        for i = 1:nCells
            
            energyTotal += UVᵢ.(cellAreas[i], edgeLengths[i], l₀, cellPressures[i],tens_sum[i])
        end
    elseif energyModel == "ventilation_rational"
        # Ventilation energy with rational tension law

        energyTotal = 0
        for j = 1:nEdges
            energyTotal += l₀*W(edgeLengths[j]/l₀)
        end
        energyTotal -= P₀*sum(cellAreas)
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

    return energyPerCell
end

export energy, energyCells
end


