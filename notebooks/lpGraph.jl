using VertexModel
using UnPack
using LinearAlgebra
using CairoMakie
using Statistics


pressures = LinRange(-0.5,1.1,500)
stability = Dict{Float64,String}()
dus = Dict{Float64,Float64}()
edge_Lengths = Float64[]

for P in pressures
    integ,_,_,_ = vertexModel(initialSystem = "one",
        nRows = 3,
        realCycleTime = 2*86400.0,
        l₀ = 0.15,
        A₀ = 1.0,
        Pᵢ = P,
        boundaryToggle = 0,
        edgeCellsToggle = 0,
        outputToggle = 0,
        frameDataToggle = 0,
        frameImageToggle = 0,
        printToggle = 0,
        videoToggle = 0,
        energyModel = "ventilation_rational")

    (params, matrices) = integ.p
    u_prev = copy(integ.u)
    

    du = similar(integ.u)
    integ.f(du,integ.u, integ.p, integ.t)
    @unpack edgeLengths = matrices
    push!(edge_Lengths, mean(edgeLengths))
    dus[P] = norm(du)
    println("Pᵢ = $P, mean edge length = $(mean(edgeLengths)), ||du|| = $(norm(du))")
    if norm(du) < 1e-8
        stability[P] = "STABLE"
    else
        stability[P] = "UNSTABLE"
        if norm(du) > 1e+3
            for P in pressures
                if haskey(stability, P)
                    continue
                else
                    stability[P] = "UNSTABLE"
                    push!(edge_Lengths, mean(edgeLengths))
                end
                
            end
            break
        end
    end
end

stable_P = Float64[]
stable_L = Float64[]

unstable_P = Float64[]
unstable_L = Float64[]

for (i,P) in enumerate(pressures)
    if stability[P] == "STABLE"
        push!(stable_P, P)
        push!(stable_L,  edge_Lengths[i])
    else
    
        push!(unstable_P, P)
        push!(unstable_L,  edge_Lengths[i])
    
    end
end    

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Cell Pressure Pᵢ",
    ylabel = "Mean Edge Length l",
    title = "Mean Edge Length vs Cell Pressure"
)

scatter!(ax, stable_P, stable_L; color = :blue, label = "Stable")
scatter!(ax, unstable_P, unstable_L; color = :red, label = "Unstable")

axislegend(ax; position = :rt)
fig
