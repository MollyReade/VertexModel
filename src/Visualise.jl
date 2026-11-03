#
#  Visualise.jl
#  VertexModel
#
#

module Visualise

# Julia packages
using Printf
using LinearAlgebra
using ColorSchemes
using Colors
using UnPack
using GeometryBasics
using Random
using Makie
using CairoMakie
using StaticArrays
using SparseArrays
using CircularArrays
using FromFile
using DrWatson
using Statistics

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "AnalysisFunctions.jl" using AnalysisFunctions

function visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks, strainLimits, pressureLimits; limits = (-1.0, 1.0, -1.0, 1.0))

    #TODO Change this to be post processing so correct limits can be used for the whole simulation, by loading each FramaData which contains params and matrices, so limits can be calculated.


    @unpack cellEdgeCount,
        cellαᵢs,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        cellPressures,
        edgeMidpoints,
        edgeLengths,
        edgel₀s,
        F,
        FEdges,
        edgeMidpointLinks,
        μ = matrices
    @unpack nEdges,
        nVerts,
        nCells,
        Pᵢ,
        Amp,
        ω,
        boundaryCondition,
        P₀,
        ipModel,
        l₀ = params

    empty!(ax)

    # Modify to fix frame size, so that it doesn't change as the tissue expands and contracts
    ax.limits = limits

    if boundaryCondition == "displacement"
        ax.title = "t = $(@sprintf("%.3f", t)), s(t) = $(@sprintf("%.3f",(Amp/ω)*(1-cos(ω*t)))))" #(Amp/ω)*(1-cos(ω*(t)))
    elseif boundaryCondition == "force" && ipModel == "sinusoidal"
        ax.title = "t = $(@sprintf("%.3f", t)), Pₚ = $(@sprintf("%.3f",P₀ + Amp*sin(t*ω))))" #P₀ + Amp * sin(t)
    else
        ax.title = "t = $(@sprintf("%.3f", t)), Pₚ = $(@sprintf("%.3f",P₀)))"
    end

    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R, params, matrices)
        minPressure = minimum(cellPressures)
        maxPressure = maximum(cellPressures)
        println(cellPressures)
        for i = 1:nCells
            color = pressureToColour(cellPressures[i], pressureLimits[1], pressureLimits[2])
            #strokeColour = stressToColour(edgeLengths[i]/l₀, 0, 7) 
            # Problem: colours all edges on the polygon the same colour, but we want to colour each edge based on its own strain. So we need to plot each edge separately, rather than the whole polygon at once.

            poly!(ax, cellPolygons[i], color=color, strokecolor=(:black,1.0), strokewidth=2)
        end
        # Calculate edges
        edges = matrices.A * R  
        # give directions not start and end points, so actual edge start and end points are edgeMidpoints - 0.5*edges and edgeMidpoints + 0.5*edges
        edges2 = [Point{2,Float64}.(edgeMidpoints .- 0.5*edges) Point{2,Float64}.(edgeMidpoints .+ 0.5*edges)]
        #println(edges2)
        for j = 1:nEdges
            edgeColour = stressToColour(edgeLengths[j]/edgel₀s[j], 0, strainLimits[2])
            # Plot each edge separately with its own colour based on its strain
            lines!(ax, edges2[j,:], color=edgeColour, linewidth=3)
        end
        Colorbar(fig[1,2], limits = (pressureLimits[1],pressureLimits[2]), colormap = :viridis, label = "Pressure", height= Relative(0.5))
        # Need Labelling of colorbar to match usage - 0 is minimum, 1 in middle and 7 is maximum. So need to adjust position of ticks on colorbar
        Colorbar(fig[1,3], limits = (0,strainLimits[2]), colormap = :seismic, label = "Edge Strain", height= Relative(0.5), ticks = ([0, strainLimits[2]/2, strainLimits[2]], ["0", "1", string(strainLimits[2])]))
    end

    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax, Point{2,Float64}.(R), color=:green)
        #annotations!(ax, string.(collect(1:length(R))), Point{2,Float64}.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax, Point{2,Float64}.(edgeMidpoints), color=:blue)
        # Label each edge with its unstressed length
        #annotation!(ax,Point{2,Float64}.(edgeMidpoints), text .= string.(round.(edgel₀s, sigdigits=3)),  color=:blue)
    end

    # Scatter cell positions
    if scatterCells == 1
        scatter!(ax, Point{2,Float64}.(cellPositions), color=:red)
        # Label each cell with its resistance
        
        #annotation!(ax,Point{2,Float64}.(cellPositions), text .= string.(cellαᵢs),  color=:red)
    end

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    # FIXME Forces out of date
    if plotForces == 1
        arrows2d!(ax, Point{2,Float64}.(R), Vec2f.(sum(F.+FEdges, dims=2)), color=:green)
    end

    if plotEdgeMidpointLinks == 1
        for i = 1:nCells
            for j = 1:cellEdgeCount[i]
                lines!(ax,
                    Point{2,Float64}.([edgeMidpoints[cellEdgeOrders[i][j]],(edgeMidpoints[cellEdgeOrders[i][j]] .+ edgeMidpointLinks[i, cellVertexOrders[i][j]])]),
                    linestyle=:dot,
                    color=:black)
            end
        end
    end

    # Set limits
    #add_ruler!(ax,matrices,params,R)
    reset_limits!(ax)

    # Add frame to movie 
    recordframe!(mov)

    return nothing

end

function add_ruler!(ax, matrices, params,R; xpos_frac=0.1, ypos_frac=0.01)

    @unpack edgeLengths = matrices
    @unpack nRows = params
    
    xrange = extrema(r[1] for r in R)

    xstart =  xpos_frac 
    xend = xstart + xrange[end]
    ypos = -Statistics.mean(edgeLengths)*(nRows÷2 + 1.5)


    #lines!(ax,[xstart, xend], [ypos - 0.01 * xrange[end], ypos + 0.01 * xrange[end]], color=:red, linewidth=5)
    lines!(ax,[xstart, xstart + 7.00], [ypos - 0.01, ypos + 0.01], color=:red, linewidth=5)
    #TODO find a better way to calculate width of ruler used for whole simulation
    lines!(ax,[xstart - 7.00, xstart], [ypos - 0.01, ypos + 0.01], color=:white)

    text!(ax, string(round(xstart+7, sigdigits=3)), position=((xstart+xend)/2, ypos + 0.02 * xrange[end]), color=:red, align = (:center, :bottom), fontsize=20)
end

function pressureToColour(pressure, minPressure, maxPressure)
    range = maxPressure - minPressure

    if !isfinite(range) || abs(range) < 1e-12
        return get(ColorSchemes.viridis, 0.5)  # neutral middle colour
    end
    normPressure = clamp((pressure - minPressure) / (maxPressure - minPressure), 0, 1)
    return get(ColorSchemes.viridis, normPressure)
end

function stressToColour(stress, minStress, maxStress)
    range = maxStress - minStress

    if !isfinite(range) || abs(range) < 1e-12
        return get(ColorSchemes.seismic, 0.5)  # neutral middle colour
    end
    # This maps evenly accross the range, but I want stress=1 to be at the middle of the colour scale, so I need to adjust the mapping accordingly, with stress < 1 being mapped to the lower half of the colour scale and stress > 1 being mapped to the upper half of the colour scale.

    if stress < 1
        normStress = clamp((stress - minStress) / (1 - minStress), 0, 0.5)
    else
        normStress = clamp(0.5 + (((stress - 1) / (maxStress - 1)) * 0.5), 0.5, 1)
    end
    #normStress = clamp((stress - minStress) / (maxStress - minStress), 0, 1)
    return get(ColorSchemes.seismic, normStress)
end

export visualise

end
