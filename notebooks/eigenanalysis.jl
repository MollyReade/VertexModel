# There are as many Lₚ eigenmodes as cells, so ran using a small 5 row simulation, with sinusoidal pleural pressure and edge cells, κ = 0.6, P₀ = 0.7, amplitude = 0.3

using FromFile
using CairoMakie
set_theme!(theme_latexfonts())
using Printf
using JLD2
using ColorSchemes
using GeometryBasics

using VertexModel

dataFolder = "C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\data\\sims\\Full-2D-Experiment\\ventilation_rational\\26-09-29-09-51-46_Pᵢ=0.0_l₀=1.0_nCells=19_realTimetMax=50.0\\frameData"

@from "C:\\Users\\mmr20\\OneDrive\\Documents\\GitHub\\VertexModel\\src\\AnalysisFunctions.jl" using AnalysisFunctions 


function plot_mode_grid(vecs,vals; n_modes=19)
    rows = ceil(n_modes/3)
    fig = Figure(size=(1200, rows*400))

    for k in 1:n_modes
        val = vals[k]
        ax = Axis(fig[(k-1) ÷ 3,(k-1) % 3],
            title = @sprintf("Mode %d (Eigenvalue %.3f)",k, val ),xlabel = "Index, i",ylabel="Eigenvector value i", titlesize=17, xlabelsize = 16, ylabelsize = 16
            ,xticks = 1:(n_modes), limits = (0.9,(n_modes)+0.1,-1.1,1.1)
        )

        mode = vecs[:,k]
        mode ./= maximum(abs.(mode))

        x = 1:length(mode)
        y = zeros(length(mode))

        scatter!(ax, x, y, color=:black)

        colors = [m > 0 ? :red : :blue for m in mode]

        arrows2d!(ax, x, y, zeros(length(mode)),mode,color = colors)
    end

    return fig
end

function plot_spatial_mode(ax, R, matrices, mode, mode_idx, eigenvalue; colorrange)


    for cell_idx in axes(matrices.B, 1)
        polygon = Point2f.(R[matrices.cellVertexOrders[cell_idx]])
        poly!(ax, polygon;
            color=mode[cell_idx],
            colormap=:bluesreds,
            colorrange=colorrange,
            strokecolor=:black,
            strokewidth=1
        )
    end
    scatter!(ax, Point2f.(matrices.cellPositions), color=:black, markersize=5)

    ax.title = @sprintf("Mode %d (Eigenvalue %.3f)", mode_idx, eigenvalue)
    ax.titlesize = 20
    ax.aspect = DataAspect()
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidedecorations!(ax)
    return ax
end

function plot_spectrum(eigenvalues)
    n_modes = length(eigenvalues)

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
        xlabel="Mode index",
        ylabel="Eigenvalue",
        title="Lₚ eigenvalue spectrum"
    )
    mode_indices = 1:n_modes
    scatter!(ax, mode_indices, real.(eigenvalues), color=:black, markersize=9)
    return fig
end

matrices, params, R = load("$(dataFolder)\\systemData$(lpad(string(99),3,'0')).jld2", "matrices", "params", "R")

eigenvectors, eigenvalues = eigenmodesLp(params,matrices)
println(eigenvalues)
println(eigenvectors)

# Normalise eigenvectors so that the maximum absolute value of each eigenvector is 1
for v in eachcol(eigenvectors)
    if v[1] < 0
        v .= -v
    end
    v ./= maximum(abs.(v))
end

# Want spatial plots to be in one figure
nModes = length(eigenvalues)
nRows = 5
nCols = 4

spatialModes = real.(eigenvectors[:,1:nModes])

spatialScale = maximum(abs.(spatialModes))
spatialScale = spatialScale > 0 ? spatialScale : 1.0

spatialFigure = Figure(resolution = (1400, 1800))

for mode_idx in 1:nModes
    row = (mode_idx - 1) ÷ nCols + 1
    col = (mode_idx - 1) % nCols + 1

    ax = Axis(spatialFigure[row, col])
    plot_spatial_mode(ax, R, matrices, spatialModes[:, mode_idx], mode_idx, eigenvalues[mode_idx]; colorrange=(-spatialScale, spatialScale))
end

Colorbar(spatialFigure[1:nRows, nCols+1], colormap=:bluesreds, limits=(-spatialScale, spatialScale), label="Eigenvector value", labelsize=20, ticklabelsize=16)

# for i in 1:length(eigenvalues)
#     mode_idx = i
#     fig_spatial = plot_spatial_mode(R, matrices, eigenvectors, eigenvalues; mode_idx=mode_idx)
    
#     display(fig_spatial)
# end

fig_spectrum = plot_spectrum(eigenvalues)
display(fig_spectrum)
display(spatialFigure)
# mode_idx = 2
# fig_spatial = plot_spatial_mode(R, matrices, eigenvectors, eigenvalues; mode_idx=mode_idx)
# fig_spectrum = plot_spectrum(eigenvalues; mode_idx=mode_idx)
# display(fig_spatial)
# display(fig_spectrum)

# for v in eachcol(eigenvectors)
#     if v[1] < 0
#         v .= -v
#     end
# end 
# figGrid = plot_mode_grid(eigenvectors,eigenvalues)
# figGrid[begin-1,:] = Label(figGrid,"Lₚ",fontsize=25, font = "Serif")
# display(figGrid)

# eigenvectors2, eigenvalues2 = eigenmodesLa(params,matrices)
# println(eigenvalues2)

# for v in eachcol(eigenvectors2)
#     if v[1] < 0
#         v .= -v
#     end
# end 
# figGrid2 = plot_mode_grid(eigenvectors2,eigenvalues2)
# figGrid2[begin-1,:] = Label(figGrid2,"Lₐ",fontsize=25, font = "Serif")
# display(figGrid2)