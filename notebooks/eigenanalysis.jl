# There are as many eigenmodes as vertices, so ran using a small 5 row simulation, with sinusoidal pleural pressure and edge cells, κ = 0.6, P₀ = 0.7, amplitude = 0.3

using FromFile
using CairoMakie
set_theme!(theme_latexfonts())
using Printf
using JLD2

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

matrices, params = load("$(dataFolder)\\systemData$(lpad(string(99),3,'0')).jld2", "matrices", "params")

eigenvectors, eigenvalues = eigenmodesLp(params,matrices)
println(eigenvalues)

for v in eachcol(eigenvectors)
    if v[1] < 0
        v .= -v
    end
end 
figGrid = plot_mode_grid(eigenvectors,eigenvalues)
figGrid[begin-1,:] = Label(figGrid,"Lₚ",fontsize=25, font = "Serif")
display(figGrid)

eigenvectors2, eigenvalues2 = eigenmodesLa(params,matrices)
println(eigenvalues2)

for v in eachcol(eigenvectors2)
    if v[1] < 0
        v .= -v
    end
end 
figGrid2 = plot_mode_grid(eigenvectors2,eigenvalues2)
figGrid2[begin-1,:] = Label(figGrid2,"Lₐ",fontsize=25, font = "Serif")
display(figGrid2)