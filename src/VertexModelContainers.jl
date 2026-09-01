#
#  VertexModelContainers.jl
#  VertexModel
#
#  Structs to hold system parameters and matrices, allowing us avoid passing dozens of arguments to each function

module VertexModelContainers

using SparseArrays
using StaticArrays
using Random
using Distributions
using CircularArrays

@kwdef mutable struct ParametersContainer
    initialSystem      ::String             # System used for initialising simulations
    nCells             ::Int64              # Number of cells
    nEdges             ::Int64              # Number of edges
    nVerts             ::Int64              # Number of vertices
    γ                  ::Float64            # Parameter in energy relaxation
    λ                  ::Float64            # Parameter in energy relaxation
    Λ                  ::Float64            # Parameter in ventilation
    κ                  ::Float64            # Parameter in vertex dissipation
    L₀                 ::Float64            # Cell preferred perimeter length L₀ = -λ/(2*γ)
    A₀                 ::Float64            # Cell preferred area
    l₀                 ::Float64            # Edge preferred length
    Pᵢ                 ::Float64            # Internal pressure per cell
    Pₘ                 ::Float64            # Pressure at mouth
    P₀                 ::Float64            # Pressure at system periphery
    Amp                ::Float64            # Amplitude of sinusoidal internal pressure 
    boundaryToggle     ::Bool              # Flag controlling whether boundary conditions are applied
    edgeCellsToggle    ::Bool              # Flag controlling whether edge cells are included in the system
    outputTotal        ::Int64              # Total number of data outputs
    outputInterval     ::Float64            # Non dimensionalised data output interval
    viscousTimeScale   ::Float64            # Relaxation rate, approx from Sarah's data.
    realTimetMax       ::Float64            # Dimensionalised run time in seconds
    tMax               ::Float64            # Non dimensionalised run time
    realCycleTime      ::Float64            # Cell cycle time in seconds
    nCycles            ::Float64            # Number of cell cycles in simulation 
    nRows              ::Int64              # Number of rows of cells in initial system
    nonDimCycleTime    ::Float64            # Non dimensionalised cell cycle time
    t1Threshold        ::Float64            # Length of edge below which a T1 transition occurs
    peripheralTension  ::Float64            # Tension at system periphery
    seed               ::Int64              # Random number seed 
    rng                ::MersenneTwister    # Random number sequence
    distLogNormal      ::LogNormal{Float64} # Log normal distribution 
    energyModel        ::String             # Choice of energy function (eg log, quadratic)
    resistanceModel    ::String             # Choice of resistance function (eg linear, quadratic)
    folderName         ::String             # Name of directory in which simulation results are stored
    currentTime        ::Float64            # Current simulation time
    ipModel            ::String             # Choice of internal pressure model (eg sinusoidal, constant)
end

@kwdef mutable struct MatricesContainer
    A                ::SparseMatrixCSC{Int64, Int64}                # Incidence matrix mapping edges to vertices. Rows => edges; columns => vertices.
    B                ::SparseMatrixCSC{Int64, Int64}                # Incidence matrix mapping cells to edges. Rows => cells; columns => edges. (Values +/-1 for orientation)
    Aᵀ               ::SparseMatrixCSC{Int64, Int64}                # Transpose of incidence matrix A
    Ā                ::SparseMatrixCSC{Int64, Int64}                # Undirected adjacency matrix mapping edges to vertices (abs.(A))
    Āᵀ               ::SparseMatrixCSC{Int64, Int64}                # Transpose of adjacency matric Ā
    Bᵀ               ::SparseMatrixCSC{Int64, Int64}                # Transpose of incidence matrix B
    B̄                ::SparseMatrixCSC{Int64, Int64}                # Undirected adjacency matrix mapping cells to edges (abs.(B))
    B̄ᵀ               ::SparseMatrixCSC{Int64, Int64}                # Transpose of adjacency matric B̄
    C                ::SparseMatrixCSC{Int64, Int64}                # Adjacency matrix mapping cells to vertices. Rows => cells; Columns => vertices. (C=0.5*B̄*Ā)
    areaJacobian     ::SparseMatrixCSC{SVector{2, Float64}, Int64}  # Sparse matrix of Jacobian determinants for each cell, indexed by cell and vertex
    cellEdgeCount    ::Vector{Int64}                                # Vector of number of edges around each cell
    cellVertexOrders ::Vector{CircularVector{Int64, Vector{Int64}}} # Storing the ordering of vertices around the cell
    cellEdgeOrders   ::Vector{CircularVector{Int64, Vector{Int64}}} # Storing the ordering of edges around the cell
    boundaryVertices ::Vector{Int64}                                # Vector of 1s and 0s denoting vertices that lie on the system boundary
    boundaryEdges    ::Vector{Int64}                                # Vector of 1s and 0s denoting edges that lie on the system boundary
    cellPositions    ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors for each cell centre of mass
    cellPerimeters   ::Vector{Float64}                              # Vector of scalar cell perimeter lengths
    cellOrientedAreas::Vector{SMatrix{2, 2, Float64}}               # Vector of oriented cell areas. Each row is a 2x2 antisymmetric static matrix of the form [0 A / -A 0] where A is the scalar cell area
    cellAreas        ::Vector{Float64}                              # Vector of scalar cell areas
    cellA₀s          ::Vector{Float64}                              # Vector of scalar cell preferred areas
    cellL₀s          ::Vector{Float64}                              # Vector of scalar cell preferred perimeters
    cellαᵢs          ::Vector{Float64}                              # Vector of scalar cell shape factors
    cellTensions     ::Vector{Float64}                              # Vector of boundary tensions for each cell
    cellPressures    ::Vector{Float64}                              # Vector of time left until division for each cell
    μ                ::Vector{Float64}                              # Vector of cell stiffness factors 
    Γ                ::Vector{Float64}                              # Vector of factors determinind relative strength of cell tension and internal pressure per cell 
    edgeCellNormals  ::SparseMatrixCSC{SVector{2, Float64}, Int64}  # Sparse array of vectors normal to each edge, indexed by adjacent cell and edge for each normal vector
    edgeLengths      ::Vector{Float64}                              # Vector of lengths for each edge in the system
    edgeTangents     ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors containing edge length and direction as a 2D vector
    edgeTensions    ::Vector{Float64}                              # Vector of scalar edge tensions
    edgeMidpoints    ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors containing edge midpoints as (x,y) positions
    edgel₀s          ::Vector{Float64}                              # Vector of scalar edge preferred lengths
    edgeMidpointLinks::SparseMatrixCSC{SVector{2, Float64}, Int64}  # Sparse array of vectors connecting adjacent edge midpoints, indexed by adjacent vertex and cell for each midpoint link 
    timeSinceT1      ::Vector{Float64}                              # Vector of times since each edge last underwent a T1 transition
    normEdges        ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors containing unit vectors along each edge
    vertexAreas      ::Vector{Float64}                              # Vector of areas of triangles surrounding vertices    
    F                ::SparseMatrixCSC{SVector{2, Float64}, Int64}  # Matrix of 2D static vectors containing force vectors acting on each vertex and cell
    FEdges           ::SparseMatrixCSC{SVector{2, Float64}, Int64}  # Matrix of 2D static vectors containing force vectors acting on each vertex from edges
    externalF        ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors containing total force applied to each vertex by external pressure
    totalF           ::Vector{SVector{2, Float64}}                  # Vector of 2D static vectors containing resultant force vectors acting on each vertex
    ϵ                ::SMatrix{2, 2, Float64, 4}                    # Antisymmetric rotation matrix
    cellShapeTensor  ::Vector{SMatrix{2, 2, Float64}}               # Shape tensor of a cell
end

export ParametersContainer,MatricesContainer

end
