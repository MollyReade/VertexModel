#
#  InitialSystemLayout.jl
#  VertexModel
#
#  Function to create a hexagonal grid of cells. 
#  Given number of rows nRows, central row has length nRows, each adjacent row has length nRows-1 etc. 
#  Number of cells is then nRows*(nRows-1) - (floor(Int64, nRows/2)+1)*(floor(Int64, nRows/2)+2) + nRows

module InitialSystemLayout

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using FromFile
using Random

@from "SenseCheck.jl" using SenseCheck

function initialSystemLayout(;
        nRows = 3,
        boundaryToggle = 0,
        initialEdgeLength = 5.0*0.75/6, # Need to find a better value for this than 5*L₀/6
        edgeCellsToggle = 0
    )

    # nRows = 9 # Must be an odd number
    cellPoints = [SVector(x, 0.0) for x = 1:nRows]
    for j = 1:(floor(Int64,nRows/2))
        for i = 1:nRows-j
            # Need to add a small amount of randomness to prevent errors in voronoi tessellation 
            push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, j * sqrt(1 - 0.5^2) + rand() * 0.001 - 0.0005))
            push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, -j * sqrt(1 - 0.5^2) + rand() * 0.001 - 0.0005))
        end
    end
    xs = [x[1] for x in cellPoints]
    ys = [x[2] for x in cellPoints]

    ptsArray = zeros(Float64, (2, length(cellPoints)))
    for (i, point) in enumerate(cellPoints)
        ptsArray[1, i] = point[1]
        ptsArray[2, i] = point[2]
    end

    
    triangulation_unconstrained = triangulate(ptsArray)
    tessellation_constrained = voronoi(triangulation_unconstrained, clip=true)
    
    #Exclude points outside constraining boundary
    usableVertices = Int64[]
    for a in values(tessellation_constrained.polygons)
        push!(usableVertices, a...)
    end
    sort!(unique!(usableVertices))
    outerVertices = setdiff(collect(1:num_polygon_vertices(tessellation_constrained)), usableVertices)
    
    #println(usableVertices)
    
    # Map vertex indices in tessellation to vertex indices in incidence matrices (after excluding outer vertices)
    vertexIndexingMap = Dict(usableVertices .=> collect(1:length(usableVertices)))

    Rtmp = SVector.(tessellation_constrained.polygon_points[usableVertices])

    # Find pairs of vertices connected by edges in tessellation 
    # Use incidence matrix indexing for vertices, and exclude outer vertices 
    pairs = [(vertexIndexingMap[p[1]], vertexIndexingMap[p[2]]) for p in keys(tessellation_constrained.adjacent.adjacent) if p[1] ∈ usableVertices && p[2] ∈ usableVertices]
    # Ensure lowest index is first in tuple, and remove duplicates 
    orderedPairs = unique([(min(p...), max(p...)) for p in pairs])

    nVerts = length(Rtmp)
    nEdges = length(orderedPairs)
    nCells = length(cellPoints)

    # Construct A matrix mapping tessellation edges to tessellation vertices 
    A = spzeros(Int64, nEdges, nVerts)
    for (edgeIndex, vertices) in enumerate(orderedPairs)
        A[edgeIndex, vertices[1]] = 1
        A[edgeIndex, vertices[2]] = -1
    end

    # NB get_polygon(tessellation_constrained,x) or tessellation_constrained.polygons[x] return indices of vertices around cell x ordered anti-clockwise, with first and last element the same

    # Construct B matrix mapping voronoi cell around each fibril to surrounding edges between vertices in tessellation
    # NB assume ϵᵢ is a clockwise rotation so cell orientation is into page. 
    B = spzeros(Int64, nCells, nEdges)
    for c = 1:nCells
        for i = 2:length(tessellation_constrained.polygons[c])
            vertexLeading = vertexIndexingMap[tessellation_constrained.polygons[c][i-1]]  # Leading with respect to *clockwise* direction around cell
            vertexTrailing = vertexIndexingMap[tessellation_constrained.polygons[c][i]]
            # Find index of edge connecting these vertices 
            edge = (findall(x -> x != 0, @view A[:, vertexLeading])∩findall(x -> x != 0, @view A[:, vertexTrailing]))[1]
            if A[edge, vertexLeading] > 0
                B[c, edge] = 1
            else
                B[c, edge] = -1
            end
        end
    end

    # Prune peripheral vertices with 2 edges that both belong to the same cell
    # Making the assumption that there will never be two such vertices adjacent to each other
    
    verticesToRemove = Int64[]
    edgesToRemove = Int64[]
    for i = 1:nVerts
        edges = findall(x -> x != 0, @view A[:, i])
        cells1 = findall(x -> x != 0, @view B[:, edges[1]])
        cells2 = findall(x -> x != 0, @view B[:, edges[2]])
        if cells1 == cells2
            # If the lists of cells to which both edges of vertex i belong are identical, this implies that the edges are peripheral and only belong to one cell, so edge i should be removed.
            push!(verticesToRemove, i)
            push!(edgesToRemove, edges[1])
        end
    end
    for i in verticesToRemove
        edges = findall(x -> x != 0, @view A[:, i])
        otherVertexOnEdge1 = setdiff(findall(x -> x != 0, @view A[edges[1], :]), [i])[1]
        A[edges[2], otherVertexOnEdge1] = A[edges[2], i]
        A[edges[1], otherVertexOnEdge1] = 0
    end
    A = A[setdiff(1:size(A, 1), edgesToRemove), setdiff(1:size(A, 2), verticesToRemove)]
    B = B[:, setdiff(1:size(B, 2), edgesToRemove)]
    Rtmp = Rtmp[setdiff(1:size(Rtmp, 1), verticesToRemove)]

    

    senseCheck(A, B; marker="Removing peripheral vertices")


    R = SVector{2, Float64}[]
    for r in Rtmp 
        push!(R, SVector(initialEdgeLength*(r[1] - (nRows-1)/2 - 1.0 ), initialEdgeLength*r[2]))
    end

    if edgeCellsToggle == 0
        A, B, R = removeBoundaryCells(A,B,R)
    end

    return A, B, R

end

function removeBoundaryCells(A,B,R)

    nVerts = size(A)[2]
    nEdges = size(A)[1]
    nCells = size(B)[1]

    #Find how many cells each edge belongs to
    edgesP = vec((abs.(ones(nCells)'*B)))
    #Find which vertices belong to peripheral edges
    vertsP = vec((0.5.*edgesP'*abs.(A))')
    #Edges connected to boundary vertices
    edgesB = vec(abs.(A*(ones(nVerts).-vertsP)))
    #Edges without peripheral edges or edges connected to boundary vertices
    edgesI = vec(ones(nEdges).-edgesB.-edgesP)
    #1 for interior vertices, 0 for boundary vertices
    vertsI = vec(ones(nVerts) - vertsP)
    #1 if cell has a boundary edge, 0 otherwise
    cellsB = vec(abs.(abs.(B)*edgesP))
    cellsB[findall(x -> x != 0, cellsB)] .= 1
    #1 if interior cell, 0 otherwise
    cellsI=ones(nCells)-cellsB

    #A,B,R with only interior cells, edges and vertices
    Bint = B[cellsI.>0,edgesI.>0]
    Aint = A[edgesI.>0,vertsI.>0]
    Rint = R[vertsI.>0]
    return Aint,Bint,Rint

end

function initialSmallConfig(n)
    if n=="one"
        ATmp= [-1.0 1.0 0.0 0.0 0.0 0.0
             0.0 -1.0 1.0 0.0 0.0 0.0
             0.0 0.0 -1.0 1.0 0.0 0.0
             0.0 0.0 0.0 -1.0 1.0 0.0
             0.0 0.0 0.0 0.0 -1.0 1.0
             1.0 0.0 0.0 0.0 0.0 -1.0]
        A = sparse(ATmp)
        B = -1.0.*ones(1,6)
        R = Array{SVector{2,Float64}}(undef,6)
        for k = 1:6
            R[k] = SVector{2}([cos((k*π)/3.0),sin((k*π)/3.0)])
        end

        R .*= 1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))
    end
    return A,B,R
end

export initialSystemLayout, initialSmallConfig 

end