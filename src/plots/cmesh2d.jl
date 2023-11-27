"""
    Pt

Cartesian 2D point.

# Fields
- `x::Float64`: x coordinate
- `y::Float64`: y coordinate
"""
@with_kw mutable struct Pt
    x::Float64 = 0.0
    y::Float64 = 0.0
end


"""
    Facet

Polygon facet.

# Fields
- `vertices::Vector{Pt}`: vertices of the polygon
- `index::Int`: array index number in the mesh
- `col::Int`: column number in the mesh
- `row::Int`: row number in the mesh
- `label::String`: text, if any, appears for this facet
- `anchor::Pt`: point, if any, to which this facet belongs
- `occupants::Vector{Pt}`: one or more points occupy this facet
"""
@with_kw mutable struct Facet
    vertices::Vector{Pt} = []
    index::Int = 0
    col::Int = 0
    row::Int = 0
    label::String = ""
    anchor::Pt = Pt()
    occupants::Vector{Pt} = []
end


"""
    CMesh2D

2D mesh (or grid) of facets.

# Fields
- `origin::Pt`: origin of the mesh (lower left corner)
- `Ncol::Int`: number of columns
- `Nrow::Int`: number of rows
- `Δw::Float64`: width of each column
- `Δh::Float64`: height of each row
- `facets::Vector{Facet}`: mesh of facets (see Notes)

# Notes
- A mesh comprise several non-ovelapping facets, where each facet
  may contain zero or more points (`occupants`).
- Each facet may have only zero or one text label. A label
  describes the `anchor` point.

  (row no. ↑)
       +---+---+---+---+---+
    4  | 16| 17| 18| 19| 20|
       +---+---+---+---+---+
    3  | 11| 12| 13| 14| 15|
       +---+---+---+---+---+
    2  | 6 | 7 | 8 | 9 | 10|
       +---+---+---+---+---+
    1  | 1 | 2 | 3 | 4 | 5 |
       +---+---+---+---+---+
         1   2   3   4   5    (col no. →)

    where +---+
          |   |
          +---+  is one facet.

"""
@with_kw mutable struct CMesh2D
    origin::Pt = Pt()
    Ncol::Int = 0
    Nrow::Int = 0
    Δw::Float64 = 0.0
    Δh::Float64 = 0.0
    facets::Vector{Facet} = []
end


"""
    CMesh2D(dim, oxy, Δ)

Constructor of CMesh2D.

# Arguments
- `dim`: mesh dimensions: no. of columns and rows
- `oxy`: mesh origin: x- and y-coordinates of the lower left corner of the mesh
- `Δ`: mesh width and height of each column and row, respectively

# Returns
- `CMesh2D`
"""
function CMesh2D(dim, oxy, Δ)
    origin = Pt(oxy[1], oxy[2])
    Ncol, Nrow = dim[1], dim[2]
    Δw, Δh = Δ[1], Δ[2]

    v::Vector{Pt} = []
    for r ∈ 1:Nrow+1, c ∈ 1:Ncol+1
        y = origin.y + (r - 1) * Δh
        x = origin.x + (c - 1) * Δw
        push!(v, Pt(x, y))
    end

    facets = map(1:(Ncol * Nrow)) do i
        q, r = divrem(i, Ncol)
        row = q + (r > 0)
        col = Ncol - (Ncol - r) % Ncol
        Δ = Ncol + 1
        i1 = (row - 1) * Δ + col
        i2 = i1 + Δ
        i3 = i2 + 1
        i4 = i1 + 1
        vertices = map(i -> Pt(v[i].x, v[i].y), [i1, i2, i3, i4])
        Facet(vertices=vertices, index=i, col=col, row=row)
    end

    CMesh2D(
        origin = origin,
        Ncol = Ncol,
        Nrow = Nrow,
        Δw = Δw,
        Δh = Δh,
        facets = facets
    )
end


"""
    distance(pt1::Pt, pt2::Pt)

Distance between two given points.

    x <------> x
    pt1        pt2

# Arguments
- `pt1::Point`: first point
- `pt2::Point`: second point

# Returns
- `Float64`: distance between both points
"""
distance(pt1::Pt, pt2::Pt) = sqrt((pt2.x - pt1.x)^2 + (pt2.y - pt1.y)^2)


"""
    mean_dist(pt::Pt, vert::AbstractVector{Pt})

Mean distance between a point and a polygon/facet.

            b    c
   x   ⇚    +---+
       ⇚    |   |
       ⇚    +---+
            a   d

where ⇚ is the mean distance between all points (a, b, c, and d)
of the vertices and the point.

# Arguments
- `pt::Point`: point
- `vert::AbstractVector{Pt}`: polygon vertices

# Returns
- `Float64`: mean distance between the point and polygon
"""
function mean_dist(pt::Pt, vert::AbstractVector{Pt})
    sum(distance.(Ref(pt), vert)) / length(vert)
end


"""
    min_dist(pt::Pt, vert::AbstractVector{Pt})

Shortest distance between a point and a polygon/facet.

             b   c
   x <-----> +---+
             |   |
             +---+
            a    d

where the distance between x and point b is returned, as point b
is closest to the point.

# Arguments
- `pt::Point`: point
- `vert::AbstractVector{Pt}`: polygon vertices

# Returns
- `Float64`: shortest (minimum) distance between the point and polygon
"""
function min_dist(pt::Pt, vert::AbstractVector{Pt})
    minimum(distance.(Ref(pt), vert))
end


"""
    mindex(mesh::CMesh2D, col::Int, row::Int)

Facet number, given its column and row number.

# Arguments
- `mesh::CMesh2D`: mesh
- `col::Int`: column number of facet
- `row::Int`: row number of facet

# Notes
- Facet number begins from lower left corner, following this flow:

  (row no. ↑)
       +---+---+---+---+---+
    4  | 16| 17| 18| 19| 20|
       +---+---+---+---+---+
    3  | 11| 12| 13| 14| 15|
       +---+---+---+---+---+
    2  | 6 | 7 | 8 | 9 | 10|
       +---+---+---+---+---+
    1  | 1 | 2 | 3 | 4 | 5 |
       +---+---+---+---+---+
         1   2   3   4   5    (col no. →)

  where column 2 and row 4 returns facet number 17.

# Returns
- `Int`: facet number
"""
mindex(mesh::CMesh2D, col::Int, row::Int) = (row - 1) * mesh.Ncol + col


"""
    mpos(mesh::CMesh2D, pt::Pt)

Facet column and row number, given a point in the mesh.

  (row no. ↑)
       +---+---+---+---+---+
    4  |   |   |   |   |   |
       +---+---+---+---+---+
    3  |   |   |   |   |   |
       +---+---+---+---+---+
    2  |   |   | x |   |   |
       +---+---+---+---+---+
    1  |   |   |   |   |   |
       +---+---+---+---+---+
         1   2   3   4   5    (col no. →)

  where point x is the mesh is located in facet at column 3 and row 2.

# Arguments
- `mesh::CMesh2D`: mesh
- `pt::Pt`: point

# Returns
- `NamedTuple(:col, :row){Int, Int}`: facet column and row number
"""
function mpos(mesh::CMesh2D, pt::Pt)
    col = Int(floor((pt.x - mesh.origin.x) / mesh.Δw)) + 1
    row = Int(floor((pt.y - mesh.origin.y) / mesh.Δh)) + 1
    col = min(mesh.Ncol, col)
    row = min(mesh.Nrow, row)
    (; col, row)
end


"""
    mpos(mesh::CMesh2D, idx::Int)

Facet column and row number, given its facet number.

  (row no. ↑)
       +---+---+---+---+---+
    4  | 16| 17| 18| 19| 20|
       +---+---+---+---+---+
    3  | 11| 12| 13| 14| 15|
       +---+---+---+---+---+
    2  | 6 | 7 | 8 | 9 | 10|
       +---+---+---+---+---+
    1  | 1 | 2 | 3 | 4 | 5 |
       +---+---+---+---+---+
         1   2   3   4   5    (col no. →)

  where facet number 17 is at column 2 and row 4.

# Arguments
- `mesh::CMesh2D`: mesh
- `idx::Int`: facet number

# Returns
- `NamedTuple(:col, :row){Int, Int}`: facet column and row number
"""
function mpos(mesh::CMesh2D, idx::Int)
    Ncol = mesh.Ncol
    q, r = divrem(idx, Ncol)
    row = q + (r > 0)
    col = Ncol - (Ncol - r) % Ncol
    (; col, row)
end


"""
    Base.getindex(mesh::CMesh2D, idx::Int)

Return the facet given its facet number.

# Arguments
- `mesh::CMesh2D`: mesh
- `idx::Int`: facet number

# Returns
- `Facet`: `Facet` object
"""
Base.getindex(mesh::CMesh2D, idx::Int) = mesh.facets[idx]


"""
    Base.getindex(mesh::CMesh2D, pt::Pt)

Return the `Facet` object given a point in the mesh.

# Arguments
- `mesh::CMesh2D`: mesh
- `pt::Pt`: point in mesh

# Returns
- `Facet`: `Facet` object
"""
Base.getindex(mesh::CMesh2D, pt::Pt) = mesh.facets[mindex(mesh, mpos(mesh, pt)...)]


"""
    Base.getindex(mesh::CMesh2D, col::Int, row::Int)

Return the `Facet` object given a facet column and row number in the mesh.

# Arguments
- `mesh::CMesh2D`: mesh
- `col::Int`: column number of facet
- `row::Int`: row number of facet

# Returns
- `Facet`: `Facet` object
"""
Base.getindex(mesh::CMesh2D, col::Int, row::Int) = mesh.facets[mindex(mesh, col, row)]


"""
    offsets(mesh::CMesh2D, order::Int, col::Int, row::Int)

Find all pairs of column and row numbers of neighboring facets (INTERNAL USE).

# Arguments
- `mesh::CMesh2D`: mesh
- `order::Int`: how far away from reference to find neighboring facets
- `col::Int`: column number of reference facet
- `row::Int`: row number of reference facet

# Returns
- `Vector{Tuple{Int, Int}}`: column and row numbers of neigboring facets
"""
function offsets(mesh::CMesh2D, order::Int, col::Int, row::Int)
    ar = []
    for c ∈ -order:order, r ∈ -order:order
        if (abs(c) == order) || (abs(r) == order)
            Δc, Δr = col + c, row + r
            if (1 <= Δc <= mesh.Ncol) && (1 <= Δr <= mesh.Nrow)
                push!(ar, (Δc, Δr))
            end
        end
    end
    ar
end


"""
    get_neighbors(order::Int, mesh::CMesh2D, col::Int, row::Int)

Return all neighboring facets given a facet column and row number in the mesh.

# Arguments
- `order::Int`: number of facets away from reference, where order 1 is one
                facet away from reference, 2 is two facets away, 3 is three
                facets away, and so on.
- `mesh::CMesh2D`: mesh
- `col::Int`: column number of reference facet
- `row::Int`: row number of reference facet

# Notes
- Relative to reference facet X, the facet order is as follows:

       +---+---+---+---+---+---+---+
       | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 2 | 2 | 2 | 2 | 2 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 2 | 1 | 1 | 1 | 2 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 2 | 1 | X | 1 | 2 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 2 | 1 | 1 | 1 | 2 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 2 | 2 | 2 | 2 | 2 | 3 |
       +---+---+---+---+---+---+---+
       | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
       +---+---+---+---+---+---+---+

   where the numbers are the orders, with reference to facet X. Orders 1, 2,
   and 3 indicate their facets are located 1, 2, and 3 facets away, respectively,
   from the reference facet X.

# Returns
- `Vector{Facet}`: nearest/neighboring `Facet` objects
"""
function get_neighbors(order::Int, mesh::CMesh2D, col::Int, row::Int)
    ofs = offsets(mesh, order, col, row)
    map(o -> getindex(mesh, o[1], o[2]), ofs)
end


"""
    get_neighbors(order::Int, mesh::CMesh2D, pt::Pt)

Return all neighboring facets given a point in the mesh.

# Arguments
- `order::Int`: number of facets away from reference, where order 1 is one
                facet away from reference, 2 is two facets away, 3 is three
                facets away, and so on.
- `mesh::CMesh2D`: mesh
- `pt::Pt`: reference point in the mesh

# Returns
- `Vector{Facet}`: nearest/neighboring `Facet` objects
"""
function get_neighbors(order::Int, mesh::CMesh2D, pt::Pt)
    get_neighbors(order, mesh, mpos(mesh, pt)...)
end


"""
    centroid(vert::AbstractVector{Pt})

Center of polygon/facet.

# Arguments
- `vert::AbstractVector{Pt}`: polygon vertices

# Returns
- `Pt`: center of polygon
"""
function centroid(vert::AbstractVector{Pt})
    n = length(vert)
    a = 0.0
    i1 = 2
    for i ∈ 1:n
        a += vert[i].x * vert[i1].y - vert[i1].x * vert[i].y
        i1 = (i1 % n) + 1
    end

    a *= 0.5
    cx = cy = 0.0
    i1 = 2
    for i ∈ 1:n
        t = vert[i].x * vert[i1].y - vert[i1].x * vert[i].y
        cx += (vert[i].x + vert[i1].x) * t
        cy += (vert[i].y + vert[i1].y) * t
        i1 = (i1 % n) + 1
    end

    cx /= (6.0 * a)
    cy /= (6.0 * a)
    Pt(cx, cy)
end
