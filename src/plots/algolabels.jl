"""
    is_available(facet::Facet)

Return `true` if a given facet has no points and no label.

# Arguments
- `facet::Facet`: facet to test

# Returns
- `Bool`: `true` if facet has no points and no label, else `false`.
"""
is_available(facet::Facet) = (length(facet.occupants) == 0) && isempty(facet.label)


"""
    has_label(facet::Facet)

Return `true` if a given facet has label.

# Arguments
- `facet::Facet`: facet to test

# Returns
- `Bool`: `true` if facet has label, else `false` for no label.
"""
has_label(facet::Facet) = !isempty(facet.label)


"""
    reset_place!(mesh::CMesh2D)

Clear a mesh so that all its facet labels, anchor points, and scores are reset.
Does not clear other information in `CMesh2D`. (INTERNAL USE)

# Arguments
- `mesh::CMesh2D`: mesh

# Returns
- nothing
"""
function reset_place!(mesh::CMesh2D)
    facets = [f for f ∈ mesh.facets if has_label(f)]
    foreach(facets) do f
        f.label = ""
        f.anchor = Pt()
    end
end


"""
    save_places(mesh::CMesh2D)

Return a copy of all facets that have labels (INTERNAL USE).

# Arguments
- `mesh::CMesh2D`: mesh

# Returns
- `Vector{Facet}`: list of saved facets that have labels
"""
function save_places(mesh::CMesh2D)
    map([f for f ∈ mesh.facets if has_label(f)]) do f
        Facet(index = f.index,
              col = f.col,
              row = f.row,
              label = f.label,
              anchor = f.anchor)
    end
end


"""
    load_places(mesh::CMesh2D, saved_places::AbstractVector{Facet})

Load a copy of all facets that have labels into current mesh (INTERNAL USE).

# Arguments
- `mesh::CMesh2D`: mesh
- `saved_places`: copy of all `Facets` to load into mesh

# Returns
- nothing
"""
function load_places!(mesh::CMesh2D, saved_places::AbstractVector{Facet})
    reset_place!(mesh)
    foreach(saved_places) do p
        facet = mesh[p.index]
        facet.label = p.label
        facet.anchor = p.anchor
    end
end


"""
    swap_place!(lhs::Facet, rhs::Facet)

Swap label and anchor information between two given facets.

# Arguments
- `lhs::Facet`: facet to swap (left facet <-> right facet)
- `rhs::Facet`: facet to swap (left facet <-> right facet)

# Returns
- nothing
"""
function swap_place!(lhs::Facet, rhs::Facet)
    lhs.label, rhs.label = rhs.label, lhs.label
    lhs.anchor, rhs.anchor = rhs.anchor, lhs.anchor
end


"""
    place!(mesh::CMesh2D, pts::AbstractVector{Pt},
           labels::AbstractVector{<:AbstractString}, max_order::Int)

Attempt to find the facets to place labels, so that the distance
between the facets and their respective anchor points are as close
to another as possible (INTERNAL USE).

# Arguments
- `mesh::CMesh2D`: mesh
- `pts::AbstractVector{Pt}`:
- `labels::AbstractVector{<:AbstractString}`: list of labels.
- `max_order::Int`: maximum number of facets away a label can be
                    placed from its anchor point

# Returns
- `Int`: total fit score
"""
function place!(mesh::CMesh2D, pts::AbstractVector{Pt},
                labels::AbstractVector{<:AbstractString}, max_order::Int)
    reset_place!(mesh)
    score = 0
    for i ∈ randperm(length(pts))
        pt = pts[i]
        order = 1
        idx = -1

        while (idx < 0) && (order <= max_order)
            nf = [f for f ∈ get_neighbors(order, mesh, pt) if is_available(f)]
            d = map(f -> mean_dist(pt, f.vertices), nf)
            if length(d) > 0
                idx = findmin(d)[2]
                n = nf[idx]
                n.label = labels[i]
                n.anchor = pt
                score += order - 1
                break
            end
            order += 1
        end
    end
    score
end


"""
    find_place!(mesh::CMesh2D, xy, labels::AbstractVector{<:AbstractString},
                max_order::Int, niter::Int)

Find where in the mesh to place all labels such that they do not
overlap one another, yet are as close as needed to their corresponding
anchor points.

# Arguments
- `mesh::CMesh2D`: mesh
- `xy`: list of points
- `labels::AbstractVector{<:AbstractString}`: list of labels for the points
- `max_order::Int`: maximum number of facets away a label can be placed
                    from its anchor point
- `niter::Int`: maximum attempts to place labels

# Notes
- Labels are placed in such a way to minimize the total (overall) score.
- Best `score` is zero. This indicates all facets are immediately
  next to their corresponding `anchor` points. In practice, this means every
  point in the chart has its corresponding label just beside it.
- The higher the score value, the farther away the labels are from their
  `anchor` points.

# Returns
- `Tuple{Int, Int}`: no. of iterations to complete label placement and best score.
"""
function find_places!(mesh::CMesh2D, xy, labels::AbstractVector{<:AbstractString},
                      max_order::Int, niter::Int)
    x, y = xy[1], xy[2]
    n = length(x)
    pts = map(i -> Pt(x[i], y[i]), 1:n)

    best_score = 999999999
    saved_places = []
    i = 1
    while (i <= niter) && (best_score > 0)
        score = place!(mesh, pts, labels, max_order)
        if score < best_score
            best_score = score
            saved_places = save_places(mesh)
        end
        i += 1
    end

    length(saved_places) > 0 && load_places!(mesh, saved_places)
    (iterations=i, score=best_score)
end


"""
    labelsize(txt::AbstractString, fontpt::Int, resolution, rgx, rgy; inflate=1.0)

Text dimensions (width and height). This dimension will be used to create the
mesh.

# Arguments
- `txt::AbstractString`: label text
- `fontpt::Int`: font size (points)
- `resolution`: figure size (pixels)
- `rgx`: extrema (range) of the x values
- `rgy`: extrema (range) of the y values

# Keywords
- `inflate`: by how much to increase/decrease the text dimension (default: 1.0 for
             no inflation)

# Returns
- `NamedTuple{(:width, :height), {Float64, Float64)}`: width and height of text label
"""
function labelsize(txt::AbstractString, fontpt::Int, resolution, rgx, rgy; inflate=1.0)
    wpx = resolution[1]
    hpx = resolution[2]

    sz = []
    let
        # create a dummy figure just to get the text dimensions:
        fig = Figure()
        axs = Axis(fig[1, 1], width=wpx, height=hpx)
        t = text!(axs, position=(rgx[1], rgy[1]), txt, fontsize=fontpt)
        sz = boundingbox(t).widths      # width and height of text
    end

    Δx = (rgx[2] - rgx[1])
    Δy = (rgy[2] - rgy[1])
    width = sz[1] / wpx * Δx * inflate
    height = sz[2] / hpx * Δy * inflate
    (; width, height)
end


"""
    fill_mesh!(mesh::CMesh2D, x::AbstractVector, y::AbstractVector,
               labels::AbstractVector{<:AbstractString};
               niter::Int,
               max_order::Int)

Place non-overlapping labels as close to their anchor points.

# Arguments
- `x::AbstractVector`: list of x points
- `y::AbstractVector`: list of y points
- `labels::AbstractVector{<:AbstractString}`: list of labels for every (x, y) point

# Keywords
- `niter::Int`: maximum attempts to place labels
- `max_order::Int`: maximum number of facets away a label can be placed
                    from its anchor point

# Returns
- `CMesh2D`: mesh
"""
function fill_mesh!(mesh::CMesh2D, x::AbstractVector, y::AbstractVector,
                    labels::AbstractVector{<:AbstractString};
                    niter::Int,
                    max_order::Int)

    for z ∈ zip(x, y)
        pt = Pt(z[1], z[2])
        push!(mesh[pt].occupants, pt)
    end

    find_places!(mesh, (x, y), labels, niter, max_order)
end


"""
    create_mesh(rgx, rgy,
                labels::AbstractVector{<:AbstractString};
                resolution,
                flabel,
                fontpt::Int)

Construct a blank 2D mesh/grid.

# Arguments
- `rgx`: range of x (smallest and largest x values)
- `rgy`: range of y (smallest and largest y values)
- `labels::AbstractVector{<:AbstractString}`: list of labels

# Keywords
- `resolution`: figure size (width x height in pixels)
- `flabel`: by how much to increase/decrease the text dimension for labels
- `fontpt::Int`: font size (points) of text labels

# Returns
- `CMesh2D`: mesh
"""
function create_mesh(rgx, rgy,
                     labels::AbstractVector{<:AbstractString};
                     resolution,
                     flabel,
                     fontpt::Int)

    nmax = findmax([length(lbl) for lbl ∈ labels])[2]   # maximum label width
    Δw, Δh = labelsize(labels[nmax], fontpt, resolution, rgx, rgy; inflate=flabel)

    ncols = Int(ceil((rgx[2] - rgx[1]) / Δw))
    nrows = Int(ceil((rgy[2] - rgy[1]) / Δh))
    CMesh2D((ncols, nrows), (rgx[1], rgy[1]), (Δw, Δh))
end


"""
    draw_mesh!(axs, mesh::CMesh2D)

Draw a mesh.

# Arguments
- `axs`: chart axis
- `mesh::CMesh2D`: mesh

# Returns
- nothing
"""
function draw_mesh!(axs, mesh::CMesh2D)
    color = :red
    lw = 0.5
    ls = :dash

    Ncol, Nrow = mesh.Ncol, mesh.Nrow
    v0 = mesh[1, 1].vertices[1]
    v1 = mesh[Ncol, 1].vertices[4]
    lines!(axs, [v0.x, v1.x], [v0.y, v1.y], color=color, linewidth=lw, linestyle=ls)
    v1 = mesh[1, Nrow].vertices[2]
    lines!(axs, [v0.x, v0.x], [v0.y, v1.y], color=color, linewidth=lw, linestyle=ls)

    for r ∈ 1:Nrow
        v0 = mesh[1, r].vertices[2]
        v1 = mesh[Ncol, r].vertices[3]
        lines!(axs, [v0.x, v1.x], [v0.y, v1.y],
                color=color, linewidth=lw, linestyle=ls)
    end

    for c ∈ 1:Ncol
        v0 = mesh[c, 1].vertices[4]
        v1 = mesh[c, Nrow].vertices[3]
        lines!(axs, [v0.x, v0.x], [v0.y, v1.y],
                color=color, linewidth=lw, linestyle=ls)
    end
end


"""
    segments_intersect(ls1, ls2)

Determine whether two line segments intersect each other.

# Arguments
- `ls1`: first line segment, comprising a pair of `Pt` end points
- `ls1`: second line segment, comprising a pair of `Pt` end points

# Returns
- `Bool`: `true` if both line segments intersect each other, else `false`
"""
function segments_intersect(ls1, ls2)
    s10_x = ls1[2].x - ls1[1].x
    s10_y = ls1[2].y - ls1[1].y
    s32_x = ls2[2].x - ls2[1].x
    s32_y = ls2[2].y - ls2[1].y

    denom = s10_x * s32_y - s32_x * s10_y
    isapprox(denom, 0) && return false    # collinear
    denom_p = (denom > 0)

    s02_x = ls1[1].x - ls2[1].x
    s02_y = ls1[1].y - ls2[1].y
    s_numer = s10_x * s02_y - s10_y * s02_x
    (s_numer < 0) == denom_p && return false   # no collision

    t_numer = s32_x * s02_y - s32_y * s02_x
    (t_numer < 0) == denom_p && return false   # no collision

    !(((s_numer > denom) == denom_p) || ((t_numer > denom) == denom_p))
end


"""
    draw_labels!(axs, mesh::CMesh2D, fdist, fontpt::Int)

Plot or draw the labels on the chart.

# Arguments
- `axs`: chart axis
- `mesh::CMesh2D`: mesh
- `fdist`: inflate the max distance between facet and its anchor point
           before a line will be drawn between them
- `fontpt::Int`: font size (points) for text labels

# Notes
- No lines will be drawn to connect the text labels to their anchor points
  unless labels have a certain distance from their anchor points.
- This distance threshold is taken as the smaller of either the facet width
  or height. Argument `fdist` can extend or shorten this distance.
- Normally, `fdist` should be < 1.0.
- This method attempts to draw non-intersecting lines between labels and
  their corresponding anchor points.

# Returns
- nothing
"""
function draw_labels!(axs, mesh::CMesh2D, fdist, fontpt::Int)
    color = :black
    lw = 0.5

    length_thld = min(mesh.Δw, mesh.Δh) * fdist
    facets = [f for f ∈ mesh.facets if has_label(f)]

    lfacets = []
    for f ∈ facets
        pt = centroid(f.vertices)
        text!(axs, position=(pt.x, pt.y), f.label,
              fontsize=fontpt, align=(:center, :center))

        (min_dist(f.anchor, f.vertices) > length_thld) && push!(lfacets, f)
    end

    n = length(lfacets)
    for i ∈ 1:n
        f1 = lfacets[i]
        ls1 = (centroid(f1.vertices), f1.anchor)
        for f2 ∈ lfacets[[j for j ∈ 1:n if j != i]]
            ls2 = (centroid(f2.vertices), f2.anchor)
            segments_intersect(ls1, ls2) && swap_place!(f1, f2)
        end
    end

    r = 0.5 * min(mesh.Δw, mesh.Δh)
    for f ∈ lfacets
        pt = centroid(f.vertices)
        dt = distance(pt, f.anchor)
        t = (dt - r) / dt
        xt = (1 - t) * f.anchor.x + t * pt.x
        yt = (1 - t) * f.anchor.y + t * pt.y
        lines!(axs, [f.anchor.x, xt], [f.anchor.y, yt],
               color=color, linewidth=lw)
    end
end


"""
    viz(x::AbstractVector, y::AbstractVector,
        labels::AbstractVector{<:AbstractString};
        resolution=(600, 600),
        flabel=1.0,
        fdist=1.0,
        fontpt::Int=13,
        niter::Int=100,
        max_order::Int=3,
        mesh_on::Bool=false)

Plot a scatter chart and non-overlapping labels for each marker.

# Arguments
- `x::AbstractVector`: list of x points
- `y::AbstractVector`: list of y points
- `labels::AbstractVector{<:AbstractString}`: list of labels for every (x, y) point

# Keywords
- `resolution`: figure size (default: 600 x 600 pixels)
- `flabel`: by how much to increase/decrease the text dimension (default: 1.0)
- `fdist`: inflate the max distance between facet and its anchor point
           before a line will be drawn between them (default: 1.0)
- `fontpt::Int`: font size (points) for text labels (default: 13 points)
- `niter::Int`: maximum attempts to place labels (default: 100)
- `max_order::Int`: maximum number of facets away a label can be placed
                    from its anchor point (default: 3)
- `mesh_on::Bool`: whether to draw an overlying mesh over the chart
                   Used for debugging (default: `false` for no mesh drawing)

# Returns
- nothing
"""
function viz(x::AbstractVector, y::AbstractVector,
             labels::AbstractVector{<:AbstractString};
             resolution=(600, 600),
             flabel=1.0,
             fdist=1.0,
             fontpt::Int=13,
             niter::Int=100,
             max_order::Int=3,
             mesh_on::Bool=false)

    fig = Figure(resolution=resolution)
    axs = Axis(fig[1, 1])

    scatter!(axs, x, y)

    mesh = create_mesh(extrema(x), extrema(y), labels;
                       resolution=resolution, flabel=flabel, fontpt=fontpt)

    fill_mesh!(mesh, x, y, labels; niter=niter, max_order=max_order)

    mesh_on && draw_mesh!(axs, mesh)

    draw_labels!(axs, mesh, fdist, fontpt)

    display(fig)
end
