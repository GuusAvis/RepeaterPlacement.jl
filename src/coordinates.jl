"""
Numerical types that are allowed within this module.
`StochasticTriple` is needed for compatibility with automatic differentiation.
"""
const AllowedElTypes = Union{Float64, T} where T <: StochasticAD.StochasticTriple

"""
Vector types that are allowed within this module.
"""
const AllowedVecTypes = Vector{T} where T <: AllowedElTypes

"""
Matrix types that are allowed within this module.
"""
const AllowedMatTypes = Matrix{T} where T <: AllowedElTypes

"""
    Coordinates(end_node_coordinate_matrix, repeater_coordinate_matrix)

A collection of end-node and repeater coordinates.

Can be either initialized by passing coordinate matrices for the end nodes and repeaters
explicitly, or by creating an empty object first (by ommitting arguments from the
constructor) and then adding them using `add_end_node` and `add_repeater`.

If the matrices only contain the coordinates of a single node, they can also be passed
as a vector.
"""
struct Coordinates{T <: AllowedMatTypes}
    end_node_coordinate_matrix::Matrix{Float64}
    repeater_coordinate_matrix::T
end

empty_mat = Matrix{Float64}(undef, 0, 0)

function Coordinates(end_node_coordinate_matrix=empty_mat,
        repeater_coordinate_matrix=empty_mat
    )
    if ! (end_node_coordinate_matrix isa AllowedMatTypes)
        end_node_coordinate_matrix = convert(Matrix{Float64}, end_node_coordinate_matrix)
    end
    if ! (repeater_coordinate_matrix isa AllowedMatTypes)
        repeater_coordinate_matrix = convert(Matrix{Float64}, repeater_coordinate_matrix)
    end
    Coordinates(end_node_coordinate_matrix, repeater_coordinate_matrix)
end

vec_to_mat(x) = reshape([i for i in Iterators.flatten(x)], :, 1)
function Coordinates(end_node_coordinate_matrix::Vector,
        repeater_coordinate_matrix=empty_mat)
    Coordinates(vec_to_mat(end_node_coordinate_matrix), repeater_coordinate_matrix)
end
function Coordinates(end_node_coordinate_matrix, repeater_coordinate_matrix::Vector)
    Coordinates(end_node_coordinate_matrix, vec_to_mat(repeater_coordinate_matrix))
end
function Coordinates(end_node_coordinate_matrix::Vector, repeater_coordinate_matrix::Vector)
    Coordinates(vec_to_mat(end_node_coordinate_matrix),
        vec_to_mat(repeater_coordinate_matrix))
end


"""
    end_nodes(x::Coordinates)

Return coordinate matrix with end-node coordinates only.
"""
end_nodes(x::Coordinates) = x.end_node_coordinate_matrix

"""
    repeaters(x::Coordinates)

Return coordinate matrix with repeater coordinates only.
"""
repeaters(x::Coordinates) = x.repeater_coordinate_matrix

Base.copy(x::Coordinates) = Coordinates(copy(end_nodes(x)), copy(repeaters(x)))

"""
    nodes(x::Coordinates)

Return the coordinate matrix of all nodes (end nodes and repeaters).
"""
function nodes(x::Coordinates)
    isempty(repeaters(x)) && return end_nodes(x)
    isempty(end_nodes(x)) && return repeaters(x)
    [end_nodes(x) repeaters(x)]
end

"""
    num_end_nodes(x::Coordinates)

Return the number of end nodes represented by a `Coordinates` object.
"""
num_end_nodes(x::Coordinates) = size(end_nodes(x), 2)

"""
    num_repeaters(x::Coordinates)

Return the number of repeaters represented by a `Coordinates` object.
"""
num_repeaters(x::Coordinates) = size(repeaters(x), 2)

"""
    num_nodes(x::Coordinates)

Return total number of nodes (end nodes and repeaters) represented by `Coordinates` object.
"""
num_nodes(x::Coordinates) = num_end_nodes(x) + num_repeaters(x)

"""
    end_node(x::Coordinates, i)

Return the coordinates of the ith end node as a vector.
"""
end_node(x::Coordinates, i) = end_nodes(x)[:, i]

"""
    repeater(x::Coordinates, i)

Return the coordinates of the ith repeater as a vector.
"""
repeater(x::Coordinates, i) = repeaters(x)[:, i]

"""
    node(x::Coordinates, i)

Return the coordinates of the ith node as a vector. First indices are end nodes.
"""
node(x::Coordinates, i) = nodes(x)[:, i]

"""
    add_end_node(x::Coordinates, coord)

Add an end node to a `Coordinates` object with coordinate vector `coord`.
"""
function add_end_node(x::Coordinates, coord)
    new = isempty(end_nodes(x)) ? [coord;;] : [end_nodes(x) coord]
    Coordinates(new, repeaters(x))
end

"""
    add_repeater(x::Coordinates, coord)

Add a repeater to a `Coordinates` object with coordinate vector `coord`.
"""
function add_repeater(x::Coordinates, coord)
    new = isempty(repeaters(x)) ? [coord;;] : [repeaters(x) coord]
    Coordinates(end_nodes(x), new)
end

"""
    distance(p1, p2)

Calculate the Euclidean distance between coordinate vectors `p1` and `p2`.
"""
distance(p1, p2) = sqrt(sum((p1 .- p2) .^ 2))

"""
    nearest_neighbour(coord, coords)

Find the nearest neighbour to a coordinate in a list of coordinates.

# Arguments
- `coord`: The coordinate to find the nearest neighbour for.
- `coords`: A collection of coordinates to find the nearest neighbour from;
    should be a matrix where each column is a coordinate vector.
    Can also be a `Coordinates` object; in that case, the coordinate matrix is automatically
    extracted.

# Returns
- `nearest_neighbour_index`: The index of the nearest neighbour
    (the index of the column in the coordiante matrix, such that when `coords` is a
    `Coordinates` its coordinates can be retrieved using
    `node(coords, nearest_neighbour_index)`).
- `dist`: The Euclidean distance between `coord` and the nearest neighbour.

# Examples
```julia-repl
julia> coords = [[0, 0] [10, 0]]
2×2 Matrix{Int64}:
 0  10
 0   0

julia> coord = [3, 0]
2-element Vector{Int64}:
 3
 0

julia> nearest_neighbour(coord, coords)
(1, 3.0)

julia> coordinates_nearest_neighbour = coords[:, 1]
2-element Vector{Int64}:
 0
 0
```
"""
function nearest_neighbour(coord, coords)
    nearest_neighbour_index = nothing
    dist_to_nearest_neighbour = Inf
    for (neighbour_index, neighbour_coords) in enumerate(eachcol(coords))
        dist = distance(coord, neighbour_coords)
        if dist < dist_to_nearest_neighbour
            nearest_neighbour_index = neighbour_index
            dist_to_nearest_neighbour = dist
        end
    end
    nearest_neighbour_index, dist_to_nearest_neighbour
end

nearest_neighbour(coord, coords::Coordinates) =
    nearest_neighbour(coord, nodes(coords))

"""
    max_dist(coords)

Find the largest distance between any two nearest neighbours in a collection of coordinates.

# Arguments
- `coords`: a matrix where every column is a coordinate vector.
    Can also be a `Coordinates` object; in that case, the coordinate matrix is automatically
    extracted.
"""
function max_dist(coords)
    sol = 0
    for i in axes(coords, 2)
        # we change the matrix every step, this may be inefficient
        # is there a more efficient way to exclude a node from its own set of neighbors?
        coord = coords[:, i]
        other_coords = [coords[:, begin:i-1] coords[:, i+1:end]]
        _, dist = nearest_neighbour(coord, other_coords)
        dist > sol && (sol = dist)
    end
    sol
end

max_dist(x::Coordinates) = max_dist(nodes(x))

"""
    adjacency_matrix(coords, radius=Inf)

Build the adjacency matrix for a geometric graph with a given radius.

# Arguments
- `coords`: a matrix where every column is a coordinate vector.
    If `coords` is a `Coordinates` object, the coordinate matrix is automatically extracted.
- `radius`: the maximum distance between nodes for them to be connected.

The matrix is symmetric where `matrix[i, j]` is the distance between node i and node j
if they are connected, and 0 if they are not (i.e., if the distance would have exceeded
the radius).
"""
function adjacency_matrix(coords::Matrix, radius=Inf)
    num_coords = size(coords, 2)
    mat = zeros(eltype(coords), num_coords, num_coords)
    for i = 1:num_coords
        for j = 1:(i - 1)
            dist = distance(coords[:, i], coords[:, j])
            dist <= radius && (mat[i, j] = mat[j, i] = dist)
        end
    end
    mat
end

adjacency_matrix(coords::Coordinates, radius=Inf) =
    adjacency_matrix(nodes(coords), radius)

"""
    build_graph(coords, radius=Inf)

Return a graph with Euclidean distances as weights.

# Arguments
- `coords`: a matrix where every column is a coordinate vector.
    If `coords` is a `Coordinates` object, the coordinate matrix is automatically extracted.
- `radius`: the maximum distance between nodes for them to be connected.
"""
function build_graph(coords, radius=Inf)
    SimpleWeightedGraph(adjacency_matrix(coords, radius))
end

"""
    build_waxman_graph(coords::Coordinates, beta=0.4, alpha=0.1, L=nothing)

Create a graph from the node coordinates probabilistically using the Waxman formula.

A Waxman graph is a graph where nodes are assigned coordinates uniformly at random,
and an edge is added probabilistically between each pair of nodes with probability
```math
p = \beta \exp(-\frac{d}{\alpha L})
```
where `d` is the distance between the nodes, `L` is a scaling factor (default is the maximum
distance between any two nodes), and `\beta` and `\alpha` are parameters that control the
edge addition probability.

This function builds a `SimpleWeightedGraph` where edges are added probabilistically using
the above formula according to the distances between nodes in `coords`.

This function has been written to be consistent with the [implementation in the NetworkX
Python library](https://networkx.org/documentation/stable/reference/generated/networkx.generators.geometric.waxman_graph.html).

# References
- Waxman, S. J. (1990). "Routing of multipoint connections in a packet-switched network."
  IEEE Journal on Selected Areas in Communications, 8(9), 1558-1567.

"""
function build_waxman_graph(coords::Coordinates, beta=0.4, alpha=0.1, L=nothing)
    g = SimpleWeightedGraph(num_nodes(coords))
    node_pairs = [(i, j) for i in 1:num_nodes(coords) for j in 1:i - 1]
    distances = [distance(node(coords, i), node(coords, j)) for (i, j) in node_pairs]
    isnothing(L) && (L = maximum(distances))
    for ((i, j), d) in zip(node_pairs, distances)
        p = beta * exp(-d / alpha / L)
        rand() < p && add_edge!(g, i, j, d)
    end
    g
end

"""
    initialize_line(num_reps, dist=100)

Initialize `num_reps` repeaters placed randomly on the line between two end nodes separated
by a distance `dist` km.
"""
function initialize_line(num_reps, dist=100, rng=Random.default_rng())
    end_node_coords = Float64[[0, 0] [dist, 0]]
    initial_sols = [rand(rng) * dist for _ in 1:num_reps]
    initial_sol = [[i for i in initial_sols] [0 for _ in initial_sols]]'
    Coordinates(end_node_coords, initial_sol)
end

"""
    initialize_square(num_reps, dist=100)

Initialize `num_reps` repeaters placed randomly in a square of end nodes, where the edges
of the square are `dist` km long.
"""
function initialize_square(num_reps, dist=100, regular=false, rng=Random.default_rng())
    end_node_coords = Float64[[0, 0] [dist, dist] [0, dist] [dist, 0]]
    if regular
        num_reps_on_grid = floor(sqrt(num_reps)) ^ 2
        num_reps_left_out = num_reps - num_reps_on_grid
        spacing = dist / (sqrt(num_reps_on_grid) + 1)
        initial_sols =  []
        for i in 1:sqrt(num_reps_on_grid), j in 1:sqrt(num_reps_on_grid)
            append!(initial_sols, [spacing * i, spacing * j])
        end
        append!(initial_sols, [rand(rng) * dist for _ in 1:(2num_reps_left_out)])
    else
        initial_sols = [rand(rng) * dist for _ in 1:(2num_reps)]
    end
    initial_sol = reshape(initial_sols, 2, :)
    Coordinates(end_node_coords, initial_sol)
end

"""
    initialize_random(num_end_nodes, num_reps, scale=100)

Initialize `num_end_nodes` end nodes and `num_reps` repeaters, both placed randomly.

All coordinates are picked uniformly at random within the square of side length `scale`.
"""
function initialize_random(num_end_nodes, num_reps, scale=100, rng=Random.default_rng())
    coords = Coordinates()
    for _ in 1:num_end_nodes
        coords = add_end_node(coords, [rand(rng) * scale, rand(rng) * scale])
    end
    for _ in 1:num_reps
        coords = add_repeater(coords, [rand(rng) * scale, rand(rng) * scale])
    end
    coords
end

"""
    waxman_graph(num_end_nodes, num_reps, beta=0.4, alpha=0.1, L=1.)

Create a Waxman graph and a corresponding `Coordinates` object.

A Waxman graph is a graph where nodes are assigned coordinates uniformly at random,
and an edge is added probabilistically between each pair of nodes with probability
```math
p = \beta \exp(-\frac{d}{\alpha L})
```
where `d` is the distance between the nodes, `L` is a scaling factor (default is the maximum
distance between any two nodes), and `\beta` and `\alpha` are parameters that control the
edge addition probability.

This function creates a Waxman graph by first creating a `Coordinates` object with
`num_end_nodes` end nodes and `num_reps` repeaters, each placed uniformly at random
using `initialize_random`.
Then a corresponding `SimpleWeightedGraph` is built using `build_waxman_graph`.

This function has been written to be consistent with the [implementation in the NetworkX
Python library](https://networkx.org/documentation/stable/reference/generated/networkx.generators.geometric.waxman_graph.html).

# References
- Waxman, S. J. (1990). "Routing of multipoint connections in a packet-switched network."
  IEEE Journal on Selected Areas in Communications, 8(9), 1558-1567.

"""
function waxman_graph(num_end_nodes, num_reps, beta=0.4, alpha=0.1, L=1.)
    coords = initialize_random(num_end_nodes, num_reps, L)
    g = build_waxman_graph(coords, beta, alpha, L)
    g, coords
end
