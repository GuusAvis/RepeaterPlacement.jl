"""
Path through a graph. Specifies both the traversed nodes and the segment lengths.

Paths are non-directional; if a path `p1` is the reverse of a path `p2`, then `p1 == p2`.
"""
struct Path{T <: AllowedVecTypes}
    nodes::Vector{Int}
    lengths::T
    function Path(nodes::Vector{Int}, lengths::T) where T <: AllowedVecTypes
        length(nodes) - 1 == length(lengths) ||
            length(nodes) == length(lengths) == 0 ||
                throw(DimensionMismatch("nodes and lengths don't match"))
        new{T}(nodes, lengths)
    end
end

"""
    edges(p::Path)

Return the edges of the path as a list.
"""
function Graphs.edges(p::Path)
    edges = []
    last_node = p.nodes[1]
    for n in p.nodes[begin+1:end]
        edge = (last_node, n)
        push!(edges, edge)
        last_node = n
    end
    edges
end

"""
    nodes(p::Path)

Return the nodes of the path as a list.
"""
nodes(p::Path) = p.nodes

#TODO define Graphs.vertices?

Path() = Path(Int[], Float64[])

function Path(nodes, graph::AbstractSimpleWeightedGraph)
    for i in 1:length(nodes) - 1
        has_edge(graph, nodes[i], nodes[i + 1]) ||
            throw(ArgumentError("Edge $((nodes[i], nodes[i + 1])) not in graph."))
    end
    weights = [get_weight(graph, nodes[i], nodes[i + 1]) for i in 1:length(nodes) - 1]
    Path(nodes, weights)
end

Path(nodes, coords::Coordinates) = Path(nodes, build_graph(coords))

function Base.:(+)(p1::Path, p2::Path)
    p1.nodes[end] == p2.nodes[begin] || ArgumentError("p2 must start at end node in p1")
    Path([p1.nodes; p2.nodes[begin + 1:end]], [p1.lengths; p2.lengths])
end

Base.isempty(p::Path) = isempty(p.nodes) && isempty(p.lengths)

"""
    hash(p::Path)

Standard hash function for `Path`.
This implementation is consistent with `Base.:(==)(p1::Path, p2::Path)`, which
is important for algorithms which rely on this property (e.g. `Dict`, `Set`).

In particular, since every `Path` is equal to its own reverse, every `Path`
hashes to the same value as its own reverse.
"""
function Base.hash(p::Path)
    """
        lexical_lt(a::Vector, b::Vector)
    
    Lexical ordering function for Vectors.
    
    This function is specific for this use case: we are only handling the case
    where both vectors are equally sized.
    """
    function lexical_lt(a::Vector, b::Vector)
        for i in 1:length(a)
            if a[i] < b[i]
                return true
            end
            if b[i] < a[i]
                return false
            end
        end
        return false
    end

    # Each `Path` is equal to its own reverse, so given `p::Path`, we need to
    # decide if we will hash `p` or `reverse(p)`.

    p_reversed = Path(reverse(p.nodes), reverse(p.lengths))
    nodes_reversion_symmetry = p_reversed.nodes == p.nodes
    need_to_reverse = nodes_reversion_symmetry ? lexical_lt(p_reversed.lengths, p.lengths) : lexical_lt(p_reversed.nodes, p.nodes)
    p_ordered = need_to_reverse ? p_reversed : p

    return hash(p_ordered.nodes, hash(p_ordered.lengths))::UInt
end

Base.:(==)(p1::Path, p2::Path) = (p1.nodes == p2.nodes && p1.lengths == p2.lengths) ||
    (p1.nodes == reverse(p2.nodes) && p1.lengths == reverse(p2.lengths))

"""
    length(p::Path)

Number of segments in a path.
"""
Base.length(p::Path) = length(p.lengths)

"""
    path_length(p::Path)

Total length of a path, i.e., sum of lengths of the edges.
"""
path_length(p::Path) = sum(p.lengths)

"""
    find_all_paths(graph, source, destination, visited=Int[], path=Int[])

Find all paths between `source` and `destination` in `graph`.

Note that a path, by definition, cannot visit any node more than once.
Each path is a vector of node indices (ordered from `source` to `destination`).
This function returns a vector holding all the paths.
"""
function find_all_paths(graph, source, destination, path=Path())
    paths = Vector{Path}[]
    if isempty(path)  # this is the first function call
        path = Path([source], Float64[])
    else
        new_edge = Path([path.nodes[end], source], graph)
        path = path + new_edge
    end
    source == destination && return [path]
    for neighbor in neighbors(graph, source)
        !(neighbor in path.nodes) &&
            (paths = [paths; find_all_paths(graph, neighbor, destination, path)])
    end
    paths
end

"""
    find_k_shortest_paths(graph, source, destination, k)

Find the k shortest `Path`s between `source` and `destination` in `graph`.
"""
function find_k_shortest_paths(graph, source, destination, k)
    yenstate = yen_k_shortest_paths(graph, source, destination, Graphs.weights(graph), k)
    [Path(path, graph) for path in yenstate.paths]
end
"""
    find_paths(graph, source, destination, max_num_paths=nothing)

Find either all `Path`s or the `max_num_paths` shortest `Path`s between points in a graph.
"""
function enumerate_paths(graph, source, destination, max_num_paths=nothing)
    max_num_paths === nothing && return find_all_paths(graph, source, destination)
    find_k_shortest_paths(graph, source, destination, max_num_paths)
end

"""
    sample_path(paths, f, temperature=0., boltzmann_constant=0.)

Return a path and its cost, sampled using the Boltzmann factor.

When the temperature is zero, this function deterministically returns the path with the
lowest cost. Otherwise, the relative weight for sampling path p is given by 
``e^{-cost(p) / temperature}``.

`f` should be a function that takes a `Path` object as its argument and returns
a tuple `(cost, extra_cache)`, where `cost` is the cost of the path and `extra_cache` is a
dictionary that can be used to store additional information about the path.

"""
function sample_path(paths, f; temperature=0., boltzmann_constant=0., kwargs...)
    isempty(paths) && throw(ArgumentError("paths cannot be empty"))
    effective_temperature = boltzmann_constant * temperature
    cost_fct_results = [f(p; temperature=temperature, kwargs...)
        for p in paths]
    costs = [c[1] for c in cost_fct_results]
    extra_caches = [c[2] for c in cost_fct_results]
    if ! iszero(effective_temperature)
        # normalize smallest cost to zero to avoid numerical problems
        costs .-= minimum(costs)
        w = Weights(exp.(-costs / effective_temperature))
        if ! iszero(w)
            path = sample(paths, w)
            cost, extra_cache = f(path)
            path, cost, extra_cache
        end
    end
    cost, index = findmin(costs)
    path = paths[index]
    extra_cache = extra_caches[index]
    path, cost, extra_cache
end

"""
    path_finding_through_enumeration(path_cost_fct, graph, s, t;
        max_num_paths=nothing, prune_fct=identity, kwargs...)

Find the best path by first enumerating paths and then calculating the cost of each.

# Arguments
- `cost_function`: Function that takes a `Path` object and returns two objects:
    1. a `Float64` that represents the cost of the path;
    2. a `Dict{Symbol, Any}` that may contain extra information.
- `graph`: the `SimpleWeightedGraph` through which to find the shortest path.
- `source`: the source node.
- `destination`: the destination node.

# Keyword arguments
- `max_num_paths`: if nothing, all paths are enumerated else, only this many shortest paths.
- `prune_fct`: a function that takes a vector of paths and returns a pruned vector.
    This function is used to select a subset of enumerated paths to calculate the cost for.

# Returns
- `path`: the shortest path from `source` to `destination`, as a `Path` object.
- `cost`: the cost of the path.
- `cache`: a `Dict{Symbol, Any}` containing information about the solution returned
    by the cost function.
"""
function path_finding_through_enumeration(cost_fct, graph, source, destination;
        max_num_paths=nothing, prune_fct=identity, kwargs...)
    paths = prune_fct(enumerate_paths(graph, source, destination, max_num_paths))
    isempty(paths) && return nothing, nothing, nothing
    # isempty(paths) && throw(ErrorException)
    sample_path(paths, cost_fct; kwargs...)
end

"""
    extended_dijkstra(cost_function, graph, source, destination)

Extended Dijkstra algorithm to determine the path that minimizes the cost function.

This implementation is based on the extended Dijkstra algorithm as included in the paper
"Concurrent Entanglement Routing for Quantum Networks: Model and Designs" by
Shouqian Shi and Chen Qian.
While the regular Dijkstra algorithm only works for additive cost functions, this one
works for any cost function that is both monotonic (meaning that the cost of path p + e is
greater than the cost of p for any path p and edge e) and isotonic (meaning that if the
cost of p is better than that of p', then the cost of p + e will be better than that of
p' + e).
In the paper it is claimed that the algorithm works for any monotonic cost function, but
this is not true if the cost function is not isotonic.

# Arguments
- `cost_function`: Function that takes a `Path` object and returns two objects:
    1. a `Float64` that represents the cost of the path;
    2. a `Dict{Symbol, Any}` that may contain extra information.
- `graph`: the `SimpleWeightedGraph` through which to find the shortest path.
- `source`: the source node.
- `destination`: the destination node.

# Returns
- `path`: the shortest path from `source` to `destination`, as a `Path` object.
- `cost`: the cost of the path.
- `cache`: a `Dict{Symbol, Any}` that may contain extra information provided by
    `cost_function` for the shortest path.

# Note
The regular Dijkstra algorithm can be obtained by setting the cost function to
`path_length`. It will be less efficient than a dedicated implementation of the regular
Dijkstra algorithm though, as the total length of the path will be recalculated from
scratch each time, instead of by just adding an edge length to the distance to the
previous node.
"""
function extended_dijkstra(cost_function, graph, source, destination; kwargs...)
    costs = Dict{eltype(graph), Float64}(i => Inf for i in vertices(graph))
    caches = Dict{eltype(graph), Dict{Symbol, Any}}(i => Dict{Symbol, Any}()
        for i in vertices(graph))
    previous_nodes = Dict{eltype(graph), Union{Nothing, eltype(graph)}}(
        i => nothing for i in vertices(graph))
    visited = Dict{eltype(graph), Bool}(i => false for i in vertices(graph))
    costs[source] = 0.
    q = PriorityQueue{eltype(graph), Float64}()
    q[source] = 0.

    while !isempty(q)
        current_node = dequeue!(q)
        if current_node == destination
            path = _path_from_previous_nodes(previous_nodes, destination, graph)
            cost = costs[current_node]
            cache = caches[current_node]
            return path, cost, cache
        end
        visited[current_node] && continue
        visited[current_node] = true
        path_to_current = _path_from_previous_nodes(previous_nodes, current_node, graph)
        for v in neighbors(graph, current_node)
            visited[v] && continue
            path = path_to_current + Path([current_node, v], graph)
            cost, cache = cost_function(path; kwargs...)
            if cost < costs[v]
                costs[v] = cost
                caches[v] = cache
                previous_nodes[v] = current_node
                q[v] = cost
            end
        end
    end
end

function _path_from_previous_nodes(previous_nodes::Dict{U, V}, destination, graph
        ) where {U, V}
    node = destination
    nodes = U[]
    while !isnothing(node)
        pushfirst!(nodes, node)
        node = previous_nodes[node]
    end
    Path(nodes, graph)
end

"""
    PathEnumerationMethod

Abstract type used to define different methods for enumerating paths in a graph,
with its primary use being pathfinding algorithms.
"""
abstract type PathEnumerationMethod end

"""
    DepthFirstSearch

Enumerate longer paths first. If multiple paths are the longest, pick the best path.

This can be used to perform pathfinding in a graph such that some path between the source
and destination is found as soon as possible, by constantly expanding upon long paths
such that a large part of the graph is traversed early on in the search.
The path may be long and have a bad cost.
"""
struct DepthFirstSearch <: PathEnumerationMethod end

"""
    BreadthFirstSearch

Enumerate shorter paths first. If multiple paths are the shortest, pick the best path.

This can be used to perform pathfinding in a graph such that a broad range of paths is
explored at the same time.
"""
struct BreadthFirstSearch <: PathEnumerationMethod end

"""
    WaveFrontSearch

Enumerate paths based on their cost.

This can be used to perform a Dijkstra-like search through the space of all paths between
a source and a destination, where the next path that is expanded is always the one with the
best cost.
"""
struct WavefrontSearch <: PathEnumerationMethod end

"""
    path_finding_for_monotonic_costs(cost_function, graph, source, destination;
        method::PathEnumerationMethod=DepthFirstSearch(),
        upper_bound=Inf, node_sequence_guesses=Vector{Vector{eltype(graph)}}[],
        prune_edges=false, kwargs...)

Find the best path through a graph assuming only monotonicity, not isotoniciy.

Typical path finding algorithms, like (extended) Dijsktra, assume a cost function that is
both monotonic (cost(p + e) > cost(p)) and isotinic (cost(p) > cost(p') implies that
cost(p + e) > cost(p' + e)). This algorithm drops the isotonicity assumption, and as a
consequence it is unable to take some of the shortcuts that are included in e.g. Dijkstra.
Basically the only thing you can do is enumerate all the paths and choose the best one.
This algorithm does exactly that, but it enumerates in a smarter way based on monotonicity.
Each time a subpath is enumerated for which the cost is higher than the current best path,
it is discarded, as all paths based on that subpath are be worse than the current best path.

# Arguments
- `cost_function`: Function that takes a `Path` object and returns two objects:
    1. a `Float64` that represents the cost of the path;
    2. a `Dict{Symbol, Any}` that may contain extra information.
- `graph`: the `SimpleWeightedGraph` through which to find the shortest path.
- `source`: the source node.
- `destination`: the destination node.

# Returns
- `path`: the shortest path from `source` to `destination`, as a `Path` object.
- `cost`: the cost of the path.
- `cache`: a `Dict{Symbol, Any}` that may contain extra information provided by
    `cost_function` for the shortest path.

The order in which paths are enumerated is determined by the `method` keyword argument.
The default is `DepthFirstSearch()`, which extends the longest paths first, hoping to reach
the destination node as soon as possible and allow for early pruning of bad paths.

If an upper bound on the cost of the best path is known a priori, it can be passed as the
`upper_bound` keyword argument. This will allow the algorithm to discard paths early on.
However, if the upper bound is set lower than the actual cost of the best path,
the algorithm will fail to find a path.
In that case, it returns `nothing, nothing`.

An upper bound can also be provided implicitly by providing guesses for good paths.
These should be passed as a vector of `Vector{eltype(graph)}`, where each inner vector is a
different guess represented as a sequence of nodes to be traversed by the path.
For each of the guesses, the cost function is automatically calculated.
Then, the most stringent upper bound provided by the `upper_bound` keyword and these guesses
is used as the real upper bound in path finding.

If the keyword argument `prune_edges` is set to `true`, the graph is pruned before path
finding commences. This means that all edges with a cost higher than the upper bound are
removed from the graph. Moreover, each time the upper bound is updated, the graph is pruned
again. This can result in a significant speedup, but the method assumes that the cost Of
a path consisting of a single edge grows monotonically with the weight of that edge.
This may not be true for every cost function, and hence the default value of `prune_edges`
is `false`.
"""
function path_finding_for_monotonic_costs(cost_function, graph, source, destination;
        method::PathEnumerationMethod=DepthFirstSearch(),
        upper_bound=Inf, node_sequence_guesses=Vector{Vector{eltype(graph)}}[],
        prune_edges=false, kwargs...)
    valid_guesses = filter(node_sequence_guesses) do n
        (first(n) == source && last(n) == destination) ||
            (first(n) == destination && last(n) == source)
    end
    path_guesses = Path[]
    for guess in valid_guesses
        try push!(path_guesses, Path(guess, graph))
        catch e
            e isa ArgumentError && continue
            rethrow(e)
        end
    end
    path_guess_results = cost_function.(path_guesses; kwargs...)
    path_guess_costs = first.(path_guess_results)
    upper_bounds = vcat(path_guess_costs, upper_bound)
    upper_bound = minimum(upper_bounds)
    prune_edges && (graph = prune_graph(graph, cost_function, upper_bound; kwargs...))
    paths = PathQueue(method)
    initial_path = Path([source], graph)
    initial_cost = -Inf
    enqueue!(paths, initial_path, initial_cost)
    final_paths = PriorityQueue{Path, Float64}()
    caches = Dict{Path, Dict{Symbol, Any}}()
    while !isempty(paths)
        p, c = dequeue_pair!(paths)
        c > upper_bound && continue
        for path in grow_path(p, graph)
            cost, cache = cost_function(path; kwargs...)
            cost > upper_bound && continue
            if path.nodes[end] == destination
                upper_bound = cost
                prune_edges && (
                    graph = prune_graph(graph, cost_function, upper_bound; kwargs...))
                enqueue!(final_paths, path, cost)
                caches[path] = cache
                continue
            end
            enqueue!(paths, path, cost)
        end
    end
    if isempty(final_paths)
        isempty(path_guesses) && return nothing, nothing, nothing
        index_best_guess = argmin(path_guess_costs)
        best_guess_path = path_guesses[index_best_guess]
        best_guess_cost = path_guess_costs[index_best_guess]
        best_guess_cache = path_guess_results[index_best_guess][2]
        return best_guess_path, best_guess_cost, best_guess_cache
    end
    # isempty(final_paths) && throw(ErrorException)
    path, cost = dequeue_pair!(final_paths)
    cache = caches[path]
    path, cost, cache
end

"""
    PathQueue{M <: PathEnumerationMethod, T}
    PathQueue{M} where M <: PathEnumerationMethod
    PathQueue(method::PathEnumerationMethod)

Priority queue that orders paths based on the method used to enumerate them.

`PathQueue{M, T}` is a parametric type, with `M` indicating the used enumeration method
and `T` the data type used to order paths correctly in the inner priority queue on which
`PathQueue` is based. The value of `T` is considered an implementation detail and should
not be set or accessed by the user. Methods defined on `PathQueue` should dispatch on `M`.

Paths can be added or removed from the queue using the `DataTypes.enqueue!`,
`DataTypes.dequeue!` and `DataTypes.dequeue_pair!` functions.
Methods for `PathQueue{M}` have been defined such that these functions automatically
lead to ordering of paths as determined by the value of the `M <: PathEnumerationMethod`.

The simplest way to construct an instance of this type is using the constructor
`PathQueue(method::PathEnumerationMethod)`, where `method` defines the desired enumeration
method.
"""
struct PathQueue{M <: PathEnumerationMethod, T}
    priority_queue::PriorityQueue{Path, T}
end
PathQueue{M}() where M <: Union{DepthFirstSearch, BreadthFirstSearch} =
    PathQueue{M, Tuple{Int, Float64}}(PriorityQueue{Path, Tuple{Int, Float64}}())
PathQueue{WavefrontSearch}() =
    PathQueue{WavefrontSearch, Float64}(PriorityQueue{Path, Float64}())
PathQueue(::M) where M <: PathEnumerationMethod = PathQueue{M}()

function DataStructures.enqueue!(q::PathQueue{DepthFirstSearch}, p::Path, cost::Float64)
    enqueue!(q.priority_queue, p, (-length(p), cost))
end
function DataStructures.enqueue!(q::PathQueue{BreadthFirstSearch}, p::Path, cost::Float64)
    enqueue!(q.priority_queue, p, (length(p), cost))
end
function DataStructures.enqueue!(q::PathQueue{WavefrontSearch}, p::Path, cost::Float64)
    enqueue!(q.priority_queue, p, cost)
end
function DataStructures.dequeue_pair!(q::PathQueue{M}) where
        M <: Union{DepthFirstSearch, BreadthFirstSearch}
    path, (_, cost) = dequeue_pair!(q.priority_queue)
    path => cost
end
DataStructures.dequeue_pair!(q::PathQueue{WavefrontSearch}) =
    dequeue_pair!(q.priority_queue)
DataStructures.dequeue!(q::PathQueue) = dequeue!(q.priority_queue)
Base.isempty(q::PathQueue) = isempty(q.priority_queue)

function grow_path(path, graph)
    paths = typeof(path)[]
    for neighbor in neighbors(graph, path.nodes[end])
        !(neighbor in path.nodes) &&
            push!(paths, path + Path([path.nodes[end], neighbor], graph))
    end
    paths
end

"""
    prune_graph(graph, cost_function, upper_bound; kwargs...)

Prune a graph by removing edges for which the cost is worse than some upper bound.
"""
function prune_graph(graph, cost_function, upper_bound; kwargs...)
    g = copy(graph)
    edges_to_search = [e for e in edges(g)]
    while !isempty(edges_to_search)
        med = median(weight.(edges_to_search))
        med_cost, _ = cost_function(Path([1, 2], [med]); kwargs...)
        if med_cost > upper_bound
            for _ in 1:length(edges_to_search)
                edge = popfirst!(edges_to_search)
                weight(edge) ≥ med && (rem_edge!(g, edge); continue)
                push!(edges_to_search, edge)
            end
        else
            for _ in 1:length(edges_to_search)
                edge = popfirst!(edges_to_search)
                weight(edge) > med && push!(edges_to_search, edge)
            end
        end
    end
    g
end

"""
    path_from_linear_interpolation(coords::Coordinates, start_node, end_node, num_points=10)

Create a guess for a best path between two nodes by linearly interpolating between them.

This is a technique that only works for geometric graphs, where the nodes have coordinates
associated with them instead of just edges. The guess is created by drawing a straight line
between the start and end node, and at regular intervals along that line, choosing the
nearest node to that point as the next node in the path.

# Arguments
- `coords`: a `Coordinates` object that holds the node locations.
- `start_node`: ther index of the start node of the path.
- `end_node`: the index of the end node of the path.
- `num_points`: the number of points considered on the line between the start and end node.

# Returns
- `path`: a `Path` object that represents the guess.
"""
function path_from_linear_interpolation(coords::Coordinates, start_node, end_node,
        num_points=10)
    node_list = [start_node]
    start_coord = node(coords, start_node)
    end_coord = node(coords, end_node)
    for i in 1:num_points
        next_coord_guess = start_coord + i * (end_coord - start_coord) / (num_points + 1)
        next_node_index, _ = nearest_neighbour(next_coord_guess, coords)
        next_node_index != node_list[end] && push!(node_list, next_node_index)
    end
    node_list[end] != end_node && push!(node_list, end_node)
    Path(node_list, coords)
end