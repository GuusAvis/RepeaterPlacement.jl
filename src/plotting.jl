"""
    plot_graph(sol::CoordinateSolution, kwargs...)
    plot_graph(coords::Coordinates, radius=Inf; path=nothing, kwargs...)
    plot_graph(g::SimpleWeightedGraph, locs_x=nothing, locs_y=nothing; path=nothing,
        kwargs...)
    plot_graph(g::SimpleWeightedGraph, coords::Coordinates; path=nothing,
        kwargs...)

Make and plot graph with nodes at correct locations and edges labeled with their lengths.

If a `CoordinateSolution` is provided, a `Coordinates` and `Path` are extracted from it.
If a `Coordinates` is provided, this function uses `build_graph` to create the graph.
If a graph is provided, the X and Y coordinates of the nodes must be provided separately
as `locs_x` and `locs_y` in the form of vectors, or using a `Coordinates` object.
If a `Path` is provided, the edges in the path will be colored differently from the others.
Plotting is handled by `GraphPlot.gplot()`.
"""
plot_graph(sol::CoordinateSolution) = 
    plot_graph(Coordinates(sol), sol.X.kwargs[:radius];
        path=sol.X.cache[:path])

function plot_graph(coords::Coordinates, radius=Inf; path=nothing, kwargs...)
    g = build_graph(coords, radius)
    locs_x = [node(coords, i)[1] for i in 1:num_nodes(coords)]
    locs_y = [-node(coords, i)[2] for i in 1:num_nodes(coords)]
    plot_graph(g, locs_x, locs_y; path=path, kwargs...)
end

function plot_graph(g::SimpleWeightedGraph, locs_x=nothing, locs_y=nothing;
        path=nothing, kwargs...)
    edge_weights = [get_weight(g, e.src, e.dst) for e in edges(g)]
    edge_weights = round.(edge_weights, digits=1)
    if path === nothing
        locs_x === nothing && locs_y === nothing &&
            return gplot(g, edgelabel=edge_weights, nodelabel=1:nv(g))
        locs_x !== nothing && locs_y !== nothing &&
            return gplot(g, locs_x, locs_y, edgelabel=edge_weights, nodelabel=1:nv(g))
        throw(ArgumentError("locs_x and locs_y must be both specified or both omitted"))
    end
    edge_colors = []
    for e in edges(g)
        index_src = findfirst(x->x==e.src, path.nodes)
        index_dst = findfirst(x->x==e.dst, path.nodes)
        if index_src !== nothing && index_dst !== nothing &&
                abs(index_src - index_dst) == 1
            push!(edge_colors, colorant"orange")
        else
            push!(edge_colors, colorant"lightgray")
        end
    end
    locs_x === nothing && locs_y === nothing &&
        return gplot(g, edgelabel=edge_weights, nodelabel=1:nv(g), edgestrokec=edge_colors)
    gplot(g, locs_x, locs_y; edgelabel=edge_weights, nodelabel=1:nv(g),
        edgestrokec=edge_colors)#, kwargs...)
end

function plot_graph(g::SimpleWeightedGraph, coords::Coordinates, path=nothing, kwargs...)
    locs_x = [node(coords, i)[1] for i in 1:num_nodes(coords)]
    locs_y = [-node(coords, i)[2] for i in 1:num_nodes(coords)]
    plot_graph(g, locs_x, locs_y; path=path, kwargs...)
end

extract_edges(e::Pair{Int, Int}) = [e]
extract_edges(p::Path) = edges(p)
extract_edges(g::R) where R <: AbstractGraph = [(e.src, e.dst) for e in Graphs.edges(g)]

"""
    plot_node_locations(coords::Coordinates, edges=[], special_edges=[];
        draw_lengths=true, size=(600, 600), legend=false, padding=20, markersize=17,
        tickfontsize=12, guidefontsize=15, annotation_offset=15, annotation_fontsize=15,
        kwargs...)

Create a scatter plot showing repeater locations and end-node locations.

If `edges` is provided as a list of pairs of node indices, these edges are drawn in the
figure. If `special_edges` is provided, these edges are drawn with a thicker line.
Both lists are also allowed to contain `Path` or `SimpleGraph` objects, in which case the
edges of the paths or graphs are extracted and drawn in the figure.

If `draw_lengths` is set to `true`, the length of each connection is annotated in the figure
(rounded to a whole number), and the distance of the annotation from its corresponding edge
and its font size can be controlled using the keyword arguments `annotation_offset` and
`annotation_fontsize`.
"""
function plot_node_locations(coords::Coordinates, edges=[], special_edges=[];
        draw_lengths=true, size=(600, 600), legend=false, padding=20, markersize=17,
        tickfontsize=12, guidefontsize=15, annotation_offset=15, annotation_fontsize=15,
        kwargs...)

    # draw the node locations
    x_min = minimum([nodes(coords)[1, i] for i in 1:num_nodes(coords)])
    x_max = maximum([nodes(coords)[1, i] for i in 1:num_nodes(coords)])
    y_min = minimum([nodes(coords)[2, i] for i in 1:num_nodes(coords)])
    y_max = maximum([nodes(coords)[2, i] for i in 1:num_nodes(coords)])
    p = plot(;tickfontsize=tickfontsize, guidefontsize=guidefontsize, size=size, kwargs...)
    if ! isempty(end_nodes(coords))
        scatter!(p, end_nodes(coords)[1, :], end_nodes(coords)[2, :], label="end nodes";
            legend=legend, markershape=:square, markersize=markersize,
            aspect_ratio=:equal, kwargs...)
    end
    if ! isempty(repeaters(coords))
    scatter!(p, repeaters(coords)[1, :], repeaters(coords)[2, :], label="repeater nodes";
        legend=legend, markershape=:circle, markersize=markersize, kwargs...)
    end
    xlabel!("x coordinate [km]")
    ylabel!("y coordinate [km]")
    xlims!(x_min - padding, x_max + padding)
    ylims!(y_min - padding, y_max + padding)

    edges = reduce(vcat, extract_edges.(edges); init=[])
    edges = Set(Set.(edges))
    special_edges = reduce(vcat, extract_edges.(special_edges); init=[])
    special_edges = Set(Set.(special_edges))

    # draw all edges
    # for e in Iterators.flatten([edges, special_edges])
    for e in union(edges, special_edges)
        node1, node2 = [node(coords, i) for i in e]
        lw = e in special_edges ? 2 : 1
        col = e in special_edges ? :black : :grey
        plot!(p, [node1[1], node2[1]], [node1[2], node2[2]], color=col, linewidth=lw,
            z_order=:back, label=nothing)
        if draw_lengths
            dist = distance(node1, node2)
            dist = text(Int(round(dist, digits=0)), annotation_fontsize, color=:black)
            mid_x = mean([node1[1], node2[1]])
            mid_y = mean([node1[2], node2[2]])
            line_direction = [node2[1] - node1[1], node2[2] - node1[2]]
            normal = [-line_direction[2], line_direction[1]]
            normal = normal / (normal[1]^2 + normal[2]^2)^0.5  # normalize normal vector
            pos = [mid_x, mid_y] + annotation_offset * normal
            annotate!(p, pos[1], pos[2], dist)
        end
    end
    p
end

"""
    plot_node_locations(sol::CoordinateSolution{F, S}; kwargs...) where
        {F<:PathwiseCostFct, S}

Create a scatter plot showing repeater and end-node locations, as well as the best paths.

The best paths to draw in the figure are obtained from the cache of the solution.
"""
plot_node_locations(sol::CoordinateSolution{F, S}; kwargs...) where
        {F<:PathwiseCostFct, S} = 
    plot_node_locations(Coordinates(sol), sol.X.cache[:all_paths]; kwargs...)

"""
    plot_node_locations(sol::CoordinateSolution{LargestPathwiseCostFct, S}; kwargs...)
        where S

Create scatter plot showing repeaters, end nodes and best paths, with worst best path thick.

The best paths and worst best path are obtained from the cache of the solution.
"""
plot_node_locations(sol::CoordinateSolution{F, S}; kwargs...) where
        {F<:LargestPathwiseCostFct, S} =
    plot_node_locations(Coordinates(sol), sol.X.cache[:all_paths], [sol.X.cache[:path]];
        kwargs...)

"""
    plot_node_locations(trace::CoordinateTrace, stepsize=1; show_temperature=true)

Create a gif showing the evolution of the node placement over the iterations of the
optimization.

Choose a larger stepsize to skip some iterations in the gif, making it shorter.
"""
function plot_node_locations(trace::CoordinateTrace, stepsize=1; show_temperature=true,
        show_paths=false, kwargs...)
    initial_sol = trace[begin]
    coords = Coordinates(initial_sol)
    dist = 0
    for n1 in end_nodes(coords), n2 in end_nodes(coords)
        dist = max(dist, distance(n1, n2))
    end
    padding = 0.1 * dist
    @animate for i in 1:stepsize:length(trace)
        if show_paths
            paths = trace[i].X.cache[:all_paths]
            special_paths = [trace[i].X.cache[:path]]
        else
            paths = []; special_paths = []
        end
        plot_node_locations(Coordinates(trace[i]), paths, special_paths; legend=false)
        if show_temperature && :temperature in keys(trace[i].X.kwargs)
            temp = round(trace[i].X.kwargs[:temperature], sigdigits=2)
            annotate!(2 * padding, dist - padding, "T = $temp", fontsize=15)
        end
        xlims!(-padding, dist + padding)
        ylims!(-padding, dist + padding)
    end
end

"""
    plot_cache_value(trace::CoordinateTrace, key; ylabel=string(key),
        show_temperature=true, kwargs...)

Plot the value stored for `key` in the cache of each solution in the trace.
"""
function plot_cache_value(trace::CoordinateTrace, key, error_key=nothing;
    ylabel=string(key), show_temperature=true, kwargs...)
    vals = [StochasticAD.value(sol.X.cache[key]) for sol in trace]
    if isnothing(error_key)
        p = plot(vals; label=ylabel, xlabel="iteration", ylabel=ylabel, kwargs...)
    else
        errors = [sol.X.cache[error_key] for sol in trace]
        p = plot(vals; ribbon=errors, fillalpha=0.4, label=ylabel, xlabel="iteration",
            ylabel=ylabel, kwargs...)
    end
    if show_temperature && :temperature in keys(trace[begin].X.kwargs)
        temps = [sol.X.kwargs[:temperature] for sol in trace]
        pt = twinx()
        plot!(pt, temps; color=:red, linestyle=:dash,
        label="temperature", ylabel="temperature") # , yaxis=:log)
    end
    p
end

plot_cost_trace(trace::CoordinateTrace; show_temperature=true, kwargs...) =
    plot_cache_value(trace, :cost_function_value; ylabel="cost",
    show_temperature=show_temperature, kwargs...)

"""
    create_figures(sol::CoordinateSolution; skip_evaluation=false)
    create_figures(sol::CoordinateSolution, evaluated_sol::CoordinateSolution)
    create_figures(trace::CoordinateTrace; skip_evaluation=false)
    create_figures(trace::CoordinateTrace, evaluated_trace::CoordinateTrace;
        extra_plot_functions=Function[])

Use data and potentially evaluated data to make a number of informative figures.

If `skip_evaluation=false`, the provided data is processed using `evaluate` before
making figures. Most figures can be made if both unevaluated and evaluated data are provided
together.

The most figures can be made when both a trace and evaluated trace are provided, and hence
that is the recommended way to use this function. In this case, the keyword argument
`extra_plot_functions` can be used to pass a list of additional functions that
take `trace` and `evaluated_trace` and return a figure.
"""
function create_figures(sol::CoordinateSolution; skip_evaluation=false, kwargs...)
        skip_evaluation || (sol = evaluate(sol; kwargs...))
    p1 = plot_node_locations(sol; kwargs...)
    # p2 = plot_graph(s)  # not sure how to save graph created with GraphPlot.gplot().
    # [p1, p2]
    Any[p1]
end

function create_figures(sol::CoordinateSolution, evaluated_sol::CoordinateSolution;
        kwargs...)
    create_figures(evaluated_sol; skip_evaluation=true, kwargs...)
end

function create_figures(trace::CoordinateTrace; kwargs...)
    evaluated_trace = evaluate_complete(trace; kwargs...)
    create_figures(trace, evaluated_trace; kwargs...)
end

function create_figures(trace::CoordinateTrace, evaluated_trace::CoordinateTrace;
        extra_plot_functions=Function[], kwargs...)
    ps = create_figures(evaluated_trace[end]; skip_evaluation=true, kwargs...)
    push!(ps, plot_node_locations(evaluated_trace; kwargs...))
    push!(ps, plot_cost_trace(trace; kwargs...))
    for f in extra_plot_functions
        p = f(trace, evaluated_trace; kwargs...)
        push!(ps, p)
    end
    ps
end