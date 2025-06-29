module RepeaterPlacement
using StochasticAD, Statistics, Optimisers, Plots, JLD2,
    ParameterSchedulers, StatsBase, Graphs, SimpleWeightedGraphs, GraphPlot, UUIDs,
    DataFrames, CSV, DataStructures, Dates, Random

export
    Coordinates, nodes, num_end_nodes, num_repeaters, num_nodes, end_node, repeater, node,
    add_end_node, add_repeater,
    build_graph, build_waxmann_graph, initialize_line, initialize_square, initialize_random,
    waxmann_graph,
    Path, path_length, enumerate_paths, path_finding_through_enumeration, extended_dijkstra,
    path_finding_for_monotonic_costs,
    PathEnumerationMethod, DepthFirstSearch, BreadthFirstSearch, WavefrontSearch,
    create_prune_fct_from_dominance, directly_dominates, completely_dominates,
    incompletely_dominates,
    CoordinateCostFunction,
    PathwiseCostFct, path_cost_fct, LargestPathwiseCostFct, MeanPathwiseCostFct,
    CoordinateSolution, CoordinateTrace, total_time, populate_cache!, populate_cache,
    evaluate_complete,
    plot_graph, plot_node_locations, plot_cache_value, plot_cost_trace

include("coordinates.jl")
include("path_finding.jl")
include("path_pruning.jl")
include("coordinate_cost_function.jl")
include("path_cost_function.jl")
include("optimization.jl")
include("plotting.jl")
include("data_handling.jl")

end