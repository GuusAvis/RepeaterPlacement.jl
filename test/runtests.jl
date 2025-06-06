using Test
using RepeaterPlacement, Graphs, SimpleWeightedGraphs, StochasticAD, ParameterSchedulers,
    DataStructures, Plots, Random

@testset "Coordinates" begin
    include("test_coordinates.jl")
end
@testset "Paths and path finding" begin
    include("test_path_finding.jl")
end
@testset "Path pruning" begin
    include("test_path_pruning.jl")
end
@testset "Coordinate cost function" begin
    include("test_coordinate_cost_function.jl")
end
@testset "Path cost function" begin
    include("test_path_cost_function.jl")
end
@testset "Optimization" begin
    include("test_optimization.jl")
end
@testset "Integration testing" begin
    include("test_integration.jl")
end