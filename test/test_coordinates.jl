@testset "Coordinates class" begin
    @testset "properties" begin
        coords = Coordinates([[0, 0] [5, 5] [10, 0]], [[3, 3] [3, 7]])
        @test RepeaterPlacement.end_nodes(coords) == [[0, 0] [5, 5] [10, 0]]
        @test RepeaterPlacement.repeaters(coords) == [[3, 3] [3, 7]]
        @test RepeaterPlacement.nodes(coords) == [0 5 10 3 3; 0 5 0 3 7]
        @test num_end_nodes(coords) == 3; @test num_repeaters(coords) == 2
        @test num_nodes(coords) == 5
        @test end_node(coords, 2) == [5, 5]
        @test node(coords, 2) == [5, 5]
        @test repeater(coords, 1) == [3, 3]
        @test node(coords, 4) == [3, 3]
        coords = Coordinates()
        @test num_repeaters(coords) == num_end_nodes(coords) == 0
        coords = Coordinates([[1, 1];;])
        @test num_end_nodes(coords) == 1; @test num_repeaters(coords) == 0
        coords = Coordinates([[1, 1]])  # pass as vector instead
        @test num_end_nodes(coords) == 1; @test num_repeaters(coords) == 0
        coords = Coordinates([1, 1], [2, 2])
        @test num_end_nodes(coords) == 1; @test num_repeaters(coords) == 1
    end
    @testset "add nodes" begin
        coords = Coordinates([[0, 1];;], [[10, 11];;])
        coords = add_end_node(coords, [1, 2])
        coords = add_repeater(coords, [50, 100])
        @test RepeaterPlacement.end_nodes(coords) == [[0, 1] [1, 2]]
        @test RepeaterPlacement.repeaters(coords) == [[10, 11] [50, 100]]
        coords = Coordinates()
        coords = add_end_node(coords, [1, 2])
        coords = add_repeater(coords, [50, 100])
        @test RepeaterPlacement.end_nodes(coords) == [[1, 2];;]
        @test RepeaterPlacement.repeaters(coords) == [[50, 100];;]
    end
    @testset "distance functions" begin
        # distance
        @test RepeaterPlacement.distance([0, 0], [1, 0]) == 1
        @test RepeaterPlacement.distance([0, 0], [1, 1]) == sqrt(2)

        # nearest neighbour
        nn, dist = RepeaterPlacement.nearest_neighbour([0, 0], [[1, 1] [2, 2] [3, 3]])
        @test nn == 1
        @test dist == sqrt(2)
        coords = Coordinates([[1.0, 1.0] [2.0, 2.0]], [[3.0, 3.0] [4.0, 4.0]])
        nn, dist = RepeaterPlacement.nearest_neighbour([0, 0], coords)
        @test nn == 1
        @test node(coords, nn) == [1., 1.]
        @test dist == sqrt(2)
        nn, dist = RepeaterPlacement.nearest_neighbour([3.3, 3.3], coords)
        @test nn == 3
        @test node(coords, nn) == [3., 3.]
        @test dist ≈ sqrt(.3^2 + .3^2)

        # max distance
        @test RepeaterPlacement.max_dist([[0, 0] [10, 0] [15, 0]]) == 10
        coords = Coordinates([[0, 0] [10, 0]], [[15, 0] [50, 0]])
        @test RepeaterPlacement.max_dist(coords) == 35
        coords = Coordinates([[0, 0] [50, 0]], [[55, 0] [60, 0]])
        @test RepeaterPlacement.max_dist(coords) == 50
        coords = Coordinates([[55, 0] [0, 0]], [[40, 0] [60, 0]])
        @test RepeaterPlacement.max_dist(coords) == 40
    end
    @testset "build graph" begin
        # adjacency matrix
        mat = RepeaterPlacement.adjacency_matrix([[0, 0] [0, 3] [10, 0]], 5)
        @test mat[1, 2] == mat[2, 1] == 3
        @test mat[1, 1] == mat[2, 2] == mat[3, 3] == 0
        @test mat[1, 3] == mat[3, 1] == mat[2, 3] == mat[3, 2] == 0
        coords = Coordinates([[0, 0] [0, 3]], [[10, 0];;])
        mat = RepeaterPlacement.adjacency_matrix(coords, 5)
        @test mat[1, 2] == mat[2, 1] == 3
        @test mat[1, 1] == mat[2, 2] == mat[3, 3] == 0
        @test mat[1, 3] == mat[3, 1] == mat[2, 3] == mat[3, 2] == 0
        mat = RepeaterPlacement.adjacency_matrix(coords)
        @test mat[1, 2] == mat[2, 1] == 3
        @test mat[1, 1] == mat[2, 2] == mat[3, 3] == 0
        @test mat[1, 3] == mat[3, 1] == 10
        @test mat[2, 3] == mat[3, 2] == sqrt(3^2 + 10^2)

        # graph
        g = build_graph(coords)
        @test nv(g) == 3
        @test ne(g) == 3
        @test get_weight(g, 1, 2) == 3
        @test get_weight(g, 1, 3) == 10
        @test get_weight(g, 2, 3) == sqrt(3^2 + 10^2)
        g = build_graph(coords, 5)
        @test nv(g) == 3
        @test ne(g) == 1
        @test get_weight(g, 1, 2) == 3
        @test !has_edge(g, 1, 3)
        @test !has_edge(g, 2, 3)

        # only end nodes
        coords = Coordinates([[0, 0] [0, 3] [10, 0]])
        g = build_graph(coords, 5)
        @test nv(g) == 3
        @test ne(g) == 1
        @test get_weight(g, 1, 2) == 3
        @test !has_edge(g, 1, 3)
        @test !has_edge(g, 2, 3)
    end

    @testset "build_waxmann_graph" begin
        coords = Coordinates([[0, 0] [0, 3] [10, 0]])
        g = build_waxmann_graph(coords)
        @test nv(g) == 3
        @test 0 <= ne(g) <= 3

        # probability 1 for each edge
        g = build_waxmann_graph(coords, 1, Inf)
        @test nv(g) == 3
        @test ne(g) == 3
        @test get_weight(g, 1, 2) == 3

        # probability 0 for each edge
        g = build_waxmann_graph(coords, 0)
        @test nv(g) == 3
        @test ne(g) == 0
    end

    @testset "random initialization" begin
        # line
        line = RepeaterPlacement.initialize_line(3, 10)
        @test num_end_nodes(line) == 2
        @test num_repeaters(line) == 3
        @test num_nodes(line) == 5
        @test Set(end_node(line, i) for i in 1:2) == Set([[0, 0], [10, 0]])
        for i in 1:3
            @test repeater(line, i)[2] == 0
            @test 0 < repeater(line, i)[1] < 10
        end

        # square
        square = RepeaterPlacement.initialize_square(3, 10)
        @test num_end_nodes(square) == 4
        @test num_repeaters(square) == 3
        @test num_nodes(square) == 7
        @test Set(end_node(square, i) for i in 1:4) ==
            Set([[0, 0], [10, 10], [0, 10], [10, 0]])
        for i in 1:3
            @test 0 < repeater(square, i)[1] < 10
            @test 0 < repeater(square, i)[2] < 10
        end

        # random end nodes
        coords = RepeaterPlacement.initialize_random(3, 2, 10)
        @test num_end_nodes(coords) == 3
        @test num_repeaters(coords) == 2
        @test num_nodes(coords) == 5
        for i in 1:3, j in 1:2
            @test 0 < end_node(coords, i)[j] < 10
        end
        for i in 1:2, j in 1:2
            @test 0 < repeater(coords, i)[j] < 10
        end

        # consistent seeding
        coords_1 = RepeaterPlacement.initialize_random(3, 2, 10, Xoshiro(1234))
        coords_2 = RepeaterPlacement.initialize_random(3, 2, 10, Xoshiro(1234))
        @test RepeaterPlacement.repeaters(coords_1) == RepeaterPlacement.repeaters(coords_2)
        @test RepeaterPlacement.end_nodes(coords_1) == RepeaterPlacement.end_nodes(coords_2)
    end

    @testset "waxmann_graph" begin
        g, coords = waxmann_graph(3, 2, 10.)
        @test num_end_nodes(coords) == 3
        @test num_repeaters(coords) == 2
        @test num_nodes(coords) == 5
        for i in 1:3, j in 1:2
            @test 0 < end_node(coords, i)[j] < 10
        end
        for i in 1:2, j in 1:2
            @test 0 < repeater(coords, i)[j] < 10
        end
        @test nv(g) == 5
        @test 0 <= ne(g) <= 10

        # deterministic edges
        g, coords = waxmann_graph(3, 2, 1., Inf)
        @test nv(g) == 5
        @test ne(g) == 5 * 4 / 2

        # no edges
        g, coords = waxmann_graph(3, 2, 0.)
        @test nv(g) == 5
        @test ne(g) == 0

    end
end