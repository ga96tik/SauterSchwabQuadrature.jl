@testset "reorder!" begin
    quadPW = SauterSchwabQuadrature._legendre(4, 0, 1)
    p(x, y, z) = SVector(x, y, z)

    @testset "Triangle" begin
        t = SVector(
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
        )

        cases = (
            (
                CommonVertex(quadPW),
                SVector(
                    p(0.0, 0.0, 0.0),
                    p(2.0, 0.0, 0.0),
                    p(0.0, 2.0, 0.0),
                ),
            ),
            (
                CommonEdge(quadPW),
                SVector(
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                    p(1.0, 1.0, 0.0),
                ),
            ),
            (
                CommonFace(quadPW),
                SVector(
                    p(0.0, 1.0, 0.0),
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                ),
            ),
        )

        for (strat, s) in cases
            I, J, K, L = reorder(t, s, strat)

            Ip = Vector{Int}(undef, 3)
            Jp = Vector{Int}(undef, 3)
            Kp = Vector{Int}(undef, 3)
            Lp = Vector{Int}(undef, 3)
            Ip2, Jp2, Kp2, Lp2 = reorder!(Ip, Jp, Kp, Lp, t, s, strat)

            @test Ip2 === Ip
            @test Jp2 === Jp
            @test Kp2 === Kp
            @test Lp2 === Lp
            @test Ip == I
            @test Jp == J
            @test Kp == K
            @test Lp == L

            reorder!(Ip, Jp, Kp, Lp, t, s, strat)
            @test (@allocated reorder!(Ip, Jp, Kp, Lp, t, s, strat)) == 0
        end
    end

    @testset "Quadrilateral" begin
        t = SVector(
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
        )

        cases = (
            (
                CommonVertexQuad(quadPW),
                SVector(
                    p(0.0, -1.0, 0.0),
                    p(-1.0, -1.0, 0.0),
                    p(-1.0, 0.0, 0.0),
                    p(0.0, 0.0, 0.0),
                ),
            ),
            (
                CommonEdgeQuad(quadPW),
                SVector(
                    p(0.0, -1.0, 0.0),
                    p(1.0, -1.0, 0.0),
                    p(1.0, 0.0, 0.0),
                    p(0.0, 0.0, 0.0),
                ),
            ),
            (
                CommonFaceQuad(quadPW),
                SVector(
                    p(1.0, 1.0, 0.0),
                    p(0.0, 1.0, 0.0),
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                ),
            ),
        )

        for (strat, s) in cases
            I, J, K, L = reorder(t, s, strat)

            Ip = Vector{Int}(undef, 4)
            Jp = Vector{Int}(undef, 4)
            Ip2, Jp2, Kp2, Lp2 = reorder!(Ip, Jp, nothing, nothing, t, s, strat)

            @test Ip2 === Ip
            @test Jp2 === Jp
            @test Ip == I
            @test Jp == J
            @test K === nothing
            @test L === nothing
            @test Kp2 === nothing
            @test Lp2 === nothing

            reorder!(Ip, Jp, nothing, nothing, t, s, strat)
            @test (@allocated reorder!(Ip, Jp, nothing, nothing, t, s, strat)) == 0
        end
    end
end
