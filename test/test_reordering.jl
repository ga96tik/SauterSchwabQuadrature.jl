
@testset "Reordering" begin
    
    @testset "Common Face" begin 
        
        p11 = point(1,0,0)
        p12 = point(3,0,1)
        p13 = point(3,2,1)

        Sourcechart = SVector(p11,p12,p13)
        Testchart   = SVector(p13,p12,p11)

        I, J, K, L = reorder(Testchart, Sourcechart, CommonFace(1))

        @test I == [1, 2, 3]
        @test J == [3, 2, 1]
        @test K == [1, 2, 3]
        @test L == [3, 2, 1]
    end

    @testset "Common Edge" begin 

        p11 = point(1,0,0)
        p12 = point(3,0,0)
        p13 = point(3,2,0)

        p21 = point(0,0,0) 
        p22 = point(1,0,0) 
        p23 = point(3,2,0)

        Sourcechart = SVector(p11,p12,p13)
        Testchart   = SVector(p21,p22,p23)

        I, J, K, L = reorder(Testchart, Sourcechart, CommonEdge(1))

        @test I == [3, 1, 2]
        @test J == [3, 2, 1]
        @test K == [2, 3, 1]
        @test L == [3, 2, 1]
    end

    @testset "Common Vertex" begin 

        p11 = point(1,0,0)
        p12 = point(3,0,0)
        p13 = point(3,2,0)

        p21 = point(0,0,0) 
        p22 = point(1,0,0) # common point
        p23 = point(0,1,2)

        Sourcechart = SVector(p11,p12,p13)
        Testchart   = SVector(p21,p22,p23)

        I, J, K, L = reorder(Testchart, Sourcechart, CommonVertex(1))

        @test I == [2, 1, 3]
        @test J == [1, 2, 3]
        @test K == [2, 1, 3]
        @test L == [1, 2, 3]
    end
end
