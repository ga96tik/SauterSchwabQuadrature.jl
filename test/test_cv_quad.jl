

@testset "Common Vertex" begin 

    # quadrilateral defined by four points
    p1 = SVector(0.,0.,0.) # the same for both quads
    p2 = SVector(2.,0.,0.)
    p3 = SVector(2.,2.,0.)
    p4 = SVector(0.,2.,0.)

    q1 = Quadrilateral(p1,p2,p3,p4)

    p2 = SVector(-2., 0.,0.) # choose points of second quad such that parametrizations align as required in [1]
    p3 = SVector(-2.,-2.,0.) 
    p4 = SVector( 0.,-2.,0.)

    q2 = Quadrilateral(p1,p2,p3,p4)

    @testset "Regular Kernel" begin

        # --- Kernel
        K = regularKernel(q1,q2)

        # --- Sauter Schwab strategy
        quadPW = SauterSchwabQuadrature._legendre(10, 0, 1) # quadpoints
        intSauSw = sauterschwab_parameterized(K, CommonVertexQuad(quadPW)) # compute

        # --- Double quadrature
        qOut = SauterSchwabQuadrature._legendre(25, 0, 1)
        qOin = SauterSchwabQuadrature._legendre(25, 0, 1) # quadpoints

        intDQ = sum(w1*w2*w3*w4*K((η1, η2), (η3, ξ))
                    for (η1, w1) in qOut, (η2, w2) in qOut, (η3, w3) in qOin, (ξ, w4) in qOin) # compute
        
        # --- difference
        @test abs(intSauSw - intDQ) / intSauSw < 1e-13
    end

    @testset "Singular Kernel" begin
        
        # --- Kernel
        K = singularKernel(q1,q2)

        # --- Sauter Schwab strategy
        quadPW = SauterSchwabQuadrature._legendre(10, 0, 1) # quadpoints
        intSauSw = sauterschwab_parameterized(K, CommonVertexQuad(quadPW)) # compute

        # --- Double quadrature
        qOut = SauterSchwabQuadrature._legendre(15, 0, 1)
        qOin = SauterSchwabQuadrature._legendre(15, 0, 1) # quadpoints

        intDQ = sum(w1*w2*w3*w4*K((η1, η2), (η3, ξ))
                    for (η1, w1) in qOut, (η2, w2) in qOut, (η3, w3) in qOin, (ξ, w4) in qOin) # compute

        # --- difference
        @test abs(intSauSw - intDQ) / intSauSw < 1e-7
    end
end
