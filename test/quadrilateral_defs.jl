

# ------------------------------------- Geometry definition with parametrization
struct Quadrilateral
    p1::SVector{3}
    p2::SVector{3}
    p3::SVector{3}
    p4::SVector{3}
end


# parametrize planar quadrilateral with u, v ∈ [0,1]
function (quad::Quadrilateral)(u, v)    

    return quad.p1 + u * (quad.p2 - quad.p1) + v * (quad.p4 - quad.p1) + u * v * (quad.p3 - quad.p4 + quad.p1 - quad.p2) # see, e.g., [1] page 187
end


# Jacobi determinant of parametrization
function jacobianDet(quad::Quadrilateral, u)

    aux = quad.p3 - quad.p4 + quad.p1 - quad.p2

    ∂ru = quad.p2 - quad.p1 + u[2] * aux
    ∂rv = quad.p4 - quad.p1 + u[1] * aux

    D = (∂ru[2] * ∂rv[3] - ∂ru[3] * ∂rv[2])^2 + (∂ru[3] * ∂rv[1] - ∂ru[1] * ∂rv[3])^2 + (∂ru[1] * ∂rv[2] - ∂ru[2] * ∂rv[1])^2 

    return sqrt(D)
end




# ------------------------------------- Kernel definitions
struct singularKernel
    quad1::Quadrilateral
    quad2::Quadrilateral
end


function (sK::singularKernel)(u, v)

    x = sK.quad1(u...)
    y = sK.quad2(v...)

    return jacobianDet(sK.quad1, u) * jacobianDet(sK.quad2, v) / norm(x - y)
end




struct regularKernel
    quad1::Quadrilateral
    quad2::Quadrilateral
end

function (rK::regularKernel)(u, v)

    x = rK.quad1(u...)
    y = rK.quad2(v...)

    return jacobianDet(rK.quad1, u) * jacobianDet(rK.quad2, v) * norm(x - y)
end