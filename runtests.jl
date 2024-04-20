using Test

using LinearAlgebra

include("scripts.jl")

@testset "Testi za funkcijo interpoliraj" begin
    # Primer poračunan na roke
    x = [0.9 1.3 1.9 2.1]
    y = [1.3 1.5 1.85 2.1]
    knots = [0.0 -5/19.0 195/76 0.0] 
    koef = [
        1.3 59/114 0.0 -25/228;
        1.5 173/456 -5/38.0 1075/1368;
        1.85 41/38 195/152 -325/152;
    ]
    Z = interpoliraj(x, y)
    @test isapprox(Z.koef, koef) # Test 1
    # Primer poračunan z uporabo CubicSplines
end
