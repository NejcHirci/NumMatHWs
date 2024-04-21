using Test

using LinearAlgebra

include("spline.jl")

@testset verbose=true "Testiranje naravnega kubicnega zlepka" begin
@testset "Testi funkcije interpoliraj" begin
    # Primer poračunan na roke
    x = [0.9 1.3 1.9 2.1]
    y = [1.3 1.5 1.85 2.1]
    koef = [
        -25/228 0.0 859/228 13/4;
        325/456 -25/342 1289/456 48/19;
        0.0 325/152 21/2 1393/152; 
    ]
    Z = interpoliraj(x, y)
    @test isapprox(Z.koef, koef) # Test 1
end
@testset "Testi funkcije vrednost" begin
    # Primerjamo še razliko vrednosti v dani točki, če interpoliramo točke znane funkcije
    x = [a for a in range(0, stop=2pi, length=400)]
    y = sin.(x)
    Z = interpoliraj(x, y)
    @test abs(vrednost(Z, 3.14) - sin(3.14)) < 1e-6 # Test 2
    @test abs(vrednost(Z, 1.57) - sin(1.57)) < 1e-6 # Test 3
    @test abs(vrednost(Z, 4.5) - sin(4.5)) < 1e-6 # Test 4
end
end
