using Test, Distributions

include("../src/main.jl");

tol = 1e-10

@testset verbose=true "Testiranje porazdelitvene funkcije normalne porazdelitve" begin

    @testset verbose=true "s tabeliranimi vrednostmi" begin
        @test abs(gaussian_cdf(0.0, 20, 1e-12) - 0.5) / 0.5 < tol
        @test abs(gaussian_cdf(1.0, 20, 1e-12) - 0.8413447460685429) / 0.8413447460685429 < tol
        @test abs(gaussian_cdf(-100.0, 20, 1e-12) - 0.0) < tol
        @test abs(gaussian_cdf(10.0, 20, 1e-12) - 1.0) / 1.0 < tol
    end

    @testset verbose=true "z uporabo paketa Distributions.jl" begin
        normal_dist = Normal()
        for val in range(-5.0, 5.0, length=100)
            @test abs(gaussian_cdf(val, 20, 1e-12) - cdf(normal_dist, val)) < tol
        end
    end
end 