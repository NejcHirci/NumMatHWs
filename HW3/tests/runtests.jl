using Test

include("../src/main.jl");

tol = 1e-10

@testset verbose=true "Testiranje Runge-Kutta metode 4. reda" begin

    @testset verbose=true "globalna napaka na primeru dy=-y+1, y(0) = 2" begin
        rk_f = (x,y) -> -y .+ 1.0
        t = 1.0
        y0 = [2.0]
        n = 100 # korak velikosti 0.01
        h = t/n
        xs = range(0, stop=t, length=n+1)
        y_rk4 = rk4(rk_f, y0, t, n)
        y_true = 1.0 .+ exp.(-xs)
        error = maximum(abs.(y_rk4 .- y_true))
        @test error < h^4
    end

    @testset verbose=true "ocena lokalne napake na primeru nihala" begin
        l = 1.0
        theta0 = 0.1
        dtheta0 = 0.1
        t = 1.0
        n = 100
        h = t/n
        local_err = maximum(pendulum_error(l, t, theta0, dtheta0, 100))
        @test local_err < 100.0 * h^5
    end    
end 