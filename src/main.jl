using Plots, LinearAlgebra


"""
    yn = rocket_rk4(y0, t, n, mass_moon=7.348e22, mass_earth=5.974e24)

Funkcija izračuna položaj sonde v omejenem sistemu treh masnih
teles. Implementacija uporablja Runge-Kutta metodo 4. reda.

# Argumenti
- `yn::Vector{Float64}`: začetni položaj in hitrost sonde
- `t::Float64`: čas simulacije
- `n::Int`: število korakov
- `mass_moon::Float64`: masa Lune
- `mass_earth::Float64`: masa Zemlje
"""
function rocket_rk4(y0, t, n, mass_moon=7.348e22, mass_earth=5.974e24)
    function rocket_dF(xn, yn, μ=mass_moon / (mass_earth + mass_moon))

        x, y, z, xdot, ydot, zdot = yn
        out = zeros(6)
        out[1:3] = [xdot, ydot, zdot]
    
        # Razdalja med sondo in Zemljo
        r_M = sqrt((x+μ)^2 + y^2 + z^2)
    
        # Razdalja med sondo in Luno
        r_m = sqrt((x-1.0+μ)^2 + y^2 + z^2)

        out[4] = 2 * ydot + x - (1.0-μ) * (x + μ) / r_M^3 - μ / r_m^3 * (x - 1.0 + μ)
        out[5] =-2 * xdot + y - (1.0-μ) * y / r_M^3 - μ * y / r_m^3
        out[6] = -(1.0-μ) / r_M^3 * z - μ / r_m^3 * z
        return out
    end
    return rk4(rocket_dF, y0, t, n) 
end


"""

    find_moon_flyby(init_low, init_high, rk4_steps=10000, secant_steps=100, secant_tol=1e-6)

Funkcija, ki z uporabo strelske metode izračuna optimalne začetne pogoje
za obhod Lune. Za optimizacijo uporabljam sekantno metodo.

# Argumenti
- `init_low::Float64`: spodnja meja začetnega pogoja
- `init_high::Float64`: zgornja meja začetnega pogoja
- `rk4_steps::Int`: število korakov Runge-Kutta metode
- `secant_steps::Int`: število korakov sekantne metode
- `secant_tol::Float64`: toleranca za konvergenco
"""
function find_moon_flyby(init_low, init_high, rk4_steps=10000, secant_steps=100, secant_tol=1e-6)
    tend = 2.0
    rk4_h = tend / rk4_steps
    println("h: ", rk4_h)

    # Objektivna funkcija za strelsko metodo
    function objective(init_cond, target)
        rk_yn = rocket_rk4(init_cond, tend, rk4_steps)
        mid_pos = rk_yn[div(end, 2), 1:3]
        target_mid_pos = [1.0, 0.0, 0.0]
        final_pos = rk_yn[end, 1:3]
        target_pos = target[1:3]
        error = norm(final_pos - target_pos) + norm(mid_pos - target_mid_pos)
        return error
    end

    x0 = init_low
    x1 = init_high

    for i=1:secant_steps
        g0 = objective(x0, x0)
        g1 = objective(x1, x1)
        if abs(g0 - g1) < secant_tol
            println("g0: ", g0, " g1: ", g1)
            println("Converged in ", i, " steps")
            break
        end
        x2 = x1 - g1 * (x0 - x1) ./ (g0 - g1)
        if norm(x2 - x1) < secant_tol
            println("Last steps: ", x2, " ", x1)
            println("Converged in ", i, " steps")
            break
        end
        x0 = x1
        x1 = x2
    end
    return x1
end


"""
    yn = rk4(f, y0, t, n)

Funkcija, ki izračuna rešitev diferencialne enačbe y' = f(y) z uporabo
Runge-Kutta metode 4. reda. Funkcija vrne vektor y_n, ki predstavlja
vrednost y ob času t.

# Argumenti
- `f::Function`: funkcija, ki opisuje desno stran sistema diferencialnih enačb
- `y0::Vector{Float64}`: začetni pogoji
- `t::Float64`: čas simulacije
- `n::Int`: število korakov

# Rezultat
- `yn::Vector{Float64}`: vrednost y ob času t
"""
function rk4(f, y0, t, n)
    h = t/n
    x = 0.0
    y = y0
    out = zeros(n+1, length(y0))
    out[1, :] = y
    for i=2:n+1
        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y .+ k1/2.0)
        k3 = h * f(x + h/2, y .+ k2/2.0)
        k4 = h * f(x + h, y .+ k3)
        y = y .+ (k1 .+ 2.0 * k2 .+ 2.0 * k3 .+ k4) ./ 6.0
        x = x + h
        out[i, :] = y
    end
    return out
end


"""
    odmik = nihalo(l,t,theta0,dtheta0,n)

Funkcija simulira nihanje matematičnega nihala z uporabo 
Runge-Kutta 4. reda. Funkcija vrne odmik nihala v odvisnosti od časa.

# Argumenti
- `l::Float64`: dolžina nihala
- `t::Float64`: čas simulacije
- `theta0::Float64`: začetni odmik v radianih
- `dtheta0::Float64`: začetna hitrost
- `n::Int`: število korakov
"""
function nihalo(l,t,theta0,dtheta0,n)
    # Definiramo desno stran sistema diferencialnih enačb
    function pendulum_dF(y, l, g=9.80665)
        θ, ω = y
        dθ = ω
        dω = -g/l * sin(θ)
        return [dθ, dω]
    end
    y0 = [theta0, dtheta0]
    return rk4((x,y) -> pendulum_dF(y, l), y0, t, n)
end


"""

    err = pendulum_error(l, t, theta0, dtheta0, n)

Funkcija vrne oceno za lokalno napako simulacije nihala.

# Argumenti
- `l::Float64`: dolžina nihala
- `t::Float64`: čas simulacije
- `theta0::Float64`: začetni odmik v radianih
- `dtheta0::Float64`: začetna hitrost
- `n::Int`: število korakov
"""
function pendulum_error(l, t, theta0, dtheta0, n)
    # Ocenimo napako korakov kot razliko med rezultati 
    # za koraka h in h/2.

    out_h = nihalo(l, t, theta0, dtheta0, n)
    out_h2 = nihalo(l, t, theta0, dtheta0, 2*n)

    # Poračunamo lokalno napako tako, da vzamemo vsak drugi korak
    # in izračunamo razliko med rezultati za koraka h in h/2.
    local_err = abs.(out_h2[1:2:end, 1] - out_h[:, 1])
 
    return local_err
end


"""
    odmik = harmonicno_nihalo(l,t,theta0,dtheta0,n)

Funkcija simulira nihanje harmoničnega nihala, kjer uporabimo
aproksimacijo sinθ ≈ θ.

# Argumenti
- `l::Float64`: dolžina nihala
- `t::Float64`: čas simulacije
- `theta0::Float64`: začetni odmik v radianih
- `dtheta0::Float64`: začetna hitrost
- `n::Int`: število korakov
"""
function harmonicno_nihalo(l,t,theta0,dtheta0, n)
    function harmonic_dF(y, l, g=9.80665)
        θ, ω = y
        dθ = ω
        dω = -g/l * θ
        return [dθ, dω]
    end
    y0 = [theta0, dtheta0]
    return rk4((x,y) -> harmonic_dF(y, l), y0, t, n)
end


"""

    plot_pendulum(odmik, t)

Funkcija iz rezultatov simulacije nihala izriše grafe
odmika, kotne hitrosti in energije sistema.

# Argumenti
- `odmik`: tabela odmikov in kotnih hitrosti nihala
- `t::Float64`: čas simulacije
"""
function plot_pendulum(odmik, t, l)
    theta = odmik[:, 1]
    dtheta = odmik[:, 2]

    n, _ = size(odmik)
    xs = LinRange(0, t, n)
    xodmik = l .* sin.(theta)
    
    display(plot(xs, xodmik, label="Odmik", xlabel="Čas", ylabel="Odmik"))
    display(plot(xs, dtheta, label="Kotna hitrost", xlabel="Čas", ylabel="Kotna hitrost"))

    # Narišemo še energijo nihala
    g = 9.80665
    m = 1.0
    true_E = 0.5 * m * (l * dtheta[1])^2 + m * g * l * (1.0 - cos(theta[1]))
    E = 0.5 * m * (l * dtheta).^2 .+ m * g * l * (1.0 .- cos.(theta))
    ylims = (minimum(E), maximum(E))
    display(plot(xs, E, label="Energija", xlabel="Čas", ylims=ylims, ylabel="Energija"))
    println("Relativna napaka ohranitve energije: ", (maximum(E) - minimum(E)) / true_E)
end

"""

    t = nihajni_cas(l, theta0, dtheta0, n)

Funkcija izračuna nihajni čas matematičnega nihala, kjer najprej
uporabi Runge-Kutta metodo 4. reda za izračun odmikov in nato 
poišče časovni interval med dvema zaporednima prehodoma. Za izboljšanje
natančnosti je prehod izračun kot interpolacija med zaporednima točkama
odmika z nasprotnim predznakom.

# Argumenti
- `l::Float64`: dolžina nihala
- `theta0::Float64`: začetni odmik
- `dtheta0::Float64`: začetna hitrost
- `n::Int`: število korakov
"""
function nihajni_cas(l, theta0, dtheta0, n, tol=1e-6)
    # Definiramo čas simulacije kot 2x aproksimacije periode
    t = 2.0 * 2.0 * pi * sqrt(l / 9.80665)
    odmik = nihalo(l, t, theta0, dtheta0, n)

    # Poiščemo zaporedne čase, ko ima odmik nasproten predznak
    # in izračunamo čas med njima.
    n, _ = size(odmik)
    xs = LinRange(0, t, n)
    xodmik = l .* sin.(odmik[:, 1])
    sign_change = sign.(xodmik)
    crossings = findall(diff(sign_change) .!= 0)

    # Če ni dovolj prehodov, podvojimo čas simulacije
    while length(crossings) < 3
        t *= 2
        odmik = nihalo(l, t, theta0, dtheta0, n)
        xodmik = l .* sin.(odmik[:, 1])
        sign_change = sign.(xodmik)
        crossings = findall(diff(sign_change) .!= 0)
    end

    # Poiščemo časovni interval med dvema zaporednima prehodoma
    x0, x1 = xs[crossings[1]], xs[crossings[1]+1]
    y0, y1 = xodmik[crossings[1]], xodmik[crossings[1]+1]
    xmid1 = x1 - y1 * (x1 - x0) / (y1 - y0)

    x2, x3 = xs[crossings[3]], xs[crossings[3]+1]
    y2, y3 = xodmik[crossings[3]], xodmik[crossings[3]+1]
    xmid2 = x3 - y3 * (x3 - x2) / (y3 - y2)

    return xmid2 - xmid1
end


"""

    plot_compare(odmik_p, odmik_h, t, l)

Funkcija izriše primerjavo harmoničnega in matematičnega nihala.

# Argumenti
- `odmik_p`: odmik in hitrost matematičnega nihala
- `odmik_h`: odmik in hitrost harmoničnega nihala
- `t::Float64`: čas simulacije
- `l::Float64`: dolžina nihala
- `title::String`: naslov grafa
"""
function plot_compare(odmik_p, odmik_h, t, l, title="")
    n, _ = size(odmik_p)
    xs = LinRange(0, t, n)
    xodmik_p = l .* sin.(odmik_p[:, 1])
    xodmik_h = l .* sin.(odmik_h[:, 1])
    
    plot(xs, xodmik_p, label="Matematično")
    display(plot!(xs, xodmik_h, label="Harmonično", xlabel="Čas", ylabel="Odmik", title=title))
end
