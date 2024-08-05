using Plots

"""
    yn = dF(xn, m, M)

Funkcija vrne rezultat desne strani sistema diferencialnih enačb
glede na začetna pogoja `m` in `M`.

# Argumenti
- `xn::Vector{Float64}`: trenutni položaj in hitrost sonde

"""
function dF(xn, rat=0.1)
    # [dxdt] = [v_x]
    # [dydt] = [v_y]
    # [dzdt] = [v_z]
    # [dvxdt] = [dv_x]
    # [dvydt] = [dv_y]
    x, y, z, v_x, v_y, v_z = xn

    # Razdalja med sondo in Zemljo
    r_M = sqrt((x+rat)^2 + y^2 + z^2)

    # Razdalja med sondo in Luno
    r_m = sqrt((x-1.0+rat)^2 + y^2 + z^2)

    dv_x = 
        2 * v_y
        + x 
        - (1.0-rat) * (x + rat) / r_M^3 
        - rat / r_m^3 * (x - 1.0 + rat)
    dv_y = 
        -2 * v_x 
        + y
        - (1.0-rat) * y / r_M^3 
        - rat * y / r_m^3
    dv_z = 
        -(1.0-rat) / r_M^3 * z 
        - rat / r_m^3 * z
    return Vector{Float64}([v_x, v_y, v_z, dv_x, dv_y, dv_z])
end


"""
    yn = rocket_rk4(yn, mass_moon=1, mass_earth=1, h=0.01, steps=100)

Funkcija izračuna položaj sonde v omejenem sistemu treh masnih
teles, glede na podan začetni položaj `yn`, korak `h` in funkcijo
položaja `f`. Implementacija uporablja Runge-Kutta metodo 4. reda.
"""
function rocket_rk4(yn, mass_moon=1, mass_earth=1, h=0.01, steps=100)
    f = x -> dF(x, mass_moon/(mass_earth + mass_moon))

    out_yn = zeros(steps, 6);
    out_tn = zeros(steps, 1);

    out_yn[1, :] = yn

    tn = 0.0
    for i=2:steps
        # Update the speed and position of everything
        k1 = f(yn)
        k2 = f(yn + h*k1/2.0)
        k3 = f(yn + h*k2/2.0)
        k4 = f(yn + h*k3)
        yn = yn + h/6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        tn = tn + h
        out_yn[i, :] = yn
        out_tn[i] = tn
    end
    return out_yn, out_tn
end


"""
    find_moon_flyby()

Funkcija, ki z uporabo strelske metode izračuna optimalne začetne pogoje
za obhod Lune. Za optimizacijo uporabljam bisekcijo
"""
function find_moon_flyby(init_low, init_high, rk4_steps=10000, bisect_steps=100, bisect_tol=1e-6)
    # t0=1.0 predstavlja en obhod Lune v nasem koordinatnem sistemu
    # 10 dni v enotah obhoda Lune
    tend = 10 / 27.3 

    mass_earth = 5.974e24;
    mass_moon = 7.348e22;
    rat = mass_moon / (mass_earth + mass_moon)

    rk4_h = tend / rk4_steps
    println("h: ", rk4_h)
    rk4_F = x -> dF(x, rat)

    # Objektivna funkcija za strelsko metodo
    function objective(init_cond, target)
        rk_yn, _ = rocket_rk4(rk4_F, init_cond, rk4_h, rk4_steps)
        final_pos = rk_yn[end, 1:3]
        println("final_pos: ", final_pos)
        target_pos = target[1:3]
        error = sum(target_pos - final_pos)
        return error
    end

    x0 = init_low
    xmid = init_low
    x1 = init_high

    for i=1:bisect_steps
        f0 = objective(x0, x1)
        f1 = objective(x1, x1)
        xmid = (x0 + x1) / 2.0
        fmid = objective(xmid, xmid)
        println("f0: ", f0, "f1: ", f1, " fmid: ", fmid)

        if abs(fmid) < bisect_tol
            println("Konvergirano po ", i, " korakih")
            break
        end
        if sign(f1) == sign(fmid)
            x1 = xmid
        else
            x0 = xmid
        end
    end
    return xmid
end


"""

    dy = pendulum_dF(y, l, g)

Funkcija vrne rezultat desne strani sistema diferencialnih enačb
za matematično nihalo.
"""
function pendulum_dF(y, l, g=9.80665)
    θ, ω = y
    dθ = ω
    dω = -g/l * sin(θ)
    return [dθ, dω]
end

"""
    odmik = nihalo(l,t,theta0,dtheta0,n)

Funkcija simulira nihanje matematičnega nihala z uporabo 
Runge-Kutta 4. reda. Funkcija vrne odmik nihala v odvisnosti od časa.
"""
function nihalo(l,t,theta0,dtheta0,n)
    h = t/n
    y = [theta0, dtheta0]
    out = zeros(n, 2)
    out[1, :] = y
    for i=2:n
        k1 = pendulum_dF(y, l)
        k2 = pendulum_dF(y + h*k1/2.0, l)
        k3 = pendulum_dF(y + h*k2/2.0, l)
        k4 = pendulum_dF(y + h*k3, l)
        y = y + h/6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i, :] = y
    end
    return out
end

"""

    dy = pendulum_dF(y, l, g)

Funkcija vrne rezultat desne strani sistema diferencialnih enačb
za matematično nihalo.
"""
function harmonic_dF(y, l, g=9.80665)
    θ, ω = y
    dθ = ω
    dω = -g/l * θ
    return [dθ, dω]
end


"""
    odmik = harmonicno_nihalo(l,t,theta0,dtheta0,n)

Funkcija simulira nihanje harmoničnega nihala, kjer uporabimo
aproksimacijo sinθ = θ, kjer je začetni odmik dovolj majhen.
"""
function harmonicno_nihalo(l,t,theta0,dtheta0, n)
    h = t/n
    y = [theta0, dtheta0]
    out = zeros(n, 2)
    out[1, :] = y
    for i=2:n
        k1 = harmonic_dF(y, l)
        k2 = harmonic_dF(y + h*k1/2.0, l)
        k3 = harmonic_dF(y + h*k2/2.0, l)
        k4 = harmonic_dF(y + h*k3, l)
        y = y + h/6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i, :] = y
    end
    return out
end

"""

    plot_pendulum(odmik, t)

Funkcija iz rezultatov simulacije izriše grafe
odmika, hitrosti, pospeška in prikaz začetnega stanja.
"""
function plot_pendulum(odmik, t, l)
    n, _ = size(odmik)
    xs = LinRange(0, t, n)
    xodmik = l .* sin.(odmik[:, 1])
    
    display(plot(xs, xodmik, label="Odmik", xlabel="Čas", ylabel="Odmik"))
    display(plot(xs, odmik[:,2], label="Kotna hitrost", xlabel="Čas", ylabel="Kotna hitrost"))
end

"""
    
    plot_compare(odmik_p, odmik_h, t, l)

Funkcija izriše primerjavo harmoničnega in matematičnega nihala.
"""
function plot_compare(odmik_p, odmik_h, t, l)
    n, _ = size(odmik_p)
    xs = LinRange(0, t, n)
    xodmik_p = l .* sin.(odmik_p[:, 1])
    xodmik_h = l .* sin.(odmik_h[:, 1])
    
    plot(xs, xodmik_p, label="Matematično")
    display(plot!(xs, xodmik_h, label="Harmonično", xlabel="Čas", ylabel="Odmik"))

    plot(xs, odmik_p[:,2], label="Matematično")
    display(plot!(xs, odmik_h[:,2], label="Harmonično", xlabel="Čas", ylabel="Kotna hitrost"))
end
