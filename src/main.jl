using Plots

"""
    yn = dF(xn, m, M)

Funkcija vrne rezultat desne strani sistema diferencialnih enačb
glede na začetna pogoja `m` in `M`.

# Argumenti
- `xn::Vector{Float64}`: trenutni položaj in hitrost sonde

"""
function dF(xn, rat=0.0123)
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

    dv_x = x + 2 * v_y - (1.0-rat) / r_M^3 * (x + rat) - rat / r_m^3 * (x - 1.0 + rat)
    dv_y = y - 2 * v_x - (1.0-rat) / r_M^3 * y - rat / r_m^3 * y
    dv_z = - (1.0-rat) / r_M^3 * z - rat / r_m^3 * z

    return Vector{Float64}([v_x, v_y, v_z, dv_x, dv_y, dv_z])
end


"""
    yn = rocket_rk4(f, yn, h, steps=100)

Funkcija izračuna položaj sonde v omejenem sistemu treh masnih
teles, glede na podan začetni položaj `yn`, korak `h` in funkcijo
položaja `f`. Implementacija uporablja Runge-Kutta metodo 4. reda.
"""
function rocket_rk4(f, yn, h, steps=100)
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
za obhod Lune. Za optimizacijo uporabljam sekantno metodo.
"""
function find_moon_flyby(init_guess, rk4_steps=1000, secant_steps=100, secant_tol=1e-6)
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
        target_pos = target[1:3]
        error = sqrt(sum((final_pos - target_pos).^2))
        return error
    end

    x0 = init_guess
    x1 = init_guess + [-1.0, 0.0, 0.0, 0.0, 20.0, 0.0]

    for i=1:secant_steps
        f0 = objective(x0, x0)
        f1 = objective(x1, x1)

        if abs(f1 - f0) < secant_tol
            println("Konvergirano po ", i, " korakih")
            break
        end

        if sqrt(sum((x0 - x1).^2)) < secant_tol
            println("Konvergirano po ", i, " korakih")
            break
        end

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, x1 = x1, x2
    end
    println("Optimalni zacetni pogoji: ", x1)
    return x1
end
