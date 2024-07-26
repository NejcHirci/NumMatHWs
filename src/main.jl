using Plots

"""
    yn = dF(xn, m, M)

Funkcija vrne rezultat desne strani sistema diferencialnih enačb
glede na začetna pogoja `m` in `M`.

# Argumenti
- `xn::Vector{Float64}`: trenutni položaj in hitrost sonde

"""
function dF(tn, xn, m=7.3 * 10^22, M=5.97 * 10^24)
    mean = m / (M+m)
    r_M = (x,y,z) -> sqrt((x+mean)^2 + y^2 + z^2)
    r_m = (x,y,z) -> sqrt((x-mean)^2 + y^2 + z^2)

    # [x]   = [x_prev + v_x * h]
    # [y]   = [y_prev + v_y * h]
    # [z]   = [z_prev + v_z * h]
    # [v_x] = [v_x_prev + dv_x]
    # [v_y] = [v_y_prev + dv_x]
    # [v_z] = [v_z_prev + dv_x]

    x, y, z, v_x, v_y, v_z = xn
    x = x + v_x * tn
    y = y + v_y * tn
    z = z + v_z * tn

    curR_M = r_M(x,y,z)
    curR_m = r_m(x,y,z)

    dv_x = x + 2 * v_y - (1-mean)/curR_M^3 * (x + mean) - mean/curR_m^3 * (x-mean)
    dv_y = y - 2 * v_x - (1-mean)/curR_M^3 * y - mean/curR_m^3 * y
    dv_z = - (1-mean)/curR_M^3 * z - mean/curR_m^3 * z

    return Vector{Float64}([x, y, z, v_x + dv_x * tn, v_y + dv_y * tn, v_z + dv_z * tn])
end


"""
    yn = rk4(f, yn, h, steps=100)

Funkcija izračuna položaj sonde v omejenem sistemu treh masnih
teles, glede na podan začetni položaj `yn`, korak `h` in funkcijo
položaja `f`. Implementacija uporablja Runge-Kutta metodo 4. reda.
"""
function rk4(f, yn, h, steps=100)
    out_yn = zeros(steps, 6)
    out_tn = zeros(steps, 1)

    tn = 0.0
    for i=1:steps
        k1 = f(tn, yn)
        k2 = f(tn + h/2, yn + h * k1/2)
        k3 = f(tn + h/2, yn + h * k2/2)
        k4 = f(tn + h, yn + h * k3)
        yn = yn + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        tn = tn + h
        out_yn[i, :] = yn
        out_tn[i] = tn
    end
    println("Last position: ", yn)
    return out_yn, out_tn
end


"""
    plot_2d(yn, tn)

Funkcija izriše 2D graf položaja sonde v omejenem sistemu treh masnih
teles, glede na podane podatke `yn` in `tn`.
"""
function plot_2d(yn, tn)
    x = [i[1] for i in yn]
    y = [i[2] for i in yn]
    plot(x, y, label="Sonda", xlabel="x", ylabel="y", title="2D graf položaja sonde")
end