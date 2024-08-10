"""
    y = gaussian_CDF(x)

Izračuna vrednost kumulativne porazdelitve normalne porazdelitve v točki `x` z uporabo Rombergove metode 
preko integrala funkcije napake `erf`.

# Argumenti
- `x::Float64`: vrednost, pri kateri računamo kumulativno porazdelitev
- `n::Int64`: število korakov Rombergove metode
- `tol::Float64`: zahtevana natančnost
- `ignore_tol::Bool`: če je `true`, funkcija ne vrne napake, če Rombergova metoda ne skonvergira
"""
function gaussian_cdf(x::Float64, n=15, tol=1e-6, ignore_tol=false)
    err_f = t -> exp(-t^2)
    return 0.5 + 1 / sqrt(pi) * romberg_method(err_f, 0, x / sqrt(2.0), n, tol, ignore_tol)
end

"""
    y = romberg_method(f, a, b, n, tol=1e-12)

Izračuna integral funkcije `f` na intervalu `[a, b]` z uporabo Rombergove metode v maksimalno `n` korakih 
ali dokler ni dosežena zahtevana natančnost `tol`.

# Argumenti
- `f::Function`: funkcija, katere integral računamo
- `a::Float64`: spodnja meja intervala
- `b::Float64`: zgornja meja intervala
- `n::Int64`: število korakov Rombergove metode
- `tol::Float64`: zahtevana natančnost
- `ignore_tol::Bool`: če je `true`, funkcija ne vrne napake, če Rombergova metoda ne skonvergira
"""
function romberg_method(f, a, b, n::Int64, tol, ignore_tol)
    R_prev = zeros(n)    # Prejsnja vrstica
    R = zeros(n)         # Trenutna vrstica
    h = b - a
    R_prev[1] = 0.5 * h * (f(a) + f(b))
    
    for i=2:n
        # Najprej izračunamo vsoto za naslednji nivo v tabeli po sestavljeni trapezni metodi
        h = h / 2.0
        f_sum = 0.0
        for k=1:2^(i-2)
            f_sum += f(a + (2*k - 1)*h)
        end
        R[1] = h * f_sum + 0.5 * R_prev[1]

        # Ekstrapoliramo približek integrala
        for j in 2:i
            factor = 4.0^(j-1.0)
            R[j] = (factor * R[j-1] - R_prev[j-1]) / (factor - 1.0)
        end

        # Če je razlika med zadnjima dvema vrednostima v vrstici manjša od zahtevane natančnosti, končaj
        if !ignore_tol && abs(R_prev[i-1] - R[i]) < tol && i > 2
            return R[i]
        end
        R_prev, R = R, R_prev
    end
    if !ignore_tol
        throw("Rombergova metoda ni skonvergirala po $(n) korakih s toleranco $(tol).")
    else
        return R_prev[n]
    end
end


"""
    point = hypotrochoid(t, a, b)

Izračuna točko na hipotrohoidi s parametri `a`, `b` in `t`.

# Argumenti
- `t::Float64`: vrednost parametra t
- `a::Float64`: prvi parameter hipotrohoide
- `b::Float64`: drugi parameter hipotrohoide
"""
function hypotrochoid(t, a=1.0, b=-11.0/7.0)
    x = (a + b) * cos(t) + b * cos((a + b) / b * t)
    y = (a + b) * sin(t) + b * sin((a + b) / b * t)
    return [x, y]
end

"""
    point = hypotrochoid_derivative(t, a, b)

Izračuna odvod hipotrohoide s parametri `a`, `b` in `t`.

# Argumenti
- `t::Float64`: vrednost parametra t
- `a::Float64`: prvi parameter hipotrohoide
- `b::Float64`: drugi parameter hipotrohoide
"""
function hypotrochoid_derivative(t, a=1.0, b=-11.0/7.0)
    x = - (a + b) * sin(t) - b * sin((a + b) / b * t) * (a + b) / b
    y = (a + b) * cos(t) + b * cos((a + b) / b * t) * (a + b) / b
    return [x, y]
end

"""
    presecisci = hypotrochoid_intersect(tol)

Izračuna prvi dve samopresečišči hipotrohoide s toleranco `tol`.

# Argumenti
- `tol::Float64`: zahtevana natančnost
"""
function hypotrochoid_intersect(tol=1e-20)
    hip_f = t -> hypotrochoid(t, 1.0, -11.0/7.0)
    
    # Vemo, da je prvo presečišče pri vrednosti y = 0, poiščemo vrednost t, med [8, 10] z bisekcijo
    # za dane parametre hipotrohoide
    a = 5.0
    t0 = 5.0
    b = 10.0
    maxiter = 200
    for i = 1:maxiter
        t0 = (a + b) / 2.0
        if hip_f(t0)[2] == 0.0 || (b - a) / 2.0 < tol
            break
        end
        if sign(hip_f(a)[2]) == sign(hip_f(t0)[2])
            a = t0
        else
            b = t0
        end
    end
    #println("Vrednost parametra t0 za prvo presečišče: $t0")

    # Za dane parametre vemo, da hipotrohoida sestavlja 7 simetričnih krakov, torej 
    # se naslednje samopresečišče pojavi pod kotom 2pi/7 glede na x os.
    a = 10.0
    b = 11.0
    t1 = 10.0

    # Iščemo presečišče hipotrohoide s premico y = tan(2pi/7) * x.
    dif_y = t -> hip_f(t)[2] - tan(2.0 * pi / 7.0) * hip_f(t)[1]
    for i=1:maxiter
        t1 = (a + b) / 2.0
        if dif_y(t1) == 0.0 || (b - a) / 2.0 < tol
            break
        end
        if sign(dif_y(a)) == sign(dif_y(t1))
            a = t1
        else
            b = t1
        end
    end
    #println("Vrednost parametra t1 za drugo presečišče: $t1")
    return [t0, t1]
end


"""
    area = hypotrochoid_area(tol)

Izračuna ploščino omejeno s hipotrohoido s toleranco `tol`.

# Argumenti
- `n::Int64`: število korakov Rombergove metode
- `tol::Float64`: zahtevana natančnost
- `ignore_tol::Bool`: če je `true`, funkcija ne vrne napake, če Rombergova metoda ne skonvergira
"""
function hypotrochoid_area(n::Int64, tol=1e-20, ignore_tol=false)
    # Najprej poiščemo presečišči hipotrohoide 
    t0, t1 = hypotrochoid_intersect(tol)

    # Uporabimo formulo za ploščino krivočrtnega trikotnika
    function integrand(t)
        x, y = hypotrochoid(t)
        dx, dy = hypotrochoid_derivative(t)
        return x * dy - y * dx
    end

    # Izračunamo ploščino enega kraka hipotrohoide
    area0 = 0.5 * romberg_method(integrand, t0, t1, n, tol, ignore_tol)

    # Vrnemo ploščino vseh sedmih krakov
    return 7 * area0
end



    