"""
    y = gaussian_CDF(x)

Izračuna vrednost kumulativne porazdelitve normalne porazdelitve v točki `x` z uporabo Rombergove metode 
preko integrala funkcije napake `erf`.

# Argumenti
- `x::Float64`: vrednost, pri kateri računamo kumulativno porazdelitev
- `n::Int64`: število korakov Rombergove metode
- `tol::Float64`: zahtevana natančnost
"""
function gaussian_cdf(x::Float64, n=50, tol=1e-20)
    err_f = t -> exp(-t^2)
    return 0.5 + 1 / sqrt(pi) * romberg_method(err_f, 0, x / sqrt(2.0), n, tol)
end

"""
    y = romberg_method(f, a, b, n, tol=1e-12)

Izračuna integral funkcije `f` na intervalu `[a, b]` z uporabo Rombergove metode v maksimalno `n` korakih 
ali dokler ni dosežena zahtevana natančnost `tol`.
"""
function romberg_method(f, a, b, n, tol)
    R = zeros(n, n)  # Rombergova tabela
    R[1,1] = (b-a)/2 * (f(a) + f(b))

    for i=2:n
        # Najprej izračunamo vsoto za naslednji nivo v tabeli po sestavljeni trapezni metodi
        h = (b - a) / 2.0^(i-1)
        trapezoidal_sum = sum(f(a + (2*k - 1)*h) for k in 1:2^(i-2.0))
        R[i, 1] = 0.5 * R[i-1, 1] + h * trapezoidal_sum

        # Ekstrapoliramo približek integrala
        for j in 2:i
            R[i, j] = R[i, j-1] + 1.0/(4.0^(j-1.0) - 1.0) * (R[i, j-1] - R[i-1, j-1])
        end

        # Če je razlika med zadnjima dvema vrednostima v vrstici manjša od zahtevane natančnosti, končaj
        if abs(R[i, i-1] - R[i, i]) < tol
            return R[i, i]
        end
    end
    throw(Error("Rombergova metoda ni skonvergirala po $(n) korakih."))
end


"""
    point = hipotrohoida(a, b, t)

Izračuna točko na hipotrohoidi s parametri `a`, `b` in `t`.
"""
function hipotrohoida(a, b, t)
    x = (a + b) * cos(t) + b * cos((a + b) / b * t)
    y = (a + b) * sin(t) + b * sin((a + b) / b * t)
    return [x, y]
end