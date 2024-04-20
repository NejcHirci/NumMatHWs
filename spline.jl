using LinearAlgebra

"""
    Zlepek

Struktura, ki predstavlja naravni kubični zlepek.
Vsebuje `točke``, ki jih zlepek interpolira, in koeficiente zlepka v `koef`.
"""
struct Zlepek
    tocke::Matrix{Float64}
    koef::Matrix{Float64}
end

"""
    interpoliraj(x, y)

Funkcija interpolira dane točke (x, y) z naravnim kubičnim zlepkom.

Argument `x` je vektor x koordinat točk

Argument `y`  je vektor y koordinat točk.
"""
function interpoliraj(x::Array{Float64}, y::Array{Float64})
    n = length(x)
    if n < 2
        throw(ArgumentError("Potrebno je vsaj dve točki za interpolacijo."))
    end

    if n == 2
        println("Dve točki")
        # Če imamo samo dve točki, vrnemo linearni zlepek
        return Zlepek([(x[1], y[1]), (x[2], y[2])], [y[1] 0 0 0; y[2] 0 0 0])
    end

    h = zeros(n)
    b = zeros(n)
    v = zeros(n)
    u = zeros(n)
    
    # Izvedemo vnaprejsnji izracun vrednosti h, b, v, u za konstrukcijo linearnega sistema
    for i in 1:n-1
        h[i] = x[i+1] - x[i]
        b[i] = (y[i+1] - y[i]) / h[i]
        if i > 1
            v[i] = 2 * (h[i-1] + h[i])
            u[i] = 6 * (b[i] - b[i-1])
        end
    end

    A = Tridiagonal(h[1:end-3], v[2:end-1], h[1:end-3])
    knots = A \ u[2:end-1]

    # Dodaj robne vozle za naravni zlepek z0 = zn = 0
    knots = [0; knots; 0]

    koef = zeros(n-1, 4)
    tocke = zeros(n, 2)
    tocke[:, 1] = x
    tocke[:, 2] = y
    # Izračunaj koeficiente iz vozlov in vrednosti
    for i in 1:n-1
        koef[i, 1] = knots[i+1] / (6 * h[i])
        koef[i, 2] = knots[i] / (6 * h[i])
        koef[i, 3] = y[i+1] / h[i] - knots[i+1] * h[i] / 6
        koef[i, 4] = y[i] / h[i] - knots[i] * h[i] / 6
    end
    return Zlepek(tocke, koef)
end

"""
    vrednost(zlepek, x)

Funkcija vrne vrednost zlepka v točki x.

Argument `zlepek` je objekt strukture `Zlepek`.

Argument `x` je vrednost, v kateri želimo izračunati vrednost zlepka.
"""
function vrednost(zlepek::Zlepek, x::Float64)
    n = size(zlepek.tocke, 1)
    i = 1
    if x <= zlepek.tocke[2,1]
        i = 1
    elseif x >= zlepek.tocke[n-1,1]
        i = n-1
    else
        for j = 1:n-1
            if zlepek.tocke[j,1] <= x && x <= zlepek.tocke[j+1,1]
                i = j
                break
            end
        end
    end
    
    a = zlepek.koef[i, 1]
    b = zlepek.koef[i, 2]
    c = zlepek.koef[i, 3]
    d = zlepek.koef[i, 4]

    xi0 = zlepek.tocke[i, 1]
    xi1 = zlepek.tocke[i+1, 1]
    
    return a * (x - xi0)^3 + b * (xi1 - x)^3 + c * (x - xi0) + d * (xi1 - x)
end

using Plots
"""
    plotZlepek(zlepek)

Funkcija nariše graf zlepka, tako da na vsakem segmentu zlepka nariše želeno število točk.

Argument `zlepek` je objekt strukture `Zlepek`.

Argument `tocke_na_seg` je število točk, ki jih želimo narisati na vsakem segmentu zlepka.
"""
function plotZlepek(Z::Zlepek, tocke_na_seg::Integer=10)
    n = size(Z.tocke, 1)
    x = []
    y = []
    colors = []
    for i in 1:n-1
        x_new = range(Z.tocke[i, 1], Z.tocke[i+1, 1], length=tocke_na_seg)
        y_new = [vrednost(Z, xi) for xi in x_new]
        x = [x; x_new]
        y = [y; y_new]
        col = i % 2 == 0 ? "red" : "blue"
        colors = [colors; fill(col, length(x_new))]
    end

    Plots.plot(Z.tocke[:,1], Z.tocke[:,2], seriestype = :scatter, color = :black, label = "Točke", aspect_ratio=1)
    Plots.plot!(x, y, seriestype = :line, color = colors, label = "Zlepek", aspect_ratio=1)
end