inlcude("spline.jl")

# Primer uporabe

# Definiramo točke skozi katere želimo, da gre krivulja
x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
y = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# Definiramo zlepek
zlepek = interpoliraj(x, y)

# Narišemo zlepek
plotZlepek(zlepek)

# Izračunamo vrednost zlepka v točki 2.5
println("Vrednost v točki x=2.5:" vrednost(zlepek, 2.5))