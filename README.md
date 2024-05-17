# 2. domača naloga

Avtor: **Nejc Hirci**

Izvorna koda za nalogi implementacije porazdelitvene funkcije normalne porazdelitve in izračun ploščine hipotrohoide.

## Porazdelitvena funkcija normalne slučajne spremenljivke

### Opis naloge

Napišite učinkovito funkcijo, ki izračuna vrednosti porazdelitvene funkcije za standardno normalno porazdeljeno slučajno spremenljivko $X \sim N(0, 1)$.

### Uporaba implementacije

Rešitev je implementirana v Julia programskem jeziku, pri čemer se žele funkcije nahajo v datoteki `src/main.jl`. Funkcija `gaussian_cdf(x, tol=1e-10)` z uporabo Rombergove metode numerično izračuna vrednost porazdelitvene funkcije normalne slučajne spremenljivke za vrednost `x`. Uporabnik lahko določi tudi natančnost izračuna z izbirnim argumentom `tol`. Primer uporabe je prikazan v Jupyter zvezku `docs/porocilo.ipynb` in pa v datoteki `docs/demo.jl`.

### Zagon testov

Testi so implementirani v datoteki `test/runtests_part1.jl` in se zaženejo z ukazom:

```julia
include("test/runtests_part1.jl")
```

## Ploščina hipotrohoide

### Opis naloge

Izračunajte ploščino območja, ki ga omejuje [hipotrohoida](https://en.wikipedia.org/wiki/Hypotrochoid) podana s parametričnimi enačbami.

### Uporaba implementacije

### Zagon testov

Testi so implementirani v datoteki `test/runtests_part2.jl` in se zaženejo z ukazom:

```julia
include("test/runtests_part2.jl")
```

## Ustvarjanje poročila

Za analizo rezultatov in testiranje rešitve sem uporabil Jupyter zvezek. Zato je potrebno za ustvarjanje poročila namestiti Jupyter:

```julia
using Pkg
Pkg.add("IJulia")
Pkg.actiavte(".")
```

Nato lahko z uporabo ukaza `jupyter notebook` poženemo Jupyter in odpremo `docs/porocilo.ipynb`. V njem lahko nato izvozimo poročilo v obliki HTML ali PDF datoteke.