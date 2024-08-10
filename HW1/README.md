# NM-Spline

Avtor: **Nejc Hirci**

Izvorna koda za nalogo naravnega kubičnega zlepka pri predmetu Numerične metode na FRI.

## Opis Naloge

Naloga implementira algoritem za izračun naravnega kubičnega zlepka, ki je interpolacijski zlepek sestavljen iz kubičnih polinomov, tako vrednosti polinomov ustrezajo podanim vrednostim $(x_i, f_i)$ in ustreza pogojem podrobneje opisanim v poročilu.

## Navodila za uporabo

Rešitev je implementirana v Julia programskem jeziku, pri čemer se želene funkcije nahajajo v datoteki `spline.jl`. Primer uporabe skript je na voljo v Jupyter zvezku `report.ipynb` in pa v datoteki `demo.jl`. Za pravilno delovanje je potrebno namestiti paket `Plots`.

## Zagon testov

Testi so implementirani v datoteki `runtests.jl` in se zaženejo z ukazom:

```julia
include("runtests.jl")
```

## Ustvarjanje poročila

Za analizo rezultatov in testiranje rešitve sem uporabil Jupyter zvezek. Zato je potrebno za
ustvarjanje poročila namestiti Jupyter:

```julia
using Pkg
Pkg.add("IJulia")
```

Nato lahko z uporabo ukaza `jupyter notebook` poženemo Jupyter in odpremo `report.ipynb`. V njem lahko nato izvozimo poročilo v obliki HTML ali PDF datoteke.