# 3. domača naloga

Avtor: **Nejc Hirci**

Izvorna koda za reševanje sistema diferencialnih enačb matematičnega nihala.

## Opis naloge

Naloga opiše gibanje nedušenega nitnega nihala podanega z diferencialno enačbo drugega reda. Nalogo bomo reševali z metodo Runge-Kutta 4. reda. Implementirali bomo tudi funkcijo za izračun periode nihanja nihala in jo preverili z analitično rešitvijo.

## Uporaba implementacije

Rešitev je implementirana v Julia programskem jeziku, pri čemer se želene funkcije nahajo v datoteki `src/main.jl`. Primer uporabe je prikazan v Jupyter zvezku `docs/porocilo.ipynb`.

## Zagon testov

Testi so implementirani v datoteki `test/runtests.jl` in se zaženejo z ukazom:

```julia
include("test/runtests_part1.jl")
```

## Ustvarjanje poročila

Za analizo rezultatov in testiranje rešitve sem uporabil Jupyter zvezek in programski jezik Julia. Za pravilno uporabo je potrebno aktivirati projekt in instancirati okolje z zahtevanimi paketi. To storimo z ukazom:

```julia
# V paketnem načinu naložimo okolje
activate .
instantiate
```

Za generiranje poročila poženemo Jupyter zvezek z ukazi:
```julia
using IJulia;
notebook(dir=".")
```

v brskalniku se bo odprl Jupyter vmesnik, kjer lahko odpremo zvezek `docs/porocilo.ipynb` in ga izvedemo. Poročilo je bilo nato generirano z 
izvozom v obliki LaTeX dokumenta (potrebno je imeti nameščen [Inkscape](https://inkscape.org/) in dodan v `PATH`). Nato pa je bil dokument še prevede v PDF
format z uporabo `pdflatex`.
