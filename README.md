# 3. domača naloga

Avtor: **Nejc Hirci**

Izvorna koda za reševanje sistema diferencialnih enačb omejenega krožnega problema treh teles.

### Opis naloge

Naloga opiše gibanje sonde v sistemu treh teles, kjer se sonda giblje v gravitacijskem polju Zemlje in Lune.


### Uporaba implementacije

Rešitev je implementirana v Julia programskem jeziku, pri čemer se žele funkcije nahajo v datoteki `src/main.jl`. Primer uporabe je prikazan v Jupyter zvezku `docs/porocilo.ipynb`.

### Zagon testov

Testi so implementirani v datoteki `test/runtests_part1.jl` in se zaženejo z ukazom:

```julia
include("test/runtests_part1.jl")
```

## Ploščina hipotrohoide

### Opis naloge

Izračunajte ploščino območja, ki ga omejuje [hipotrohoida](https://en.wikipedia.org/wiki/Hypotrochoid) podana s parametričnimi enačbami.

### Uporaba implementacije

Implementirana je funkcija `hypotrochoid_area(tol=1e-10)`, ki najde dve samo-presečišči v krivulji hipotrohoide (začetek in konec enega kraka) in preko
Rombergove metode izračuna površino, ki jo omejuje krivulja.

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
