#!/usr/bin/env python3
"""Paso 5. Reconstruye la Tabla 4 del reporte desde cero: toma las coordenadas de
CDS de cada guia en sgRNAs_y_oligos.csv, las mapea a coordenadas genomicas de
IRGSP-1.0 con la API REST de Ensembl Plants y comprueba cuatro cosas:

  1. que cada 23-mero (protospacer + PAM) devuelve un UNICO bloque de mapeo,
     que es la prueba de que no cruza una frontera exon-intron;
  2. a que exon del transcrito cae cada uno;
  3. que la distancia genomica entre los cortes de sg1 y sg2 es la delecion
     corta anunciada (104 pb en OsKSL4, 119 pb en OsCPS4);
  4. que la distancia entre los cortes de sg1 y sg3 es la delecion larga
     (3.964 y 3.053 pb).

Sale con codigo 1 si alguna comprobacion falla.

Uso:  python3 13_coordenadas_ensembl.py
      python3 13_coordenadas_ensembl.py --autotest   # sin red: valida la aritmetica
                                                     # contra los bloques registrados
Requiere: conexion a rest.ensembl.org. Sin dependencias externas."""
import os as _os, sys, csv, json, time
import urllib.request, urllib.error

_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
REST  = _os.environ.get("ENSEMBL_REST", "https://rest.ensembl.org")

TRANSCRITO = {"OsKSL4": "Os04t0179700-01", "OsCPS4": "Os04t0178300-03"}
ESPERADO_CORTA = {"OsKSL4": 104, "OsCPS4": 119}
ESPERADO_LARGA = {"OsKSL4": 3964, "OsCPS4": 3053}

def get(ruta, intentos=4):
    url = REST + ruta + ("&" if "?" in ruta else "?") + "content-type=application/json"
    for n in range(intentos):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 503) and n < intentos - 1:
                time.sleep(2 * (n + 1)); continue
            raise
        except Exception:
            if n < intentos - 1:
                time.sleep(2 * (n + 1)); continue
            raise

def mapear(transcrito, ini, fin):
    """CDS [ini,fin] -> bloques genomicos. Devuelve la lista de bloques."""
    d = get(f"/map/cds/{transcrito}/{ini}..{fin}")
    return d.get("mappings", [])

def exones(transcrito):
    d = get(f"/lookup/id/{transcrito}?expand=1")
    return [(e["start"], e["end"]) for e in d.get("Exon", [])]

def exon_de(pos, lista):
    for n, (a, b) in enumerate(lista, 1):
        if a <= pos <= b:
            return n
    return None


# --- autocomprobacion sin red -------------------------------------------------
# Bloques de mapeo registrados el 9-oct-2026 contra rest.ensembl.org (IRGSP-1.0).
# Permiten validar la aritmetica CDS<->genoma de este script sin conexion.
BLOQUES_REGISTRADOS = {
    "KSL4-sg1": dict(a=98,   b=120,  start=5433845, end=5433867, strand=-1, corte=103,  genomico=5433862),
    "KSL4-sg2": dict(a=191,  b=213,  start=5433752, end=5433774, strand=-1, corte=207,  genomico=5433758),
    "KSL4-sg3": dict(a=1769, b=1791, start=5429881, end=5429903, strand=-1, corte=1774, genomico=5429898),
    "CPS4-sg1": dict(a=268,  b=290,  start=5319075, end=5319097, strand=1,  corte=284,  genomico=5319091),
    "CPS4-sg2": dict(a=387,  b=409,  start=5319194, end=5319216, strand=1,  corte=403,  genomico=5319210),
    "CPS4-sg3": dict(a=1086, b=1108, start=5322128, end=5322150, strand=1,  corte=1102, genomico=5322144),
}

def autotest():
    fallos = 0
    pos = {}
    for nom, d in BLOQUES_REGISTRADOS.items():
        off = d["start"] - d["a"] if d["strand"] == 1 else d["end"] + d["a"]
        g = (d["corte"] + off) if d["strand"] == 1 else (off - d["corte"])
        ok = g == d["genomico"]
        fallos += not ok
        pos[nom] = g
        print(f"{'PASA' if ok else 'FALLA'}  {nom}: corte CDS {d['corte']} -> {g:,}".replace(",", "."))
    for pre, corta, larga in (("KSL4", 104, 3964), ("CPS4", 119, 3053)):
        c = abs(pos[f"{pre}-sg1"] - pos[f"{pre}-sg2"])
        l = abs(pos[f"{pre}-sg1"] - pos[f"{pre}-sg3"])
        for et, obt, esp in (("corta", c, corta), ("larga", l, larga)):
            ok = obt == esp
            fallos += not ok
            print(f"{'PASA' if ok else 'FALLA'}  {pre} delecion {et}: {obt} pb (esperado {esp})")
    print("\nFALLOS:", fallos)
    return fallos

def main():
    filas = list(csv.DictReader(open(OUT + "sgRNAs_y_oligos.csv")))
    if not filas:
        sys.exit("sgRNAs_y_oligos.csv vacio o no encontrado en " + OUT)

    ex = {g: exones(t) for g, t in TRANSCRITO.items()}
    for g, lista in ex.items():
        print(f"{g}: {len(lista)} exones en {TRANSCRITO[g]}")

    cortes, fallos = {}, 0
    print(f"\n{'sgRNA':10} {'23-mero en IRGSP-1.0':28} {'corte':>10} {'exon':>5}  bloques")
    for f in filas:
        gen, nom = f["gen"], f["sgRNA"]
        ini, fin = int(f["inicio_CDS"]), int(f["fin_CDS"])
        corte_cds = int(f["corte_CDS"])
        # el 23-mero incluye el PAM: aguas abajo en hebra +, aguas arriba en hebra -
        a, b = (ini, fin + 3) if f["hebra"] == "+" else (ini - 3, fin)
        bloques = mapear(TRANSCRITO[gen], a, b)
        n = len(bloques)
        if n != 1:
            print(f"  FALLO {nom}: {n} bloques de mapeo -> cruza frontera exon-intron")
            fallos += 1
            continue
        m = bloques[0]
        gi, gf, hebra = m["start"], m["end"], m["strand"]
        # relacion lineal CDS<->genoma dentro del bloque
        if hebra == 1:
            off = gi - a
            g_de = lambda p: p + off
        else:
            off = gf + a
            g_de = lambda p: off - p
        corte = g_de(corte_cds)
        e = exon_de(corte, ex[gen])
        cortes[nom] = corte
        print(f"{nom:10} {gi:,}-{gf:,}".replace(",", ".").ljust(40)
              + f"{corte:>10,}".replace(",", ".") + f"{e:>6}  {n}")

    print()
    for gen, pre in (("OsKSL4", "KSL4"), ("OsCPS4", "CPS4")):
        try:
            corta = abs(cortes[f"{pre}-sg1"] - cortes[f"{pre}-sg2"])
            larga = abs(cortes[f"{pre}-sg1"] - cortes[f"{pre}-sg3"])
        except KeyError:
            fallos += 1; continue
        for etiqueta, obt, esp in (("corta sg1-sg2", corta, ESPERADO_CORTA[gen]),
                                   ("larga sg1-sg3", larga, ESPERADO_LARGA[gen])):
            ok = obt == esp
            fallos += not ok
            print(f"{gen} delecion {etiqueta}: {obt} pb "
                  f"({'coincide' if ok else 'NO COINCIDE'} con los {esp} pb del reporte)")

    print("\nFALLOS:", fallos)
    sys.exit(1 if fallos else 0)

if __name__ == "__main__":
    if "--autotest" in sys.argv:
        sys.exit(1 if autotest() else 0)
    main()
