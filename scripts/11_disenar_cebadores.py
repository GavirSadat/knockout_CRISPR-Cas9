#!/usr/bin/env python3
"""Paso 4. Disena los cebadores de genotipado sobre el locus genomico, que es
lo unico valido: los pares de corte pueden atravesar intrones y un cebador
disenado sobre CDS no existiria en el ADN de la planta.

Para cada gen produce dos amplicones: uno que abarca los cortes de sg1 y sg2,
cuya delecion se lee por desplazamiento de tamano, y otro centrado en sg3.

Uso: python3 11_disenar_cebadores.py OsKSL4_locus_genomico.gb OsCPS4_locus_genomico.gb
Requiere: biopython."""
import sys, csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

MIN_L, MAX_L = 18, 25
TM_MIN, TM_MAX = 57.0, 63.0
GC_MIN, GC_MAX = 40.0, 60.0
MARGEN = 60          # distancia minima entre el cebador y el corte mas cercano

def valido(p):
    if not (GC_MIN <= gc_fraction(p) * 100 <= GC_MAX): return None
    if any(b * 4 in p for b in "ACGT"): return None
    if p[-1] not in "GC": return None                      # pinza 3'
    if p[-5:].count("G") + p[-5:].count("C") > 3: return None   # evita 3' demasiado GC
    tm = mt.Tm_NN(Seq(p))
    return tm if TM_MIN <= tm <= TM_MAX else None

def buscar(seq, ini, fin, sentido):
    """Mejor cebador dentro de la ventana [ini, fin) en el sentido indicado."""
    mejores = []
    for L in range(MIN_L, MAX_L + 1):
        for i in range(ini, fin - L):
            p = seq[i:i + L] if sentido == "F" else str(Seq(seq[i:i + L]).reverse_complement())
            tm = valido(p)
            if tm is None: continue
            if seq.count(p) + seq.count(str(Seq(p).reverse_complement())) != 1: continue
            mejores.append((abs(tm - 60), p, i, i + L, round(tm, 1),
                            round(gc_fraction(p) * 100)))
    return sorted(mejores)[0] if mejores else None

def main(rutas):
    filas = []
    for ruta in rutas:
        rec = SeqIO.read(ruta, "genbank")
        s = str(rec.seq).upper()
        guias = sorted([(int(f.location.start), int(f.location.end),
                         (f.qualifiers.get("label", [""])[0]))
                        for f in rec.features if f.type == "misc_binding"])
        if len(guias) != 3:
            sys.exit(f"{ruta}: se esperaban 3 guias anotadas y hay {len(guias)}")
        gen = rec.name
        print(f"\n== {gen}  ({len(s)} pb)  guias en {[g[2] for g in guias]}")
        pares = [("corta (sg1+sg2)", guias[0][0], guias[1][1]),
                 ("sg3", guias[2][0], guias[2][1])]
        for etiqueta, a, b in pares:
            F = buscar(s, max(0, a - 450), a - MARGEN, "F")
            R = buscar(s, b + MARGEN, min(len(s), b + 450), "R")
            if not (F and R):
                print(f"  {etiqueta}: no se encontro par valido; amplia la ventana")
                continue
            amp = R[3] - F[2]
            print(f"  amplicon {etiqueta}: {amp} pb")
            print(f"     F {F[1]}  Tm {F[4]} GC {F[5]}%  pos {F[2]+1}-{F[3]}")
            print(f"     R {R[1]}  Tm {R[4]} GC {R[5]}%  pos {R[2]+1}-{R[3]}")
            filas.append(dict(gen=gen, amplicon=etiqueta, tamano_wt=amp,
                              cebador_F=F[1], tm_F=F[4], gc_F=F[5],
                              cebador_R=R[1], tm_R=R[4], gc_R=R[5]))
    with open("cebadores_genotipado.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas[0].keys())); w.writeheader(); w.writerows(filas)
    print("\nEscrito: cebadores_genotipado.csv")

if __name__ == "__main__":
    main(sys.argv[1:] or sys.exit("Indica los .gb de los loci genomicos"))
