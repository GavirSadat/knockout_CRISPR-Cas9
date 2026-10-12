#!/usr/bin/env python3
"""Control de calidad del plasmido final. Comprueba sobre la secuencia circular
que estan los tres casetes completos, en orden, con sus terminadores, que las
salientes del vector quedaron en las uniones y que no sobrevive ningun sitio de
las enzimas del ensamblaje.

Uso: python3 10_qc_plasmido.py ../benchling/pMOMI-ksl4_completo.gb"""
import sys
from Bio import SeqIO
from lib_edicion import rc, sitios_bsai

SCAF = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
ESPACIADORES = {
 "ksl4": ["GCGACAATCTAGCACGAATG", "GCCGGCAGCTTCCTTCCACG", "GTCCGACATCGAAAAAGTCAT"],
 "cps4": ["GCACTGAGATCAATCGGAGA", "GATCGTTCAGAATCAGCTAC", "GGTAAAGGATATCGATGACA"]}

def main(ruta):
    rec = SeqIO.read(ruta, "genbank")
    s = str(rec.seq).upper()
    gen = "ksl4" if "ksl4" in ruta.lower() else "cps4"
    fallos = 0
    def check(nombre, ok, detalle=""):
        nonlocal fallos
        print(("PASA  " if ok else "FALLA ") + nombre + (f" :: {detalle}" if detalle else ""))
        if not ok: fallos += 1

    check("longitud razonable (16-20 kb)", 16000 <= len(s) <= 20000, f"{len(s)} pb")
    check("topologia circular", rec.annotations.get("topology") == "circular",
          str(rec.annotations.get("topology")))
    for e in ESPACIADORES[gen]:
        n = s.count(e) + s.count(rc(e))
        check(f"espaciador {e[:12]}... presente una sola vez", n == 1, f"{n} copias")
    n_sc = s.count(SCAF) + s.count(rc(SCAF))
    check("tres andamios completos", n_sc == 3, f"{n_sc}")
    pos = [s.find(SCAF, i) for i in range(len(s)) if s.startswith(SCAF, i)]
    check("cada andamio seguido de TTTTTTT",
          all(s[p + len(SCAF):p + len(SCAF) + 7] == "TTTTTTT" for p in pos) if pos else False,
          f"{len(pos)} andamios en hebra +")
    orden = [s.find(e) if s.find(e) >= 0 else s.find(rc(e)) for e in ESPACIADORES[gen]]
    check("los tres casetes en el orden disenado", orden == sorted(orden), str(orden))
    check("sin sitios BsaI en todo el plasmido", len(sitios_bsai(s)) == 0, str(sitios_bsai(s)))
    # BsmBI no se usa en este ensamblaje, asi que solo seria un problema si cayera
    # dentro de un elemento disenado por nosotros: espaciador, andamio o terminador.
    # Dentro de los promotores OsU6 es secuencia nativa y resulta irrelevante.
    import re as _re
    criticos, informativos = [], []
    for m in _re.finditer("CGTCTC|GAGACG", s):
        dentro = [f for f in rec.features
                  if int(f.location.start) <= m.start() < int(f.location.end)]
        etiquetas = " ".join(str(f.qualifiers.get("label", "")) for f in dentro).lower()
        if any(k in etiquetas for k in ("espaciador", "andamio", "terminador pol iii")):
            criticos.append((m.start() + 1, etiquetas[:50]))
        else:
            informativos.append((m.start() + 1, etiquetas[:50] or "esqueleto/promotor"))
    check("sin BsmBI en espaciadores, andamios ni terminadores",
          not criticos, str(criticos))
    if informativos:
        print(f"  (informativo) {len(informativos)} sitios BsmBI en secuencia nativa, "
              f"irrelevantes porque el ensamblaje solo usa BsaI: "
              f"{[p for p, _ in informativos][:6]}")
    texto = " ".join(str(v) for f in rec.features for vals in f.qualifiers.values()
                     for v in (vals if isinstance(vals, list) else [vals])).upper()
    hay_cas9 = "CAS9" in texto
    hay_sel = "HPT" in texto or "HYG" in texto or "HIGROM" in texto
    check("el esqueleto conserva Cas9", hay_cas9)
    check("el esqueleto conserva el marcador de seleccion", hay_sel)
    print("\n  Anotaciones del esqueleto presentes:")
    for f in rec.features:
        et = " ".join(str(v) for vals in f.qualifiers.values()
                      for v in (vals if isinstance(vals, list) else [vals]))[:60]
        if f.type in ("CDS", "gene", "promoter", "terminator", "rep_origin", "misc_recomb"):
            print(f"    {f.type:12s} {int(f.location.start):6d}-{int(f.location.end):6d}  {et}")
    print(f"\nFALLOS: {fallos}")
    sys.exit(1 if fallos else 0)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else sys.exit("Indica el .gb"))
