#!/usr/bin/env python3
"""Paso 3. Off-target. Dos modos:

  a) exhaustivo sobre un FASTA dado (paralogos, un cromosoma, un locus):
     python3 03_offtarget.py --fasta paralogos.fa --max-mm 3
  b) genoma completo delegando en bowtie, que es lo que usa CHOPCHOP por dentro:
     python3 03_offtarget.py --bowtie-index IRGSP-1.0 --max-mm 3 \
         --genoma-fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
     Escribe guias.fa, lanza bowtie -v 3 -a --best --strata y, con --genoma-fa,
     recupera del genoma los 3 nt contiguos a cada alineamiento y descarta los
     que no llevan PAM NGG: bowtie alinea los 20 nt de la guia y es ciego al PAM.

Salida: offtarget.tsv con gen, guia, sitio, desapareamientos totales y en la
semilla de 12 nt, y si el sitio lleva PAM NGG funcional."""
import argparse, csv, subprocess, sys
from lib_edicion import buscar_offtarget, desapareamientos, rc

GUIAS = {
 "KSL4-sg1": "GCGACAATCTAGCACGAATG", "KSL4-sg2": "GCCGGCAGCTTCCTTCCACG",
 "KSL4-sg3": "TCCGACATCGAAAAAGTCAT",  "CPS4-sg1": "GCACTGAGATCAATCGGAGA",
 "CPS4-sg2": "GATCGTTCAGAATCAGCTAC", "CPS4-sg3": "GGTAAAGGATATCGATGACA"}
PARALOGOS = {"LOC_Os02g36210": "OsCPS2", "LOC_Os02g36140": "OsKSL7",
             "LOC_Os11g28530": "OsKSL8", "LOC_Os12g30824": "OsKSL10"}

def leer_fasta(ruta):
    nombre, trozos = None, []
    for linea in open(ruta):
        if linea.startswith(">"):
            if nombre: yield nombre, "".join(trozos)
            nombre, trozos = linea[1:].strip(), []
        else:
            trozos.append(linea.strip())
    if nombre: yield nombre, "".join(trozos)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta"); ap.add_argument("--bowtie-index")
    ap.add_argument("--genoma-fa", help="FASTA del genoma, para comprobar el PAM de cada hit")
    ap.add_argument("--guias-fa", help="FASTA con guias a cribar en vez del set del diseno")
    ap.add_argument("--max-mm", type=int, default=3)
    ap.add_argument("--salida", default="offtarget.tsv")
    args = ap.parse_args()
    filas = []

    if args.fasta:
        for cab, seq in leer_fasta(args.fasta):
            for nombre, g in GUIAS.items():
                for pos, hebra, sitio, pam, mm, mm_semilla in buscar_offtarget(g, seq, args.max_mm):
                    filas.append(dict(guia=nombre, referencia=cab, posicion=pos, hebra=hebra,
                                      sitio=sitio, PAM=pam, mm=mm, mm_semilla=mm_semilla,
                                      nota="EXACTO" if mm == 0 else
                                           ("CRITICO: semilla intacta" if mm_semilla == 0 else "")))
    elif args.bowtie_index:
        entrada = args.guias_fa or "guias.fa"
        if not args.guias_fa:
            with open("guias.fa", "w") as fh:
                for n, g in GUIAS.items():
                    fh.write(f">{n}\n{g}\n")
        else:
            GUIAS.clear()
            for linea in open(entrada):
                if linea.startswith(">"): nom = linea[1:].strip()
                elif linea.strip(): GUIAS[nom] = linea.strip().upper()
        # sin --strata: con esa bandera bowtie reporta solo el mejor estrato y,
        # como cada guia tiene su sitio perfecto, ocultaria los de 1-3
        # desapareamientos, que son precisamente los off-target que se buscan
        cmd = ["bowtie", "-v", str(args.max_mm), "-a", "--best",
               "-f", "-x", args.bowtie_index, entrada, "--sam", "offtarget.sam"]
        print("Ejecutando:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        hits = []
        for linea in open("offtarget.sam"):
            if linea.startswith("@"): continue
            c = linea.rstrip("\n").split("\t")
            if c[2] == "*": continue
            nm = [x for x in c[11:] if x.startswith("NM:i:")]
            hits.append(dict(guia=c[0], ref=c[2], pos=int(c[3]),
                             hebra="-" if int(c[1]) & 16 else "+", sitio=c[9],
                             mm=int(nm[0][5:]) if nm else 0))
        print(f"bowtie devolvio {len(hits)} alineamientos; comprobando PAM...")
        pams = {}
        if args.genoma_fa:
            # el PAM esta inmediatamente 3' del protospacer en la hebra del hit
            necesarios = {}
            for h in hits:
                L = len(h["sitio"])
                a, b = ((h["pos"] + L, h["pos"] + L + 2) if h["hebra"] == "+"
                        else (h["pos"] - 3, h["pos"] - 1))     # 1-based inclusivo
                necesarios.setdefault(h["ref"], []).append((a, b, id(h)))
            nombre, trozos = None, []
            def resolver(nombre, seq):
                for a, b, k in necesarios.get(nombre, []):
                    if a < 1 or b > len(seq): continue
                    p = seq[a - 1:b]
                    pams[k] = p if len(p) == 3 else ""
            for linea in open(args.genoma_fa):
                if linea.startswith(">"):
                    if nombre: resolver(nombre, "".join(trozos))
                    nombre, trozos = linea[1:].split()[0], []
                elif nombre in necesarios:
                    trozos.append(linea.strip())
            if nombre: resolver(nombre, "".join(trozos))
        for h in hits:
            p = pams.get(id(h), "")
            if args.genoma_fa and h["hebra"] == "-" and p:
                p = rc(p)
            if args.genoma_fa and p and p[1:].upper() != "GG":
                continue                                   # sin PAM NGG no hay corte
            filas.append(dict(guia=h["guia"], referencia=h["ref"], posicion=h["pos"],
                              hebra=h["hebra"], sitio=h["sitio"],
                              PAM=p.upper() if p else "sin verificar",
                              mm=h["mm"], mm_semilla="",
                              nota="EXACTO" if h["mm"] == 0 else
                                   ("CRITICO: <=2 desapareamientos" if h["mm"] <= 2 else "")))
    else:
        sys.exit("Indica --fasta o --bowtie-index")

    with open(args.salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(filas)
    criticos = [f for f in filas if f["mm"] and f["mm"] <= 2]
    print(f"{len(filas)} sitios con PAM NGG tras el filtro")
    for n in sorted(GUIAS):
        prop = [f for f in filas if f["guia"] == n]
        ex = sum(1 for f in prop if f["mm"] == 0)
        print(f"  {n}: {ex} exacto(s), "
              + ", ".join(f"{sum(1 for f in prop if f['mm']==k)} a {k} mm" for k in (1, 2, 3)))
    print(f"{len(criticos)} sitios con <=2 desapareamientos: si alguno no es el propio, "
          "esa guia se descarta")
    print("Parologos a revisar nominalmente:", ", ".join(PARALOGOS.values()))

if __name__ == "__main__":
    main()
