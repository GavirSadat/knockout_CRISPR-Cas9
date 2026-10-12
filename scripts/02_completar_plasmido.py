#!/usr/bin/env python3
"""Paso 2. Descarga el vector binario de NCBI, lo digiere in silico con BsaI,
ensambla el arreglo de sgRNA y escribe el plasmido completo en GenBank, sin
tramos de N.

Uso:  python3 02_completar_plasmido.py --email tu@correo.edu \\
          --vector KR029109 --array ../benchling/sgRNA_array_ksl4.gb \\
          --salida pMOMI-ksl4_completo.gb
Requiere: biopython, conexion a NCBI."""
import argparse, sys
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from lib_edicion import sitios_bsai, digerir_vector, ensamblar

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True)
    ap.add_argument("--vector", default="KR029109")
    ap.add_argument("--array", required=True)
    ap.add_argument("--salida", required=True)
    ap.add_argument("--permitir-N", action="store_true",
                    help="ensambla aunque el arreglo conserve marcadores de posicion")
    args = ap.parse_args()
    Entrez.email = args.email

    with Entrez.efetch(db="nuccore", id=args.vector, rettype="gb", retmode="text") as h:
        vec = SeqIO.read(h, "genbank")
    vseq = str(vec.seq).upper()
    print(f"Vector {args.vector}: {len(vseq)} pb, {len(vec.features)} features")

    sitios = sitios_bsai(vseq)
    print("Sitios BsaI:", [(p, oh, o) for p, oh, o in sitios])
    esqueleto, oh5, oh3 = digerir_vector(vseq)
    print(f"Esqueleto {len(esqueleto)} pb; el inserto debe llevar {oh5} en 5' y {oh3} en 3'")

    arr = SeqIO.read(args.array, "genbank")
    ins = str(arr.seq).upper()
    if "N" in ins[4:-4] and not args.permitir_N:
        sys.exit("El arreglo aun tiene marcadores de posicion (N): sustituye los "
                 "promotores OsU6 por su secuencia real o usa --permitir-N")
    ins = oh5 + ins[4:-4] + oh3                       # fija las salientes del vector
    plasmido = ensamblar(vseq, ins)
    print(f"Plasmido ensamblado: {len(plasmido)} pb; sitios BsaI restantes: "
          f"{len(sitios_bsai(plasmido))}")

    rec = SeqRecord(Seq(plasmido), id=args.salida.split('.')[0], name="pMOMI",
                    description=f"{args.vector} con el arreglo de {args.array}")
    rec.annotations.update(molecule_type="DNA", topology="circular",
                           organism="synthetic construct")
    desp = len(esqueleto)
    rec.features.append(SeqFeature(FeatureLocation(desp, desp + len(ins)),
                        type="misc_feature", qualifiers={"label": "arreglo sgRNA insertado"}))
    # traslada las anotaciones del esqueleto: el vector se abrio y se rote, asi que
    # cada feature cambia de coordenada segun el tramo en el que estaba
    sitios = sitios_bsai(vseq); p1 = sitios[0][0]; p2 = sitios[1][0]
    cola = len(vseq) - (p2 + 4)          # longitud del tramo que ahora va primero
    trasladadas = 0
    for f in vec.features:
        a, b = int(f.location.start), int(f.location.end)
        if a >= p2 + 4 and b <= len(vseq):
            na, nb = a - (p2 + 4), b - (p2 + 4)
        elif b <= p1:
            na, nb = a + cola, b + cola
        else:
            continue                      # cruza el relleno extraido: se descarta
        rec.features.append(SeqFeature(FeatureLocation(na, nb, strand=f.location.strand),
                                       type=f.type, qualifiers=dict(f.qualifiers)))
        trasladadas += 1
    print(f"Anotaciones del esqueleto trasladadas: {trasladadas} de {len(vec.features)}")
    for f in arr.features:                             # traslada la anotacion del inserto
        if int(f.location.start) >= 4 and int(f.location.end) <= len(ins):
            rec.features.append(SeqFeature(
                FeatureLocation(desp + int(f.location.start), desp + int(f.location.end),
                                strand=f.location.strand),
                type=f.type, qualifiers=dict(f.qualifiers)))
    SeqIO.write(rec, args.salida, "genbank")
    print("Escrito:", args.salida)

if __name__ == "__main__":
    main()
