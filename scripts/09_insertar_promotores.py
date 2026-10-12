#!/usr/bin/env python3
"""Paso 2a. Sustituye los marcadores de posicion (N) de los promotores OsU6 en
el arreglo de sgRNA por la secuencia real tomada de los plasmidos pYLsgRNA.

Uso:
  python3 09_insertar_promotores.py --email tu@correo.edu \\
      --accesiones KR0291xx KR0291yy KR0291zz \\
      --array ../benchling/sgRNA_array_ksl4.gb \\
      --salida ../benchling/sgRNA_array_ksl4_completo.gb

Las accesiones van en el orden OsU6a, OsU6b, OsU6c, que es el de los casetes.
Se leen de la linea 'Accession number' de cada pagina de Addgene (#66194,
#66196, #66197). Alternativa sin NCBI: --fasta promotores.fa con tres
secuencias en el mismo orden."""
import argparse, re, sys
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from lib_edicion import rc

SCAF = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

def promotor_de_plasmido(rec, largo=400):
    """En pYLsgRNA la estructura es: promotor U6 - sitios BsaI con relleno -
    andamio - TTTTTTT. El plasmido es CIRCULAR y el registro suele empezar
    dentro del sitio de clonacion, de modo que el promotor puede quedar al
    final de la secuencia lineal: por eso se busca sobre la secuencia
    duplicada. Si el registro anota el promotor, ese feature manda."""
    s = str(rec.seq).upper()
    for f in rec.features:
        txt = " ".join(str(v) for vals in f.qualifiers.values()
                       for v in (vals if isinstance(vals, list) else [vals])).upper()
        if ("U6" in txt or "U3" in txt) and ("PROMOT" in txt or f.type in ("promoter", "regulatory")):
            sub = str(f.extract(rec.seq)).upper()
            if 150 <= len(sub) <= 800:
                return sub, f"feature anotado '{txt[:40]}'"
    for hebra, base in (("+", s), ("-", rc(s))):
        ext = base + base
        i = ext.find(SCAF[:30])
        if i < 0:
            continue
        ventana = ext[max(0, i - 2000):i]
        sitios = [m.start() + max(0, i - 2000) for m in re.finditer("GGTCTC|GAGACC", ventana)]
        if not sitios:
            continue
        corte = sitios[-2] if len(sitios) >= 2 else sitios[-1]
        if corte - largo < 0:
            continue
        prom = ext[corte - largo:corte]
        if len(prom) == largo:
            return prom, f"{largo} pb antes del primer BsaI (hebra {hebra}, pos {corte % len(base)})"
    raise ValueError("no se pudo localizar el promotor; revisa el registro a mano")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email"); ap.add_argument("--accesiones", nargs=3)
    ap.add_argument("--fasta"); ap.add_argument("--array", required=True)
    ap.add_argument("--salida", required=True); ap.add_argument("--largo", type=int, default=400)
    args = ap.parse_args()

    proms = []
    if args.fasta:
        proms = [str(r.seq).upper() for r in SeqIO.parse(args.fasta, "fasta")]
        origen = ["fasta"] * len(proms)
    elif args.accesiones:
        if not args.email: sys.exit("--email es obligatorio para consultar NCBI")
        Entrez.email = args.email
        origen = []
        for acc in args.accesiones:
            with Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as h:
                rec = SeqIO.read(h, "genbank")
            p, como = promotor_de_plasmido(rec, args.largo)
            print(f"{acc}: plasmido {len(rec.seq)} pb -> promotor de {len(p)} pb ({como})")
            if len(p) < 200:
                sys.exit(f"ABORTADO: {len(p)} pb es demasiado corto para un promotor U6. "
                         "Revisa el registro y pasa las secuencias con --fasta.")
            proms.append(p); origen.append(acc)
    else:
        sys.exit("Indica --accesiones o --fasta")
    if len(proms) != 3:
        sys.exit(f"Se esperaban 3 promotores y hay {len(proms)}")

    arr = SeqIO.read(args.array, "genbank")
    s = str(arr.seq).upper()
    tramos = list(re.finditer("N{20,}", s))
    if len(tramos) != 3:
        sys.exit(f"El arreglo tiene {len(tramos)} tramos de N, se esperaban 3")
    nuevo, fin, cambios = "", 0, []
    for t, p in zip(tramos, proms):
        nuevo += s[fin:t.start()] + p
        cambios.append((t.start(), t.end(), len(p) - (t.end() - t.start()), len(p)))
        fin = t.end()
    nuevo += s[fin:]

    def mapear(pos):
        """Traslada una coordenada del arreglo original al arreglo con promotores."""
        return pos + sum(d for ini, fin_, d, _ in cambios if fin_ <= pos)

    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = SeqRecord(Seq(nuevo), id=arr.id, name=arr.name)
    rec.annotations.update(molecule_type="DNA", topology="linear",
                           organism="synthetic construct", source="synthetic construct")
    for f in arr.features:
        a, b = int(f.location.start), int(f.location.end)
        solapa = [c for c in cambios if a < c[1] and c[0] < b]
        if solapa:                      # es el feature del promotor: toma el tramo nuevo
            ini, fin_, _, largo = solapa[0]
            na, nb = mapear(ini), mapear(ini) + largo
            q = dict(f.qualifiers)
            et = q.get("label", [""])
            et = et[0] if isinstance(et, list) else et
            q["label"] = et.replace(" (marcador de posicion)", "")
            q["note"] = "secuencia real tomada del plasmido pYLsgRNA"
        else:
            na, nb, q = mapear(a), mapear(b), dict(f.qualifiers)
        rec.features.append(SeqFeature(FeatureLocation(na, nb, strand=f.location.strand),
                                       type=f.type, qualifiers=q))
    rec.description = arr.description + " (promotores reales insertados)"
    SeqIO.write(rec, args.salida, "genbank")
    print(f"Arreglo completo: {len(s)} -> {len(nuevo)} pb; promotores sin N: "
          f"{'N' not in nuevo[4:-4]} (los NNNN de los extremos son las salientes "
          f"del vector y los fija 02_completar_plasmido.py)")
    print("Siguiente paso: 02_completar_plasmido.py --array", args.salida)

if __name__ == "__main__":
    main()
