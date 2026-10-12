#!/usr/bin/env python3
"""Paso 1. Descarga de NCBI las secuencias del diseno, comprueba que las CDS
usadas son identicas a las del registro, baja la region genomica del locus,
mapea los exones y verifica que cada protospacer (con su PAM) cae dentro de un
unico exon. Escribe un informe TSV y un GenBank del locus anotado.

Uso:  python3 01_verificar_ncbi.py --email tu@correo.edu
Requiere: biopython, conexion a NCBI."""
import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
import argparse, csv, sys
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from lib_edicion import mapear_exones, protospacer_en_exon, rc

GENES = {
    "OsKSL4": dict(mrna="AB126934", locus="Os04g0179700", gene_id="4335094",
                   cds_local=OUT + "OsKSL4_CDS_sgRNAs.gb"),
    "OsCPS4": dict(mrna="AB066270", locus="Os04g0178300", gene_id=None,
                   cds_local=OUT + "OsCPS4_CDS_sgRNAs.gb"),
}

def q(f, clave, por_defecto=""):
    """Biopython devuelve los qualifiers como listas al leer un GenBank y como
    cadenas al construirlo en memoria; esta funcion normaliza ambos casos."""
    v = f.qualifiers.get(clave, por_defecto)
    return v[0] if isinstance(v, list) and v else (v if isinstance(v, str) else por_defecto)

def comparar_cds(local, ncbi):
    """Devuelve (veredicto, posiciones_1based_que_difieren)."""
    if local == ncbi:
        return "IDENTICAS", []
    if local == ncbi + local[-3:] and local[-3:] in ("TAA", "TAG", "TGA"):
        return "IDENTICAS (el registro anota la CDS sin el codon de parada)", []
    if len(local) != len(ncbi):
        return f"LONGITUD DISTINTA: {len(local)} vs {len(ncbi)}", []
    dif = [i + 1 for i, (a, b) in enumerate(zip(local, ncbi)) if a != b]
    return f"{len(dif)} sustituciones", dif

def cds_de_registro(acc):
    """Devuelve la CDS anotada en el registro de GenBank."""
    with Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as h:
        rec = SeqIO.read(h, "genbank")
    for f in rec.features:
        if f.type == "CDS":
            return str(f.extract(rec.seq)).upper(), rec
    raise ValueError(f"{acc} no tiene feature CDS")

def region_genomica(gene_id, locus, margen=3000):
    """Coordenadas del locus en el cromosoma y descarga de la region."""
    if gene_id is None:
        with Entrez.esearch(db="gene", term=f"{locus}[All Fields] AND Oryza sativa[Organism]") as h:
            gene_id = Entrez.read(h)["IdList"][0]
    with Entrez.esummary(db="gene", id=gene_id) as h:
        info = Entrez.read(h)["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
    acc = info["ChrAccVer"]
    ini, fin = int(info["ChrStart"]) + 1, int(info["ChrStop"]) + 1
    hebra = 1 if ini <= fin else 2
    a, b = (min(ini, fin) - margen, max(ini, fin) + margen)
    with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text",
                       seq_start=a, seq_stop=b, strand=hebra) as h:
        rec = SeqIO.read(h, "fasta")
    return str(rec.seq).upper(), acc, a, b, hebra

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True)
    ap.add_argument("--salida", default="verificacion_exones.tsv")
    args = ap.parse_args()
    Entrez.email = args.email
    filas, problemas = [], 0

    for nombre, meta in GENES.items():
        local = SeqIO.read(meta["cds_local"], "genbank")
        cds_local = str(local.seq).upper()
        cds_ncbi, _ = cds_de_registro(meta["mrna"])
        veredicto, difs = comparar_cds(cds_local, cds_ncbi)
        identica = veredicto.startswith("IDENTICAS")
        print(f"[{nombre}] CDS local {len(cds_local)} pb vs {meta['mrna']} "
              f"{len(cds_ncbi)} pb -> {veredicto}")
        if difs:
            print(f"   posiciones que difieren (nt de la CDS): {difs[:40]}"
                  + (" ..." if len(difs) > 40 else ""))

        gen, acc, a, b, hebra = region_genomica(meta["gene_id"], meta["locus"])
        bloques = mapear_exones(cds_local, gen)
        print(f"[{nombre}] {acc}:{a}-{b} hebra {hebra}; {len(bloques)} exones; "
              f"intrones: {[bloques[i+1][2]-bloques[i][3] for i in range(len(bloques)-1)]}")
        # prueba decisiva: la CDS reconstruida desde los exones del genoma debe
        # ser identica a la del diseno. El genoma manda sobre el clon de cDNA.
        reconstruida = "".join(gen[g0:g1] for _, _, g0, g1 in bloques)
        print(f"   CDS reconstruida desde el genoma {acc}: "
              f"{'IDENTICA a la del diseno' if reconstruida == cds_local else 'DIFIERE -> revisar'}")
        if reconstruida != cds_local:
            problemas += 1

        glob = SeqRecord(Seq(gen), id=f"{nombre}_locus", name=nombre,
                         description=f"{nombre} locus genomico {acc}:{a}-{b}")
        glob.annotations.update(molecule_type="DNA", topology="linear",
                                organism="Oryza sativa Japonica Group")
        for n, (c0, c1, g0, g1) in enumerate(bloques, 1):
            glob.features.append(SeqFeature(FeatureLocation(g0, g1), type="exon",
                                 qualifiers={"label": f"exon {n}", "note": f"CDS {c0+1}-{c1}"}))

        guias = [f for f in local.features
                 if f.type == "misc_binding" and "protospacer" in q(f, "label")]
        if not guias:
            print(f"   AVISO: {meta['cds_local']} no tiene protospacers anotados")
        for f in guias:
            etiqueta = q(f, "label").split()[0]
            ini, fin = int(f.location.start) + 1, int(f.location.end)
            # el PAM va inmediatamente 3' del protospacer en la hebra de la guia
            ini_t, fin_t = (ini, fin + 3) if f.location.strand == 1 else (ini - 3, fin)
            ok, coord = protospacer_en_exon(bloques, ini_t, fin_t)
            # una sustitucion del clon de cDNA dentro de la guia o de su PAM la invalida
            choques = [d for d in difs if ini_t <= d <= fin_t]
            print(f"   {etiqueta}: {'OK, exonico' if ok else 'RECHAZADA, cruza un intron'}"
                  + (f"  genomico {coord[0]}-{coord[1]}" if ok else "")
                  + (f"  OJO: difiere del registro en {choques}" if choques else ""))
            if not ok or choques:
                problemas += 1
            else:
                glob.features.append(SeqFeature(
                    FeatureLocation(coord[0] - 1, coord[1], strand=f.location.strand),
                    type="misc_binding", qualifiers={"label": etiqueta}))
            filas.append(dict(gen=nombre, sgRNA=etiqueta, cds_ini=ini, cds_fin=fin,
                              exonico="si" if ok else "NO",
                              gen_ini=coord[0] if ok else "", gen_fin=coord[1] if ok else "",
                              cds_vs_ncbi=veredicto,
                              difiere_dentro_de_la_guia=",".join(map(str, choques)) or "no"))
        SeqIO.write(glob, f"{nombre}_locus_genomico.gb", "genbank")

    if not filas:
        sys.exit("No se anoto ninguna guia: revisa que los .gb esten en benchling/")
    with open(args.salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(filas)
    print(f"\nInforme: {args.salida}. Problemas detectados: {problemas}")
    sys.exit(1 if problemas else 0)

if __name__ == "__main__":
    main()
