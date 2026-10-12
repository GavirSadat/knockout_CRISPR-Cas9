"""Sustituye las NLS nativas del vector (SV40 en N, nucleoplasmina en C) por la
bipartita de VirD2 con codones de arroz, y repara la anotacion del promotor
OsU6a que se habia perdido al trasladar features."""
import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
IN    = _os.environ.get("MOMI_IN", OUT)   # donde 02_completar_plasmido.py dejo los *_completo.gb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
RICE = {'A':'GCC','R':'AGG','N':'AAC','D':'GAC','C':'TGC','Q':'CAG','E':'GAG','G':'GGC',
        'H':'CAC','I':'ATC','L':'CTG','K':'AAG','M':'ATG','F':'TTC','P':'CCG','S':'AGC',
        'T':'ACC','W':'TGG','Y':'TAC','V':'GTG'}
opt = lambda p: "".join(RICE[a] for a in p)
VIRD2, LINK, FLAG = "KRPRDRHDGELGGRKRAR", "SGGSGGS", "DYKDHDGDYKDHDIDYKDDDDK"
SV40_NT  = "ATGGCTCCTAAGAAGAAGCGGAAGGTT"                        # M A P K K K R K V
NUCL_NT  = "AAGCGTCCTGCTGCCACCAAAAAGGCCGGACAGGCTAAGAAAAAGAAG"  # KRPAATKKAGQAKKKK
NUEVO_N  = "ATG" + opt(VIRD2) + opt(LINK)
NUEVO_C  = opt(LINK) + opt(VIRD2) + opt(FLAG)
BAD = ["GGTCTC","GAGACC","CGTCTC","GAGACG"]

def mover(f, cortes):
    """Desplaza un feature segun las sustituciones hechas aguas arriba."""
    d = lambda p: p + sum(delta for ini, delta in cortes if ini < p)
    return SeqFeature(FeatureLocation(d(int(f.location.start)), d(int(f.location.end)),
                                      strand=f.location.strand),
                      type=f.type, qualifiers=dict(f.qualifiers))

for src, gen in [("pMOMI-ksl4_completo.gb","ksl4"), ("pMOMI-cps4_completo.gb","cps4")]:
    ruta = _os.path.join(IN, src)
    if not _os.path.exists(ruta):
        raise SystemExit(f"Falta {ruta}. Ejecuta antes 02_completar_plasmido.py.")
    rec = SeqIO.read(ruta, "genbank")
    s = str(rec.seq).upper()
    i_n = s.find(SV40_NT); i_c = s.find(NUCL_NT)
    assert i_n >= 0 and i_c > i_n, (gen, i_n, i_c)
    nuevo = s[:i_n] + NUEVO_N + s[i_n+len(SV40_NT):i_c] + NUEVO_C + s[i_c+len(NUCL_NT):]
    cortes = [(i_n, len(NUEVO_N)-len(SV40_NT)), (i_c, len(NUEVO_C)-len(NUCL_NT))]
    out = rec[:0]; out.seq = Seq(nuevo); out.id = f"pMOMI-{gen}-VirD2"; out.name = f"pMOMI_{gen}_VirD2"
    out.annotations.update(molecule_type="DNA", topology="circular", organism="synthetic construct")
    out.description = (f"pMOMI-{gen} con las NLS del vector sustituidas por la bipartita "
                       "de VirD2 en ambos extremos, codones de arroz")
    for f in rec.features:
        g = mover(f, cortes)
        if "Cas9" in str(f.qualifiers.get("product","")):
            q = dict(g.qualifiers); q.pop("translation", None)
            q["label"] = "NLS VirD2 - Cas9p - NLS VirD2 - 3xFLAG"
            q["note"] = "Cas9 de codones de planta del vector, con ambas NLS sustituidas"
            g = SeqFeature(g.location, type="CDS", qualifiers=q)
        out.features.append(g)
    # anotaciones nuevas
    ini_n = i_n; ini_c = i_c + cortes[0][1]
    out.features.append(SeqFeature(FeatureLocation(ini_n, ini_n+len(NUEVO_N)), type="misc_feature",
        qualifiers={"label":"NLS bipartita VirD2 N-terminal + enlazador","note":"codones de arroz"}))
    out.features.append(SeqFeature(FeatureLocation(ini_c, ini_c+len(NUEVO_C)), type="misc_feature",
        qualifiers={"label":"enlazador + NLS VirD2 C-terminal + 3xFLAG","note":"codones de arroz"}))
    # promotor OsU6a, que se habia perdido
    arr = [f for f in out.features if "arreglo sgRNA insertado" in str(f.qualifiers.get("label",""))][0]
    a0 = int(arr.location.start)
    out.features.append(SeqFeature(FeatureLocation(a0, a0+447), type="promoter",
        qualifiers={"label":"promotor OsU6a","note":"secuencia real tomada del plasmido pYLsgRNA"}))
    # controles
    cds = [f for f in out.features if f.type=="CDS" and "VirD2" in str(f.qualifiers.get("label",""))][0]
    prot = str(Seq(str(cds.extract(out.seq))).translate())
    ok_frame = prot.count("*")==1 and prot.endswith("*")
    print(f"pMOMI-{gen}: {len(s)} -> {len(nuevo)} pb")
    print(f"   N-term: {prot[:30]}")
    print(f"   C-term: ...{prot[-32:]}")
    print(f"   marco correcto: {ok_frame} | sin BsaI/BsmBI nuevos: "
          f"{sum(nuevo.count(b) for b in BAD)==sum(s.count(b) for b in BAD)}")
    SeqIO.write(out, OUT + f"pMOMI-{gen}_VirD2.gb", "genbank")
