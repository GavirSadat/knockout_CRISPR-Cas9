import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
_SEQS = _os.path.join(_AQUI, "seqs.py")
_EFF  = _os.environ.get("CRISPOR_EFF", _os.path.join(_RAIZ, "externo", "crisporEffScores.py"))
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
BAD=["GGTCTC","GAGACC","CGTCTC","GAGACG"]

# ---- codones de alto GC3 tipicos de arroz (japonica) ----
RICE={'A':'GCC','R':'AGG','N':'AAC','D':'GAC','C':'TGC','Q':'CAG','E':'GAG','G':'GGC','H':'CAC',
      'I':'ATC','L':'CTG','K':'AAG','M':'ATG','F':'TTC','P':'CCG','S':'AGC','T':'ACC','W':'TGG',
      'Y':'TAC','V':'GTG'}
def opt(pep): return "".join(RICE[a] for a in pep)

LINK="SGGSGGS"; VIRD2="KRPRDRHDGELGGRKRAR"; FLAG3="DYKDHDGDYKDHDIDYKDDDDK"
nterm="ATG"+opt(VIRD2)+opt(LINK)
cterm=opt(LINK)+opt(VIRD2)+opt(FLAG3)+"TGA"
for n,s in (("N-term NLS",nterm),("C-term NLS+FLAG",cterm)):
    gc3=100*sum(1 for i in range(2,len(s.rstrip('TGA')),3) if s[i] in "GC")/max(1,len(s)//3)
    print(n,len(s),"pb  BsaI/BsmBI libre:",not any(b in s or b in rc(s) for b in BAD),
          " traduccion:",Seq(s).translate())

SCAF="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
TERM="TTTTTTT"
spacers={"ksl4":[("KSL4-sg1","GCGACAATCTAGCACGAATG"),("KSL4-sg2","GCCGGCAGCTTCCTTCCACG"),
                 ("KSL4-sg3","GTCCGACATCGAAAAAGTCAT")],
         "cps4":[("CPS4-sg1","GCACTGAGATCAATCGGAGA"),("CPS4-sg2","GATCGTTCAGAATCAGCTAC"),
                 ("CPS4-sg3","GGTAAAGGATATCGATGACA")]}
PROM=[("OsU6a","pYLsgRNA-OsU6a, Addgene 66194"),("OsU6b","pYLsgRNA-OsU6b, Addgene 66196"),
      ("OsU6c","pYLsgRNA-OsU6c, sistema pYL (Ma et al. 2015)")]
OH=["AGGA","CTAC"]        # salientes BsaI internas verificadas ortogonales
PLEN=350                  # longitud nominal del promotor OsU6 (placeholder)

def feat(rec,a,b,typ,label,note=None,strand=1):
    q={"label":label}
    if note: q["note"]=note
    rec.features.append(SeqFeature(FeatureLocation(a,b,strand=strand),type=typ,qualifiers=q))

def array(gen):
    seq=""; recfeats=[]
    rec=SeqRecord(Seq(""),id=f"array_{gen}",name=f"sgRNA_array_{gen}")
    parts=[]
    seq+="NNNN"; parts.append((0,4,"misc_feature","saliente 5' del vector",
        "4 nt generados por BsaI en pYLCRISPR/Cas9Pubi-H (KR029109): leer del vector e insertar aqui"))
    for i,(name,sp) in enumerate(spacers[gen]):
        pn,src=PROM[i]
        a=len(seq); seq+="N"*PLEN
        parts.append((a,len(seq),"promoter",f"promotor {pn} (marcador de posicion)",
            f"~{PLEN} pb; copiar la secuencia real de {src}. Pol III, inicia en G"))
        a=len(seq); seq+=sp
        parts.append((a,len(seq),"misc_RNA",f"{name} espaciador",
            f"{len(sp)} nt; {'incluye G anadida en 5 prima (no apareada con el genoma)' if len(sp)==21 else 'coincide con el genoma'}"))
        a=len(seq); seq+=SCAF
        parts.append((a,len(seq),"misc_RNA","andamio sgRNA (SpCas9)","76 nt, andamio canonico"))
        a=len(seq); seq+=TERM
        parts.append((a,len(seq),"terminator",f"terminador Pol III {name}","TTTTTTT"))
        if i<2:
            a=len(seq); seq+=OH[i]
            parts.append((a,len(seq),"misc_feature",f"saliente BsaI {OH[i]}",
                "cicatriz de 4 nt entre casetes; verificada ortogonal (>=2 desapareamientos frente a las demas y a sus complementos reversos)"))
    a=len(seq); seq+="NNNN"
    parts.append((a,len(seq),"misc_feature","saliente 3' del vector","segunda saliente BsaI del vector"))
    rec.seq=Seq(seq)
    for a,b,t,l,n in parts: feat(rec,a,b,t,l,n)
    rec.description=f"Arreglo de 3 casetes sgRNA para {gen.upper()} (inserto Golden Gate BsaI en pYLCRISPR/Cas9Pubi-H)"
    rec.annotations.update(molecule_type="DNA",topology="linear",date="02-OCT-2026",
        organism="synthetic construct",source="synthetic construct",
        comment="DISENO. Las regiones N son marcadores de posicion de secuencias que se toman del vector parental "
                "(Addgene 66187 / GenBank KR029109) y de los plasmidos pYLsgRNA; todo lo demas es secuencia definitiva. "
                "En Benchling: importar el vector, digerir con BsaI y ensamblar este inserto.")
    SeqIO.write(rec,OUT+f"sgRNA_array_{gen}.gb","genbank")
    return len(seq)

for g in spacers: print("array",g,array(g),"pb")

# ---- mapa de T-DNA (diseno con marcadores de posicion) ----
def tdna(gen):
    arr=str(SeqIO.read(OUT+f"sgRNA_array_{gen}.gb","genbank").seq)
    blocks=[("LB",25,"misc_recomb","borde izquierdo del T-DNA (pCAMBIA1300)",None),
            ("2x35S",800,"promoter","promotor 2x CaMV 35S",None),
            ("HPT",1026,"CDS","higromicina fosfotransferasa (seleccion en planta)",None),
            ("T35S",250,"terminator","terminador CaMV 35S",None)]
    seq=""; parts=[]
    for name,ln,typ,lab,_ in blocks:
        a=len(seq); seq+="N"*ln
        parts.append((a,len(seq),typ,lab+" (marcador de posicion)",f"{ln} pb nominales; tomar de KR029109"))
    a=len(seq); seq+=arr
    parts.append((a,len(seq),"misc_feature",f"arreglo sgRNA {gen.upper()} (3 casetes)","secuencia definitiva salvo los promotores OsU6"))
    a=len(seq); seq+="N"*1990
    parts.append((a,len(seq),"promoter","promotor Pubi de maiz (marcador de posicion)","~1990 pb; tomar de KR029109"))
    a=len(seq); seq+=nterm
    parts.append((a,len(seq),"CDS","NLS bipartita VirD2 N-terminal + enlazador","M-KRPRDRHDGELGGRKRAR-SGGSGGS, codones de arroz"))
    a=len(seq); seq+="N"*4101
    parts.append((a,len(seq),"CDS","ORF Cas9p (marcador de posicion)","4101 pb sin codon de parada; Cas9 optimizada para codones de planta, KR029109"))
    a=len(seq); seq+=cterm
    parts.append((a,len(seq),"CDS","enlazador + NLS bipartita VirD2 + 3xFLAG + parada","SGGSGGS-KRPRDRHDGELGGRKRAR-DYKDHDGDYKDHDIDYKDDDDK, codones de arroz; segunda copia de la NLS"))
    a=len(seq); seq+="N"*253
    parts.append((a,len(seq),"terminator","terminador NOS (marcador de posicion)","253 pb nominales"))
    a=len(seq); seq+="N"*25
    parts.append((a,len(seq),"misc_recomb","borde derecho del T-DNA",None))
    rec=SeqRecord(Seq(seq),id=f"pMOMI-{gen}",name=f"pMOMI_{gen}",
        description=f"Diseno del T-DNA pMOMI-{gen}: 3 sgRNA + Cas9p con NLS + HPT")
    for a,b,t,l,n in parts: feat(rec,a,b,t,l,n)
    rec.annotations.update(molecule_type="DNA",topology="linear",date="02-OCT-2026",
        organism="synthetic construct",source="synthetic construct",
        comment="ARCHIVO DE DISENO, no una secuencia verificada. Los tramos de N son piezas del vector parental "
                "pYLCRISPR/Cas9Pubi-H (Addgene 66187 / GenBank KR029109), con longitud nominal: importar el vector "
                "en Benchling y sustituirlas. Secuencia definitiva: espaciadores, andamios, terminadores Pol III, "
                "salientes BsaI y las fusiones NLS.")
    SeqIO.write(rec,OUT+f"pMOMI-{gen}_TDNA_diseno.gb","genbank")
    print("T-DNA",gen,len(seq),"pb")
for g in spacers: tdna(g)
