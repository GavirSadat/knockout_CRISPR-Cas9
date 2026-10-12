import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
_SEQS = _os.path.join(_AQUI, "seqs.py")
_EFF  = _os.environ.get("CRISPOR_EFF", _os.path.join(_RAIZ, "externo", "crisporEffScores.py"))
import re, math, csv
if not _os.path.exists(_EFF):
    raise SystemExit("Falta crisporEffScores.py de CRISPOR. Ver externo/README.md, o exporta CRISPOR_EFF=/ruta/crisporEffScores.py")
src=open(_EFF).read()
ns={'math':math}
# extraer solo las tablas de parametros y las dos funciones puras
def grab(start_pat,end_pat):
    i=src.index(start_pat); j=src.index(end_pat,i); return src[i:j]
exec(grab('paramsCRISPRscan = [','def listToSvml'),ns)
exec(grab('doenchParams = [','def calcSscScores'),ns)
# doctests de referencia del propio codigo
print("control CRISPRscan (esperado 77):",ns['calcCrisprScanScores'](["TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA"]))

exec(open(_SEQS).read().split('import re')[0])
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
cds={"OsKSL4":ksl4,"OsCPS4":cps4}
guides=[("OsKSL4","KSL4-sg1","GCGACAATCTAGCACGAATG","GGG","-"),
        ("OsKSL4","KSL4-sg2","GCCGGCAGCTTCCTTCCACG","GGG","+"),
        ("OsKSL4","KSL4-sg3","TCCGACATCGAAAAAGTCAT","CGG","-"),
        ("OsCPS4","CPS4-sg1","GCACTGAGATCAATCGGAGA","GGG","+"),
        ("OsCPS4","CPS4-sg2","GATCGTTCAGAATCAGCTAC","CGG","+"),
        ("OsCPS4","CPS4-sg3","GGTAAAGGATATCGATGACA","CGG","+")]
SCAF="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

def hairpin(a,b,minlen=4):
    """cuenta tramos complementarios >= minlen entre a y b (b se busca como complemento reverso)"""
    hits=[]
    for L in range(len(a),minlen-1,-1):
        for i in range(len(a)-L+1):
            sub=a[i:i+L]
            if rc(sub) in b:
                hits.append(sub)
        if hits: break
    return hits

rows=[]
for gen,name,proto,pam,strand in guides:
    seq = cds[gen] if strand=="+" else rc(cds[gen])
    i = seq.find(proto+pam)
    c30 = seq[i-4:i+26]; c35 = seq[i-6:i+29]
    d14 = ns['calcDoenchScores']([c30])[0]
    cscan = ns['calcCrisprScanScores']([c35])[0]
    gc = 100*(proto.count("G")+proto.count("C"))/20
    polyT = "TTTT" in proto
    self_hp = hairpin(proto,proto)
    scaf_hp = hairpin(proto,SCAF[:24])
    mm0 = sum(x.count(proto)+x.count(rc(proto)) for x in (ksl4,cps4))
    rows.append(dict(sgRNA=name,contexto30=c30,Doench2014=d14,CRISPRscan=cscan,GC=int(gc),
                     polyT=("si" if polyT else "no"),
                     autocompl=(self_hp[0] if self_hp else "-"),
                     compl_andamio=(scaf_hp[0] if scaf_hp else "-"),
                     MM0_en_CDS=mm0))
    print(f"{name}  30mer={c30}  Doench2014={d14}  CRISPRscan={cscan}  GC={int(gc)}%  polyT={polyT}  "
          f"autocompl={self_hp or '-'}  vs_andamio={scaf_hp or '-'}  MM0={mm0}")
with open(OUT + 'chopchop_scores.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
