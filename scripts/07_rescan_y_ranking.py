import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
_SEQS = _os.path.join(_AQUI, "seqs.py")
_EFF  = _os.environ.get("CRISPOR_EFF", _os.path.join(_RAIZ, "externo", "crisporEffScores.py"))
import re, math
if not _os.path.exists(_EFF):
    raise SystemExit("Falta crisporEffScores.py de CRISPOR. Ver externo/README.md, o exporta CRISPOR_EFF=/ruta/crisporEffScores.py")
src=open(_EFF).read()
ns={'math':math}
def grab(a,b):
    i=src.index(a); return src[i:src.index(b,i)]
exec(grab('paramsCRISPRscan = [','def listToSvml'),ns)
exec(grab('doenchParams = [','def calcSscScores'),ns)
exec(open(_SEQS).read().split('import re')[0])
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
BAD=["GGTCTC","GAGACC","CGTCTC","GAGACG"]
SCAF="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
def maxhp(a,b,minlen=4):
    best=0
    for L in range(len(a),minlen-1,-1):
        for i in range(len(a)-L+1):
            if rc(a[i:i+L]) in b: return L
    return best
def scan(cdsseq,other,lo,hi,requireG=True):
    out=[];L=len(cdsseq)
    for strand in "+-":
        s = cdsseq if strand=="+" else rc(cdsseq)
        for m in re.finditer("(?=([ACGT]{21}GG))", s):
            i=m.start(); p=m.group(1)[:20]; pam=m.group(1)[20:23]
            if i<6 or i+29>len(s): continue
            cut = i+17 if strand=="+" else L-(i+17)
            if not lo<=cut<=hi: continue
            if requireG and p[0]!="G": continue
            gc=100*(p.count("G")+p.count("C"))/20
            if not 40<=gc<=70: continue
            if "TTTT" in p: continue
            if re.search(r"(A{5}|C{5}|G{5}|T{5})",p): continue
            if any(b in p or b in rc(p) for b in BAD): continue
            if sum(x.count(p[-12:])+x.count(rc(p[-12:])) for x in (ksl4,cps4))>1: continue
            d14=ns['calcDoenchScores']([s[i-4:i+26]])[0]
            cs =ns['calcCrisprScanScores']([s[i-6:i+29]])[0]
            self_hp=maxhp(p,p); scaf_hp=maxhp(p,SCAF[:24])
            out.append(dict(p=p,pam=pam,st=strand,cut=cut,gc=int(gc),d14=d14,cs=cs,hp=self_hp,shp=scaf_hp,
                            cod=(cut+2)//3,tot=d14+cs))
    return out
for gen,seq,other,win,mot in (("OsKSL4",ksl4,cps4,(100,700),(1730,1830)),
                              ("OsCPS4",cps4,ksl4,(100,700),(1050,1140))):
    print("="*8,gen,"ventana temprana")
    for c in sorted(scan(seq,other,*win),key=lambda d:-d["tot"])[:8]:
        print(f'  {c["p"]} {c["pam"]} {c["st"]} corte {c["cut"]:5d} cod {c["cod"]:4d} GC{c["gc"]} D14={c["d14"]:3d} CS={c["cs"]:3d} hp={c["hp"]} shp={c["shp"]}')
    print("-"*8,gen,"ventana del motivo")
    cands=scan(seq,other,*mot) or scan(seq,other,*mot,requireG=False)
    for c in sorted(cands,key=lambda d:-d["tot"])[:6]:
        print(f'  {c["p"]} {c["pam"]} {c["st"]} corte {c["cut"]:5d} cod {c["cod"]:4d} GC{c["gc"]} D14={c["d14"]:3d} CS={c["cs"]:3d} hp={c["hp"]} shp={c["shp"]}')
