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
# exones verificados sobre NC_089038.1 (coordenadas 1-based de la CDS)
EX_KSL4=[(1,241),(242,291),(292,574),(575,755),(756,976),(977,1097),(1098,1195),
         (1196,1410),(1411,1523),(1524,1619),(1620,1835),(1836,1977),(1978,2223),(2224,2529)]
EX_CPS4=[(1,133),(134,219),(220,523),(524,707),(708,922),(923,1043),(1044,1156),
         (1157,1389),(1390,1505),(1506,1604),(1605,1784),(1785,1872),(1873,2082),(2083,2304)]
def exon_de(exones,a,b,margen=0):
    for i,(e0,e1) in enumerate(exones,1):
        if e0<=a-margen and b+margen<=e1: return i
    return None
def maxhp(a,b,minlen=4):
    for L in range(len(a),minlen-1,-1):
        for i in range(len(a)-L+1):
            if rc(a[i:i+L]) in b: return L
    return 0
def scan(seq,exones,lo,hi):
    out=[];L=len(seq)
    for strand in "+-":
        s = seq if strand=="+" else rc(seq)
        for m in re.finditer("(?=([ACGT]{21}GG))", s):
            i=m.start(); p=m.group(1)[:20]; pam=m.group(1)[20:23]
            if i<6 or i+29>len(s): continue
            # coordenadas 1-based del 23-mero completo en la CDS
            if strand=="+": a23,b23=i+1,i+23
            else: a23,b23=L-(i+22),L-i
            cut = i+17 if strand=="+" else L-(i+17)
            if not lo<=cut<=hi: continue
            ex = exon_de(exones,a23,b23,margen=6)   # 23-mero + contexto de los scores
            if ex is None: continue
            if p[0]!="G": continue
            gc=100*(p.count("G")+p.count("C"))/20
            if not 40<=gc<=70: continue
            if "TTTT" in p: continue
            if re.search(r"(A{5}|C{5}|G{5}|T{5})",p): continue
            if any(b in p or b in rc(p) for b in BAD): continue
            if sum(x.count(p[-12:])+x.count(rc(p[-12:])) for x in (ksl4,cps4))>1: continue
            d14=ns['calcDoenchScores']([s[i-4:i+26]])[0]
            cs =ns['calcCrisprScanScores']([s[i-6:i+29]])[0]
            out.append(dict(p=p,pam=pam,st=strand,cut=cut,a=a23,b=b23,ex=ex,gc=int(gc),
                            d14=d14,cs=cs,hp=maxhp(p,p),shp=maxhp(p,SCAF[:24]),tot=d14+cs,
                            cod=(cut+2)//3))
    return out
print("### OsKSL4, candidatos intraexonicos en el primer tercio (corte 60-800)")
for c in sorted(scan(ksl4,EX_KSL4,60,800),key=lambda d:-d["tot"])[:12]:
    print(f'  {c["p"]} {c["pam"]} {c["st"]} corte {c["cut"]:4d} cod {c["cod"]:4d} exon {c["ex"]:2d} '
          f'({c["a"]}-{c["b"]}) GC{c["gc"]} D14={c["d14"]:3d} CS={c["cs"]:3d} hp={c["hp"]} shp={c["shp"]}')
print("\n### control: las seis guias actuales, con su exon")
for n,p,pam,st,seq,ex in [("KSL4-sg1","GTCCCGGATATGTTGAGGGG","AGG","+",ksl4,EX_KSL4),
                          ("KSL4-sg2","GGTTCGGAGCAAATAAGGAG","AGG","+",ksl4,EX_KSL4),
                          ("KSL4-sg3","TCCGACATCGAAAAAGTCAT","CGG","-",ksl4,EX_KSL4),
                          ("CPS4-sg1","GCACTGAGATCAATCGGAGA","GGG","+",cps4,EX_CPS4),
                          ("CPS4-sg2","GATCGTTCAGAATCAGCTAC","CGG","+",cps4,EX_CPS4),
                          ("CPS4-sg3","GGTAAAGGATATCGATGACA","CGG","+",cps4,EX_CPS4)]:
    if st=="+": i=seq.find(p+pam); a,b=i+1,i+23
    else: i=seq.find(rc(p+pam)); a,b=i+1,i+23
    print(f'  {n}: 23-mero {a}-{b} -> exon {exon_de(ex,a,b)} (con margen 6: {exon_de(ex,a,b,6)})')
