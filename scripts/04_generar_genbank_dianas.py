import os as _os
_AQUI = _os.path.dirname(_os.path.abspath(__file__))
_RAIZ = _os.path.dirname(_AQUI)
OUT   = _os.environ.get("MOMI_OUT", _os.path.join(_RAIZ, "benchling")) + _os.sep
_SEQS = _os.path.join(_AQUI, "seqs.py")
_EFF  = _os.environ.get("CRISPOR_EFF", _os.path.join(_RAIZ, "externo", "crisporEffScores.py"))
exec(open(_SEQS).read().split('import re')[0])
import re, csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]

SCAFFOLD="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
TERM="TTTTTTT"
BSA="GGTCTC"

guides={
 "KSL4":[("KSL4-sg1","GCGACAATCTAGCACGAATG","GGG","-"),
         ("KSL4-sg2","GCCGGCAGCTTCCTTCCACG","GGG","+"),
         ("KSL4-sg3","TCCGACATCGAAAAAGTCAT","CGG","-")],
 "CPS4":[("CPS4-sg1","GCACTGAGATCAATCGGAGA","GGG","+"),
         ("CPS4-sg2","GATCGTTCAGAATCAGCTAC","CGG","+"),
         ("CPS4-sg3","GGTAAAGGATATCGATGACA","CGG","+")]}
cds={"KSL4":ksl4,"CPS4":cps4}
info={"KSL4":dict(locus="LOC_Os04g10060",rap="Os04g0179700",gene="OsKSL4",
        prod="syn-pimara-7,15-diene synthase (OsDTS2)",acc="AB126934 / AY616862",motif=(1771,1785,"DDxxD (DDFFD)")),
      "CPS4":dict(locus="LOC_Os04g09900",rap="Os04g0178300",gene="OsCPS4",
        prod="syn-copalyl diphosphate synthase (OsCyc1)",acc="AB066270 / AK100631",motif=(1093,1104,"DxDD (DIDD)"))}

def locate(seq,proto,pam,strand):
    """devuelve (inicio1,fin1) del protospacer y posicion del corte en coords 1-based del CDS"""
    if strand=="+":
        i=seq.find(proto+pam)
        assert i>=0, proto
        return i+1, i+len(proto), i+len(proto)-3   # corte 3 nt antes del PAM
    else:
        t=rc(proto+pam); i=seq.find(t)
        assert i>=0, proto
        start=i+len(pam)+1; end=i+len(t)
        return start, end, start+2                 # corte 3 nt desde el PAM en la hebra -

rows=[]
for g,seq in cds.items():
    meta=info[g]
    rec=SeqRecord(Seq(seq), id=meta["locus"], name=meta["gene"],
        description=f'{meta["gene"]} CDS (Oryza sativa Japonica, IRGSP-1.0/MSU7) con sitios diana CRISPR-Cas9')
    rec.annotations.update(molecule_type="DNA", topology="linear", organism="Oryza sativa Japonica Group",
        source="Oryza sativa Japonica Group", date="02-OCT-2026",
        comment=f'Referencia: {meta["rap"]} / {meta["locus"]}; GenBank {meta["acc"]}. '
                'Coordenadas en CDS. Verificar en el locus genomico que cada protospacer caiga dentro de un unico exon.')
    rec.features.append(SeqFeature(FeatureLocation(0,len(seq)),type="CDS",
        qualifiers={"gene":meta["gene"],"product":meta["prod"],"codon_start":"1","transl_table":"1","label":meta["gene"]+" CDS"}))
    m0,m1,mn=meta["motif"]
    rec.features.append(SeqFeature(FeatureLocation(m0-1,m1),type="misc_feature",
        qualifiers={"label":"motivo catalitico "+mn,"note":"sitio activo; una indel aqui inactiva la enzima aunque sea en marco"}))
    for name,proto,pam,strand in guides[g]:
        s,e,cut=locate(seq,proto,pam,strand)
        st=1 if strand=="+" else -1
        rec.features.append(SeqFeature(FeatureLocation(s-1,e,strand=st),type="misc_binding",
            qualifiers={"label":name+" protospacer","note":f"{proto} ({len(proto)} nt), PAM {pam}, hebra {strand}"}))
        pam_s = e+1 if strand=="+" else s-3
        rec.features.append(SeqFeature(FeatureLocation(pam_s-1,pam_s+2,strand=st),type="misc_feature",
            qualifiers={"label":name+" PAM","note":pam}))
        rec.features.append(SeqFeature(FeatureLocation(cut-1,cut,strand=st),type="misc_feature",
            qualifiers={"label":name+" corte Cas9","note":f"corte romo tras nt {cut} del CDS (codon {(cut+2)//3})"}))
        rows.append(dict(gen=meta["gene"],sgRNA=name,protospacer=proto,PAM=pam,hebra=strand,
            inicio_CDS=s,fin_CDS=e,corte_CDS=cut,codon=(cut+2)//3,
            GC=round(100*(proto.count("G")+proto.count("C"))/len(proto)),
            espaciador_transcrito=("G"+proto if proto[0]!="G" else proto),
            oligo_F="GGCA"+("G"+proto if proto[0]!="G" else proto),
            oligo_R="AAAC"+rc("G"+proto if proto[0]!="G" else proto),
            sgRNA_completo=("G"+proto if proto[0]!="G" else proto)+SCAFFOLD))
    cuts=sorted(r["corte_CDS"] for r in rows if r["gen"]==meta["gene"])
    rec.features.append(SeqFeature(FeatureLocation(cuts[0],cuts[-1]),type="misc_feature",
        qualifiers={"label":"delecion maxima esperada","note":f"{cuts[-1]-cuts[0]} pb entre sg1 y sg3"}))
    SeqIO.write(rec, OUT+f'{meta["gene"]}_CDS_sgRNAs.gb',"genbank")
    print(meta["gene"],"cortes:",cuts,"delec sg1-sg2:",cuts[1]-cuts[0],"sg2-sg3:",cuts[2]-cuts[1])

with open(OUT+"sgRNAs_y_oligos.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("filas:",len(rows))
