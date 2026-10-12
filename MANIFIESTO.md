# Taller 3 — Knockout CRISPR-Cas9 de momilactonas en *Oryza sativa*
## Manifiesto del paquete 

Diseño dirigido a *OsKSL4* (LOC_Os04g10060) y *OsCPS4* (LOC_Os04g09900), tres sgRNA
por gen, entregados en dos binarios independientes construidos sobre
pYLCRISPR/Cas9Pubi-H.

## Lo que contiene este paquete

### `benchling/` — secuencias del diseño

| Archivo | Qué es |
| --- | --- |
| `OsKSL4_CDS_sgRNAs.gb` | CDS de 2.529 pb con el motivo DDFFD, los tres protospacers, sus PAM y los cortes |
| `OsCPS4_CDS_sgRNAs.gb` | CDS de 2.304 pb con el motivo DIDD y las tres guías |
| `sgRNA_array_ksl4.gb` · `sgRNA_array_cps4.gb` | Inserto de tres casetes; los promotores OsU6 son tramos de N que rellena el paso 09 |
| `pMOMI-ksl4_VirD2.gb` · `pMOMI-cps4_VirD2.gb` | **Los plásmidos del reporte**, 17.714 y 17.713 pb, con la NLS bipartita de VirD2 en ambos extremos de la Cas9 |
| `pMOMI-ksl4_TDNA_diseno.gb` · `pMOMI-cps4_TDNA_diseno.gb` | Mapa del T-DNA con las fusiones NLS; versión de diseño, con tramos de N |
| `sgRNAs_y_oligos.csv` | Las seis guías con coordenadas, %GC, sgRNA completos y oligos de clonación |
| `chopchop_scores.csv` | Contexto de 30 nt, Doench 2014, CRISPRscan, autocomplementariedad |

### `scripts/` — código, reproducible

| Archivo | Función |
| --- | --- |
| `lib_edicion.py` | Mapeo exón-intrón con regla GT-AG, verificación exónica, digestión BsaI, Golden Gate, off-target |
| `tests_offline.py` | 15 pruebas de esa lógica sobre secuencias sintéticas; todas pasan |
| `01_verificar_ncbi.py` | Descarga de NCBI, identidad de CDS, exones, verificación de las guías |
| `09_insertar_promotores.py` | Promotores OsU6 reales desde los pYLsgRNA |
| `02_completar_plasmido.py` | Ensamblaje del plásmido final con traslado de anotaciones |
| `10_qc_plasmido.py` | Control de calidad de doce puntos sobre el plásmido |
| `12_sustituir_nls_virD2.py` | Sustituye SV40 y nucleoplasmina por la NLS de VirD2 con codones de arroz |
| `03_offtarget.py` | Off-target sobre FASTA o genoma completo vía bowtie, con filtro de PAM |
| `11_disenar_cebadores.py` | Cebadores de genotipado sobre el locus genómico |
| `04`–`08` | Generadores de los GenBank, del arreglo y del T-DNA, scores y reescaneo de candidatos |
| `seqs.py` | Las dos CDS, insumo de `04`–`08` |
| `candidatos_ksl4_sg2.fa` | Los diez candidatos cribados contra el genoma para reemplazar la guía descartada |
| `offtarget_dianas.tsv` | Las seis guías contra las dos CDS diana |
| `README.md` | Orden de ejecución y dependencias |

## Lo que no está versionado y se regenera localmente

Estos archivos los generaron las corridas con NCBI y con el genoma, que no se
pueden ejecutar sin red. Se regeneran con los pasos `01`, `09`, `02`, `12`, `03`
y `11`:

| Archivo | Lo produce |
| --- | --- |
| `pMOMI-ksl4_completo.gb` · `pMOMI-cps4_completo.gb` | `02` — el ensamblaje intermedio, 17.570 y 17.569 pb, todavía con las NLS del vector parental |
| `sgRNA_array_ksl4_completo.gb` · `sgRNA_array_cps4_completo.gb` | `09` — inserto con los promotores reales |
| `OsKSL4_locus_genomico.gb` · `OsCPS4_locus_genomico.gb` | `01` — loci de NC_089038.1 con 14 exones y las guías |
| `verificacion_exones.tsv` | `01` — verificación exónica guía por guía |
| `offtarget.tsv` | `03` — cribado genómico con bowtie |
| `cebadores_genotipado.csv` | `11` — los cuatro pares de la Tabla 7 |

**Ojo con la NLS:** los `_completo.gb` conservan la Cas9p del vector parental,
con SV40 en el extremo N y nucleoplasmina en el C. Los plásmidos finales que
describe el reporte son los `_VirD2.gb` (17.714 y 17.713 pb), que produce `12`
a partir de aquellos y que **ya están versionados** en `benchling/`. Para reconstruirlos desde cero basta seguir el orden del README de la raíz.

## Las seis guías finales

| sgRNA | Espaciador 5'→3' | PAM | Corte (nt CDS) | Codón | Exón |
| --- | --- | --- | --- | --- | --- |
| KSL4-sg1 | GCGACAATCTAGCACGAATG | GGG | 103 | 35 | 1 |
| KSL4-sg2 | GCCGGCAGCTTCCTTCCACG | GGG | 207 | 69 | 1 |
| KSL4-sg3 | GTCCGACATCGAAAAAGTCAT | CGG | 1774 | 592 | 11 |
| CPS4-sg1 | GCACTGAGATCAATCGGAGA | GGG | 284 | 95 | 3 |
| CPS4-sg2 | GATCGTTCAGAATCAGCTAC | CGG | 403 | 135 | 3 |
| CPS4-sg3 | GGTAAAGGATATCGATGACA | CGG | 1102 | 368 | 7 |

Las tres verificaciones ejecutadas: seis de seis intraexónicas con su PAM sobre
NC_089038.1; una coincidencia exacta por guía en el genoma completo y ningún
sitio a uno o dos desapareamientos; y los dos plásmidos con cero sitios BsaI y
doce de doce controles superados.
