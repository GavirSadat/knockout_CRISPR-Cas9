# `scripts/` — qué hace cada archivo

El flujo de trabajo completo, con el orden de ejecución y las variables de
entorno, está en el [README de la raíz](../README.md). Esta tabla es solo el
inventario.

| Archivo | Qué hace |
| --- | --- |
| `lib_edicion.py` | Biblioteca: mapeo exón-intrón con regla GT-AG, verificación exónica de protospacers, digestión BsaI, ensamblaje Golden Gate y búsqueda de off-target |
| `seqs.py` | Las dos CDS de referencia, como texto |
| `tests_offline.py` | 15 pruebas sobre un gen sintético de 3 exones y un vector sintético con dos sitios BsaI enfrentados. **Es el primer paso de cualquier ejecución.** |
| `01_verificar_ncbi.py` | Descarga de NCBI, identidad de las CDS, mapeo de exones y verificación exónica de las seis guías. Sale con código 1 si alguna cruza un intrón |
| `02_completar_plasmido.py` | Digiere el vector real con BsaI y ensambla el plásmido completo, sin tramos de N |
| `03_offtarget.py` | Off-target exhaustivo sobre FASTA, o genoma completo vía bowtie con filtro posterior de PAM |
| `04_generar_genbank_dianas.py` | GenBank de las CDS diana con guías, PAM y cortes anotados, más `sgRNAs_y_oligos.csv` |
| `05_generar_array_y_tdna.py` | Arreglo de sgRNA, fusiones de NLS y mapa lineal del T-DNA |
| `06_scores_eficiencia.py` | Doench 2014 y CRISPRscan de las seis guías (requiere `CRISPOR_EFF`) |
| `07_rescan_y_ranking.py` | Reescaneo de todos los candidatos con los filtros y ranking |
| `08_rescan_intraexonico.py` | El reescaneo que descartó dos candidatos por cruzar frontera de exón; conserva las guías superadas como registro de por qué se cambiaron |
| `09_insertar_promotores.py` | Sustituye los marcadores de posición de los promotores OsU6 por su secuencia real de los pYLsgRNA |
| `10_qc_plasmido.py` | Doce controles sobre el plásmido circular final |
| `11_disenar_cebadores.py` | Cebadores de genotipado sobre el locus genómico (nunca sobre la CDS) |
| `13_coordenadas_ensembl.py` | Coordenadas genómicas sobre IRGSP-1.0 por la API de Ensembl Plants, verificación de bloque único por 23-mero y de las cuatro distancias de deleción; `--autotest` lo valida sin red |
| `12_sustituir_nls_virD2.py` | Sustituye las NLS del vector parental (SV40 en N, nucleoplasmina en C) por la bipartita de VirD2 con codones de arroz |
| `candidatos_ksl4_sg2.fa` | Los candidatos evaluados para reemplazar la KSL4-sg2 descartada |
| `offtarget_dianas.tsv` | Resultado ya ejecutado: las seis guías contra las dos CDS diana |

`04`, `05` y `06` son deterministas: regeneran byte a byte los archivos
versionados en `benchling/`.

`04`, `06`, `07` y `08` leen `seqs.py`; `06`, `07` y `08` necesitan además
`crisporEffScores.py` (ver [`../externo/README.md`](../externo/README.md)).
