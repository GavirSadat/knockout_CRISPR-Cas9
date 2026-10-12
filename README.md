
# Knockout CRISPR-Cas9 de la ruta de momilactonas en *Oryza sativa*
**Taller 3 — Edición Genética | Maestría en Biociencias, Universidad EAFIT**  
**Autores:** Gavir Sadat Carmona López, Juan Diego Rojas Saldarriaga, Yair Antonio Vega

## 📌 Resumen del Proyecto
Diseño *in silico* reproducible y verificable de una estrategia CRISPR-Cas9 para generar alelos nulos en los dos primeros pasos de la biosíntesis de momilactonas en arroz cultivar japónica: **OsCPS4** (LOC_Os04g09900) y **OsKSL4** (LOC_Os04g10060). 

La estrategia emplea un enfoque multicapa (tres sgRNA por gen) que garantiza:
1. Corrimiento de marco temprano.
2. Deleción intraexónica diagnosticable por PCR.
3. Destrucción del sitio activo catalítico (motivos DDFFD y DIDD).

El sistema utiliza **SpCas9 optimizada para codones de planta**, flanqueada por secuencias de localización nuclear (NLS) bipartitas de **VirD2** con sesgo de codones para arroz. El ensamblaje de los tres casetes de sgRNA se realiza mediante Golden Gate (BsaI) sobre el vector binario `pYLCRISPR/Cas9Pubi-H`.

## 📂 Estructura del Repositorio

* `benchling/`: Archivos `.gb` (GenBank) con las secuencias plasmídicas completas, los mapas del T-DNA anotados, insertos y archivos `.csv` con los oligos de clonación y puntajes de eficiencia. Listos para importar en Benchling.
* `scripts/`: Código fuente en Python (Biopython) que genera las secuencias, mapea intrones/exones, ensambla los plásmidos *in silico*, evalúa el *off-target* (vía bowtie) y ejecuta el control de calidad.
* `MANIFIESTO.md`: Inventario detallado del paquete, estado de las verificaciones computacionales y tabla definitiva de las 6 guías (sgRNAs) seleccionadas.

*(Nota: Los reportes técnicos en PDF se añadirán en la versión final de entrega).*

## ⚙️ Requisitos y Dependencias
Para ejecutar los scripts de diseño y verificación computacional desde cero, el entorno debe contar con:
* **Python 3.x**
* **Biopython** (`pip install biopython`)
* **Bowtie 1.3.x** (Requerido por `03_offtarget.py` para el rastreo genómico)
* Conexión a internet (para las consultas a NCBI y Ensembl Plants en los pasos `01` y `13`).

## 🚀 Flujo de Trabajo y Orden de Ejecución
Toda la lógica computacional es reproducible. Para regenerar los plásmidos y las verificaciones desde cero, ejecuta los scripts de la carpeta `scripts/` en el siguiente orden estricto:

1. `01_verificar_ncbi.py`: Descarga genómica desde NCBI, mapeo de exones y validación intraexónica de guías.
2. `03_offtarget.py`: Cribado exhaustivo a escala de genoma vía Bowtie (filtro PAM NGG).
3. `04` al `08`: Generación de archivos GenBank base, arreglos, T-DNA y reescaneo intraexónico.
4. `09_insertar_promotores.py`: Inserción de promotores reales OsU6 (a, b, c) en el arreglo.
5. `02_completar_plasmido.py`: Ensamblaje Golden Gate del plásmido intermedio.
6. `12_sustituir_nls_virD2.py`: Sustitución de señales SV40/Nucleoplasmina por la NLS VirD2 optimizada para arroz. **(Genera los plásmidos finales)**.
7. `10_qc_plasmido.py`: Control de calidad riguroso (12 puntos) sobre el plásmido ensamblado.
8. `11_disenar_cebadores.py`: Diseño termodinámico de cebadores para genotipado.
9. `13_coordenadas_ensembl.py`: Reconstrucción y mapeo de coordenadas CDS sobre el genoma IRGSP-1.0.
