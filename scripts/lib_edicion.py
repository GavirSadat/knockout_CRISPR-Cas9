"""Funciones nucleo del diseno: mapeo exon-intron, verificacion de protospacers
y ensamblaje Golden Gate in silico. Sin dependencias de red: se validan con
tests_offline.py sobre secuencias sinteticas."""
from __future__ import annotations

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def rc(s: str) -> str:
    return s.translate(COMP)[::-1]

# ---------------------------------------------------------------- exones ----
def mapear_exones(cds: str, genomico: str, min_ancla: int = 25,
                  intron_max: int = 50000):
    """Alinea una CDS contra su region genomica y devuelve los bloques exonicos
    como [(ini_cds, fin_cds, ini_gen, fin_gen)] en coordenadas 0-based
    semiabiertas. Ajusta cada frontera para que el intron empiece en GT y
    termine en AG cuando sea posible (regla GT-AG)."""
    cds, genomico = cds.upper(), genomico.upper()
    inicio = genomico.find(cds[:min_ancla])
    if inicio < 0:
        raise ValueError("El inicio de la CDS no aparece en la region genomica")
    bloques, ci, gi = [], 0, inicio
    while ci < len(cds):
        n = 0
        while ci + n < len(cds) and gi + n < len(genomico) and cds[ci + n] == genomico[gi + n]:
            n += 1
        if ci + n >= len(cds):                      # ultimo exon
            bloques.append((ci, len(cds), gi, gi + (len(cds) - ci)))
            break
        k = min(min_ancla, len(cds) - (ci + n))
        siguiente = genomico.find(cds[ci + n:ci + n + k], gi + n, gi + n + intron_max)
        if siguiente < 0:
            raise ValueError(f"No se encontro el exon siguiente tras el nt {ci+n} de la CDS")
        # correccion GT-AG: retroceder hasta 6 nt para respetar el sitio de corte
        ajuste = 0
        for d in range(0, 7):
            if (genomico[gi + n - d:gi + n - d + 2] == "GT"
                    and genomico[siguiente - d - 2:siguiente - d] == "AG"):
                ajuste = d
                break
        n -= ajuste
        siguiente -= ajuste
        bloques.append((ci, ci + n, gi, gi + n))
        ci, gi = ci + n, siguiente
    return bloques

def protospacer_en_exon(bloques, ini_cds: int, fin_cds: int):
    """ini/fin en coordenadas 1-based de la CDS, incluyendo el PAM.
    Devuelve (True, (ini_gen, fin_gen)) si el tramo cae entero en un solo exon."""
    a, b = ini_cds - 1, fin_cds
    for c0, c1, g0, g1 in bloques:
        if c0 <= a and b <= c1:
            return True, (g0 + (a - c0) + 1, g0 + (b - c0))
    return False, None

# ------------------------------------------------------------ Golden Gate ---
def sitios_bsai(seq: str, circular: bool = True):
    """Devuelve [(posicion_corte_hebra_superior, saliente, orientacion)] para
    cada sitio BsaI. GGTCTC(N1) corta dejando 4 nt salientes 5'."""
    s = seq.upper()
    ext = s + s[:10] if circular else s
    sitios = []
    i = ext.find("GGTCTC")
    while i != -1 and i < len(s):
        sitios.append((i + 7, ext[i + 7:i + 11], "+"))
        i = ext.find("GGTCTC", i + 1)
    j = ext.find("GAGACC")
    while j != -1 and j < len(s):
        sitios.append((j - 5, ext[j - 5:j - 1], "-"))
        j = ext.find("GAGACC", j + 1)
    return sorted(sitios)

def digerir_vector(vector: str):
    """Vector circular con exactamente dos sitios BsaI en orientacion opuesta.
    Devuelve (esqueleto, saliente_5p_inserto, saliente_3p_inserto): tras el
    esqueleto va la saliente 5' del inserto y antes de el la 3'."""
    sitios = sitios_bsai(vector)
    if len(sitios) != 2:
        raise ValueError(f"Se esperaban 2 sitios BsaI y hay {len(sitios)}")
    (p1, oh1, o1), (p2, oh2, o2) = sitios
    if o1 == o2:
        raise ValueError("Los dos sitios BsaI deben apuntar en sentidos opuestos")
    # el relleno es lo que queda entre los dos cortes; el esqueleto es el resto
    esqueleto = vector[p2 + 4:] + vector[:p1]
    return esqueleto, oh1, oh2

def ensamblar(vector: str, inserto: str):
    """Ensambla el inserto (con sus 4 nt salientes en cada extremo) en el
    vector digerido con BsaI y devuelve el plasmido circular resultante."""
    esqueleto, oh_izq, oh_der = digerir_vector(vector)
    if inserto[:4].upper() not in (oh_izq, "NNNN"):
        raise ValueError(f"Saliente 5' del inserto {inserto[:4]} != {oh_izq} del vector")
    if inserto[-4:].upper() not in (oh_der, "NNNN"):
        raise ValueError(f"Saliente 3' del inserto {inserto[-4:]} != {oh_der} del vector")
    cuerpo = oh_izq + inserto[4:-4] + oh_der
    return esqueleto + cuerpo

# ------------------------------------------------------------- off-target ---
def desapareamientos(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a.upper(), b.upper()))

def buscar_offtarget(guia: str, secuencia: str, max_mm: int = 3, pam: str = "NGG"):
    """Busca sitios con PAM NGG y hasta max_mm desapareamientos frente a la guia.
    Devuelve [(posicion_1based, hebra, sitio, pam, n_mm, mm_en_semilla)]."""
    g = guia.upper()
    L = len(g)
    salidas = []
    for hebra, s in (("+", secuencia.upper()), ("-", rc(secuencia.upper()))):
        for i in range(len(s) - L - 2):
            p = s[i + L:i + L + 3]
            if p[1:] != "GG":
                continue
            sitio = s[i:i + L]
            mm = desapareamientos(g, sitio)
            if mm <= max_mm:
                salidas.append((i + 1, hebra, sitio, p, mm,
                                desapareamientos(g[-12:], sitio[-12:])))
    return salidas
