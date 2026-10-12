"""Validacion de la logica sin red: se construyen un gen sintetico con intrones
y un vector sintetico con dos sitios BsaI, y se comprueba que el mapeo exonico,
la verificacion de protospacers y el ensamblaje Golden Gate dan el resultado
exacto esperado."""
import random, sys
from lib_edicion import (mapear_exones, protospacer_en_exon, sitios_bsai,
                         digerir_vector, ensamblar, buscar_offtarget, rc)

random.seed(7)
az = lambda n: "".join(random.choice("ACGT") for _ in range(n))
fallos = 0

def check(nombre, cond, detalle=""):
    global fallos
    print(("PASA  " if cond else "FALLA ") + nombre + (" :: " + detalle if detalle else ""))
    if not cond: fallos += 1

# --- 1. gen sintetico: 3 exones (300, 240, 180) y 2 intrones GT..AG ---------
ex = [az(300), az(240), az(180)]
intr = ["GT" + az(96) + "AG", "GT" + az(146) + "AG"]
genomico = az(500) + ex[0] + intr[0] + ex[1] + intr[1] + ex[2] + az(400)
cds = "".join(ex)
bloques = mapear_exones(cds, genomico)
check("mapeo exonico: 3 bloques", len(bloques) == 3, str(len(bloques)))
check("mapeo exonico: longitudes 300/240/180",
      [b[1]-b[0] for b in bloques] == [300, 240, 180], str([b[1]-b[0] for b in bloques]))
check("mapeo exonico: inicio genomico correcto", bloques[0][2] == 500, str(bloques[0][2]))

# --- 2. protospacer interno a un exon vs protospacer que cruza el intron ----
ok_int, coord = protospacer_en_exon(bloques, 100, 122)      # dentro del exon 1
check("protospacer interno se acepta", ok_int and coord == (600, 622), str(coord))
ok_cruz, _ = protospacer_en_exon(bloques, 290, 312)         # cruza exon1/exon2
check("protospacer que cruza un intron se rechaza", not ok_cruz)
ok_lim, _ = protospacer_en_exon(bloques, 279, 301)          # 23 nt justo al borde
check("protospacer que se pasa 1 nt del exon se rechaza", not ok_lim)

# --- 3. vector sintetico con dos sitios BsaI enfrentados --------------------
relleno = az(200)
# los dos sitios apuntan hacia afuera del relleno: asi se eliminan con el
# relleno y el plasmido ensamblado ya no es sustrato de BsaI
vector = (az(1000) + "ATTG" + "A" + "GAGACC" + relleno + "GGTCTC" + "A" + "CTGA" + az(800))
sitios = sitios_bsai(vector)
check("se detectan 2 sitios BsaI", len(sitios) == 2, str([s[1:] for s in sitios]))
esq, oh_i, oh_d = digerir_vector(vector)
check("salientes del vector leidas", (oh_i, oh_d) == ("ATTG", "CTGA"), f"{oh_i}/{oh_d}")
check("el relleno desaparece del esqueleto", relleno not in esq)

inserto = "ATTG" + az(1200) + "CTGA"
plasmido = ensamblar(vector, inserto)
check("ensamblaje: longitud esperada",
      len(plasmido) == len(esq) + len(inserto), f"{len(plasmido)}")
check("ensamblaje: el inserto queda integro", inserto[4:-4] in plasmido)
check("ensamblaje: no quedan sitios BsaI", len(sitios_bsai(plasmido)) == 0,
      str(sitios_bsai(plasmido)))

# --- 4. off-target: sitio exacto, 2 desapareamientos y uno sin PAM ----------
guia = "GTCCCGGATATGTTGAGGGG"
diana = az(300) + guia + "AGG" + az(300)
variante = guia[:5] + ("A" if guia[5] != "A" else "C") + guia[6:14] + ("T" if guia[14] != "T" else "G") + guia[15:]
diana += variante + "TGG" + az(200) + guia + "AAA" + az(100)   # ultimo sin PAM
hits = buscar_offtarget(guia, diana, max_mm=3)
mms = sorted(h[4] for h in hits)
check("off-target: encuentra el exacto y el de 2 mm", mms[:2] == [0, 2], str(mms))
check("off-target: ignora el sitio sin PAM NGG", len(hits) == 2, str(len(hits)))
check("off-target: revisa tambien la hebra complementaria",
      len(buscar_offtarget(guia, rc(diana), max_mm=3)) == 2)

print()
print("FALLOS:", fallos)
sys.exit(1 if fallos else 0)
