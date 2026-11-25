#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_topo.py — schreibt template/site/topo-data.js (Topoisomerase I mit DNA für replikation.html).

Quelle: RCSB PDB, Eintrag 1A36 (menschliche Topoisomerase I, Kern- und C-terminale Domäne, Reste 215–765, mit einem
22-bp-DNA-Doppelstrang, 2,8 Å; Stewart, Redinbo, Qiu, Hol und Champoux 1998, Science 279: 1534–1541).
Abruf 25.11.2026 von files.rcsb.org/download/1A36.pdb (Datei in results/replikation/).
Im Kristall ist das Tyrosin 723 des aktiven Zentrums gegen Phenylalanin ausgetauscht (Y723F), damit das Enzym die DNA
nicht schneidet. Die Reste 634–640 sind im Kristall nicht sichtbar; die Kette wird dort in zwei Stücke geteilt.

Enthalten sind je Rest Cα, C und O mit Sekundärstruktur aus den HELIX/SHEET-Zeilen und der Domäne nach Redinbo et al.
1998 (Science 279: 1504–1513), für die DNA je Nukleotid P und C1' und die Base.

Aufruf:  python3 scripts/replikation/export_topo.py
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PDB = os.path.join(ROOT, "results", "replikation", "1A36.pdb")
OUT = os.path.join(ROOT, "template", "site", "topo-data.js")
DREI = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
        "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def domaene(nr):
    if 233 <= nr <= 319: return "kern2"      # Kern-Subdomäne II (Kappe)
    if nr <= 433: return "kern1"             # Kern-Subdomäne I (Kappe)
    if nr <= 635: return "kern3"             # Kern-Subdomäne III (Basis, umgreift die DNA)
    if nr <= 712: return "linker"            # Linker
    return "cterm"                           # C-terminale Domäne mit dem aktiven Tyrosin 723


reste, ss, dna = {}, {}, {}
for z in open(PDB):
    if z.startswith("HELIX "):
        a, b = int(z[21:25]), int(z[33:37])
        for k in range(a, b + 1): ss[k] = "H"
    elif z.startswith("SHEET "):
        a, b = int(z[22:26]), int(z[33:37])
        for k in range(a, b + 1): ss.setdefault(k, "E")
    if not z.startswith(("ATOM", "HETATM")) or z[16] not in " A": continue
    ch, name, rn, nr = z[21], z[12:16].strip(), z[17:20].strip(), int(z[22:26])
    p = [round(float(z[30:38]), 1), round(float(z[38:46]), 1), round(float(z[46:54]), 1)]
    if ch == "A" and rn in DREI and name in ("CA", "C", "O"):
        reste.setdefault(nr, {"nr": nr, "aa": DREI[rn]})[name] = p
    if ch in "BC" and name in ("P", "C1'"):
        dna.setdefault(ch, {}).setdefault(nr, {"nr": nr, "b": rn.strip()[-1]})[name] = p

# Mittelpunkt für die Ausrichtung: Schwerpunkt der C1'-Atome (DNA in der Mitte des Bildes)
c1 = [n["C1'"] for k in dna.values() for n in k.values() if "C1'" in n]
m = [round(sum(p[i] for p in c1) / len(c1), 1) for i in range(3)]
sh = lambda p: [round(p[i] - m[i], 1) for i in range(3)]

stuecke, akt = [], []
for nr in sorted(reste):
    r = reste[nr]
    if not all(k in r for k in ("CA", "C", "O")): continue
    if akt and nr != akt[-1]["nr"] + 1: stuecke.append(akt); akt = []
    akt.append({"nr": nr, "aa": r["aa"], "ss": ss.get(nr, "C"), "d": domaene(nr), "ca": sh(r["CA"]), "c": sh(r["C"]), "o": sh(r["O"])})
stuecke.append(akt)
straenge = []
for ch in "BC":
    nts = [n for _, n in sorted(dna[ch].items())]
    straenge.append([{"nr": n["nr"], "b": n["b"], "p": sh(n["P"]) if "P" in n else None, "c1": sh(n["C1'"])} for n in nts])
for s in stuecke: print("Stück", s[0]["nr"], "bis", s[-1]["nr"], len(s), "Reste")
for s in straenge: print("Strang", "".join(n["b"] for n in s), len(s))
tyr = next(r for s in stuecke for r in s if r["nr"] == 723)
print("Rest 723:", tyr["aa"], tyr["ca"])

kopf = ("// Topoisomerase I mit DNA (PDB 1A36; Stewart et al. 1998, Science 279: 1534). Erzeugt von scripts/replikation/export_topo.py,\n"
        "// nicht von Hand ändern. TOPO.stuecke: Proteinkette in zusammenhängenden Stücken, je Rest {nr, aa, ss, d (Domäne), ca, c, o};\n"
        "// TOPO.dna: zwei Stränge, je Nukleotid {nr, b (Base), p (Phosphor, fehlt am 5'-Ende), c1 (C1')}. Koordinaten in Å, um die DNA zentriert.\n")
txt = kopf + "const TOPO = " + json.dumps({"stuecke": stuecke, "dna": straenge}, separators=(",", ":")) + ";\n"
open(OUT, "w").write(txt)
print("geschrieben:", OUT, len(txt) // 1024, "KB")
