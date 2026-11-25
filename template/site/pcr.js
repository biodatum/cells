// Das Modell der PCR: Einzelstränge, Primer, Zyklen, Zählwerk
// (Quelle: template/site/pcr.js, wird von build.py in die Seiten kopiert; Test: scripts/pcr/test_pcr.js).
//
// Die Ausgangs-DNA hat M = ziel + 2·rand Basenpaare, Spalten c = 0 … M-1. Die obere Kette ("top") steht von links nach
// rechts 5'→3', die untere ("bottom") 3'→5'. Die Zielsequenz liegt in den Spalten a … b-1 (a = rand, b = rand + ziel).
//   Primer 1  gleich der oberen Kette bei a … a+primer-1; er lagert sich an die UNTERE Kette an und wird nach rechts verlängert
//             (neue obere Kette: 5'-Ende bei a, wächst 5'→3' bis zum Ende der Vorlage).
//   Primer 2  komplementär zur oberen Kette bei b-primer … b-1; er lagert sich an die OBERE Kette an und wird nach links verlängert
//             (neue untere Kette: 5'-Ende bei b-1, wächst 5'→3' nach links bis zum Ende der Vorlage).
// Ein Strang ist {art: "top" | "bottom", von, bis}: er deckt die Spalten von … bis-1 ab.
//
//   PCR.erzeuge(zufall, opt)  zufällige Ausgangs-DNA; opt: {rand: 4, ziel: 20, primer: 5}
//   PCR.plan(oben, opt)       Plan zu einer gegebenen oberen Kette (Länge ziel + 2·rand, nur ACGT):
//        M, a, b, rand, ziel, primer, oben, unten
//        basis(art, c)        Base der Kette "top" / "bottom" in Spalte c
//        primer1, primer2     Basenfolgen der beiden Primer, jeweils 5'→3' gelesen
//   PCR.start(p)              [{top, bottom}]: ein Doppelstrang aus den beiden Ausgangssträngen
//   PCR.neuUnten(p, top)      der Strang, der an der oberen Vorlage entsteht (Primer 2): {art: "bottom", von: top.von, bis: b}
//   PCR.neuOben(p, bottom)    der Strang, der an der unteren Vorlage entsteht (Primer 1): {art: "top", von: a, bis: bottom.bis}
//   PCR.zyklus(ds, p, opt)    ein Zyklus: Jeder Doppelstrang liefert zwei neue (opt.bindet = false: die Primer finden keinen Halt;
//                             opt.verlaengert = false: die Polymerase arbeitet nicht; dann bleibt alles, wie es war).
//                             Reihenfolge: zu jedem Doppelstrang (top, bottom) erst {top, neuUnten(top)}, dann {neuOben(bottom), bottom}.
//   PCR.zaehle(ds, p)         {gesamt, ziel}: Zahl der Doppelstränge und davon die mit Zielsequenz in beiden Strängen
//   PCR.zahlen(n)             {gesamt, ziel} nach n Zyklen in geschlossener Form: 2^n und 2^n − 2n (ab n = 2, davor 0)
//   PCR.tm(laenge, gc)        Faustformel für die Schmelztemperatur eines Primers (gc als Anteil 0 … 1)
//   PCR.bewerte(T, tm)        "heiss" (T > tm), "gut" (tm − 10 ≤ T ≤ tm) oder "kalt" (T < tm − 10)
//
// Thermocycler-Lauf (Ansicht „Im Thermocycler“): ein ganzer Lauf wird vorab berechnet, die Seite spielt ihn ab.
//   PCR.lauf(prog, opt)       prog: {start: {T, s}, den: {T, s}, anl: {T, s}, ela: {T, s}, n, ende: {T, s}, halt}
//                             opt:  {N0, laenge, tm, pol: "taq" | "normal", heiz, kuehl, dt}
//        → {spur: [[t, T]], zyklen: [{n, t, N, M, E, akt, fDen, fAnl, fEla}], phasen: [{name, zyklus, t0, t1, T}], dauer, tmProdukt, nmax}
//        N: Zahl der Zielmoleküle (Doppelstränge), M: Nebenprodukte (nur bei zu kalter Anlagerung), E: Ausbeute des Zyklus (0 … 1),
//        akt: verbliebene Aktivität der Polymerase (0 … 1). Modellannahmen stehen in pcr.html unter „Modell und Grenzen“.
const PCR = (() => {
  const KOMP = {A: "T", T: "A", G: "C", C: "G"};
  const STANDARD = {rand: 4, ziel: 20, primer: 5};

  function plan(oben, opt){
    const o = Object.assign({}, STANDARD, opt || {});
    const {rand, ziel, primer} = o;
    const M = ziel + 2 * rand;
    if (!/^[ACGT]+$/.test(oben) || oben.length !== M) throw new Error("Die Vorlage passt nicht zu ziel und rand");
    if (primer < 1 || 2 * primer > ziel) throw new Error("Zwei Primer müssen in die Zielsequenz passen");
    if (rand < 1) throw new Error("Die Ausgangs-DNA braucht auf beiden Seiten Randstücke, sonst ist sie schon die Zielsequenz");
    const unten = oben.replace(/./g, b => KOMP[b]);
    const a = rand, b = rand + ziel;
    return {
      M, a, b, rand, ziel, primer, oben, unten,
      basis: (art, c) => (art === "top" ? oben : unten)[c],
      primer1: oben.slice(a, a + primer),
      primer2: unten.slice(b - primer, b).split("").reverse().join("")
    };
  }

  function erzeuge(zufall, opt){
    const o = Object.assign({}, STANDARD, opt || {});
    let oben = "";
    for (let i = 0; i < o.ziel + 2 * o.rand; i++) oben += "ACGT"[Math.floor(zufall() * 4)];
    return plan(oben, o);
  }

  const start = p => [{top: {art: "top", von: 0, bis: p.M}, bottom: {art: "bottom", von: 0, bis: p.M}}];
  const neuUnten = (p, top) => ({art: "bottom", von: top.von, bis: p.b});
  const neuOben = (p, bottom) => ({art: "top", von: p.a, bis: bottom.bis});

  function zyklus(ds, p, opt){
    const o = Object.assign({bindet: true, verlaengert: true}, opt || {});
    if (!o.bindet || !o.verlaengert) return ds.map(d => ({top: Object.assign({}, d.top), bottom: Object.assign({}, d.bottom)}));
    const aus = [];
    ds.forEach(d => {
      aus.push({top: d.top, bottom: neuUnten(p, d.top)});
      aus.push({top: neuOben(p, d.bottom), bottom: d.bottom});
    });
    return aus;
  }

  const istZiel = (p, s) => s.von === p.a && s.bis === p.b;
  const zaehle = (ds, p) => ({gesamt: ds.length, ziel: ds.filter(d => istZiel(p, d.top) && istZiel(p, d.bottom)).length});
  const zahlen = n => ({gesamt: Math.pow(2, n), ziel: n < 2 ? 0 : Math.pow(2, n) - 2 * n});

  const tm = (laenge, gc) => 64.9 + 41 * (gc * laenge - 16.4) / laenge;
  const bewerte = (T, tmWert) => (T > tmWert ? "heiss" : (T >= tmWert - 10 ? "gut" : "kalt"));

  /* ---------- Thermocycler-Lauf ---------- */
  const NMAX = 1e12;                 // Obergrenze der Kopien (Primer und Nukleotide reichen in 25 µl für etwa so viele Produkte)
  const EMAX = 0.95;
  const KHALB = 2e10;                // bei so vielen Kopien ist die Ausbeute pro Zyklus nur noch halb so groß                 // Ausbeute pro Zyklus unter guten Bedingungen
  const sig = x => 1 / (1 + Math.exp(x));
  // Schmelztemperatur des Produkts (Faustformel für lange DNA, 50 mM Na+, 50 % G+C)
  const tmProdukt = laenge => 81.5 + 16.6 * Math.log10(0.05) + 0.41 * 50 - 675 / laenge;
  // Halbwertszeit der Polymerase in Sekunden bei der Temperatur T (Taq: etwa 40 min bei 95 °C und 9 min bei 97,5 °C, Lawyer et al. 1993)
  const halbwert = (pol, T) => pol === "taq" ? (T < 90 ? Infinity : 2400 * Math.pow(2, -(T - 95) / 1.16)) : (T < 45 ? Infinity : 2);
  // Geschwindigkeit der Polymerase (Basen pro Sekunde): Faustregel 1000 Basen pro Minute bei 72 °C, langsamer bei anderen Temperaturen
  const tempo = T => 1000 / 60 * Math.exp(-Math.pow((T - 72) / 15, 2));
  const fDen = (T, s, tmP) => sig((tmP + 8 - T) / 1.5) * Math.min(1, s / 5);   // im PCR-Puffer mit Mg2+ trennen sich lange Stränge erst deutlich über der Faustformel-Tm
  const fAnl = (T, s, tmPr) => sig((T - (tmPr + 3)) / 1) * Math.min(1, s / 10);   // bis etwa Tm binden fast alle Primer, darüber kaum noch
  const neben = (T, tmPr) => T < tmPr - 10 ? 2e-4 * Math.min(1, (tmPr - 10 - T) / 10) : 0;

  function lauf(prog, opt){
    const o = Object.assign({N0: 1000, laenge: 500, tm: 52, pol: "taq", heiz: 3, kuehl: 2, dt: 1}, opt || {});
    const tmP = tmProdukt(o.laenge);
    const spur = [], phasen = [], zyklen = [];
    let t = 0, T = 25, akt = 1;
    spur.push([0, T]);
    function fahre(name, Tz, s, zyk){        // Rampe zur Zieltemperatur, dann halten; Enzym altert unterwegs
      const t0 = t;
      while (Math.abs(T - Tz) > 1e-9){
        const r = Tz > T ? o.heiz : -o.kuehl, schritt = Math.min(o.dt, Math.abs((Tz - T) / r));
        T += r * schritt; t += schritt; if (Math.abs(T - Tz) < 1e-9) T = Tz;
        akt *= Math.pow(0.5, schritt / halbwert(o.pol, T)); spur.push([t, T]);
      }
      const tHalt = t;
      for (let r = s; r > 1e-9; r -= o.dt){ const sch = Math.min(o.dt, r); t += sch; akt *= Math.pow(0.5, sch / halbwert(o.pol, T)); spur.push([t, T]); }
      phasen.push({name, zyklus: zyk, t0, tHalt, t1: t, T: Tz});
    }
    let N = o.N0, M = 0;
    zyklen.push({n: 0, t: 0, N, M, E: 0, akt: 1, fDen: 0, fAnl: 0, fEla: 0});
    if (prog.start && prog.start.s > 0) fahre("Erste Denaturierung", prog.start.T, prog.start.s, 0);
    for (let n = 1; n <= prog.n; n++){
      fahre("Denaturierung", prog.den.T, prog.den.s, n);
      const a = fDen(prog.den.T, prog.den.s, tmP);
      fahre("Primeranlagerung", prog.anl.T, prog.anl.s, n);
      const b = fAnl(prog.anl.T, prog.anl.s, o.tm), nb = neben(prog.anl.T, o.tm);
      const akt0 = akt;
      fahre("Elongation", prog.ela.T, prog.ela.s, n);
      const aktM = (akt0 + akt) / 2;
      const c = Math.min(1, prog.ela.s * tempo(prog.ela.T) * aktM / o.laenge);
      const frei = Math.max(0, 1 - (N + M) / NMAX) / (1 + (N + M) / KHALB);   // viele Produkte: die Stränge paaren sich lieber wieder miteinander als mit den Primern
      const E = EMAX * a * b * c * frei;
      const Mneu = M * (1 + EMAX * a * c * frei) + N * nb * a * c * frei;   // Nebenprodukte entstehen an falschen Stellen und vermehren sich dann wie das Ziel
      N = N * (1 + E); M = Mneu;
      zyklen.push({n, t, N, M, E, akt, fDen: a, fAnl: b, fEla: c});
    }
    if (prog.ende && prog.ende.s > 0) fahre("Letzte Elongation", prog.ende.T, prog.ende.s, prog.n);
    if (prog.halt !== undefined && prog.halt !== null) fahre("Kühlen", prog.halt, 0, prog.n);
    return {spur, zyklen, phasen, dauer: t, tmProdukt: tmP, nmax: NMAX};
  }

  return {plan, erzeuge, start, neuUnten, neuOben, zyklus, zaehle, zahlen, tm, bewerte, lauf, KOMP};
})();
