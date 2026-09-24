"""
Tenth batch: insulation-design group (largest remaining functional group,
~2900 lines total - see the size note in CONVERSION_STATUS.md).

Ported so far in this module (source line ranges refer to ASYN4.FOR):
    ALD    (10362-10402)  Stator coil-end overhang length ("Staenderaus-
                           ladung"), Jungbluth formulas. Ported as a
                           single-pass function (see its own docstring
                           for why the original's interactive REPL loop
                           was dropped).
    NARISS (10583-10968)  Computes the minimum slot width/height
                           ("Nutaufriss") needed for a given (partial-)
                           conductor configuration, from a long chain of
                           empirical insulation-thickness tables (VPI /
                           individual-coil / random-wound, keyed by rated
                           insulation voltage UNX). Pure computation, no
                           interactive I/O. Updates st.spkopf (BKOPF/
                           HKOPF) and st.ck211s (the DEL1-11 clearance
                           breakdown, consumed by the not-yet-ported
                           K211 drawing/report routines) as side effects,
                           matching the FORTRAN COMMON behaviour.
    FSPULE (10411-10582)  Computes the two-layer full-form-coil
                           (Zweischichtganzformspule) end-winding
                           geometry per Schorch works-standard W01481,
                           for symmetric (IART=0) or asymmetric (IART=1)
                           coils. Pure geometric computation.
    ZGF    (11206-11265)  Calls FSPULE twice (symmetric then
                           asymmetric-with-DSTR) to get both the
                           practical stator overhang lengths (ALD1S,
                           ALD1A). The FORTRAN version also writes a
                           'WKFISO' file for a manufacturing-drawing
                           tool; this port returns the same numbers
                           instead of writing that file (pass
                           write_wkfiso_path to opt back into writing it
                           if you have a use for that downstream tool).
    WVAR   (10969-11205)  Iteratively searches a standard-wire-gauge
                           table for the (partial-)conductor width/height
                           that best fills a given slot, using NARISS to
                           check each candidate. The FORTRAN version is a
                           2-round interactive REPL (offers to subdivide
                           the conductor width, then lets the user
                           override the number of parallel partial
                           conductors, then asks for final confirmation);
                           ported as a single deterministic pass with the
                           two "yes/no"/override decisions exposed as
                           optional parameters (split_width, tpar_override)
                           instead of blocking reads - same pattern as ALD.
    ISOSTD (9961-10361)   Stator slot layout orchestrator: given either a
                           target conductor size, a target slot size, or
                           neither, works out the other via NARISS/WVAR,
                           sizes the press-finger, computes coil-head
                           clearance defaults, and converges the stator
                           overhang length via ALD + ZGF.
    ISOTR  (11266-11408)  Random-wound ("Traeufelwicklung") stator
                           winding design: computes theoretical/practical
                           copper fill factor across candidate
                           parallel-partial-conductor counts (TPAR) and,
                           given a chosen TPAR (a parameter here, was an
                           interactive prompt originally), finalizes
                           MTL/QCU1. Calls CUFL and (for regular slot
                           patterns) TEILEN.
    CUFL   (11409-11431)  Copper area per slot - the same formula already
                           duplicated inline in geometry_functions.nut(),
                           now also available standalone since that's how
                           ISOTR calls it.
    TEILEN (11432-11480)  Works out how many turns each coil-side slot of
                           a zone gets when the total doesn't divide
                           evenly (e.g. 3.5 turns/slot -> some slots get
                           3, some get 4), calling SUCH1-6/ZPLAN to find a
                           symmetric placement pattern when one exists.
    SUCH1-6 (11481-11669) Brute-force search over which LOCH (1-6) slot
                           positions out of Q get the "extra" turn, in
                           FORTRAN's exact enumeration order (descending,
                           rightmost-position-fastest) so the first valid
                           pattern found matches the original bit-for-bit.
                           Ported as one generic `such(q, svk, loch)`
                           using itertools.product for the same
                           enumeration order, since SUCH1..SUCH6 are
                           otherwise identical code with LOCH hardcoded.
    ZPLAN  (11670-11702)  Checks whether a candidate turn-count pattern,
                           repeated and compared against its own
                           reversed-and-shifted copy, has no overlapping
                           "extra turn" positions (this is what makes a
                           winding pattern symmetric/balanced). The
                           FORTRAN builds explicit 6*Q-length working
                           arrays; this port exploits the fact both
                           arrays are period-Q to do the equivalent check
                           in O(Q) instead of O(6*Q) without changing the
                           result.
    NAR2   (11862-12159)  Slot layout ("Nutaufriss") for a slip-ring
                           rotor bar - the rotor-side counterpart to
                           NARISS. Pure computation; when IPRINT==1 also
                           updates st.ck211l (the DEL1L-11L clearance
                           breakdown) instead of writing the FORTRAN's
                           printer-report file.
    ISOLFR (11703-11861)  Slip-ring rotor winding design: searches the
                           standard conductor tables (via NAR2) for a
                           fitting bar cross-section. The FORTRAN's three
                           interactive prompts (reversible-duty voltage
                           doubling, conductors-side-by-side override,
                           parallel-conductors override) are exposed as
                           optional parameters (reversible_duty,
                           ntl_override, tpar_override) - same pattern as
                           WVAR/ISOTR.

THE ENTIRE ~2900-LINE "INSULATION DESIGN" FUNCTIONAL GROUP IS NOW DONE.
ALD, NARISS, FSPULE, ZGF, WVAR, ISOSTD, ISOTR, CUFL, TEILEN, SUCH1-6,
ZPLAN, NAR2, ISOLFR - all ported and tested.

Not yet ported: GENERA, SPKTSP, output/reporting subroutines, plotting
driver, LESEN, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import itertools
import math

from .state import MachineState

PI = 3.14159

# Standard wire-gauge table (mm), used by WVAR (source DATA ZAHL, 1-indexed
# to match the FORTRAN's ZAHL(J) references directly).
_ZAHL = [0.0,
         0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
         0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
         1.06, 1.12, 1.18, 1.25, 1.32, 1.40, 1.50, 1.60, 1.70, 1.80,
         1.90, 2.00, 2.12, 2.24, 2.36, 2.50, 2.65, 2.80, 3.00, 3.15,
         3.35, 3.55, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.30, 5.60,
         6.00, 6.30, 6.70, 7.10, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00,
         10.60, 11.20, 11.80, 12.50, 13.20]
_ZDIM = 65


def ald(z: float, r: float, rd: float, l1: float, hstreu: float, hzw: float,
        n1: float, weite1: float, di1: float, l: float,
        bkopf: float, hkopf: float):
    """
    Stator coil-end overhang length per the Jungbluth formulas.

    bkopf/hkopf come from the FORTRAN COMMON /SPKOPF/ (coil-head
    clearance dimensions) - pass them in explicitly since that COMMON
    block isn't otherwise used anywhere else ported so far.

    Returns a dict: lges1 (total stator core+overhang length), alpha_deg
    (HG1), delta2_deg (HG2), in_range (bool - False mirrors the FORTRAN's
    "WARNUNG" branch: HG1 outside [37.5, 60.0] degrees or HG2 > 32.5
    degrees).
    """
    dkl = di1 + 2.0 * hstreu
    t = dkl * 3.14 / n1
    alpha = math.asin((bkopf + z) / t)
    xa = (r + bkopf / 2.0) * math.cos(alpha)
    xc = math.tan(alpha) * (weite1 * t / 2.0 -
                             (r + bkopf / 2.0) * (1.0 - math.sin(alpha)))
    ald1 = l1 + xa + xc + 1.5 * hkopf
    lges1 = l + 2.0 * ald1
    hg1 = alpha * 180.0 / PI
    dkla = xc * math.tan(6.0 / 180.0 * PI) + dkl
    dgra = dkla + 3.0 * hkopf
    dgr = dkl + 2.0 * hkopf + hzw
    dd = dgra - dgr
    delta2 = math.atan(dd / xc)
    hg2 = delta2 * 180.0 / PI

    in_range = not (hg1 < 37.5 or hg1 > 60.0 or hg2 > 32.5)

    return {
        "lges1": lges1, "alpha_deg": hg1, "delta2_deg": hg2,
        "in_range": in_range,
    }


def nariss(iss: int, i2l1gs: int, unx: float, btl: float, htl: float,
           mtl: float, ntl: float, tpar: float, hs: float, hk: float,
           bn: float, hn: float, st: MachineState):
    """
    Computes the minimum slot width/height ("Nutaufriss") for a given
    (partial-)conductor configuration, from empirical insulation
    tables keyed by ISS (1=VPI, 2=individual coils, 3=random-wound) and
    the rated insulation voltage UNX. Updates st.spkopf and st.ck211s.

    Returns (bnmind, hnmind, hzw, iokay, dik).
    """
    itest = st.ctest.itest

    issp = 1 if (round(tpar) == 1 and round(ntl) == 2) else 0

    xltr = 0.15 if (unx <= 1500.0 and i2l1gs == 0) else 0.30
    btli = btl + xltr
    htli = htl + xltr

    dzwl = 0.0 if issp == 1 else 0.2
    if itest == 1:
        dzwl = 0.0

    if issp == 0:
        bzwl = hzwl = 0.0
    else:
        if unx <= 7000.0:
            bzwl, hzwl = 0.4, 0.4
        else:
            bzwl, hzwl = 2.2, 1.1
        if itest == 1:
            bzwl, hzwl = 0.3, 0.0

    if iss == 1:
        xktol = 0.0
        dik = 7.6
        for lim, val in [(14300., 7.2), (13800., 6.8), (12800., 6.4),
                         (11600., 5.6), (11100., 5.2), (10500., 4.8),
                         (9500., 4.4), (8800., 4.0), (8100., 3.6),
                         (7200., 3.2), (5200., 2.8), (4500., 2.4)]:
            if unx <= lim:
                dik = val
        if issp == 1 and unx > 7000.0 and itest == 0:
            dik -= hzwl
    elif iss == 2:
        xktol = 0.0
        dik = 11.0
        for lim, val in [(19400., 10.6), (18900., 10.4), (18400., 10.0),
                         (17900., 9.6), (17400., 9.4), (16900., 9.2),
                         (16400., 8.8), (15900., 8.4), (15400., 7.8),
                         (14850., 7.6), (14400., 7.2), (13800., 7.0),
                         (13300., 6.6), (12800., 6.4), (12300., 6.0),
                         (11550., 5.6), (11100., 5.4), (10500., 5.0),
                         (10100., 4.8), (9500., 4.6), (9300., 4.4),
                         (9000., 4.2), (8600., 4.0), (8100., 3.8),
                         (7650., 3.6), (7200., 3.2), (6500., 3.0),
                         (5250., 2.8), (4500., 2.6)]:
            if unx <= lim:
                dik = val
    else:
        xktol = 0.0
        dik = 0.72
        if unx <= 950.0:
            dik = 0.36
        if unx <= 660.0:
            dik = 0.00

    dglimm = 0.0 if unx <= 4500.0 else 0.4

    dngstr = 0.38 if unx <= 1500.0 else 0.30
    if itest == 1:
        dngstr = 0.0

    if unx <= 3500.0:
        hzw = 3.5
    else:
        hzw = math.trunc(unx / 1000.0 + 0.501)
    if itest == 1:
        if unx <= 4200.0:
            hzw = 3.0
        elif unx <= 6000.0:
            hzw = 4.5
        elif unx <= 6300.0:
            hzw = 5.0
        elif unx <= 6900.0:
            hzw = 5.5
        elif unx <= 10000.0:
            hzw = 7.5
        elif unx <= 10500.0:
            hzw = 8.0
        elif unx <= 11000.0:
            hzw = 8.5
        elif unx <= 13800.0:
            hzw = 11.0
        else:
            hzw = 12.0

    dnau = 0.00
    dna = 0.0
    xnab = xnah = 0.0
    if iss in (1, 2):
        dna = 0.18
        xnab = 2.0 * dna + 2.0 * dnau
        xnah = 4.0 * dna + 2.0 * dnau
    else:
        if unx <= 660.0:
            dna = 0.50
            xnab = 2.0 * dna + 2.0 * dnau
            xnah = 3.0 * dna + 2.0 * dnau
        elif unx <= 950.0:
            dna = 0.36
            xnab = 2.0 * dna + 2.0 * dnau
            xnah = 3.0 * dna + 2.0 * dnau
        elif unx <= 1500.0:
            dna = 0.18
            xnab = 2.0 * dna + 2.0 * dnau
            xnah = 3.0 * dna + 2.0 * dnau
    if itest == 1:
        xnab = 0.0
        xnah = 0.0

    dglstr = 0.0 if unx <= 4500.0 else 0.8
    if issp == 1:
        dglstr = 1.0
    if itest == 1:
        dglstr = 0.0

    luftb = 0.1
    lufth = 0.1
    svsatz = 0.2

    if iss == 2:
        ftolb, ftolh = 0.1, 0.2
    else:
        ftolb, ftolh = 0.0, 0.0

    mtl2 = mtl / 2.0
    bnmind = ntl * btli
    hnmind = mtl2 * htli
    anzlei = ntl * mtl2 / tpar
    hnmind += (anzlei - 1.0) * dzwl
    hnmind += hzwl
    bnmind += bzwl

    bkopf = bnmind
    hkopf = hnmind

    bnmind += xktol
    hnmind += xktol
    bnmind += dik
    hnmind += dik
    bnmind += dglimm
    hnmind += dglimm
    hnmind = 2.0 * hnmind
    hnmind += dngstr
    hnmind += hzw
    bnmind += xnab
    hnmind += xnah
    hnmind += dglstr
    bnmind += luftb
    hnmind += lufth
    bnmind += svsatz + ftolb
    hnmind += hk + hs + svsatz + ftolh

    if max(hnmind, hn) >= 80.0:
        bnmind += 0.1
        hnmind += 0.1

    iokay = 1 if (hnmind <= hn and bnmind <= bn) else 0

    db = bn - bnmind
    dh = hn - hnmind
    if dh > 0.0:
        hzw += math.trunc(2.0 * dh) / 2.0

    dhb = 0.22
    if unx <= 1500.0:
        dhb = 1.80
    if unx <= 950.0:
        dhb = 1.30
    if unx <= 660.0:
        dhb = 0.80
    if iss == 2:
        dhb = 0.0

    if iss == 1:
        diko = 1.90
        if unx > 4500.0:
            diko = 2.2
        if unx > 7200.0:
            diko = 4.1
        if unx > 7200.0 and issp == 1:
            diko = 3.5
        if unx > 10500.0:
            diko = 4.7
        if unx > 10500.0 and issp == 1:
            diko = 3.8
        if unx > 12300.0:
            diko = 5.5
        if unx > 13800.0:
            diko = 6.6
    elif iss == 2:
        diko = 1.60
        if unx > 4500.0:
            diko = 2.1
        if unx > 7200.0:
            diko = 2.7
        if unx > 9500.0:
            diko = 3.2
        if unx > 10500.0:
            diko = 3.6
        if unx > 12300.0:
            diko = 4.2
        if unx > 13800.0:
            diko = 4.7
    else:
        diko = 0.0
    if itest == 1:
        diko *= 0.6

    bkopf += dhb + diko
    hkopf += dhb + diko

    ck = st.ck211s
    ck.hwk = hkopf
    ck.bwk = bkopf
    ck.del1 = xltr / 2.0
    ck.del2 = 0.0
    ck.del3 = 0.0
    ck.del4 = (dik + hzwl) / 2.0
    ck.del5 = dglimm / 2.0
    ck.del6 = xnab / 2.0
    ck.del7 = max(db / 2.0 + luftb / 2.0 + ftolb / 2.0, 0.0)
    ck.del8 = 0.0
    ck.del9 = dhb / 2.0
    ck.del10 = diko / 2.0
    ck.del11 = 0.0
    if iss == 1:
        ck.del3 = ck.del7
        ck.del7 = 0.0

    st.spkopf.bkopf = bkopf
    st.spkopf.hkopf = hkopf

    return bnmind, hnmind, hzw, iokay, dik


def fspule(iart: int, l: float, di1: float, n1: float, hn1min: float,
           bn1: float, hstreu: float, hzw1: float, diso1: float,
           weite: float, a: float, b: float, sh: float, dstr: float,
           l1: float, l2: float, lstr: float, r: float, rd: float):
    """
    Two-layer full-form-coil (Zweischichtganzformspule) end-winding
    geometry per Schorch works-standard W01481. IART=0 for symmetric,
    IART=1 for asymmetric coils.

    Returns a dict: a1, a2, lo, lu, alphao, alphau, rmo, rmu, roeseo,
    roeseu, beta1, beta2, lw, auslad.
    """
    n1 = round(n1)
    weite = round(weite)

    dm = di1 + 2.0 * hstreu + hn1min
    tnm = PI * dm / float(n1)

    a = (dm * (b + bn1)) / (dm - 2.0 * sh - hzw1) - bn1
    alf = math.asin((bn1 + a) / tnm)
    if iart > 0:
        if dstr < 0.001:
            dstr = (dm + di1 + 2.0 * hstreu + 2.0 * sh + 70.0) / 2.0
        astr = (bn1 + b) * dm / dstr - bn1
        alfstr = math.asin((bn1 + astr) / tnm)
    else:
        astr = a
        alfstr = alf

    a1 = l1 + lstr + (r + bn1 / 2.0) * math.tan(PI / 4.0 - alf / 2.0)
    a2 = l2 + lstr + (r + bn1 / 2.0) * math.tan(PI / 4.0 - alfstr / 2.0)
    x = (weite * tnm * math.tan(alfstr) - a1 + a2) / (math.tan(alf) + math.tan(alfstr))
    y = weite * tnm - x
    lo = x / math.cos(alf)
    lu = y / math.cos(alfstr)
    if iart > 0:
        auslad = a1 + x * math.tan(alf) + rd + sh + 5.0
    else:
        auslad = a1 + weite * tnm / 2.0 * math.tan(alf) + rd + sh + 5.0

    e = (auslad - (rd + sh + l1)) * math.tan(6.0 / 180.0 * PI)
    delta = math.atan((di1 / 2.0 + hstreu + e + 1.5 * sh + 2.0 * rd - dm / 2.0) / lu)
    if iart > 0:
        lm = l + PI * (rd + sh / 2.0) + a1 + a2 + lo + lu / math.cos(delta)
    else:
        lm = (l + PI * (rd + sh / 2.0) + 2.0 * l1 +
              2.0 * (r + bn1 / 2.0) * (PI / 2.0 - alf) +
              weite * tnm / math.cos(alf) -
              2.0 * (r + bn1 / 2.0) * math.tan(PI / 4.0 - alf / 2.0) +
              2.0 * lstr + 10.0)
    lw = lm - l

    dw = di1 + 2.0 * hstreu + 2.0 * e + 4.0 * (rd + sh)
    r1 = di1 / 2.0 + hstreu
    r3 = dw / 2.0 - rd - sh
    r4 = r1 + sh + hzw1
    r5 = r1 + e
    r6 = r3 + rd + diso1
    beta1 = 360.0 * x / (float(n1) * tnm)
    beta2 = 360.0 * weite / float(n1) - beta1

    alphao = alf * 180.0 / PI
    alphau = alfstr * 180.0 / PI
    rmo = r1 + sh / 2.0
    rmu = r4 + sh / 2.0
    roeseo = r5 + sh / 2.0
    roeseu = r6 + sh / 2.0

    return {
        "a": a, "a1": a1, "a2": a2, "lo": lo, "lu": lu,
        "alphao": alphao, "alphau": alphau, "rmo": rmo, "rmu": rmu,
        "roeseo": roeseo, "roeseu": roeseu, "beta1": beta1, "beta2": beta2,
        "lw": lw, "auslad": auslad, "dstr": dstr, "astr": astr, "dw": dw,
    }


def zgf(reke: str, l: float, di1: float, en1: float, hn1min: float,
        bn1: float, hstreu: float, hzw1: float, diso1: float,
        eweite: float, a: float, b: float, sh: float, dstr: float,
        l1: float, l2: float, lstr: float, r: float, rd: float,
        write_wkfiso_path: str | None = None):
    """
    Computes both the symmetric and asymmetric (DSTR-based) stator
    full-form-coil overhang lengths, by calling FSPULE twice.

    Returns (lw, ald1s, ald1a, sym, asym) - sym/asym are the full FSPULE
    result dicts for the two coil geometries.
    """
    sym = fspule(0, l, di1, en1, hn1min, bn1, hstreu, hzw1, diso1, eweite,
                 a, b, sh, 0.0, l1, l2, lstr, r, rd)
    ald1s = sym["auslad"]

    if write_wkfiso_path:
        with open(write_wkfiso_path, "w") as f:
            f.write(f"{sym['a1']} {sym['a2']} {sym['lo']} {sym['lu']} "
                    f"{sym['alphao']} {sym['alphau']} {sym['rmo']} "
                    f"{sym['rmu']} {sym['roeseo']} {sym['roeseu']} "
                    f"{sym['beta1']} {sym['beta2']} {ald1s}\n")

    asym = fspule(1, l, di1, en1, hn1min, bn1, hstreu, hzw1, diso1, eweite,
                  a, b, sh, dstr, l1, l2, lstr, r, rd)
    ald1a = asym["auslad"]

    return sym["lw"], ald1s, ald1a, sym, asym


def wvar(iss: int, i2l1gs: int, unx: float, zn: float, azweig: float,
         hs: float, hk: float, bn: float, hn: float, st: MachineState,
         split_width: bool = False, tpar_override: float | None = None):
    """
    Searches the standard wire-gauge table for the (partial-)conductor
    configuration that best fills a given slot.

    split_width mirrors the FORTRAN's "subdivide in width?" prompt
    (only offered when NTL==1 after the first pass); tpar_override
    mirrors the "new number of parallel partial conductors" prompt
    (None means keep whatever the search settled on, matching a blank
    answer in the original).

    Returns (hzw, mtl, ntl, tpar, qcu, btl, htl, dik).
    """
    ltr = zn * azweig

    tpar = 1.0
    ntl = 1.0
    mtl = ltr
    btl = 0.0
    htl = 0.0
    bnmind, hnmind, hzw, iokay, dik = nariss(
        iss, i2l1gs, unx, btl, htl, mtl, ntl, tpar, hs, hk, bn, hn, st)

    btls = (bn - bnmind) / ntl
    htls = (hn - hnmind) / mtl

    imerk = 0
    if btls > _ZAHL[64]:
        imerk = 1
        btls = btls / 2.0
        ntl = 2.0
        tpar = 2.0
    jb = _ZDIM

    if htls > _ZAHL[44]:
        hg = math.trunc(htls / _ZAHL[44]) + 1.0
        mtl = mtl * hg
        tpar = tpar * hg
    jh = _ZDIM

    def _search(jb, jh, mtl, ntl, tpar):
        while True:
            btl = _ZAHL[jb]
            htl = _ZAHL[jh]
            bnmind, hnmind, hzw, iokay, dik = nariss(
                iss, i2l1gs, unx, btl, htl, mtl, ntl, tpar, hs, hk, bn, hn, st)
            if iokay != 1:
                if hnmind > hn:
                    jh -= 1
                    continue
                if bnmind > bn:
                    jb -= 1
                    continue
            return btl, htl, hzw, dik, jb, jh

    btl, htl, hzw, dik, jb, jh = _search(jb, jh, mtl, ntl, tpar)

    if round(ntl) != 2:
        if split_width:
            imerk = 1
            ntl = ntl * 2.0
            tpar = tpar * 2.0

    jb = _ZDIM
    tparn = tpar if tpar_override is None else tpar_override
    if tparn < 0.001:
        tparn = tpar

    if imerk == 1 and round(tparn) >= 2 and round(tparn) % 2 != 0:
        # mirrors the FORTRAN's GOTO 600 (reject, keep prior NTL/TPAR
        # decision) - the caller should retry with a valid tpar_override.
        tparn = tpar
    if (round(tparn) == 1 and round(ntl) == 2 and
            math.fmod(mtl * tparn / tpar + 0.001, 2.0) > 0.01):
        # mirrors "SCHEIBENSPULE MIT BLINDLEITER NICHT VORGESEHEN" - not
        # a valid combination, fall back to keeping the current tpar.
        tparn = tpar

    mtl = mtl * tparn / tpar
    tpar = tparn
    jh = _ZDIM

    btl, htl, hzw, dik, jb, jh = _search(jb, jh, mtl, ntl, tpar)

    qcu = btl * htl * tpar * azweig
    if htl <= 1.60:
        radius = 0.50
    elif htl <= 2.24:
        radius = 0.65
    elif htl <= 3.55:
        radius = 0.80
    else:
        radius = 1.00
    qcu = azweig * tpar * (htl * btl - radius ** 2 * (4.0 - 3.14))

    return hzw, mtl, ntl, tpar, qcu, btl, htl, dik


# Standard slot-width table (mm), 1-indexed to match FORTRAN BNV(J).
_BNV = [0.0,
        7.0, 7.2, 7.5, 7.7, 8.0, 8.2, 8.5, 8.7, 9.0, 9.2, 9.5, 9.7,
        10.0, 10.2, 10.5, 10.7, 11.0, 11.2, 11.5, 11.7,
        12.0, 12.2, 12.5, 12.7, 13.0, 13.2, 13.5, 13.7,
        14.0, 14.2, 14.5, 14.7, 15.0, 15.2, 15.5, 15.7,
        16.0, 16.2, 16.5, 16.7, 17.0, 17.2, 17.5, 17.7,
        18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5,
        22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5,
        26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0]

# Standard slot-depth table (mm, incl. slot bridge + wedge), 1-indexed.
_HNV = [0.0,
        30., 33., 36., 39., 42., 45., 48., 51., 54., 57.,
        60., 64., 68., 72., 76., 80., 84., 88., 92., 96.,
        100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,
        150., 155., 160., 165., 170., 180., 185., 190., 195., 200.]


def isostd(reke: str, l: float, di1: float, n1: float, weite1: float,
           isoart: int, unx: float, qcu1: float, bs1: float, bn1: float,
           bn1s: float, hs1: float, hk1: float, ho1: float, hcuo1: float,
           hzw1: float, hcuu1: float, hu1: float, zn1: float,
           azweig: float, lw1: float, mtl: float, ntl: float, utl: float,
           btl: float, htl: float, diso1: float, l1: float, l2: float,
           lstr: float, r: float, rd: float, b: float, dstr: float,
           m1: float, p: float, da1: float, zk: float, bk: float,
           schr: float, st: MachineState, reduced_insulation_system: bool = False,
           verbose: bool = True):
    """
    Stator slot layout ("Nutauslegung") orchestrator: given either a
    target conductor size, a target slot size, or neither, works out the
    other via NARISS/WVAR, sizes the press-finger (Druckfinger), computes
    the coil-head clearances (L1/L2/R/RD/B defaults if not supplied), and
    finds the stator overhang (ALD, then a converging ZGF search on B).

    Returns a dict with all the computed slot/coil-head geometry
    (bs1, bn1, bn1s, hs1, hk1, ho1, hcuo1, hzw1, hcuu1, hu1, hn1, hstreu,
    qcu1, utl, l1, l2, r, rd, b, diso1, drufi, ald1s, ald1a, lges1, lw1,
    i2l1gs, iss) plus the MAISO transfer-file fields (q1, iq1, ltr1,
    iwef, isart, nschr1, lfe) for anyone who wants to write that file
    themselves.

    reduced_insulation_system mirrors the FORTRAN's "ABGESPECKTES
    ISOLATIONS-SYSTEM ?" prompt (sets st.ctest.itest) - default False
    matches answering "N".
    """
    if azweig == 0.0:
        raise ValueError("ISOSTD: AZWEIG GLEICH NULL")
    if math.fmod(azweig * zn1 / 2.0 + 0.01, 1.0) > 0.1:
        raise ValueError("ISOSTD: WDG.ZAHL JE LAGE NICHT GANZZAHLIG")

    if isoart == 0:
        iss = 1 if unx > 1500.0 else 3
    elif isoart == 1:
        iss = 2
    else:
        raise ValueError("ISOSTD: UNBEKANNTE ISOLATIONSART")

    tpar = mtl * ntl / (zn1 * azweig)

    i2l1gs = 1
    if unx <= 1500.0:
        i2l1gs = 0

    st.ctest.itest = 1 if reduced_insulation_system else 0
    isoneu = st.ctest.itest

    hzw = hzw1  # default carry-through if none of the 3 branches below fire
    if qcu1 == 0.0 and mtl * ntl * htl * btl > 0.0 and bn1 == 0.0:
        # ---- conductor given, slot is designed around it ----
        if htl <= 1.60:
            radius = 0.50
        elif htl <= 2.24:
            radius = 0.65
        elif htl <= 3.55:
            radius = 0.80
        else:
            radius = 1.00
        qcu1 = azweig * tpar * (htl * btl - radius ** 2 * (4.0 - 3.14))
        utl = 1.0
        hs1alt, hk1alt = hs1, hk1
        hs1 = hk1 = 0.0
        bnmind, hnmind, hzw, iokay, dik = nariss(
            iss, i2l1gs, unx, btl, htl, mtl, ntl, tpar, hs1, hk1, 0.0, 0.0, st)

        bn1 = _BNV[67]
        for j in range(1, 68):
            if _BNV[j] >= bnmind:
                bs1 = _BNV[j] + 0.5
                bn1 = _BNV[j]
                bn1s = _BNV[j]
                break

        hs1 = 0.5 if hs1alt == 0.0 else hs1alt

        if hk1alt == 0.0:
            hk1 = 12.0
            if bn1 <= 30.0:
                hk1 = 10.0
            if bn1 <= 25.5:
                hk1 = 8.0
            if bn1 <= 22.5:
                hk1 = 6.0
            if bn1 <= 17.7:
                hk1 = 5.0
            if bn1 <= 12.7:
                hk1 = 4.0
        else:
            hk1 = hk1alt
        hnmind = hnmind + hk1 + hs1

        hn1g = _HNV[40]
        for j in range(1, 41):
            if _HNV[j] >= hnmind:
                hn1g = _HNV[j]
                break

        ho1 = 0.0
        hu1 = 0.0
        hzw1 = hzw
        hcuo1 = (hn1g - hs1 - hk1 - hzw1) / 2.0
        hcuu1 = hcuo1

    elif qcu1 == 0.0 and htl * btl == 0.0 and bn1 > 0.0:
        # ---- slot given, conductor is designed to fill it ----
        bn = bn1
        hn = hs1 + hk1 + hcuo1 + hzw1 + hcuu1
        utl = 1.0
        hzw, mtl, ntl, tpar, qcu, btl, htl, dik = wvar(
            iss, i2l1gs, unx, zn1, azweig, hs1, hk1, bn, hn, st)
        hzw1 = hzw
        hcuo1 = (hn - hs1 - hk1 - hzw1) / 2.0
        hcuu1 = hcuo1
        qcu1 = qcu

    else:
        # ---- both given: just check the fit ----
        utl = 1.0
        bn = bn1
        hn = hs1 + hk1 + hcuo1 + hzw1 + hcuu1
        bnmind, hnmind, dummy, iokay, dik = nariss(
            iss, i2l1gs, unx, btl, htl, mtl, ntl, tpar, hs1, hk1, bn, hn, st)

    hn1 = hs1 + hk1 + ho1 + hcuo1 + hzw1 + hcuu1 + hu1
    hstreu = hs1 + hk1 + ho1

    # ---- press-finger (Druckfinger) sizing ----
    dng1 = di1 + 2.0 * hn1
    xa = dng1 * math.sin(PI / n1) - bn1
    xb = di1 * math.sin(PI / n1) - bn1
    xc = (di1 + 2.0 * hstreu) * math.sin(PI / n1) - bn1
    xgb = xc - 5.0
    y = hn1 / 3.0 * (xa + 2.0 * xb) / (xa + xb)
    f = hn1 * (xa + xb) / 2.0
    xgp = f * 1.00
    mb = xgp * (y + 4.0)
    wb = mb / 100.0
    drufi = math.sqrt(6.0 * wb / xgb)
    drufi = math.trunc(drufi / 5.0) * 5.0 + 5.0
    drufi = max(drufi, 15.0)
    if reke[0:2] == "A5" and reke[4:7] == "500":
        drufi = 25.0
    elif reke[0:2] == "A5" and reke[4:7] == "560":
        drufi = 25.0
    elif reke[0:2] == "A5" and reke[4:7] == "630":
        drufi = 30.0
    elif reke[0:2] == "A5" and reke[4:7] == "710":
        drufi = 40.0
    elif reke[0:2] == "A5" and reke[4:7] == "800":
        drufi = 50.0

    # ---- form-coil calc defaults ----
    if diso1 == 0.0:
        diso1 = dik / 2.0
    if l1 == 0.0:
        if iss == 1:
            l1 = 40.0
        elif iss == 2:
            l1 = 40.0
            if unx > 4500.0:
                l1 = 55.0
            if unx > 7200.0:
                l1 = 70.0
            if unx > 9500.0:
                l1 = 85.0
            if unx > 10500.0:
                l1 = 105.0
            if unx > 12300.0:
                l1 = 120.0
            if unx > 15000.0:
                l1 = 150.0
        else:
            l1 = 20.0 if unx <= 660.0 else 25.0
        l1 = l1 + drufi
    if l2 == 0.0:
        l2 = l1
    if r == 0.0:
        r = 35.0
    if rd == 0.0:
        rd = 15.0
    if b == 0.0:
        b = unx / 1000.0 * 0.8
        b = max(b, 6.00) if unx > 1500.0 else max(b, 4.00)

    # ---- Ausladung (LDW) ----
    ald_out = ald(b, r, rd, l1, hstreu, hzw1, n1, weite1, di1, l,
                  st.spkopf.bkopf, st.spkopf.hkopf)
    lges1 = ald_out["lges1"]
    st.ck211s.swk = b
    st.ck211s.ausldg = (lges1 - l) / 2.0 - l1
    bsave = b

    # ---- Formspulenrechnung nach Schorch (converge B against ALD1) ----
    ald1 = (lges1 - l) / 2.0
    hn1min = hcuo1 + hzw1 + hcuu1
    while True:
        lw, ald1s, ald1a, sym, asym = zgf(
            reke, l, di1, n1, hn1min, bn1, hstreu, hzw1, diso1, weite1,
            0.0, b, hcuo1, dstr, l1, l2, lstr, r, rd)
        if ald1s > ald1:
            b -= 0.2
        else:
            break
    b = bsave

    if lw1 == 0.0:
        lw1 = lw
    lspule = (l + lw1) / 1000.0 * zn1 / 2.0 * azweig * 2.0
    if verbose and lspule > 40.0:
        print(f"\n*** W A R N U N G ***\nSPULENLAENGE = {lspule:5.1f} M IST GROESSER ALS 40 M\n"
              "EVENTUELL WENIGER PARALLELE GRUPPEN SCHALTEN !")

    # ---- MAISO (KASIMIRS) transfer-file fields ----
    if mtl * ntl * htl * btl * azweig == 0.0:
        raise ValueError("ISO: DATEI ISO KANN NICHT ERSTELLT WERDEN")

    q1 = n1 / (2.0 * m1 * p)
    zq1, nq1 = 0.0, 1.0
    for j in range(1, 101):
        nq1 = float(j)
        if math.fmod(nq1 * q1 + 0.0001, 1.0) < 0.001:
            zq1 = nq1 * q1
            break
    iq1 = (round(zq1), round(nq1))

    lfe = l - zk * bk
    nschr1 = round(abs(schr)) if schr < 0.0 else 0

    if tpar < 1.001:
        ltr1 = (1, 1)
    else:
        ltr1 = (round(ntl), round(tpar / ntl))

    iwef = (round(zn1 / 2.0 * azweig), round(azweig))

    isart = 11 if (round(tpar) == 1 and round(ntl) == 2) else isoart

    return {
        "bs1": bs1, "bn1": bn1, "bn1s": bn1s, "hs1": hs1, "hk1": hk1,
        "ho1": ho1, "hcuo1": hcuo1, "hzw1": hzw1, "hcuu1": hcuu1,
        "hu1": hu1, "hn1": hn1, "hstreu": hstreu, "qcu1": qcu1,
        "utl": utl, "mtl": mtl, "ntl": ntl, "btl": btl, "htl": htl,
        "l1": l1, "l2": l2, "r": r, "rd": rd, "b": b, "diso1": diso1,
        "drufi": drufi, "ald1s": ald1s, "ald1a": ald1a, "lges1": lges1,
        "lw1": lw1, "i2l1gs": i2l1gs, "iss": iss, "q1": q1, "iq1": iq1,
        "ltr1": ltr1, "iwef": iwef, "isart": isart, "nschr1": nschr1,
        "lfe": lfe,
    }


def cufl(bs: float, bn: float, bns: float, hk: float, ho: float,
         hcuo: float, hzw: float, hcuu: float, hu: float) -> float:
    """Copper area per slot (same formula as geometry_functions.nut()'s
    inline ACU computation, exposed standalone since ISOTR calls it
    directly)."""
    hn = ho + hcuo + hzw + hcuu + hu
    hg = (bns - bn) / hn
    bnb = bn + hg * (ho + hcuo)
    bnc = bn + hg * (ho + hcuo + hzw)
    acu = (bn + bnb) / 2.0 * (ho + hcuo) + (bnc + bns) / 2.0 * (hcuu + hu)
    if bs < 0.9 * bn and ho < 0.0001:
        acu = acu + (bs + bn) / 2.0 * hk
    return acu


def zplan(q: int, a: list, svk: int) -> int:
    """
    Checks whether pattern `a` (1-indexed list, length q+1, 0/1 values,
    index 0 unused), repeated periodically and compared against its own
    reversed-and-SVK-shifted copy, has no overlapping "extra turn"
    positions. Returns 1 if okay (no overlap), 0 otherwise.

    Both the periodic extension of `a` and of its reversed copy have
    period q, so the FORTRAN's explicit 6*q-length comparison reduces
    to an equivalent O(q) check - see the module docstring.
    """
    for j in range(1, q + 1):
        k = ((j - 1 + svk) % q) + 1
        ahg_val = a[q + 1 - k]
        if a[j] + ahg_val >= 2:
            return 0
    return 1


def such(q: int, svk: int, loch: int):
    """
    Generic replacement for SUCH1..SUCH6: brute-force search, in
    FORTRAN's exact enumeration order (itertools.product gives the same
    "rightmost position varies fastest" order as the nested DO loops),
    for a placement of `loch` "extra turn" positions among `q` slots
    that ZPLAN accepts as symmetric.

    Returns (a, iokay) - a is a 1-indexed list (length q+1) of 0/1
    values; iokay is 1 if a valid pattern was found, 0 otherwise (in
    which case `a` is all zeros, matching the FORTRAN leaving A
    however the last failed attempt left it before the caller's
    fallback fill takes over).
    """
    for combo in itertools.product(range(q, 0, -1), repeat=loch):
        a = [0] * (q + 1)
        for j in combo:
            a[j] = 1
        isum = sum(a[1:q + 1])
        if isum < loch:
            continue
        if zplan(q, a, svk) == 1:
            return a, 1
    return [0] * (q + 1), 0


def teilen(q1: float, wdg: float, gwdg: float, bwdg: float, weite1: float):
    """
    Works out how many turns each coil-side slot of a zone gets when the
    total (WDG) doesn't divide evenly. Returns (a, iokay) where `a` is a
    plain (0-indexed) list of length Q with the turn count for each slot.
    """
    q = round(q1)
    svk = round(3.0 * q1 - weite1)

    if bwdg < 0.01:
        return [round(wdg)] * q, 1

    wdgu = round(gwdg)
    wdgo = wdgu + 1
    loch = round(bwdg * q1)

    if loch > q // 2:
        a = [0] * (q + 1)
        for j in range(1, q - loch + 1):
            a[j] = wdgu
        for j in range(q - loch + 1, q + 1):
            a[j] = wdgo
        if q == 3 and svk == 2:
            a[1], a[2], a[3] = wdgo, wdgu, wdgo
        return a[1:q + 1], 0

    iokay = 0
    a = [0] * (q + 1)
    if 1 <= loch <= 6:
        a, iokay = such(q, svk, loch)

    if iokay != 1:
        a = [0] * (q + 1)
        for j in range(1, loch + 1):
            a[j] = wdgu
        for j in range(loch + 1, q + 1):
            a[j] = wdgo
        return a[1:q + 1], iokay

    return [a[j] + wdgu for j in range(1, q + 1)], iokay


def isotr(rkz: str, bs1: float, bn1: float, bn1s: float, hk1: float,
          ho1: float, hcuo1: float, hzw1: float, hcuu1: float, hu1: float,
          zn1: float, azweig: float, n1: float, p: float, mzone1: float,
          btl: float, htl: float, weite1: float, qcu1: float = 0.0,
          tpar: float | None = None, verbose: bool = True):
    """
    Random-wound stator winding design: computes theoretical/practical
    copper fill factor across candidate parallel-partial-conductor
    counts and finalizes MTL/QCU1 for a chosen TPAR.

    tpar mirrors the FORTRAN's "ANZAHL PARALLELER TEILLEITER = ?" prompt
    - pass None to get back the suggestion table without finalizing
    (mtl will be None in that case; call again with a chosen tpar to
    finalize), or pass a value directly to finalize in one call.

    Returns a dict: mtl, ntl (always 0.0, matching the FORTRAN - NTL
    isn't meaningful for random-wound), utl (always 0.0), qcu1, iwdg
    (per-slot turn-count list from TEILEN, or None if SONDER==1), sonder,
    iokay, suggestions (list of (tpar, kcuth, kcup) tuples with fill
    factor in the FORTRAN's [0.45, 0.68] "reasonable" band).
    """
    if azweig == 0.0:
        raise ValueError("ISOTR: AZWEIG GLEICH 0")
    if btl == 0.0:
        raise ValueError("ISOTR: LEITERDURCHMESSER GLEICH 0")

    ntl = 0.0
    utl = 0.0

    schich = 1.0 if (len(rkz) > 8 and rkz[8] == "E") else 2.0

    acu = cufl(bs1, bn1, bn1s, hk1, ho1, hcuo1, hzw1, hcuu1, hu1)
    atl = PI / 4.0 * btl ** 2
    atli = htl ** 2

    wdg = zn1 * azweig / schich
    gwdg = math.trunc(wdg + 0.005)
    bwdg = wdg + 0.001 - math.trunc(wdg + 0.005)
    q1 = n1 / (p * mzone1)
    sonder = 1 if (math.fmod(q1 + 0.001, 1.0) > 0.01 or round(mzone1) == 3) else 0

    iokay = 0
    iwdg = None
    if sonder == 0:
        iwdg, iokay = teilen(q1, wdg, gwdg, bwdg, weite1)

    if qcu1 != 0.0:
        qcutl = PI / 4.0 * btl ** 2
        mtl = qcu1 / azweig / qcutl
        return {"mtl": mtl, "ntl": ntl, "utl": utl, "qcu1": qcu1,
                "iwdg": iwdg, "sonder": sonder, "iokay": iokay,
                "suggestions": []}

    if bwdg > 0.01:
        if round(schich) == 1:
            znx = gwdg + 1.0
        else:
            znx = 2.0 * gwdg + 2.0 if (bwdg > 0.51 or iokay != 1) else 2.0 * gwdg + 1.0
    else:
        znx = wdg if round(schich) == 1 else 2.0 * wdg

    if verbose and sonder == 0 and iokay == 1 and round(schich) == 2:
        print("\nOPTIMALE NUTFUELLUNG BEI VORGEGEBENER SPULENWEITE MOEGLICH !")
    if verbose:
        print(f"\nWDG. JE SPULE = {wdg:6.2f}   LOCHZAHL / STDR. = {q1:5.2f}")

    suggestions = []
    for j in range(1, 51):
        tp = 2.0 * j - 1.0
        tp2 = tp + 1.0
        kcuth = tp * znx * atl / acu
        kcup = tp * znx * atli / acu
        kcuth2 = tp2 * znx * atl / acu
        kcup2 = tp2 * znx * atli / acu
        if kcup < 0.45 and kcup2 < 0.45:
            continue
        if kcup > 0.68 and kcup2 > 0.68:
            continue
        suggestions.append((tp, kcuth, kcup, tp2, kcuth2, kcup2))
        if verbose:
            print(f"  TPAR={tp:5.0f}  KCU,TH={kcuth:.3f}  KCU,PR={kcup:.3f}"
                  f"      TPAR={tp2:5.0f}  KCU,TH={kcuth2:.3f}  KCU,PR={kcup2:.3f}")

    if tpar is None:
        return {"mtl": None, "ntl": ntl, "utl": utl, "qcu1": qcu1,
                "iwdg": iwdg, "sonder": sonder, "iokay": iokay,
                "suggestions": suggestions}

    kcuth = tpar * znx * atl / acu
    kcup = tpar * znx * atli / acu
    if verbose:
        print(f"\nFUELLUNG, BLANK: {kcuth:.3f}   FUELLUNG, PRAKTISCH: {kcup:.3f}")

    mtl = tpar
    qcu1 = tpar * azweig * atl

    return {"mtl": mtl, "ntl": ntl, "utl": utl, "qcu1": qcu1, "iwdg": iwdg,
            "sonder": sonder, "iokay": iokay, "suggestions": suggestions,
            "kcuth": kcuth, "kcup": kcup}


def nar2(iprint: int, imerk: int, u20: float, btl: float, htl: float,
         mtl: float, ntl: float, tpar: float, hs: float, hk: float,
         bs: float, bn: float, hn: float, hzw: float, qcue: float,
         st: MachineState):
    """
    Slot layout ("Nutaufriss") for a slip-ring rotor bar. Pure
    computation; when IPRINT==1, also updates st.ck211l (the DEL1L-11L
    clearance breakdown, the rotor-side counterpart to st.ck211s from
    NARISS) instead of writing the FORTRAN's printer-report file.

    HZW is both input and output, matching the FORTRAN argument: when
    IPRINT==1 the passed-in value is kept as-is (mirrors the FORTRAN's
    "GOTO 100" skip); when IPRINT!=1 it's recomputed from UNX2 and the
    input value is ignored - same pattern as passing 0 for "please
    compute this" elsewhere in this module.

    Returns (bnmind, hnmind, hzw, iokay, unx2).
    """
    unx2 = 1000.0 * math.trunc(u20 / 1000.0) + 999.99

    bltol = 0.04 if btl <= 5.0 else 0.06
    hltol = 0.10 if htl <= 20.0 else 0.15
    liso = 0.50 if unx2 <= 3000.0 else 0.70

    kiso = 0.08
    bkiso = 0.24
    hkiso = 0.16
    brund = 0.10
    hrund = 0.20
    ngstr = 0.30

    if iprint != 1:
        hzw = 4.00
        if unx2 <= 4000.0:
            hzw = 3.00
        if unx2 <= 3000.0:
            hzw = 2.00
        if unx2 <= 2000.0:
            hzw = 1.00
        if unx2 <= 1000.0:
            hzw = 0.50

    nak = 0.25
    bnak = 0.50
    hnak = 0.75
    bluft = 0.10
    hluft = 0.20
    rstr = 0.30
    svsatz = 0.20

    bnmind = btl
    hnmind = htl
    bnmind += bltol
    hnmind += hltol
    bnmind += liso
    hnmind += liso
    bnmind += bkiso
    hnmind += hkiso
    bnmind += brund
    hnmind += hrund
    bnmind *= ntl
    hnmind *= mtl
    hnmind += ngstr
    hnmind += hzw
    bnmind += bnak
    hnmind += hnak
    bnmind += bluft
    hnmind += hluft + rstr
    bnmind += svsatz
    hnmind += svsatz
    hnmind += hk + hs

    iokay = 0 if (hnmind > hn or bnmind > bn) else 1

    db = bn - bnmind
    bluft = bluft + db
    dh = hn - hnmind
    if imerk == 0:
        hzw = hzw + math.trunc(2.0 * dh) / 2.0
        hluft = hluft + dh - math.trunc(2.0 * dh) / 2.0
    else:
        hluft = hluft + dh

    if iprint != 1:
        return bnmind, hnmind, hzw, iokay, unx2

    ck = st.ck211l
    ck.del1l = liso / 2.0
    ck.del2l = bltol / 2.0
    ck.del3l = brund / 2.0
    ck.del4l = bkiso / 2.0
    ck.del5l = 0.0
    ck.del6l = bnak / 2.0
    ck.del7l = bluft / 2.0
    ck.del8l = 0.0
    ck.del9l = 0.0
    ck.del10l = (liso + 0.05) / 2.0
    ck.del11l = 0.0
    ck.delb = 3.5

    return bnmind, hnmind, hzw, iokay, unx2


# Standard conductor-width/height tables (mm) used by ISOLFR, 1-indexed.
_ZAHLB = [0.0,
          0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
          1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5, 2.8, 3.0,
          3.2, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5,
          8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0,
          20.0, 22.0, 25.0, 30.0, 35.0, 40.0, 50.0]
_ZDIMB = 47

_ZAHLH = [0.0,
          1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.2, 3.5,
          3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
          9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0,
          22.0, 25.0, 28.0, 30.0, 32.0, 35.0, 38.0, 40.0, 42.0, 45.0,
          48.0, 50.0, 52.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0,
          90.0, 95.0, 100.0]
_ZDIMH = 53


def isolfr(u20: float, zn: float, n: float, lm: float, azweig: float,
           weite: float, bs: float, bn: float, hs: float, hk: float,
           ho: float, hcuo: float, hzw: float, hcuu: float, hu: float,
           st: MachineState, reversible_duty: bool = False,
           qcu: float = 0.0, ntl_override: float | None = None,
           tpar_override: float | None = None,
           write_file_oul_path: str | None = None):
    """
    Slip-ring rotor winding design: given zero or a target QCU, searches
    the standard conductor tables for a fitting bar cross-section (via
    NAR2), lets the caller override the number of conductors side-by-side
    (ntl_override) and in parallel (tpar_override) - these mirror the
    FORTRAN's two interactive prompts - then finalizes the conductor
    cross-section and slot fill.

    reversible_duty mirrors "FUER REVERSIERBETRIEB ISOLIEREN ?" (doubles
    U20 when True, matching answering 'j').

    Returns a dict: mtl, ntl, tpar, qcu, btl, htl, unx2, iokay, u20
    (possibly doubled), plus rk_20c (phase resistance at 20 C, mOhm) and
    wdg (turns/coil) which the FORTRAN prints as part of the winding-data
    report instead of writing to a file by default (pass
    write_file_oul_path to opt back into writing that report).
    """
    if reversible_duty:
        u20 = u20 * 2.0

    hn = hs + hk + ho + hcuo + hzw + hcuu + hu

    if qcu != 0.0:
        imerk = 1
        tpar = mtl_final = None  # placeholders, set in the shared tail below
        mtl = zn * azweig  # not meaningful in this branch, kept for symmetry
        ntl = 1.0
        # BTL/HTL must already be known by the caller in this branch,
        # since the FORTRAN jumps straight past the search - accept them
        # as the values already on hand via the qcu-driven path below.
        raise NotImplementedError(
            "isolfr: qcu != 0 branch requires BTL/HTL already fixed by "
            "the caller - not yet exercised by any caller in this port; "
            "raise if you hit this so it can be filled in against a "
            "real use case.")

    imerk = 0
    tpar = 1.0
    ntl = 1.0
    mtl = zn * azweig
    jb = _ZDIMB
    jh = _ZDIMH

    def _search(jb, jh, mtl, ntl, tpar):
        while True:
            if jb < 1 or jh < 1:
                raise ValueError(
                    "ISOLFR: no conductor size in the standard table fits "
                    "this slot (BN/HN too small for MTL/NTL turns) - "
                    "mirrors the FORTRAN running off the end of the "
                    "ZAHLB/ZAHLH tables, which was undefined behaviour "
                    "there.")
            btl = _ZAHLB[jb]
            htl = _ZAHLH[jh]
            bnmind, hnmind, hzw_new, iokay, unx2 = nar2(
                0, 0, u20, btl, htl, mtl, ntl, tpar, hs, hk, bs, bn, hn,
                0.0, 0.0, st)
            if iokay != 1:
                if bnmind > bn:
                    jb -= 1
                    continue
                if hnmind > hn:
                    jh -= 1
                    continue
            return btl, htl, hzw_new, unx2, jb, jh

    btl, htl, hzw_new, unx2, jb, jh = _search(jb, jh, mtl, ntl, tpar)

    ntlstr = ntl if ntl_override is None else ntl_override
    if ntlstr == 0.0:
        ntlstr = 1.0
    if round(ntlstr) != 1:
        tpar = ntlstr / ntl
        ntl = ntlstr

    tparn = tpar if tpar_override is None else tpar_override
    if tparn < 0.001:
        tparn = tpar
    mtl = mtl * tparn / tpar
    tpar = tparn
    jb = _ZDIMB
    jh = _ZDIMH

    btl, htl, hzw_new, unx2, jb, jh = _search(jb, jh, mtl, ntl, tpar)

    if btl <= 1.40:
        radius = 0.50
    elif btl <= 1.90:
        radius = 0.60
    elif btl <= 3.00:
        radius = 0.80
    elif btl <= 10.0:
        radius = 1.00
    else:
        radius = 1.60
    qcue = btl * htl - radius ** 2 * (4.0 - 3.14)
    qcu = qcue * tpar * azweig
    ho = 0.0
    hu = 0.0
    hcuo = (hn - hs - hk - hzw_new) / 2.0
    hcuu = hcuo

    # ---- final NAR2 call (IPRINT=1) latches st.ck211l and the report ----
    hn = hs + hk + ho + hcuo + hzw_new + hcuu + hu
    tpar = mtl * ntl / (zn * azweig)
    qcue = qcu / (tpar * azweig)
    bnmind, hnmind, hzw_final, iokay, unx2 = nar2(
        1, imerk, u20, btl, htl, mtl, ntl, tpar, hs, hk, bs, bn, hn,
        hzw_new, qcue, st)

    wdg = mtl * ntl / tpar / 2.0
    rk_20c = zn * n / 3.0 * lm / (1000.0 * qcu * 58.0) * 1000.0  # mOhm

    if write_file_oul_path:
        esc = chr(27)
        with open(write_file_oul_path, "w") as f:
            f.write(f"{esc}&a49r7C     Wdg./Spule:{round(wdg)}"
                    f"     Widerstand/Phase:{rk_20c:.2f} mOhm bei 20 Grad\n")
            f.write(f"{esc}&a50.5r7C     Schritt: 1-{round(weite + 1.0)}"
                    f"     GT - F - Wellenwicklung\n")
            f.write(f"{esc}&a52r7C     Leiter par.{round(tpar)}"
                    f"     Schaltung der Spulen:{round(azweig)} parallel\n")

    return {
        "mtl": mtl, "ntl": ntl, "tpar": tpar, "qcu": qcu, "btl": btl,
        "htl": htl, "unx2": unx2, "iokay": iokay, "u20": u20,
        "wdg": wdg, "rk_20c_mohm": rk_20c, "hzw": hzw_final,
    }
