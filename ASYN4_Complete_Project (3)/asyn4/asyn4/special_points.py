"""
Ninth batch of ported FORTRAN subroutines from ASYN4.FOR - special
operating points built on top of SPKT, plus input validation.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    SK     (8214-8268)  Iteratively finds the motor breakdown (pull-out)
                         slip SKIPP, starting from an ESB-derived estimate
                         and refining by scanning +/-nearby slip values
                         for the true torque maximum. Updates st.esb.skipp.
    SKGEN  (8269-8307)  Same idea for the generator-mode breakdown slip
                         (negative slip region).
    SKNY   (8308-8391)  Generates the list of slip values (st.svorga) the
                         load-curve sweep (SLAST, via SORT) should be
                         evaluated at: fine steps around the fundamental
                         breakdown slip, around each damped harmonic's own
                         breakdown slip, plus a coarse fixed 0.025 grid.
                         This is what populates st.svorga before SORT/
                         SLAST are called - discovered as a missing piece
                         during the SLAST/HOLAUF integration testing.
    PRAXIS (7934-8042)  "Praxisnah" (practical) starting/plugging-region
                         operating point for revolving-field-magnet
                         (Drehfeldmagnet) duty: like SPKT but starting
                         from a 5x-rated current estimate, with a
                         harmonic-torque correction scaled by (I1/I1TH)^2
                         from a pre-computed theoretical point.
    PRUEF  (8043-8163)  Validates the read-in machine data (st.daten /
                         st.art.rkz) for internal consistency. Returns
                         (iabort, messages) - messages is a list of
                         (is_error: bool, text: str) tuples instead of
                         FORTRAN WRITE statements to a report file.
    DFMAG  (5157-5227)  Revolving-field-magnet (Drehfeldmagnet) braking
                         behaviour: computes both the theoretical (SPKT)
                         and practical (PRAXIS) operating point at slip
                         S=1..5, for reporting.
    SPKTSP (8979-9157)  Locked-rotor (S=1) "special point" solver: like
                         SPKT but with a damping factor KKA that scales
                         each harmonic's contribution differently
                         depending on whether its pole-pair number is
                         close to the fundamental or not - used for
                         standstill/locked-rotor test-point calculations.
    GENERA (870-1171)   Generator-mode performance: iterates to 5
                         partial-load operating points (1/4, 2/4, 3/4,
                         rated, 5/4 x rated power, in that FORTRAN J
                         order) via SPKT, then finds the generator
                         breakdown (Kipp) point via SKGEN+SPKT. All the
                         original's screen-position/printer WRITE
                         statements were replaced with a returned dict
                         (points list + cfeld/efeld/pcu2nf summary arrays
                         + breakdown point), matching the "plain output
                         instead of terminal escape codes" pattern used
                         throughout. Note: COMMON /PINFO/ and /CNLAST/
                         (read by other, not-yet-ported report routines)
                         are not modelled as state.py dataclasses yet -
                         their values are included directly in GENERA's
                         return dict instead, and can be promoted to
                         proper COMMON-block dataclasses once those
                         report routines are ported and need them.

Not yet ported: output/reporting subroutines, plotting driver, LESEN,
PROGRAM ASYN4 itself. **This completes the "Performance calculations"
functional group.**
"""
from __future__ import annotations
import cmath
import math

from .state import MachineState
from .winding_functions import dxs1, lastg, sgrund, snue
from .geometry_functions import k1k2
from .operating_point import mkwert, spkt
from .harmonics import z2nys

PI = 3.14159
MY0 = 1.2566e-6


def sk(bs2i0: float, bs2in: float, ilfr: int, p: float, fn: float,
       r1w: float, xssti1: float, xsnut1: float, xssti2: float,
       xsnut2: float, r2w: float, li: float, hs2: float, lamnu2: float,
       r2s: float, r2kw: float, h42: float, br42: float, r2kwa: float,
       pn: float, u1: float, m1: float, prbg0: float, st: MachineState):
    """
    Iteratively finds the motor breakdown slip. Updates st.esb.skipp in
    place and also returns it.
    """
    esb = st.esb
    bs2i = bs2in
    xsig1 = esb.xsti1 + esb.xn1 + esb.xsd1
    hgr1 = esb.r1 / (xsig1 + esb.x1h)
    sig1 = xsig1 / esb.x1h
    r2 = esb.rk2 + esb.rr2 + esb.rs2
    xsig2 = (2.0 * esb.xsch + esb.xsd2 + esb.xring + (esb.xn2 - esb.xsteg) +
             esb.xsteg * bs2i0 / bs2i)
    sig2 = xsig2 / esb.x1h
    sig = 1.0 - 1.0 / ((1.0 + sig1) * (1.0 + sig2))
    skipp = (r2 / (sig * esb.x1h * (1.0 + sig2)) *
             math.sqrt((1.0 + hgr1 ** 2) / (1.0 + (hgr1 / sig) ** 2)))

    skanf = 0.5 * skipp
    dsk = 0.04 * skipp
    mkipp = 0.0

    for j in range(1, 51):
        s = skanf + j * dsk
        out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                   r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                   u1, m1, prbg0, 0.0, 0, st)
        if out["m"] > mkipp:
            mkipp = out["m"]
            skipp = s
        if out["m"] < mkipp - 0.01:
            break

    esb.skipp = skipp
    return skipp


def skgen(skipp: float, ilfr: int, p: float, fn: float, r1w: float,
          xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
          r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
          r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
          u1: float, m1: float, prbg0: float, st: MachineState) -> float:
    """Iteratively finds the generator-mode breakdown slip. Returns skippg."""
    skanf = -0.80 * skipp
    dsk = -0.02 * skipp
    mkipp = 0.0
    skippg = skanf

    for j in range(1, 26):
        s = skanf + j * dsk
        out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                   r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                   u1, m1, prbg0, 0.0, 0, st)
        if out["m"] < mkipp:
            mkipp = out["m"]
            skippg = s

    return skippg


def skny(ilfr: int, skipp: float, x1hp: float, hs2: float, bs2i: float,
         xssti2: float, xsnut2: float, r2w: float, li: float,
         lamnu2: float, r2s: float, r2kw: float, h42: float, br42: float,
         fn: float, p: float, x1hp0: float, st: MachineState):
    """
    Generates the slip values (st.svorga) that the load-curve sweep
    should be evaluated at: fine steps around the fundamental and each
    damped harmonic's breakdown slip, plus a coarse fixed grid.
    """
    nyinfo = st.nyinfo
    svorga = st.svorga
    hg = [1.70, 1.20, 1.00, 0.80, 0.30]

    sanz = 0
    svor = []

    for j in range(1, 9):
        sanz += 1
        svor.append(0.03125 * skipp * j)
    for j in range(1, 11):
        sanz += 1
        svor.append(0.075 * skipp * j + 0.25 * skipp)
    for j in range(1, 13):
        sanz += 1
        svor.append(0.4 * skipp * j + skipp)

    for j in range(2, nyinfo.nyanz + 1):
        x2hny = nyinfo.x1hbez[j] * x1hp / nyinfo.uer[j]
        x2hny0 = nyinfo.x1hbez[j] * x1hp0 / nyinfo.uer[j]
        xsh = (nyinfo.sd2[j] + (1.0 - nyinfo.sf[j])) * x2hny0
        if ilfr == 1:
            xs = xssti2 + xsnut2
            r2 = r2w
        elif ilfr == 2:
            xs = 2.0 * PI * fn * MY0 * li * (lamnu2 + hs2 / bs2i) + nyinfo.xsrny[j]
            r2 = r2s + r2kw + nyinfo.rrny[j]
        else:
            xs = (2.0 * PI * fn * MY0 * li *
                  (lamnu2 + h42 / br42 + 1.58 + hs2 / bs2i) + nyinfo.xsrny[j])
            r2 = r2s + r2kw + nyinfo.rrny[j]
        skj = r2 / (nyinfo.sf[j] * x2hny + xsh + xs)

        for jj in range(1, 6):
            for jjj in (1, 2):
                sny = hg[jj - 1] * (-1) ** jjj * skj
                s = sgrund(sny, nyinfo.ny[j], p)
                if s > 1.55 or s < 0.00:
                    continue
                if abs(math.fmod(s, 0.04999999)) < 0.0001:
                    continue
                sanz += 1
                svor.append(s)

    for j in range(0, 61):
        sanz += 1
        svor.append(0.025 * j)

    svorga.sanz = sanz
    for i, s in enumerate(svor, start=1):
        svorga.svor[i] = s

    return sanz


def praxis(ilfr: int, s: float, p: float, fn: float, r1w: float,
           xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
           r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
           r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
           u1: float, m1: float, x1hp0: float, mnyth: float, i1th: float,
           pc2nth: float, st: MachineState):
    """
    "Praxisnah" starting/plugging-region operating point for
    revolving-field-magnet duty: an SPKT-like iteration starting from a
    5x-rated current estimate, with the harmonic-torque contribution
    scaled from a pre-computed theoretical point (MNYTH/I1TH/PC2NTH,
    normally from a preceding SPKT call at the same slip).

    Returns a dict with i1, i2, i2o, i2u, imy, ui, pcu1, pcu2, pcu2ny, m,
    mny, n, cosph, bs2i.
    """
    nyinfo = st.nyinfo
    kka = 0.00
    jot = 1j
    om1 = 2.0 * PI * fn

    i1 = pn / (m1 * u1 * 0.9) * 5.0
    i2 = i1 / nyinfo.uei[1]
    dxs1n = dxs1(i1, st)
    lamstg = lastg(i2, st) if ilfr != 1 else None
    ui = 0.5 * u1

    bs2i = None
    xssteg = 0.0

    for _iter in range(1, 1000):
        uiv = ui
        dxs1v = dxs1n
        if ilfr != 1:
            lastgv = lamstg
            if hs2 < 1.0e-9:
                bs2i = st.daten.bs2
                xssteg = 0.0
            else:
                bs2i = hs2 / lamstg
                xssteg = om1 * MY0 * li * lamstg

        mk = mkwert(0, ui, st)
        imy = mk["imy"]
        x1hp = ui / imy

        xdvs1 = kka * nyinfo.sd1 * x1hp0
        xschr = kka * (1.0 - nyinfo.sf[1] ** 2) * x1hp0
        cz = r1w + jot * (xssti1 + xsnut1 + dxs1n + xdvs1 + xschr)
        cz2nys, cdny, i2opu, i2upu = z2nys(
            ilfr, s, p, fn, x1hp0, x1hp0, 1.0, nyinfo.uer[1], xssti2,
            xsnut2, nyinfo.sd2[1], r2w, xssteg, li, lamnu2, nyinfo.xsrny[1],
            nyinfo.rrny[1], r2s, r2kw, h42, br42, r2kwa, nyinfo.xsrnya[1],
            nyinfo.rrnya[1], 0, st)
        cz = cz + cz2nys
        ci1 = u1 / cz
        ui = abs(cz2nys * ci1)
        i1 = abs(ci1)
        i2 = abs((cdny - 1.0) * ci1) / nyinfo.uei[1]
        dxs1n = dxs1(i1, st)
        if ilfr != 1:
            lamstg = lastg(i2, st)

        if abs((ui - uiv) / ui) > 0.003:
            continue
        if abs((dxs1n - dxs1v) / (dxs1n + 1.0e-8)) > 0.05:
            continue
        if ilfr != 1:
            if lamstg > 1.0e-9 and abs(lamstg - lastgv) / lamstg > 0.05:
                continue
        break

    i1 = abs(ci1)
    i2o = i2 * i2opu
    i2u = i2 * i2upu
    m = m1 * p * x1hp0 / om1 * i1 ** 2 * (1.0 - cdny).imag
    n = fn / p * (1.0 - s)
    pcu1 = m1 * r1w * i1 ** 2
    pcu2 = 2.0 * PI * fn / p * m * s

    mny = mnyth * (i1 / i1th) ** 2
    pcu2ny = pc2nth * (i1 / i1th) ** 2
    m = m + mny
    pcu2 = pcu2 + pcu2ny
    ci1 = ci1 + mk["ife"]
    i1 = abs(ci1)
    cosph = ci1.real / abs(ci1)

    return {
        "i1": i1, "i2": i2, "i2o": i2o, "i2u": i2u, "imy": imy, "ui": ui,
        "pcu1": pcu1, "pcu2": pcu2, "pcu2ny": pcu2ny, "m": m, "mny": mny,
        "n": n, "cosph": cosph, "bs2i": bs2i,
    }


def dfmag(reke: str, ilfr: int, p: float, fn: float, r1w: float,
          xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
          r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
          r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
          u1: float, m1: float, x1h: float, st: MachineState,
          verbose: bool = True):
    """
    Revolving-field-magnet (Drehfeldmagnet) braking behaviour: computes
    both the theoretical (SPKT) and practical (PRAXIS) operating point at
    slip S=1..5. Returns a list of 5 dicts, each with 'theoretical' and
    'praxis' sub-dicts.
    """
    results = []
    if verbose:
        print(f"\n{reke}\nBETRIEB ALS DREHFELDMAGNET")

    for j in range(1, 6):
        s = float(j)
        out_th = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                      r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                      pn, u1, m1, 0.0, 0.0, 0, st)
        mnyth, i1th, pc2nth = out_th["mny"], out_th["i1"], out_th["pcu2ny"]
        k1, k2, k1r, k1ra = k1k2(s, st)

        if verbose:
            print(f"\nS={s:.2f}  THEORETISCHER WERT:")
            print(f"  M={out_th['m']:7.2f} NM  I1={out_th['i1']:7.1f} A  "
                  f"MNY={out_th['mny']:7.2f} NM  COSPH={out_th['cosph']:7.3f}  "
                  f"K1={k1:.2f} K2={k2:.2f} K1R={k1r:.2f}")

        out_pr = praxis(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                        r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                        pn, u1, m1, x1h, mnyth, i1th, pc2nth, st)
        if verbose:
            print(f"S={s:.2f}  SPECIAL-POINT:")
            print(f"  M={out_pr['m']:7.2f} NM  I1={out_pr['i1']:7.1f} A  "
                  f"MNY={out_pr['mny']:7.2f} NM  COSPH={out_pr['cosph']:7.3f}")

        results.append({"s": s, "theoretical": out_th, "praxis": out_pr})

    return results


def pruef(st: MachineState):
    """
    Validates the read-in machine data for internal consistency.
    Returns (iabort, messages) - messages is a list of (is_error, text).
    """
    d = st.daten
    rkz = st.art.rkz
    ilfr = int(d.ilfr + 0.1)
    iabort = 0
    messages = []

    def err(text):
        nonlocal iabort
        messages.append((True, text))
        iabort = 1

    def warn(text):
        messages.append((False, text))

    if ilfr == 1:
        if round(d.m1) != 3 or round(d.mzone1) != 6:
            err("M1 BZW. MZONE1 BEI SL UNZULAESSIG")

    if round(d.mzone1 / d.m1) not in (1, 2):
        err("MZONE1/M1 UNZULAESSIG")

    c9 = rkz[8] if len(rkz) > 8 else " "
    if c9 not in (" ", "E", "0", "1", "2", "3", "4"):
        err("UNZULAESSIGE 9. SPALTE DER RECHNUNGSKENNZIFFER")

    if d.lm1 < d.l:
        err("LM1 KLEINER L")

    if ilfr == 1:
        q2 = d.n2 / (6.0 * d.p)
        if math.fmod(q2 + 0.0001, 1.0) > 0.001:
            warn("WARNUNG: Q2 GEBROCHEN WIRD MIT GANZLOCHFORMELN GERECHNET")
        if d.lm2 < d.l:
            err("LM2 KLEINER L")
    else:
        if d.lstab < 0.999 * d.l:
            err("LSTAB KLEINER L")

    if ilfr == 4:
        if d.lstaba < 0.999 * d.l:
            err("LSTABA KLEINER L")

    if d.zk * d.bk > d.l / 3.0:
        warn("WARNUNG: ZK * BK GROESSER L/3")

    if d.kfe > 1.0:
        err("KFE GROESSER 1")

    if d.di1 > d.da1:
        err("DI1 GROESSER DA1")

    if d.di2 > d.di1:
        err("DI2 GROESSER DI1")

    c11 = rkz[10] if len(rkz) > 10 else " "
    if c11 != " ":
        warn("WARNUNG: SIMULATION EINES AUSGEFALLENEN SYSTEMS")

    if math.fmod(d.mtl + 0.03, 2.0) > 0.1 and round(d.izgf) == 1:
        err("TEILLEITER UEBEREINDANDER JE NUT NICHT GERADE ZAHL")

    if iabort == 1:
        messages.append((True, "*** DATENSATZ FEHLERHAFT ***"))

    return iabort, messages


def spktsp(ilfr: int, p: float, fn: float, r1w: float, xssti1: float,
           xsnut1: float, xssti2: float, xsnut2: float, r2w: float,
           li: float, hs2: float, lamnu2: float, r2s: float, r2kw: float,
           h42: float, br42: float, r2kwa: float, pn: float, u1: float,
           m1: float, kka: float, st: MachineState):
    """
    Locked-rotor (S=1) "special point" solver: like SPKT but with a
    damping factor KKA that scales each harmonic's contribution
    differently depending on whether its pole-pair number is close to
    the fundamental (|NY| < 2.01*P) or not - used for standstill/
    locked-rotor test-point calculations.

    Returns a dict with i1, i2, i2o, i2u, imy, ui, p1, pfe1, pcu1, pcu2,
    pcu2ny, m, mny, cosph, bs2i.
    """
    nyinfo = st.nyinfo
    nyanz = nyinfo.nyanz
    jot = 1j
    s = 1.0
    om1 = 2.0 * PI * fn

    i1 = pn / (m1 * u1 * 0.9)
    i2 = i1 / nyinfo.uei[1]
    dxs1n = dxs1(i1, st)
    lamstg = lastg(i2, st) if ilfr != 1 else None
    ui = 0.97 * u1

    x1hp0 = st.dp.x1hp0
    bs2i = None
    xssteg = 0.0

    for _iter in range(1, 1000):
        uiv = ui
        dxs1v = dxs1n
        if ilfr != 1:
            lastgv = lamstg
            if hs2 < 1.0e-9:
                bs2i = st.daten.bs2
                xssteg = 0.0
            else:
                bs2i = hs2 / lamstg
                xssteg = om1 * MY0 * li * lamstg

        mk = mkwert(0, ui, st)
        imy = mk["imy"]
        x1hp = ui / imy
        xschr = (1.0 - nyinfo.sf[1] ** 2) * x1hp0

        if nyinfo.ny[1] > 1.001:
            cz = r1w + jot * (xssti1 + xsnut1 + dxs1n +
                               kka * nyinfo.sd1r * x1hp0 + kka * xschr)
        else:
            cz = r1w + jot * (xssti1 + xsnut1 + dxs1n +
                               kka * nyinfo.sd1r * x1hp0 + 1.00 * xschr)

        cz2nys = cdny = None
        for j in range(nyanz, 0, -1):
            sny = snue(nyinfo.ny[j], p, s)
            x1hny = nyinfo.x1hbez[j] * x1hp
            x1hny0 = nyinfo.x1hbez[j] * x1hp0
            if j == 1:
                x1hny = x1hny0
                sd2s = nyinfo.sd2[j] * kka
            else:
                sd2s = nyinfo.sd2[j] * 0.00
            if nyinfo.ny[1] > 1.001:
                xsr = nyinfo.xsrny[j]
                xsra = nyinfo.xsrnya[j]
            else:
                xsr = nyinfo.xsrny[j] * kka
                xsra = nyinfo.xsrnya[j] * kka

            cz2nys, cdny, i2opu, i2upu = z2nys(
                ilfr, sny, nyinfo.ny[j], fn, x1hny, x1hny0, 1.00,
                nyinfo.uer[j], xssti2, xsnut2, sd2s, r2w, xssteg, li,
                lamnu2, xsr, nyinfo.rrny[j], r2s, r2kw, h42, br42, r2kwa,
                xsra, nyinfo.rrnya[j], 0, st)

            if abs(nyinfo.ny[j]) < 2.01 * p:
                cz = cz + cz2nys
            else:
                cz = cz + cz2nys / cdny * kka

        ci1 = u1 / cz
        ui = abs(cz2nys * ci1)
        i1 = abs(ci1)
        i2 = abs((cdny - 1.0) * ci1) / nyinfo.uei[1]
        dxs1n = dxs1(i1, st)
        if ilfr != 1:
            lamstg = lastg(i2, st)

        if abs((ui - uiv) / ui) > 0.01:
            continue
        if abs((dxs1n - dxs1v) / (dxs1n + 1.0e-8)) > 0.05:
            continue
        if ilfr != 1:
            if lamstg > 1.0e-9 and abs(lamstg - lastgv) / lamstg > 0.05:
                continue
        break

    i1 = abs(ci1)
    pcu2 = 0.0
    mi = 0.0
    mny = 0.0
    pcu2ny = 0.0

    for j in range(nyanz, 0, -1):
        sny = snue(nyinfo.ny[j], p, s)
        x1hny = nyinfo.x1hbez[j] * x1hp
        x1hny0 = nyinfo.x1hbez[j] * x1hp0
        if j == 1:
            x1hny = x1hny0
            sd2s = nyinfo.sd2[j] * kka
        else:
            sd2s = nyinfo.sd2[j] * 0.00
        if nyinfo.ny[1] > 1.001:
            xsr = nyinfo.xsrny[j]
            xsra = nyinfo.xsrnya[j]
        else:
            xsr = nyinfo.xsrny[j] * kka
            xsra = nyinfo.xsrnya[j] * kka

        cz2nys, cdny, i2opu, i2upu = z2nys(
            ilfr, sny, nyinfo.ny[j], fn, x1hny, x1hny0, 1.00, nyinfo.uer[j],
            xssti2, xsnut2, sd2s, r2w, xssteg, li, lamnu2, xsr,
            nyinfo.rrny[j], r2s, r2kw, h42, br42, r2kwa, xsra,
            nyinfo.rrnya[j], 0, st)

        if abs(nyinfo.ny[j]) > 2.01 * p:
            continue
        mny = m1 * nyinfo.ny[j] * x1hny / om1 * i1 ** 2 * (1.0 - cdny).imag
        pdny = om1 / nyinfo.ny[j] * mny
        pcu2ny = sny * pdny
        mi += mny
        pcu2 += pcu2ny

    ife = mk["ife"]
    ci1 = ci1 + ife
    i1 = abs(ci1)
    pcu1 = m1 * r1w * i1 ** 2
    i2o = i2opu * i2
    i2u = i2upu * i2
    pfe1 = m1 * u1 * ife
    p1 = pcu1 + pfe1 + pcu2
    cosph = p1 / (m1 * u1 * i1)
    pcu2ny = pcu2 - pcu2ny
    mny = mi - mny
    m = mi

    return {
        "i1": i1, "i2": i2, "i2o": i2o, "i2u": i2u, "imy": imy, "ui": ui,
        "p1": p1, "pfe1": pfe1, "pcu1": pcu1, "pcu2": pcu2,
        "pcu2ny": pcu2ny, "m": m, "mny": mny, "cosph": cosph, "bs2i": bs2i,
    }


def genera(reke: str, ilfr: int, p: float, fn: float, r1w: float,
           xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
           r2w: float, li: float, hs2: float, bs2: float, lamnu2: float,
           r2s: float, r2kw: float, h42: float, br42: float, r2kwa: float,
           u1: float, m1: float, prbg0: float, pn: float, pmsv: list,
           skipp: float, pfe0: float, kz: float, kt: float, qcu1: float,
           zn1: float, n1: float, di1: float, ishalt: int, ksys: float,
           st: MachineState, verbose: bool = True):
    """
    Generator-mode performance: 5 partial-load operating points (1/4,
    2/4, 3/4, rated, 5/4 x rated power - FORTRAN J=2,3,4,1,5 in that
    order) via SPKT, then the generator breakdown (Kipp) point via
    SKGEN+SPKT.

    pmsv is a 1-indexed list/array (length 6, index 1..5 used) of target
    powers, matching FORTRAN PMSV(5).

    Returns a dict: points (list of 5 per-load-point dicts, FORTRAN J
    order, each either the full operating-point data plus 'converged':
    True, or {'converged': False, 'j', 'psoll', 'pmax'} if that load
    point didn't converge in 300 iterations), cfeld/efeld/pcu2nf
    (1-indexed length-6 lists: index 1=1/4 load, 2=2/4, 3=3/4, 4=rated,
    5=5/4 load - matches the FORTRAN CFELD/EFELD/PCU2NF arrays), dtkz,
    dtkt, mngen, in_ (rated current), pzusn, breakdown (the generator
    Kipp-point dict).
    """
    from .special_points import skgen  # same module; explicit for clarity

    bs2i = bs2 if ilfr == 1 else None
    mngen = -1.0e15
    if verbose:
        print(f"\n{reke}\nG E N E R A T O R B E T R I E B")

    points = [None] * 6
    cfeld = [0.0] * 6
    efeld = [0.0] * 6
    pcu2nf = [0.0] * 6
    in_ = pzusn = dtkz = dtkt = None

    j_to_idx = {1: 4, 2: 1, 3: 2, 4: 3, 5: 5}

    for j in range(1, 6):
        p1max = 0.0
        s = -0.0005
        ds = -0.02 if abs(pn) <= 33000.0 else -0.005
        dalt = 1.0
        pms = -pmsv[j]

        converged = False
        out = None
        for _jj in range(1, 301):
            out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                       r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                       pn, u1, m1, prbg0, 0.0, 0, st)
            if out["p1"] < p1max:
                p1max = out["p1"]
            if abs((pms - out["p1"]) / pms) > 0.0040:
                dneu = (pms - out["p1"]) / pms
                if dneu * dalt < 0.0:
                    ds = -0.5 * ds
                s += ds
                dalt = dneu
                continue
            converged = True
            break

        if not converged:
            points[j] = {"converged": False, "j": j, "psoll": pmsv[j],
                         "pmax": abs(p1max)}
            if verbose:
                print(f"\n{int(4.01 * pmsv[j] / pn)}/4-LAST: KEINE KONVERGENZ "
                      f"(PSOLL={pmsv[j] / 1000.:.1f} KW, PMAX={abs(p1max) / 1000.:.1f} KW)")
            continue

        p1, i1, pcu1, pfe1 = out["p1"], out["i1"], out["pcu1"], out["pfe1"]
        if j == 1:
            in_ = i1
            pzusn = 0.005 * abs(p1)

        pdelta = p1 - pcu1 - pfe0
        pku2 = s * pdelta
        pzus = pzusn * (i1 / in_) ** 2
        eta = abs(p1) / (abs(p1) + pcu1 + pku2 + prbg0 + pfe0 + pzus)

        if j == 1:
            g1 = i1 / qcu1
            a1eff = zn1 * n1 * i1 / (PI * di1) / ksys
            dtkz = kz * g1 * a1eff * math.sqrt((pcu1 + pfe1) / pcu1)
            dtkt = kt * (2.0 * pcu1 + out["pcu2"] + pfe1)
            mngen = out["m"]
            st.laerm.nn = out["n"]
            st.laerm.in_ = i1
            st.laerm.cosfin = abs(out["cosph"])
            st.laerm.bs2in = out["bs2i"]
            st.laerm.in2 = out["i2"]
            if verbose:
                print(f"\nNENNLAST\nDTKZ ={dtkz:7.1f} K  DTKT ={dtkt:7.1f} K")
        elif verbose:
            print(f"\n{int(4.01 * pmsv[j] / pn)}/4 - LAST")

        result = dict(out)
        result.update({"converged": True, "j": j, "eta": eta, "pku2": pku2,
                        "pzus": pzus})
        points[j] = result

        idx = j_to_idx[j]
        cfeld[idx] = abs(out["cosph"])
        efeld[idx] = eta
        pcu2nf[idx] = out["pcu2ny"] / 1000.0

        if verbose:
            print(f"  P1={abs(p1) / 1000.:7.1f} KW  I1={i1:7.1f} A  "
                  f"S={s:7.4f}  COSPH={abs(out['cosph']):5.3f}\n"
                  f"  ETA={eta:7.4f}  M={abs(out['m']):7.0f} NM  "
                  f"N={out['n'] * 60.:7.1f} 1/MIN")

    skippg = skgen(skipp, ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                   r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                   u1, m1, prbg0, st)
    s = skippg
    if in_ is None or pzusn is None:
        raise RuntimeError(
            "GENERA: the rated-load point (J=1) never converged, so IN/"
            "PZUSN were never established. The FORTRAN would have silently"
            " used whatever garbage was already sitting in that memory for"
            " the breakdown-point PZUS calculation below - this port "
            "raises instead. This means PN/PMSV don't correspond to a "
            "reachable operating point for this machine's parameters; "
            "double check they're self-consistent.")
    bout = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
               li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1, m1,
               prbg0, 0.0, 0, st)
    pdelta = bout["p1"] - bout["pcu1"] - pfe0
    pku2 = s * pdelta
    pzus = pzusn * (bout["i1"] / in_) ** 2
    eta = abs(bout["p1"]) / (abs(bout["p1"]) + bout["pcu1"] + pku2 +
                              prbg0 + pfe0 + pzus)
    breakdown = dict(bout)
    breakdown.update({"s": s, "skippg": skippg, "eta": eta,
                       "m_over_mngen": bout["m"] / mngen})

    if verbose:
        print(f"\nKIPPUNKT\n  S={s * 100.:7.3f}%  M/MN(GEN)={abs(bout['m'] / mngen):5.2f}  "
              f"I1/IN={bout['i1'] / in_:5.2f}")

    return {
        "points": points[1:], "cfeld": cfeld, "efeld": efeld,
        "pcu2nf": pcu2nf, "dtkz": dtkz, "dtkt": dtkt, "mngen": mngen,
        "in_": in_, "pzusn": pzusn, "breakdown": breakdown,
    }
