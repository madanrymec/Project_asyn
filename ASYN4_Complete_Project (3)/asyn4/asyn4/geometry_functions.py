"""
Second batch of ported FORTRAN subroutines from ASYN4.FOR - core geometry
and winding pre-calculation group.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    BRUVOR (5008-5125)  Vorberechnungen fuer den Zonenfaktor
                        (fractional-slot winding pre-calc; sets COMMON /BRULO/)
    COIL   (5126-5156)  Coilbreite / max. Aussendurchmesser (frame-size lookup table)
    NUT    (7635-7663)  Nutflaeche und theoretischer Kupferfuellfaktor
    PAKET  (7664-7723)  Blechpaketaufteilung nach Beneke
    SPAK   (8732-8788)  Blechpaketaufteilung nach VAX-Programm (Roehl)
    VJ     (9518-9542)  Magnetische Jochspannung (uses RMAG from winding_functions)
    UMRI   (9495-9517)  Magnetisch wirksamer Luftspalt fuer Oberschwingungsstrom
    VLZFE  (9543-9574)  Eisenverlustziffer bei ungleich 50 Hz
    K1K2   (7233-7261)  Stromverdraengungsfaktoren (linear interpolation)
    NUE0   (7566-7615)  kleinste Polpaarzahl fuer synchrones Moment
    MNYQMP (7549-7565)  Hilfsgroesse fuer synchrone Momente (uses winding_functions)

Not yet ported from this functional area: KIKR, KIKR1, KIKRSL, KIKRSX, STD,
SLFR, KREIS, STEG, MKWERT (bigger, more interdependent subroutines - next
batch).
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState
from .winding_functions import rmag, xsiz, xsis, xsisch

PI = 3.14159
PI_BRUVOR = 3.1415927


def _rkz_char(rkz: str, pos: int) -> str:
    """FORTRAN 1-based substring RKZ(pos:pos) -> single Python character."""
    if pos - 1 < len(rkz):
        return rkz[pos - 1]
    return " "


def bruvor(rn1: float, rp: float, mzone1: float, rm1: float, weite1: float,
           st: MachineState):
    """
    Pre-calculation for the (stator) zone factor, especially for
    fractional-slot ("Bruchloch") windings. Populates st.brulo in place.

    Returns (mz1s, weite1, a1max) - weite1 may be corrected (a message is
    printed, mirroring the FORTRAN PRINT 9000, if so).
    """
    from .winding_functions import ggt as _ggt  # the tested implementation

    rkz = st.art.rkz
    b = st.brulo

    w1alt = weite1
    b.ibruch = 0
    b.zv = 0.0
    n1 = int(rn1 + 0.1)
    p = int(rp + 0.1)
    m1 = int(rm1 + 0.1)

    c9 = _rkz_char(rkz, 9)
    if c9 in ("1", "2", "3", "4"):
        mz1s = mzone1
        a1max = 2.0 * rp
        b.zv = float(int(c9))
        _maybe_print_correction(weite1, w1alt)
        return mz1s, weite1, a1max

    rq1 = rn1 / (2.0 * rp * rm1)
    if math.fmod(rq1 + 0.001, 1.0) < 0.01:
        if c9 == "E":
            weite1 = rn1 / (2.0 * rp)
            mz1s = mzone1
            a1max = rp
        else:
            mz1s = mzone1
            if rm1 + 0.1 > mzone1:
                a1max = rp / 2.0
            else:
                a1max = 2.0 * rp
        _maybe_print_correction(weite1, w1alt)
        return mz1s, weite1, a1max

    b.ibruch = 1
    if c9 == "E":
        # ---- verdoppelte Polzahl / Schleifringlaeufer-Zweig (source label 500) ----
        rt = _ggt(rn1 / 2.0, rp)
        t = int(rt + 0.1)
        n1str = n1 // t
        pstr = p // t
        if n1str % 4 == 0:
            q1 = n1str // (4 * m1)
            q2 = q1
            a1max = 2.0 * t
        else:
            q1 = (n1str // 2 + m1) // (2 * m1)
            q2 = q1 - 1
            a1max = 1.0 * t
        yk = None
        for g in range(0, 1001):
            if ((g * n1str) // 2 + 1) % pstr == 0:
                yk = ((g * n1str) // 2 + 1) // pstr
                break
        if yk is None:
            raise ValueError("BRUVOR: YK search (label 610) did not converge")
        mz1s = m1 * t / float(p)
        if q1 == q2:
            mz1s = 2.0 * mz1s
        b.bruhi1 = (PI_BRUVOR * 2 * yk * q1) / float(n1)
        b.bruhi2 = (PI_BRUVOR * yk) / float(t)
        b.bruhi3 = (PI_BRUVOR * 2 * yk * q2) / float(n1)
        b.bruhi4 = float(q1 + q2)
        b.bruhi5 = (PI_BRUVOR * 2 * yk) / float(n1)
        ihg = n1 // (2 * p)
        if ihg % 2 == 0:
            weite1 = ihg + 1.0
        else:
            weite1 = float(ihg)
        _maybe_print_correction(weite1, w1alt)
        return mz1s, weite1, a1max

    # ---- normal fractional-slot branch (source label 1000 via fallthrough) ----
    rt = _ggt(rn1, rp)
    t = int(rt + 0.1)
    n1str = n1 // t
    pstr = p // t
    if n1str % 2 == 0:
        q1 = n1str // (2 * m1)
        q2 = q1
        a1max = 2.0 * t
    else:
        q1 = (n1str + m1) // (2 * m1)
        q2 = q1 - 1
        a1max = 1.0 * t
    yk = None
    for g in range(0, 1001):
        if (g * n1str + 1) % pstr == 0:
            yk = (g * n1str + 1) // pstr
            break
    if yk is None:
        raise ValueError("BRUVOR: YK search (label 110) did not converge")
    mz1s = m1 * t / float(p)
    if q1 == q2:
        mz1s = 2.0 * mz1s
    b.bruhi1 = (PI_BRUVOR * yk * q1) / float(n1)
    b.bruhi2 = (PI_BRUVOR * yk) / float(t)
    b.bruhi3 = (PI_BRUVOR * yk * q2) / float(n1)
    b.bruhi4 = float(q1 + q2)
    b.bruhi5 = (PI_BRUVOR * yk) / float(n1)
    _maybe_print_correction(weite1, w1alt)
    return mz1s, weite1, a1max


def _maybe_print_correction(weite1: float, w1alt: float) -> None:
    if abs(weite1 - w1alt) > 0.01:
        print(f"\n WEITE1 AUF{int(weite1 + 0.001):3d} KORRIGIERT !\n")


# Frame-size (IEC "neue Reihe") coil-width / max. outer-diameter lookup table.
_COIL_TABLE = {
    "500": (800.0, 840.0),
    "560": (900.0, 950.0),
    "630": (1000.0, 1060.0),
    "710": (1120.0, 1180.0),
    "800": (1250.0, 1320.0),
}


def coil(typ: str, bcoil: float | None = None, da1max: float | None = None):
    """
    Coil width and max. outer diameter for the "neue Reihe" frame series.
    TYP(1:2) must be 'A5' and TYP(5:7) one of 500/560/630/710/800.
    Returns (bcoil, da1max), unchanged from the passed-in defaults if no
    table entry matches (mirrors the FORTRAN, which leaves BCOIL/DA1MAX
    untouched when none of the IF blocks fire).
    """
    if typ[0:2] == "A5" and typ[4:7] in _COIL_TABLE:
        bcoil, da1max = _COIL_TABLE[typ[4:7]]
    return bcoil, da1max


def nut(bs: float, bn: float, bns: float, hk: float, ho: float, hcuo: float,
        hzw: float, hcuu: float, hu: float, qcu: float, zn: float):
    """Slot cross-section area and theoretical copper fill factor."""
    anut = (bn + bns) / 2.0 * (ho + hcuo + hzw + hcuu + hu)
    if bs < 0.9 * bn:
        anut = anut + (bs + bn) / 2.0 * hk

    hn = ho + hcuo + hzw + hcuu + hu
    hg = (bns - bn) / hn
    bnb = bn + hg * (ho + hcuo)
    bnc = bn + hg * (ho + hcuo + hzw)

    # "NEU: auf Wunsch von Frau Timonski" - the second ACU assignment in the
    # FORTRAN overwrites the first; only the final formula is kept here.
    acu = (bn + bnb) / 2.0 * (ho + hcuo) + (bnc + bns) / 2.0 * (hcuu + hu)
    if bs < 0.9 * bn and ho < 0.0001:
        acu = acu + (bs + bn) / 2.0 * hk

    kcuth = zn * qcu / acu
    return anut, kcuth


def spak(le1: float, nk1: float):
    """
    Lamination-stack split per the old VAX program (Roehl).
    Returns (lpak, lepak) where lpak is a 1-indexed array (size 24, index 0 unused)
    of up to 23 slots, matching FORTRAN LPAK(23).
    """
    lpak = np.zeros(24)
    lpakm = np.zeros(24)
    lerest = np.zeros(24)
    lkeil = np.zeros(10)  # LKEIL(9), 1-indexed

    npak = nk1 + 1
    mitte = npak - (npak // 2) * 2
    npak = (npak + 1) // 2
    npak = int(npak)

    nkeil = 1
    lepak = 0.0
    n = 1
    for n in range(1, npak + 1):
        if n == npak:
            if nk1 == 1.0:
                lpakm[npak] = le1 / (nk1 + 1)
            if nk1 > 1.0:
                lpakm[npak] = lerest[n - 1] * 0.5
            lepak = lpakm[npak]
            lpak[n] = lepak
            lkeil[nkeil + 1] = int((lepak + 9.9) * 0.2) * 5
            if lkeil[nkeil] != lkeil[nkeil + 1]:
                nkeil += 1
            lerest[npak] = lkeil[nkeil] - lpakm[npak]
        elif n > 1:
            nrest = nk1 + mitte + 3 - 2 * n
            lpakm[n] = lerest[n - 1] / nrest
            if not (lpakm[n] - lpakm[1] < 5.0 or (n < 3 and n >= npak - 1)):
                nkeil += 1
                lkeil[nkeil] = lkeil[nkeil - 1] + 5
            lpak[n] = lkeil[nkeil]
            lerest[n] = lerest[n - 1] - 2.0 * lpak[n]
        else:  # n == 1
            lpakm[1] = le1 / (nk1 + 1)
            nkeil = 1
            lkeil[1] = int(lpakm[1] * 0.2) * 5
            lpak[1] = lkeil[nkeil]
            lerest[1] = le1 - (2 - mitte) * lpak[1]

    # Korrektur, wenn Teilpakete nicht monoton steigend
    hg = lepak - lpak[npak]
    for _k in range(1, 21):
        idelta = np.zeros(24)
        for j in range(1, npak):
            idelta[j] = lpak[j + 1] - lpak[j]
        corrected = False
        for j in range(npak - 1, 0, -1):
            if idelta[j] < 0:
                lpak[j + 1] += 5
                lpak[j] -= 5
                corrected = True
                break
        if not corrected:
            break
    lepak = lpak[npak] + hg
    return lpak, lepak


def paket(reke: str, l_: float, bk: float, nk1: float, nk2: float, st: MachineState):
    """Lamination-stack split (Beneke method). Populates st.cpaket.paket1/paket2."""
    cp = st.cpaket
    cp.paket1[:] = 0.0
    cp.paket2[:] = 0.0

    ink1 = round(nk1)
    ink2 = round(nk2)
    if ink1 == 0:
        return

    r = 0.0 if ink1 <= 9 else 1.0
    bms = ((1000.0 * l_ - nk1 * bk * 1000.0) / (nk1 + 1.0) + bk * 1000.0 + r) / 5.0
    bms = math.trunc(bms)
    bm = 5.0 * bms - bk * 1000.0
    be = (l_ * 1000.0 - nk1 * bk * 1000.0 - (nk1 - 1.0) * bm) / 2.0
    ihg = ink1 // 2 + 1

    for j in range(1, ihg):
        cp.paket1[j] = bm
    cp.paket1[ihg] = be

    if reke[0:2] != "A5":
        lpak, lepak = spak(1000.0 * (l_ - nk1 * bk), nk1)
        for j in range(1, ihg):
            cp.paket1[j] = float(lpak[j])
        cp.paket1[ihg] = lepak

    if ink1 == ink2:
        for j in range(1, ihg + 1):
            cp.paket2[j] = cp.paket1[j]
        return

    hg = 1000.0 * (l_ - nk2 * bk) / 2.0
    if ink2 % 2 == 0:
        hg = hg + cp.paket1[1] / 2.0
    ihg2 = ink2 // 2 + 1
    for j in range(1, ihg2):
        cp.paket2[j] = cp.paket1[j]
        hg = hg - cp.paket2[j]
    cp.paket2[ihg2] = hg


def vj(bj: float, lj: float, p: float) -> float:
    """Magnetic yoke (Joch) voltage drop."""
    if p < 1.001:
        c1, c2 = 0.500, 0.700
    else:
        c1, c2 = 0.500, 0.550
    hg1 = 100.0 * rmag(bj * 10000.0)
    hg2 = 100.0 * rmag(c1 * bj * 10000.0)
    return lj * ((1.0 - c2) * hg1 + c2 * hg2)


def umri(imue: float, st: MachineState) -> float:
    """Magnetically effective air gap for the harmonic-current system."""
    k = st.kennli
    for j in range(2, k.eanz + 1):
        if 1.5 * imue <= k.mkpool[j, 2]:
            return (k.mkpool[j - 1, 3] +
                    (k.mkpool[j, 3] - k.mkpool[j - 1, 3]) /
                    (k.mkpool[j, 2] - k.mkpool[j - 1, 2]) *
                    (1.5 * imue - k.mkpool[j - 1, 2]))
    return k.mkpool[k.eanz, 3]


def vlzfe(v10: float, fn: float) -> float:
    """Iron-loss coefficient correction for frequency != 50 Hz."""
    if v10 <= 1.25:
        kh = 0.900
    elif v10 <= 1.50:
        kh = 0.875
    elif v10 <= 1.70:
        kh = 0.850
    elif v10 <= 2.00:
        kh = 0.825
    elif v10 <= 2.60:
        kh = 0.800
    elif v10 <= 3.00:
        kh = 0.750
    elif v10 <= 3.60:
        kh = 0.667
    elif v10 <= 4.20:
        kh = 0.600
    else:
        kh = 0.500
    sigh = kh * v10
    sigw = v10 - sigh
    return sigh * (fn / 50.0) + sigw * (fn / 50.0) ** 2


def k1k2(s: float, st: MachineState):
    """Linear interpolation of the current-displacement (skin-effect) factors."""
    sv = st.kennli.svpool
    ss = abs(s)
    for j in range(2, 26):
        if ss <= sv[j, 1]:
            hg = (ss - sv[j - 1, 1]) / (sv[j, 1] - sv[j - 1, 1])
            k1 = sv[j - 1, 2] + (sv[j, 2] - sv[j - 1, 2]) * hg
            k2 = sv[j - 1, 3] + (sv[j, 3] - sv[j - 1, 3]) * hg
            k1r = sv[j - 1, 4] + (sv[j, 4] - sv[j - 1, 4]) * hg
            k1ra = sv[j - 1, 5] + (sv[j, 5] - sv[j - 1, 5]) * hg
            return k1, k2, k1r, k1ra
    hg = (ss - sv[24, 1]) / (sv[25, 1] - sv[24, 1])
    k1 = sv[24, 2] + (sv[25, 2] - sv[24, 2]) * hg
    k2 = sv[24, 3] + (sv[25, 3] - sv[24, 3]) * hg
    k1r = sv[24, 4] + (sv[25, 4] - sv[24, 4]) * hg
    k1ra = sv[24, 5] + (sv[25, 5] - sv[24, 5]) * hg
    if k2 < 1.0e-6:
        k2 = 1.0e-6
    return k1, k2, k1r, k1ra


def nue0(p: float, mz1s: float, n2: float) -> float:
    """Smallest pole-pair number NY0 satisfying the synchronous-torque condition."""
    in2 = int(n2 + 0.001)
    ip = int(p + 0.001)
    for g1 in range(1, 301):
        for gg1 in (1, 2):
            hg = p * (1.0 + mz1s * g1 * (-1) ** gg1)
            ny = int(hg + 0.001) if hg >= 0.0 else int(hg - 0.001)
            if abs(ny + ip) > 5 * in2:
                continue
            for g2 in range(1, 6):
                for gg2 in (1, 2):
                    my = ip + in2 * g2 * (-1) ** gg2
                    if ny + my == 0:
                        return float(ny)
    return 1.01e6


def mnyqmp(p: float, nyq: float, mzone1: float, n1: float, weite1: float,
           schr: float, n2: float, st: MachineState) -> float:
    """Auxiliary quantity for the synchronous-torque calculation."""
    q1 = n1 / (p * mzone1)
    return ((p / nyq) ** 2 *
            xsiz(nyq, q1, n1, st) * xsis(nyq, weite1, n1) /
            (xsiz(p, q1, n1, st) * xsis(p, weite1, n1)) *
            xsisch(nyq, schr) / xsisch(p, schr) *
            math.sin(nyq * PI / n2) / math.sin(p * PI / n2))
