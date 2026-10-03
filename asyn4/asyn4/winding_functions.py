"""
Direct ports of the small, mostly-standalone FORTRAN FUNCTIONs from ASYN4.FOR.

Each function below is annotated with the FORTRAN source line range it was
ported from, for traceability. These were chosen first because they have
few or no COMMON-block dependencies, so they can be tested in isolation
before the big interdependent subroutines (LESEN, KREIS, SLAST, HOLAUF, ...)
are ported.

Ported so far:
    XSIS    (source 9697-9709)  Sehnungswicklungsfaktor  (pitch factor)
    XSISCH  (source 9710-9726)  Schraegungsfaktor        (skew factor)
    XSIZ    (source 9727-9753)  Zonenwicklungsfaktor Staender (stator zone factor)
    XSIZB   (source 9754-9770)  Zonenfaktor bei Bruchlochwicklungen (fractional slot)
    XSIZ2   (source 9771-9783)  Zonenwicklungsfaktor Laeufer (rotor zone factor)
    ZFAK    (source 9784-9822)  Zonenfaktor bei Strangverschachtelung
    SGRUND  (source 8201-8213)  Grundfeldschlupf aus Oberfeldschlupf
    SNUE    (source 8467-8478)  Oberfeldschlupf der Polpaarzahl NY
    RMAG    (source 8164-8200)  reduzierte Rueckenfeldstaerke (B-H interpolation)
    HEISEN  (source 5421-5453)  Schnittpunkt Magnetisierungskurve / Gerade
    GGT     (source 5336-5355)  groesster gemeinsamer Teiler (GCD)
    DXS1    (source 5228-5232)  Saettigung Nutkeil
    LASTG   (source 7262-7289)  lineare Interpolation KL-Stegstreuleitwert
    MGEGEN  (source 7394-7496)  Gegenmoment-Verlauf (counter-torque curve)

Still to port (subsequent batches): LESEN, BRUVOR, KREIS, STD, SLFR, KIKR*,
POPAZA, MKWERT, SLAST, HOLAUF, KDRUCK/PDRUCK/TLAST/SYNMOM output routines,
ISOSTD/ISOTR/ISOLFR/NAR2/NARISS/WVAR insulation-fill routines, the plotting
driver (PLOT/AXIS/SYMBOL/CIRCL/...), and PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState

PI = 3.14159          # FORTRAN uses this literal (not full double precision pi) in most places
PI_ZFAK = 3.1415927    # ZFAK uses a slightly more precise literal in the original


def xsis(ny: float, weite: float, n: float) -> float:
    """Pitch (chording) winding factor."""
    val = math.sin(ny * PI * weite / n)
    if val == 0.0:
        val = 1.0e-6
    return val


def xsisch(ny: float, schr: float) -> float:
    """Skew factor."""
    if abs(schr) < 0.1:
        val = 1.0
    else:
        val = math.sin(ny * PI / schr) / (ny * PI / schr)
    if val == 0.0:
        val = 1.0e-6
    return val


def xsizb(ny: float, st: MachineState) -> float:
    """Zone factor for fractional-slot (Bruchloch) windings."""
    b = st.brulo
    jot = 1j
    cxsizb = (math.sin(b.bruhi1 * ny) - np.exp(jot * b.bruhi2 * ny) *
              math.sin(b.bruhi3 * ny)) / (b.bruhi4 * math.sin(b.bruhi5 * ny))
    return float(np.real(cxsizb))


def zfak(ny: float, q: float, n: float, zv: float) -> float:
    """Zone factor with phase-belt interleaving (Strangverschachtelung)."""
    alf = 2.0 * PI_ZFAK / n * ny
    hg1 = int(q - 2.0 * zv + 0.001)
    hg2 = int(zv + 0.001)

    # DALF is 1-indexed in FORTRAN; build with a leading dummy slot.
    total = hg1 + 2 * hg2
    dalf = np.zeros(max(total, 1) + 1)
    for j in range(1, hg2 + 1):
        dalf[j] = 2.0 * alf
    for j in range(hg2 + 1, hg2 + hg1):
        dalf[j] = alf
    for j in range(hg2 + hg1, hg2 + hg1 + hg2):
        dalf[j] = 2.0 * alf

    csum = 1.0 + 0.0j
    gam = 0.0
    for j in range(1, hg1 + 2 * hg2):
        gam += dalf[j]
        csum += np.exp(1j * gam)
    return float(abs(csum) / q)


def xsiz(ny: float, q1: float, n1: float, st: MachineState) -> float:
    """Stator zone (distribution) winding factor, including finite slot width."""
    b = st.brulo
    v = st.vorzei
    qq1 = q1 / v.ksys
    if b.ibruch == 0 and b.zv < 0.1:
        result = math.sin(ny * qq1 * PI / n1) / (qq1 * math.sin(ny * PI / n1))
    elif b.ibruch == 1:
        result = xsizb(ny, st)
    else:
        result = zfak(ny, qq1, n1, b.zv)
    result = result * math.sin(ny * v.bs1di1) / (ny * v.bs1di1)
    if abs(result) < 1.0e-8:
        result = 1.0e-8
    result = v.vzxsiz * result
    return result


def xsiz2(ny: float, q: float, n: float) -> float:
    """Rotor (slip-ring) zone winding factor."""
    val = math.sin(ny * q * PI / n) / (q * math.sin(ny * PI / n))
    if val == 0.0:
        val = 1.0e-6
    return val


def sgrund(sny: float, ny: float, p: float) -> float:
    """Fundamental-field slip corresponding to a given harmonic-field slip."""
    val = 1.0 - (1.0 - sny) * p / ny
    if val == 0.0:
        val = 1.0e-9
    return val


def snue(ny: float, p: float, s: float) -> float:
    """Harmonic-field slip for pole-pair number NY."""
    val = 1.0 - ny / p * (1.0 - s)
    if val == 0.0:
        val = 1.0e-9
    return val


_RMAG_B = np.array([0., 5000., 7000., 8000., 9000., 10000., 11000., 12000.,
                     13000., 14000., 15000., 16000., 17000., 18000., 19000.,
                     20000., 21000., 23000., 25000., 27000., 29000., 31000.])
_RMAG_H = np.array([0., 0.85, 1.3, 1.8, 2.18, 2.75, 3.45, 4.7, 6.3, 8.82,
                     12.15, 18.2, 28., 46., 70., 103.5, 155., 359., 620.,
                     980., 1700., 2900.])


def rmag(bjoch: float) -> float:
    """Reduced back-iron field strength, linearly interpolated from a B-H table
    (data from W. Nuernberg, 'Die Asynchronmaschine'). Units: Gauss, A/cm."""
    b, h = _RMAG_B, _RMAG_H
    if bjoch >= 31000.0:
        return (bjoch - 31000.0) / 1.2566 + h[21]
    for i in range(1, 22):  # FORTRAN I=2,22 over 1-indexed B/H -> here 1-indexed via numpy 0-based
        if bjoch <= b[i]:
            return (h[i] - h[i - 1]) / (b[i] - b[i - 1]) * (bjoch - b[i - 1]) + h[i - 1]
    return h[21]  # should not be reached, mirrors FORTRAN fallthrough safety


def heisen(ee: float, f: float, st: MachineState) -> float:
    """
    Intersection of the magnetization curve with the line B = EE - F*H.
    Requires st.cbsort.h (H(0:45)) and st.cbsort.my0 to be populated
    (this table is filled by ISOSTD/ISOTR-adjacent setup code, not yet
    ported - placeholder zeros will not give meaningful results until then).
    """
    h_tab = st.cbsort.h
    my0 = st.cbsort.my0
    e = abs(ee)
    ianf = int(e / 0.05)
    if e > 2.25:
        ianf = 45
    result = 0.0
    for i in range(ianf, -1, -1):
        if i >= 45:
            c = 2.25 - my0 * h_tab[45]
            d = my0
        else:
            d = 0.05 / (h_tab[i + 1] - h_tab[i])
            c = 0.05 * float(i) - d * h_tab[i]
        result = (e - c) / (f + d)
        if result >= h_tab[i]:
            break
    return math.copysign(result, ee)


def ggt(a: float, b: float) -> float:
    """Greatest common divisor of A and B (FORTRAN SUBROUTINE GGT, returned as C)."""
    ia = int(a + 0.001)
    ib = int(b + 0.001)
    imin, imax = min(ia, ib), max(ia, ib)
    ic = 1
    for j in range(2, imin + 1):
        if imax % j == 0 and imin % j == 0:
            ic = j
    return float(ic)


def dxs1(i1: float, st: MachineState) -> float:
    """Slot-wedge saturation factor."""
    k = st.keil
    return 9.0 * k.xskeil * math.exp(-i1 / k.teta0)


def lastg(i2: float, st: MachineState) -> float:
    """Linear interpolation of the KL slot-bridge leakage permeance vs. current."""
    inut = st.kennli.inut
    lamda = st.kennli.lamda
    i2d = 1.414 * abs(i2)
    if i2d <= inut[1]:
        return lamda[1]
    if i2d >= inut[19]:
        val = lamda[19] + (lamda[20] - lamda[19]) / (inut[20] - inut[19]) * (i2d - inut[19])
        if val <= 0.0:
            val = 1.0e-9
        return val
    for j in range(2, 20):
        if i2d < inut[j]:
            return lamda[j - 1] + (lamda[j] - lamda[j - 1]) / (inut[j] - inut[j - 1]) * (i2d - inut[j - 1])
    return 0.0  # unreachable if inut is monotonically increasing, mirrors FORTRAN


def mgegen(n: float, st: MachineState) -> float:
    """Counter-torque (load torque) curve, 4 selectable shapes via MG(0)."""
    mg = st.daten.mg  # 0-based, MG(0:15)
    img = int(mg[0] + 0.1)
    if img == 0:
        return 0.0

    ns = mg[3]
    if ns == 0.0:
        ns = 0.0001
    mgn = st.daten.mgn

    if img == 1:
        # konstantes Gegenmoment
        if n <= ns:
            val = (mg[1] - mg[4]) / ns ** 2 * (ns - n) ** 2 + mg[4]
        else:
            val = mg[4]
        return val * mgn

    if img == 2:
        # linear ansteigendes Gegenmoment
        if n <= ns:
            val = (mg[1] - ns * mg[4]) / ns ** 2 * (ns - n) ** 2 + ns * mg[4]
        else:
            val = n * mg[4]
        return val * mgn

    if img == 3:
        # quadratisch ansteigendes Gegenmoment
        ma, me, ms = mg[1], mg[4], mg[2]
        if n < ns:
            val = (ma - ms) / ns ** 2 * (ns - n) ** 2 + ms
        else:
            val = ms + (me - ms) * (n - ns) ** 2 / (1.0 - ns) ** 2
        return val * mgn

    if img == 4:
        # Gegenmoment in 11 Punkten (piecewise-linear table, MG(5..15))
        if n <= 0.0:
            j = 5
        elif n >= 1.0:
            j = 14
        else:
            j = int(n * 10.0) + 5
        val = (mg[j + 1] - mg[j]) / 0.1 * (n - (j - 5.0) / 10.0) + mg[j]
        return val * mgn

    raise ValueError(f"MGEGEN: unexpected MG(0) selector IMG={img}")
