"""
Fourth batch of ported FORTRAN subroutines from ASYN4.FOR - rotor-side
current-displacement (skin-effect) calculations.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    SLFR   (8400-8466)  Vorberechnungen beim Schleifringlaeufer
                         (slip-ring rotor geometry/leakage/resistance).
                         Updates st.xs2sl (XSN2/XSS2).
    KIKR1  (6012-6065)  Widerstandserhoehung bei Staenderformspulen fuer
                         Umrichter (stator form-coil skin-effect vs. drive
                         frequency table).
    KIKR   (5593-5894)  Numerische Stromverdraengung nach dem
                         Teilleiterverfahren (cage rotor bar skin-effect,
                         partial-conductor method) - the biggest subroutine
                         ported so far (~300 lines). Updates st.kennli.svpool.
    KIKRSL (5895-5991)  Stromverdraengung bei Schleifringlaeufern
                         (slip-ring rotor skin-effect). Updates
                         st.kennli.svpool (or a caller-supplied SVPOOL, see
                         note below).
    KIKRSX (5992-6011)  Hilfsroutine: sucht XSISTR mit gleichem K1-Wert wie
                         ein aequivalenter Hochstab (used by KIKRSL).

Note on SVPOOL: the FORTRAN program actually uses THREE different SVPOOL
arrays sharing the same name but bound to different storage depending on
call site (COMMON /KENNLI/ SVPOOL(25,5) for the main cage-rotor path, and
a separate local SVPOOL(25,5) passed as SVP1/SVPOOL2 for the stator
KIKR1 and slip-ring KIKRSL paths). To keep this explicit rather than
implicit (as the FORTRAN COMMON aliasing was), KIKR/UMRI/K1K2 use
st.kennli.svpool directly, while KIKR1 and KIKRSL return their own arrays
for the caller to store wherever it needs to (matching SVP1/SVPOOL in the
main program).

Not yet ported from this functional area: KREIS, MKWERT.
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState

PI = 3.14159
MUE0 = 1.2566e-6


def _rkz_char(rkz: str, pos: int) -> str:
    if pos - 1 < len(rkz):
        return rkz[pos - 1]
    return " "


def slfr(di2: float, hojoc2: float, n2: float, p: float, weite2: float,
         zn2: float, lm2: float, qcu2: float, fn: float, li: float,
         theta2: float, hs2: float, hk2: float, ho2: float, hcuo2: float,
         hzw2: float, hcuu2: float, hu2: float, bs2: float, bn2: float,
         bn2s: float, st: MachineState):
    """
    Slip-ring rotor pre-calculations: effective tooth length, tooth widths
    at 5 points, slot leakage permeance, slot/end-winding leakage
    reactances (updates st.xs2sl.xsn2/xss2), and ohmic resistance.

    Returns dict with hozah2, brzah2 (1-indexed length-6 array), lamnu2,
    r2k, r2w, xssti2, xsnut2.
    """
    hozah2 = ho2 + hcuo2 + hzw2 + hcuu2 + hu2
    brzah2 = np.zeros(6)
    for j in range(1, 6):
        brzah2[6 - j] = (PI * (di2 + 2.0 * hojoc2 + 2.0 * hozah2 * (j - 0.5) / 5.0) / n2 -
                         (bn2s + (bn2 - bn2s) * (j - 0.5) / 5.0))

    lamstg = hs2 / bs2 + hk2 / (bs2 + (bn2 - bs2) / 3.0)
    lamo = (0.33 * hcuo2 / (bn2 + (bn2s - bn2) * ho2 / hozah2) +
            ho2 / bn2 + lamstg)
    lamu = (0.33 * hcuu2 / (bn2 + (bn2s - bn2) * (ho2 + hcuo2 + hzw2) / hozah2) +
            (ho2 + hcuo2 + hzw2) / (bn2 + 0.33 * (bn2s - bn2) *
                                     (ho2 + hcuo2 + hzw2) / hozah2) + lamstg)
    lamou = (0.5 * hcuo2 / (bn2 + (bn2s - bn2) * ho2 / hozah2) +
             lamstg + ho2 / bn2)

    q2 = n2 / (6.0 * p)
    svk2 = n2 / (2.0 * p) - weite2
    lamnu2 = (lamo + lamu + lamou * (2.0 - svk2 / q2)) / 4.0

    xsnut2 = 2.0 * PI * fn * MUE0 * li * zn2 ** 2 * p / 2.0 * 4.0 * q2 * lamnu2
    xs2sl = st.xs2sl
    xs2sl.xsn2 = xsnut2 * (lamnu2 - lamstg - ho2 / bn2) / lamnu2
    xs2sl.xss2 = xsnut2 - xs2sl.xsn2
    lamst2 = 0.2

    xssti2 = 2.0 * PI * fn * MUE0 * (lm2 - li) * zn2 ** 2 * q2 * n2 / 3.0 * lamst2

    r2k = zn2 * n2 * lm2 / (qcu2 * 58.0e6 * 3.0)
    r2w = r2k * (1.0 + 0.00429 * theta2)

    return {
        "hozah2": hozah2, "brzah2": brzah2, "lamnu2": lamnu2,
        "r2k": r2k, "r2w": r2w, "xssti2": xssti2, "xsnut2": xsnut2,
    }


_FNDATA = np.array([0.0, 1., 5., 7., 11., 13., 17., 19., 23., 25., 29.,
                     31., 35., 37., 41., 43., 47., 49., 53., 55., 59.,
                     61., 65., 67., 71., 73.])


def kikr1(izgf: int, mtl: float, ntl: float, azweig: float, btl: float,
          htl: float, fn: float, li: float, lm1: float, theta1: float,
          bn1: float, zn1: float):
    """
    Resistance increase for stator form coils fed from a drive (inverter),
    across 25 harmonic frequencies. Returns svpool, a 1-indexed
    np.ndarray of shape (26, 3): column 1 = frequency, column 2 = KR
    factor (this is the caller's own SVP1(25,2) array, NOT
    st.kennli.svpool).
    """
    svpool = np.zeros((26, 3))
    for j in range(1, 26):
        svpool[j, 1] = _FNDATA[j] * fn

    if izgf != 1:
        for j in range(1, 26):
            svpool[j, 2] = 1.0
        return svpool

    ltrt = mtl * ntl
    ltr = zn1 * azweig
    tpar = ltrt / ltr
    tparue = tpar / ntl
    kappa = 57.0e6 / (1.0 + 0.00425 * theta1)
    lambda_ = (lm1 - li) / li
    ws = zn1 * azweig / 2.0

    for j in range(1, 26):
        omega = 2.0 * PI * fn * _FNDATA[j]
        alpha = math.sqrt(ntl * btl / bn1 * omega / 2.0 * MUE0 * kappa)
        xsi = alpha * htl
        krt = 1.0 + mtl ** 2 * xsi ** 4 / (9.0 * (1.0 + lambda_) * tparue ** 2) * 0.85
        krp = 1.0
        if not (tparue < 1.001 or ws < 1.9):
            alpha = math.sqrt(ntl * btl / bn1 / (1.0 + lambda_) * omega / 2.0 * MUE0 * kappa)
            hlage = mtl * htl * 1.0 / (2.0 * ws)
            xsi = alpha * hlage
            k1 = xsi * (math.sinh(2.0 * xsi) + math.sin(2.0 * xsi)) / \
                       (math.cosh(2.0 * xsi) - math.cos(2.0 * xsi))
            k3 = 2.0 * xsi * (math.sinh(xsi) - math.sin(xsi)) / \
                             (math.cosh(xsi) + math.cos(xsi))
            krp = k1 + (ws ** 2 - 1.0) / 4.0 * k3
        kr = krt + krp - 1.0
        svpool[j, 2] = kr

    return svpool


def kikr(ilfr: int,
         h22: float, br12: float, br22: float, kappa2: float,
         h32: float, br32: float, kappa3: float,
         h42: float, br42: float, kappa4: float,
         h52: float, bstr52: float, kappa5: float,
         h62: float, br52: float, br62: float, kappa6: float,
         h72: float, bstr62: float, kappa7: float,
         qring: float, kapri: float, qringa: float, kapria: float,
         fn: float, di2: float, hojoc2: float, n2: float, beta: float,
         n1: float, p: float, li: float, lstab: float,
         st: MachineState):
    """
    Numeric skin-effect (current displacement) calculation for a cage
    rotor bar using the partial-conductor method ("Teilleiterverfahren").

    ILFR selector (mirrors FORTRAN):
        0 : only the SV factors (st.kennli.svpool columns 2-5) are computed
        2 : both tooth widths (BRZAH2/HOZAH2) and SV factors are computed
        other: only tooth widths are computed (returns early, mirroring
               the FORTRAN early RETURN)

    Returns (hozah2, brzah2, kappaf, lamnu2, r2s, r2kw, qstab).
    brzah2 is a 1-indexed np.ndarray of length 6 (slots 1..5). Any output
    not computed for the given ILFR comes back as None, matching the
    fact the FORTRAN caller simply wouldn't read a stale/uninitialized
    argument in that branch.
    """
    rkz = st.art.rkz
    sv = st.kennli.svpool

    fixed_s = [0.00, 0.03, 0.10, 0.15, 0.20, 0.25, 0.30, 0.38, 0.45, 0.60,
               0.80, 1.00, 1.40]
    for j, s in enumerate(fixed_s, start=1):
        sv[j, 1] = s
    for j in range(1, 13):
        sv[13 + j, 1] = 1.4 + 2.0 * n1 / p / 12.0 * j

    maxn = 260
    h = np.zeros(maxn)
    b = np.zeros(maxn)
    kappa = np.zeros(maxn)
    ianz = 0

    if h72 > 0.0001:
        for k in range(1, 26):
            ianz += 1
            h[ianz] = h72 / 25.0
            xk = 2.0 * k
            x = -h72 + (xk - 1.0) * h[ianz] / 2.0
            b[ianz] = bstr62 * math.sqrt(max(1.0 - (x / h72) ** 2, 0.0))
            kappa[ianz] = kappa7

    if h62 > 0.0001:
        for k in range(1, 51):
            ianz += 1
            h[ianz] = h62 / 50.0
            xk = 2.0 * k
            b[ianz] = br62 + (br52 - br62) * (xk - 1.0) / 100.0
            kappa[ianz] = kappa6

    if h52 > 0.0001:
        for k in range(1, 26):
            ianz += 1
            h[ianz] = h52 / 25.0
            xk = 2.0 * k
            x = (xk - 1.0) * h[ianz] / 2.0
            b[ianz] = bstr52 * math.sqrt(max(1.0 - (x / h52) ** 2, 0.0))
            kappa[ianz] = kappa5
            if h42 > 0.0001:
                if b[ianz] < br42:
                    b[ianz] = br42
                continue
            if h22 > 0.0001:
                if b[ianz] < br22:
                    b[ianz] = br22

    if h42 > 0.0001:
        for k in range(1, 26):
            ianz += 1
            h[ianz] = h42 / 25.0
            b[ianz] = br42
            kappa[ianz] = kappa4

    if h32 > 0.0001:
        for k in range(1, 51):
            ianz += 1
            h[ianz] = h32 / 50.0
            kappa[ianz] = kappa3
            if k * h32 / 50.0 <= h32 - br32 / 2.0:
                if abs((h32 - br32) / br32) > 1.0e-6:
                    b[ianz] = br32
                else:
                    b[ianz] = 2.0 * math.sqrt(max(0.25 * br32 ** 2 -
                              (0.5 * br32 - (k - 0.5) * h32 / 50.0) ** 2, 0.0))
                    if h42 > 0.0001:
                        if b[ianz] < br42:
                            b[ianz] = br42
            else:
                b[ianz] = 2.0 * math.sqrt(max(0.25 * br32 ** 2 -
                          ((k - 0.5) * h32 / 50.0 - (h32 - br32 / 2.0)) ** 2, 0.0))
                if h22 > 0.0001:
                    if b[ianz] < br32:
                        b[ianz] = br32

    if h22 > 0.0001:
        for k in range(1, 26):
            ianz += 1
            h[ianz] = h22 / 25.0
            xk = 2.0 * k
            b[ianz] = br22 + (br12 - br22) * (xk - 1.0) / 50.0
            kappa[ianz] = kappa2

    hozah2 = None
    brzah2 = np.zeros(6)

    if ilfr != 0:
        hozah2 = h22 + h32 + h42 + h52 + h62 + h72
        for j in range(1, 6):
            hsum = 0.0
            for jj in range(1, ianz + 1):
                hsum += h[jj]
                if hsum > (j - 0.5) / 5.0 * hozah2:
                    brzah2[6 - j] = (PI * math.cos(beta) *
                                     (di2 + 2.0 * hojoc2 +
                                      2.0 * (j - 0.5) / 5.0 * hozah2 * math.cos(beta)) / n2
                                     - b[jj])
                    break
        if ilfr != 2:
            return hozah2, brzah2, None, None, None, None, None

    # ---- Teilleiterverfahren: DC values ----
    for i in range(1, ianz + 1):
        if h[i] < 1e-20:
            h[i] = 1e-20
        if b[i] < 1e-20:
            b[i] = 1e-20
        if kappa[i] < 1e-6:
            kappa[i] = 1e-6

    gbez = 0.0
    qstab = 0.0
    for i in range(1, ianz + 1):
        gbez += h[i] * b[i] * kappa[i]
        if kappa[i] > 1.0:
            qstab += h[i] * b[i]
    rbez = 1.0 / gbez
    kappaf = gbez / qstab

    ri = np.zeros(maxn)
    ris = np.zeros(maxn)
    ri[1] = 1e-20
    ris[1] = 0.0
    ris[2] = ri[1]
    for i in range(2, ianz + 1):
        ri[i] = (kappa[i] * h[i] * b[i]) / (kappa[1] * h[1] * b[1]) * ri[1]
        ris[i + 1] = ris[i] + ri[i]

    # ---- AC values, looped over the slip table ----
    for kidx in range(1, 26):
        s = sv[kidx, 1]

        rci = np.zeros(maxn)
        ici = np.zeros(maxn)
        rcis = np.zeros(maxn)
        icis = np.zeros(maxn)
        rci[1] = 1e-20
        ici[1] = 0.0
        rcis[1] = 0.0
        icis[1] = 0.0
        rcis[2] = rci[1]
        icis[2] = ici[1]
        for i in range(2, ianz + 1):
            rci[i] = (-kappa[i] * b[i] * h[i] * 2.0 * PI * s * fn * MUE0 *
                      (h[i] + h[i - 1]) / (b[i] + b[i - 1]) * icis[i] +
                      kappa[i] * b[i] * h[i] / (kappa[i - 1] * b[i - 1] * h[i - 1]) * rci[i - 1])
            ici[i] = (kappa[i] * b[i] * h[i] * 2.0 * PI * s * fn * MUE0 *
                      (h[i] + h[i - 1]) / (b[i] + b[i - 1]) * rcis[i] +
                      kappa[i] * b[i] * h[i] / (kappa[i - 1] * b[i - 1] * h[i - 1]) * ici[i - 1])
            rcis[i + 1] = rcis[i] + rci[i]
            icis[i + 1] = icis[i] + ici[i]

        zkr = zki = nkr = nki = 0.0
        for i in range(1, ianz + 1):
            zki += h[i] / b[i] * (rcis[i] ** 2 + icis[i] ** 2)
            nki += h[i] / b[i] * ris[i] ** 2
            nkr += ri[i] ** 2 / (h[i] * b[i] * kappa[i])
            zkr += (rci[i] ** 2 + ici[i] ** 2) / (h[i] * b[i] * kappa[i])

        kr = zkr / nkr * ris[ianz + 1] ** 2 / (rcis[ianz + 1] ** 2 + icis[ianz + 1] ** 2)
        ki = zki / nki * ris[ianz + 1] ** 2 / (rcis[ianz + 1] ** 2 + icis[ianz + 1] ** 2)

        if s < 0.001:
            lamnu2 = nki / ris[ianz + 1] ** 2
            k1r = 1.0
            k1ra = 1.0 if qringa > 1e-6 else 0.0
        else:
            hg = min(s, 5.0)
            if _rkz_char(rkz, 5) == "O":
                hg = 0.001
            hgr = ((0.7 + 120.0 * qring * kapri / 50.0e6) * math.sqrt(qring) *
                   math.sqrt(PI * hg * fn * MUE0 * kapri))
            k1r = hgr * (math.sinh(2.0 * hgr) + math.sin(2.0 * hgr)) / \
                        (math.cosh(2.0 * hgr) - math.cos(2.0 * hgr))
            if qringa > 1e-6:
                hgra = ((0.7 + 120.0 * qringa * kapria / 50.0e6) * math.sqrt(qringa) *
                        math.sqrt(PI * hg * fn * MUE0 * kapria))
                k1ra = hgra * (math.sinh(2.0 * hgra) + math.sin(2.0 * hgra)) / \
                             (math.cosh(2.0 * hgra) - math.cos(2.0 * hgra))
            else:
                k1ra = 1.0
            c5 = _rkz_char(rkz, 5)
            if c5 == "H":
                if s < 0.70:
                    k1r = (k1r + 1.0) / 2.0 + s / 0.7 * (k1r - 1.0) / 2.0
                    k1ra = (k1ra + 1.0) / 2.0 + s / 0.7 * (k1ra - 1.0) / 2.0
            elif c5 == "J":
                if s >= 1.0:
                    pass
                elif s >= 0.7:
                    k1r = (k1r + 1.0) / 2.0 + (k1r - 1.0) / 2.0 * (s - 0.7) / 0.3
                    k1ra = (k1ra + 1.0) / 2.0 + (k1ra - 1.0) / 2.0 * (s - 0.7) / 0.3
                else:
                    k1r = (k1r + 1.0) / 2.0
                    k1ra = (k1ra + 1.0) / 2.0

        sv[kidx, 2] = kr
        sv[kidx, 3] = ki
        sv[kidx, 4] = k1r
        sv[kidx, 5] = k1ra

    r2s = li * rbez
    r2kw = (lstab - li) * rbez

    return hozah2, brzah2, kappaf, lamnu2, r2s, r2kw, qstab


def kikrsx(xsi: float, krtstr: float):
    """
    Finds XSISTR whose K1 value matches an equivalent high-bar (Hochstab),
    for the slip-ring skin-effect calc. Returns xsistr.
    """
    xsistr = xsi
    for j in range(0, 10001):
        xsistr = xsi + 0.05 * xsi * j
        k1 = xsistr * (math.sinh(2.0 * xsistr) + math.sin(2.0 * xsistr)) / \
                      (math.cosh(2.0 * xsistr) - math.cos(2.0 * xsistr))
        if k1 >= krtstr:
            return xsistr
    return xsistr


def kikrsl(mtl2: float, ntl2: float, azwg2: float, btl2: float, htl2: float,
           fn: float, n1: float, p: float, li: float, lm2: float,
           theta2: float, bn2: float, zn2: float, st: MachineState):
    """
    Slip-ring rotor current-displacement calculation. Returns svpool,
    a 1-indexed np.ndarray of shape (26, 6) (its own array, distinct from
    st.kennli.svpool - matches the FORTRAN SVPOOL argument bound to
    SVPOOL2 at the call site in PROGRAM ASYN4).
    """
    rkz = st.art.rkz
    svpool = np.zeros((26, 6))

    fixed_s = [0.00, 0.03, 0.10, 0.15, 0.20, 0.25, 0.30, 0.38, 0.45, 0.60,
               0.80, 1.00, 1.40]
    for j, s in enumerate(fixed_s, start=1):
        svpool[j, 1] = s
    for j in range(1, 13):
        svpool[13 + j, 1] = 1.4 + 2.0 * n1 / p / 12.0 * j

    if mtl2 * ntl2 * azwg2 * btl2 * htl2 == 0.0 or _rkz_char(rkz, 5) == "O":
        for j in range(1, 26):
            svpool[j, 2] = 1.0
            svpool[j, 3] = 1.0
            svpool[j, 4] = 0.0
            svpool[j, 5] = 0.0
        return svpool

    ltrt = mtl2 * ntl2
    ltr = zn2 * azwg2
    tpar = ltrt / ltr
    tparue = tpar / ntl2
    kappa = 57.0e6 / (1.0 + 0.00425 * theta2)
    lambda_ = (lm2 - li) / li
    ws = zn2 * azwg2 / 2.0

    for j in range(1, 26):
        if j == 1:
            svpool[j, 2] = 1.0
            svpool[j, 3] = 1.0
            svpool[j, 4] = 0.0
            svpool[j, 5] = 0.0
            continue
        omega = 2.0 * PI * fn * svpool[j, 1]
        alpha = math.sqrt(ntl2 * btl2 / bn2 * omega / 2.0 * MUE0 * kappa)
        xsi = alpha * htl2
        k1 = xsi * (math.sinh(2.0 * xsi) + math.sin(2.0 * xsi)) / \
                   (math.cosh(2.0 * xsi) - math.cos(2.0 * xsi))
        k3 = 2.0 * xsi * (math.sinh(xsi) - math.sin(xsi)) / \
                         (math.cosh(xsi) + math.cos(xsi))
        krtstr = k1 + (mtl2 ** 2 - 1.0) / 3.0 * k3
        krt = (k1 + (mtl2 ** 2 - 1.0) / 3.0 * k3 + lambda_) / (1.0 + lambda_)
        krp = 1.0
        xsistr = kikrsx(xsi, krtstr)
        k2 = (1.5 / xsistr * (math.sinh(2.0 * xsistr) - math.sin(2.0 * xsistr)) /
              (math.cosh(2.0 * xsistr) - math.cos(2.0 * xsistr)))
        if not (tparue < 1.001 or ws < 1.9):
            alpha = math.sqrt(ntl2 * btl2 / bn2 / (1.0 + lambda_) * omega / 2.0 * MUE0 * kappa)
            hlage = mtl2 * htl2 * 1.0 / (2.0 * ws)
            xsi = alpha * hlage
            k1 = xsi * (math.sinh(2.0 * xsi) + math.sin(2.0 * xsi)) / \
                       (math.cosh(2.0 * xsi) - math.cos(2.0 * xsi))
            k3 = 2.0 * xsi * (math.sinh(xsi) - math.sin(xsi)) / \
                             (math.cosh(xsi) + math.cos(xsi))
            krp = k1 + (ws ** 2 - 1.0) / 4.0 * k3
        kr = krt + krp - 1.0
        svpool[j, 2] = kr
        svpool[j, 3] = k2
        svpool[j, 4] = 0.0
        svpool[j, 5] = 0.0

    return svpool
