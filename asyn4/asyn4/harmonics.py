"""
Sixth batch of ported FORTRAN subroutines from ASYN4.FOR - completes the
"harmonic-field & winding-factor engine" functional group.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    POPAZA (7789-7933)  Determines which harmonic pole-pair numbers need
                         active damping, computes the stator's doubly-
                         linked (double-coupled) leakage SD1/SD1R, and a
                         set of per-harmonic auxiliary quantities (SD2,
                         UER, UEI, SF, X1HBEZ, RRNY/RRNYA, XSRNY/XSRNYA).
                         Populates st.nyinfo in place.
    Z2NYS  (9823-9960)  Computes the damped harmonic-field impedance
                         (CZ2NYS) and complex damping factor (CDNY) for a
                         given pole-pair number NY, for all 4 rotor types
                         (slip-ring, single-cage, double-cage common-ring,
                         double-cage separate-rings). When IS1PKT==1, also
                         writes the fundamental-frequency branch values
                         back into st.esb (mirrors the FORTRAN in/out
                         COMMON /ESB/ side effects).

This completes the "Harmonic-field & winding-factor engine" functional
group (together with NUE0, MNYQMP, K1K2, ZFAK/XSIZ family from earlier
rounds).

Not yet ported: insulation design (ISOSTD/ISOTR/...), performance calcs
(GENERA/SK/...), load/starting (SLAST/HOLAUF), output/reporting, plotting
driver, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import cmath
import math
import numpy as np

from .state import MachineState
from .winding_functions import xsiz, xsis, xsisch, xsiz2
from .geometry_functions import k1k2

PI = 3.14159
MY0 = 1.2566e-6


def _rkz_char(rkz: str, pos: int) -> str:
    if pos - 1 < len(rkz):
        return rkz[pos - 1]
    return " "


def popaza(ilfr: int, p: float, fn: float, zn1: float, zn2: float,
           n1: float, n2: float, m1: float, mzone1: float, mz1s: float,
           weite1: float, weite2: float, dring: float, kapri: float,
           qring: float, dringa: float, kapria: float, qringa: float,
           schr: float, st: MachineState) -> float:
    """
    Selects the harmonic pole-pair numbers requiring active damping,
    computes the stator's double-coupled leakage (SD1/SD1R), and
    per-harmonic auxiliary quantities. Populates st.nyinfo in place.
    Returns FW2 (rotor winding factor for the fundamental, only
    meaningful for slip-ring rotors, ILFR==1).
    """
    rkz = st.art.rkz
    ny_info = st.nyinfo

    q1 = n1 / (p * mzone1)
    om1 = 2.0 * PI * fn
    nyanz = 1
    ny = np.zeros(26)
    ny[1] = p

    if _rkz_char(rkz, 8) == " ":
        # ---- select harmonics to damp ----
        g_max = int(2.0 * n1 / p / mz1s + 0.001)
        stop = False
        for g in range(1, g_max + 1):
            if stop:
                break
            for gg in (1, 2):
                hg = p * (1.0 + mz1s * g * (-1) ** gg)
                hg1 = hg * PI / n2
                mknzmn = (100.0 * p / abs(hg) *
                          (xsiz(hg, q1, n1, st) * xsis(hg, weite1, n1) /
                           xsiz(p, q1, n1, st) / xsis(p, weite1, n1) *
                           xsisch(hg, schr) * math.sin(hg1) / hg1) ** 2)
                if p < 1.001:
                    mknzmn = 2.0 * mknzmn
                if mknzmn < 0.05:
                    continue
                if nyanz >= 25:
                    print("\n ABZUDAEMPFENDE OBERFELDER TEILWEISE UNTERDRUECKT !")
                    stop = True
                    break
                nyanz += 1
                ny[nyanz] = hg

    # ---- double-coupled stator leakage ----
    hgsum = 0.0
    for g in range(1, 501):
        for gg in (1, 2):
            my = p * (1.0 + mz1s * g * (-1) ** gg)
            hgsum += (xsiz(my, q1, n1, st) * xsis(my, weite1, n1) / my) ** 2
    sd1 = hgsum / ((xsiz(p, q1, n1, st) * xsis(p, weite1, n1) / p) ** 2)

    for j in range(2, nyanz + 1):
        hgsum -= (xsiz(ny[j], q1, n1, st) * xsis(ny[j], weite1, n1) / ny[j]) ** 2
    sd1r = hgsum / ((xsiz(p, q1, n1, st) * xsis(p, weite1, n1) / p) ** 2)

    sd2 = np.zeros(26)
    uer = np.zeros(26)
    uei = np.zeros(26)
    sf = np.zeros(26)
    x1hbez = np.zeros(26)
    rrny = np.zeros(26)
    rrnya = np.zeros(26)
    xsrny = np.zeros(26)
    xsrnya = np.zeros(26)
    fw2 = 1.0

    for j in range(1, nyanz + 1):
        xsi1 = xsiz(ny[j], q1, n1, st) * xsis(ny[j], weite1, n1)
        if ilfr == 1:
            # ---- slip-ring rotor ----
            q2 = n2 / (6.0 * p)
            hgsum2 = (xsiz2(p, q2, n2) * xsis(p, weite2, n2) * PI /
                      (n2 * math.sin(p * PI / n2))) ** 2
            for g in range(1, 101):
                for gg in (1, 2):
                    my = p * (1.0 + 6.0 * g * (-1) ** gg)
                    if my > n2 / 2.0 or my <= -n2 / 2.0:
                        continue
                    hgsum2 += (xsiz2(my, q2, n2) * xsis(my, weite2, n2) *
                               PI / (n2 * math.sin(my * PI / n2))) ** 2
            xsi2 = xsiz2(ny[j], q2, n2) * xsis(ny[j], weite2, n2)
            if j == 1:
                fw2 = xsi2
            sd2[j] = hgsum2 / ((xsi2 / ny[j]) ** 2) - 1.0
            uer[j] = (zn1 * n1 * xsi1 / (zn2 * n2 * xsi2)) ** 2
            uei[j] = 1.0 / math.sqrt(uer[j])
        else:
            # ---- cage rotor ----
            hg = ny[j] * PI / n2
            fw2 = 1.0
            sd2[j] = (hg / math.sin(hg)) ** 2 - 1.0
            uer[j] = (zn1 * n1 * xsi1) ** 2 / (m1 * n2) / st.vorzei.ksys ** 2
            uei[j] = n2 / (zn1 * n1 * xsi1) * st.vorzei.ksys
            rrny[j] = PI * dring / (qring * n2 * kapri) / (2.0 * math.sin(hg) ** 2)
            xsrny[j] = om1 * 0.37 * MY0 * PI * dring / n2 / (2.0 * math.sin(hg) ** 2)
            if ilfr == 4:
                rrnya[j] = PI * dringa / (qringa * n2 * kapria) / (2.0 * math.sin(hg) ** 2)
                xsrnya[j] = om1 * 0.37 * MY0 * PI * dringa / n2 / (2.0 * math.sin(hg) ** 2)

        sf[j] = xsisch(ny[j], schr)
        x1hbez[j] = (xsiz(ny[j], q1, n1, st) * xsis(ny[j], weite1, n1) * p /
                     (xsiz(p, q1, n1, st) * xsis(p, weite1, n1) * ny[j])) ** 2

    ny_info.sd1 = sd1
    ny_info.sd1r = sd1r
    ny_info.nyanz = nyanz
    ny_info.ny[:len(ny)] = ny
    ny_info.sd2[:len(sd2)] = sd2
    ny_info.uer[:len(uer)] = uer
    ny_info.uei[:len(uei)] = uei
    ny_info.sf[:len(sf)] = sf
    ny_info.x1hbez[:len(x1hbez)] = x1hbez
    ny_info.rrny[:len(rrny)] = rrny
    ny_info.rrnya[:len(rrnya)] = rrnya
    ny_info.xsrny[:len(xsrny)] = xsrny
    ny_info.xsrnya[:len(xsrnya)] = xsrnya

    return fw2


def z2nys(ilfr: int, sny: float, ny: float, fn: float, x1hny: float,
          x1hny0: float, sf: float, uer: float, xssti2: float,
          xsnut2: float, sd2: float, r2w: float, xssteg: float, li: float,
          lamnu2: float, xsrny: float, rrny: float, r2s: float, r2kw: float,
          h42: float, br42: float, r2kwa: float, xsrnya: float,
          rrnya: float, is1pkt: int, st: MachineState):
    """
    Damped harmonic-field impedance (CZ2NYS) and complex damping factor
    (CDNY) for pole-pair number NY. When IS1PKT==1, also writes the
    fundamental-frequency branch values into st.esb (mirroring the
    FORTRAN COMMON /ESB/ side effects at S=1, NY=P).

    Returns (cz2nys, cdny, i2opu, i2upu) - cz2nys/cdny are Python complex.
    """
    esb = st.esb
    jot = 1j

    x1hn0 = x1hny0
    x1hn = x1hny

    xssd2 = sd2 * x1hn0
    xssch = (1.0 - sf) * x1hn0
    if is1pkt == 1:
        esb.x1h = x1hny
        esb.xsd2 = xssd2
        esb.xsch = xssch

    om1 = 2.0 * PI * fn

    if ilfr == 1:
        # ---- slip-ring rotor ----
        k1, k2, k1r, k1ra = k1k2(sny, st)
        czl2 = (uer * k1 * r2w / sny +
                jot * (xssch + xssd2 + uer * (xssti2 + st.xs2sl.xss2 + k2 * st.xs2sl.xsn2)))
        if is1pkt == 1:
            esb.rk2 = r2w * uer
            esb.xring = xssti2 * uer
            esb.xn2 = xsnut2 * uer
        i2opu = 0.0
        i2upu = 1.0

    elif ilfr == 2:
        # ---- single-cage rotor, arbitrary bar shape ----
        k1, k2, k1r, k1ra = k1k2(sny, st)
        hg = om1 * MY0 * li
        czl2 = (uer * (k1r * rrny + r2kw + k1 * r2s) / sny +
                jot * (xssch + xssd2 + uer * (xsrny + xssteg + hg * lamnu2 * k2)))
        if is1pkt == 1:
            esb.rk2 = r2kw * uer
            esb.rs2 = k1 * r2s * uer
            esb.rr2 = rrny * uer
            esb.xring = xsrny * uer
            esb.xsteg = xssteg * uer
            esb.xn2 = (xssteg + hg * lamnu2 * k2) * uer
        i2opu = 0.0
        i2upu = 1.0

    else:
        # ---- double-cage rotor ----
        k1, k2, k1r, k1ra = k1k2(sny, st)
        hg = om1 * MY0 * li
        xsng = hg * 0.79 + xssteg
        xsna = hg * 0.66 + xssteg
        xsnb = hg * (lamnu2 * k2 + h42 / br42 + 1.58) + xssteg
        if is1pkt == 1:
            esb.xng = xsng * uer
            esb.xna = xsna * uer
            esb.xn2 = xsnb * uer
            esb.xring = xsrny * uer
            esb.xsteg = xssteg * uer
            esb.rr2 = rrny * uer
            esb.ra2 = r2kwa * uer
            esb.rk2 = r2kw * uer
            esb.rs2 = k1 * r2s * uer

        if ilfr == 3:
            # ---- common closing ring ----
            cza = uer * (r2kwa / sny + jot * (xsna - xsng))
            czb = uer * ((r2kw + k1 * r2s) / sny + jot * (xsnb - xsng))
            czl2 = (uer * k1r * rrny / sny + jot * (xssch + xssd2 + uer * (xsrny + xsng)) +
                    cza * czb / (cza + czb))
        else:  # ilfr == 4: separate rings
            cza = uer * ((r2kwa + k1ra * rrnya) / sny + jot * (xsna - xsng + xsrnya))
            czb = uer * ((r2kw + k1 * r2s + k1r * rrny) / sny + jot * (xsnb - xsng + xsrny))
            czl2 = jot * (xssch + xssd2 + uer * xsng) + cza * czb / (cza + czb)
            if is1pkt == 1:
                esb.rra2 = rrnya * uer
                esb.xringa = xsrnya * uer

        i2opu = abs(czb / (cza + czb))
        i2upu = abs(cza / (cza + czb))

    cz2nys = jot * xssch + jot * sf * x1hn * czl2 / (jot * sf * x1hn + czl2)
    cdny = cz2nys / (jot * x1hn)

    return cz2nys, cdny, i2opu, i2upu
