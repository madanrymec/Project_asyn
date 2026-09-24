"""
Eighth batch of ported FORTRAN subroutines from ASYN4.FOR - the
slip-load and starting-curve generation group. These are the two
subroutines that actually drive SPKT/SORT across the full range of slip
values to build the machine's torque-speed characteristic.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    SLAST  (3659-3905)  Slip-load calculation: sorts the pre-generated
                         slip values (via SORT), calls SPKT once per
                         value to build the full M-N characteristic
                         (populates st.mnline), finds the saddle/pull-up
                         torque (MSAT/NSAT), and - if a starting reactor
                         (XNETZ) or reduced starting voltage (KU1) is
                         configured - repeats the whole sweep at reduced
                         voltage (populates st.mnlred).
    HOLAUF (1172-1340)  Starting-up (run-up) calculation: given the
                         already-generated M-N curve (st.mnline) and a
                         counter-torque curve (via MGEGEN), integrates
                         run-up time and adiabatic heating (stator,
                         cage-bar bottom/top, ring) using the
                         acceleration equation dt = 2*pi*J*dn /
                         (M_avg - Mg_avg), stopping once the motor
                         torque can no longer out-accelerate the load.

Both are ported as calculations returning their results in a dict/tuple,
with printing replaced by plain, human-readable console output (the
FORTRAN prints go to line-printer-formatted files with fixed-width
FORMAT statements and HP terminal escape codes for a live-progress
screen - this rewrite keeps the same information but as normal text,
per the "Open questions" note in CONVERSION_STATUS.md; pass verbose=False
to silence it entirely). Reduced-voltage numeric outputs go to
st.mnlred, matching the FORTRAN COMMON /MNLRED/.

Not yet ported: GENERA, DFMAG, PRAXIS, PRUEF, SK, SKGEN, SKNY, SPKTSP
(SPKT-adjacent special-point subroutines), output/reporting subroutines,
plotting driver, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState
from .winding_functions import mgegen
from .geometry_functions import k1k2
from .operating_point import sort, spkt

PI = 3.14159


def slast(reke: str, ilfr: int, ku1: float, xnetz: float,
          p: float, fn: float, r1w: float, xssti1: float, xsnut1: float,
          xssti2: float, xsnut2: float, r2w: float, li: float, hs2: float,
          bs2: float, lamnu2: float, r2s: float, r2kw: float, h42: float,
          br42: float, r2kwa: float, pn: float, u1: float, m1: float,
          prbg0: float, i1n: float, mn: float, qcu1: float, qcu2: float,
          qstab: float, qstaba: float, qring: float, qringa: float,
          n2: float, kma: float, kia: float, st: MachineState,
          verbose: bool = True):
    """
    Slip-load calculation: builds the full M-N characteristic across the
    pre-generated (st.svorga) slip values, at both full and (if
    configured) reduced starting voltage.

    Returns (msat, nsat): the saddle torque and the speed it occurs at.
    """
    nyinfo = st.nyinfo
    maxi = st.maxi
    mnline = st.mnline
    mnlred = st.mnlred

    bs2i = bs2 if ilfr == 1 else None
    nsyn = fn / p

    ssanz, ss = sort(st)
    mnline.nanz = ssanz
    msat = 1.0e12
    nsat = 0.0

    if verbose:
        print(f"\n{reke}")
        print("M-N-KENNLINIE (DIREKTE EINSCHALTUNG)")
        print(" N/N1     M/MN    I1/I1N   COSPH   BS2I(mm)")

    for j in range(1, ssanz + 1):
        out = spkt(ilfr, ss[j], p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                   r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                   u1, m1, prbg0, 0.0, 0, st)
        m, i1, i2, i2o, i2u, n, cosi, bs2i = (
            out["m"], out["i1"], out["i2"], out["i2o"], out["i2u"],
            out["n"], out["cosph"], out["bs2i"])

        hg = min(ss[j], 1.0)
        m = m + (kma - 1.0) * hg * m
        i1 = i1 + (kia - 1.0) * hg * i1
        i2 = i2 + (kia - 1.0) * hg * i2
        i2o = i2o + (kia - 1.0) * hg * i2o
        i2u = i2u + (kia - 1.0) * hg * i2u

        if st.esb.skipp < ss[j] <= 1.00001:
            if m < msat:
                msat = m
                nsat = n

        mnline.npu[j] = n / nsyn
        mnline.mpu[j] = m / mn
        mnline.i1pu[j] = i1 / i1n
        mnline.cosph[j] = cosi
        mnline.g[j, 1] = i1 / qcu1
        mnline.upu[j] = 1.00

        if ilfr == 1:
            mnline.g[j, 2] = i2u / qcu2
            mnline.g[j, 3] = 0.0
            mnline.g[j, 4] = 0.0
            mnline.g[j, 5] = 0.0
        else:
            mnline.g[j, 2] = i2u / qstab
            mnline.g[j, 4] = i2u / qring / (2.0 * math.sin(p * PI / n2))
            if ilfr in (3, 4):
                mnline.g[j, 3] = i2o / qstaba
                mnline.g[j, 5] = (i2o / qringa / (2.0 * math.sin(p * PI / n2))
                                   if ilfr == 4 else 0.0)
            else:
                mnline.g[j, 3] = 0.0
                mnline.g[j, 5] = 0.0

        if verbose and abs(math.fmod(ss[j], 0.04999999)) <= 0.00001 and ss[j] <= 1.00001:
            print(f" {mnline.npu[j]:6.3f}  {mnline.mpu[j]:5.2f}  {mnline.i1pu[j]:7.2f}  "
                  f"{mnline.cosph[j]:6.3f}  {bs2i * 1000.0:8.2f}")

    if verbose:
        print(f"\nSATTELMOMENT: MSAT/MN = {msat / mn:5.2f}  BEI NSAT = {nsat * 60.0:7.1f} 1/MIN")
        print("\nMAX. INNERE MOMENTE DER POLPAARZAHLEN")
        print("  NY     N/N1    MNY/MN")
        for j in range(1, nyinfo.nyanz + 1):
            print(f"{nyinfo.ny[j]:6.0f}  {1.0 - maxi.smaxi[j]:6.3f}  {maxi.mmaxi[j] / mn:7.3f}")

    if ku1 < 0.0001 and xnetz < 0.000001:
        return msat, nsat

    # ---- reduced voltage / starting reactor sweep ----
    if verbose:
        print(f"\n{reke}")
        print("M-N-KENNLINIE (RED. SPANNUNG U. VORREAKTANZ)")
        print(" N/N1     M/MN    I1/I1N   COSPH   BS2I(mm)  U1/U1N")

    u1red = ku1 * u1
    s = 1.0e-9
    spkt(ilfr, s, p, fn, r1w, xssti1 + xnetz, xsnut1, xssti2, xsnut2,
         r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
         u1red, m1, prbg0, 0.0, 1, st)

    mnlred.nanz = 0
    nanz2 = 0
    for j in range(1, ssanz + 1):
        if ss[j] > 1.00001:
            continue
        out = spkt(ilfr, ss[j], p, fn, r1w, xssti1 + xnetz, xsnut1, xssti2,
                   xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42,
                   r2kwa, pn, u1red, m1, prbg0, 0.0, 0, st)
        m, i1, i2, i2o, i2u, n, cosi = (
            out["m"], out["i1"], out["i2"], out["i2o"], out["i2u"],
            out["n"], out["cosph"])

        hg = min(ss[j], 1.0)
        kmared = 1.0 + (kma - 1.0) * math.sqrt(ku1)
        kiared = 1.0 + (kia - 1.0) * math.sqrt(ku1)
        m = m + (kmared - 1.0) * hg * m
        i1 = i1 + (kiared - 1.0) * hg * i1
        i2 = i2 + (kiared - 1.0) * hg * i2
        i2o = i2o + (kiared - 1.0) * hg * i2o
        i2u = i2u + (kiared - 1.0) * hg * i2u

        nanz2 += 1
        mnlred.npu[nanz2] = n / nsyn
        mnlred.mpu[nanz2] = m / mn
        mnlred.i1pu[nanz2] = i1 / i1n
        mnlred.cosph[nanz2] = cosi
        mnlred.g[nanz2, 1] = i1 / qcu1

        if ilfr == 1:
            mnlred.g[nanz2, 2] = i2u / qcu2
            mnlred.g[nanz2, 3] = 0.0
            mnlred.g[nanz2, 4] = 0.0
            mnlred.g[nanz2, 5] = 0.0
        else:
            mnlred.g[nanz2, 2] = i2u / qstab
            mnlred.g[nanz2, 4] = i2u / qring / (2.0 * math.sin(p * PI / n2))
            if ilfr in (3, 4):
                mnlred.g[nanz2, 3] = i2o / qstaba
                mnlred.g[nanz2, 5] = (i2o / qringa / (2.0 * math.sin(p * PI / n2))
                                       if ilfr == 4 else 0.0)
            else:
                mnlred.g[nanz2, 3] = 0.0
                mnlred.g[nanz2, 5] = 0.0

        sini = -math.sin(math.acos(cosi))
        u1str = math.sqrt((u1red + xnetz * i1 * sini) ** 2 +
                           (xnetz * i1 * cosi) ** 2)
        mnlred.upu[nanz2] = u1str / u1

        if verbose and abs(math.fmod(ss[j], 0.04999999)) <= 0.00001:
            print(f" {n / nsyn:6.3f}  {m / mn:5.2f}  {i1 / i1n:7.2f}  "
                  f"{cosi:6.3f}  {out['bs2i'] * 1000.0 if out['bs2i'] else 0.0:8.2f}  "
                  f"{u1str / u1:6.2f}")
    mnlred.nanz = nanz2

    # ---- restore original (full-voltage) saturation state ----
    s = 1.0e-9
    spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
         r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
         u1, m1, prbg0, 0.0, 1, st)

    return msat, nsat


def holauf(ilfr: int, fn: float, p: float, jsum: float, mn: float,
           k41: float, k42: float, k4a2: float, k4r2: float, k4ar2: float,
           r1w: float, m1: float, i1n: float, ski: float,
           ku1: float, xnetz: float, idruck: int, st: MachineState,
           verbose: bool = True):
    """
    Run-up (starting) calculation: integrates run-up time and adiabatic
    heating using the already-generated M-N curve in st.mnline, and the
    counter-torque (load) curve via MGEGEN.

    Returns a dict with th (total run-up time), erw1/erw2/erwa2 (stator/
    cage-bottom/cage-top adiabatic heating), erwr/erwra (ring heating),
    nends (the stationary end speed, per-unit), plus the FORTRAN-local
    (not officially returned, but harmless to expose) q1h/q2h/iqt/mmit.
    """
    mnline = st.mnline
    nanz = mnline.nanz

    k1 = k2 = k1r = k1ra = 1.0
    nsyn = fn / p
    nk = 1.0 - ski
    fak = 2.0 * PI * jsum

    nend = 0.0
    th = 0.0
    erw1 = erw2 = erwa2 = erwr = erwra = 0.0
    q1h = q2h = iqt = mmit = 0.0
    izaehl = 0
    imerk = 0
    nends = 0.0
    ihg1 = None

    if verbose and idruck == 1:
        print(f"\nERGEBNIS DER HOCHLAUFRECHNUNG\nKU1={ku1:.2f}  XNETZ={xnetz:.6f} OHM")

    for j in range(2, nanz + 1):
        if mnline.npu[j] < 0.00001:
            continue

        mmit += (mnline.mpu[j] + mnline.mpu[j - 1]) / 2.0 * (mnline.npu[j] - mnline.npu[j - 1])
        na = mnline.npu[j - 1] * nsyn
        nb = mnline.npu[j] * nsyn
        ma = mnline.mpu[j - 1] * mn
        mb = mnline.mpu[j] * mn
        mga = mgegen(mnline.npu[j - 1], st)
        mgb = mgegen(mnline.npu[j], st)
        mgpu = mga / mn
        k1, k2, k1r, k1ra = k1k2(1.0 - mnline.npu[j - 1], st)

        if mb - mgb <= 0.0 and imerk == 0:
            # ---- final interval: no further acceleration possible ----
            imerk = 1
            nend = (na * (mb - mgb) - nb * (ma - mga)) / (mb - mgb - ma + mga)
            nends = nend
            dt = fak * (nend - na) / (0.5 * (ma - mga))
            th += dt
            erw1 += (0.5 * (mnline.g[j - 1, 1] + mnline.g[j, 1])) ** 2 / k41 * dt
            erw2 += (0.5 * (mnline.g[j - 1, 2] + mnline.g[j, 2])) ** 2 / k42 * dt * k1
            erwa2 += (0.5 * (mnline.g[j - 1, 3] + mnline.g[j, 3])) ** 2 / k4a2 * dt
            erwr += (0.5 * (mnline.g[j - 1, 4] + mnline.g[j, 4])) ** 2 / k4r2 * dt * k1r
            erwra += (0.5 * (mnline.g[j - 1, 5] + mnline.g[j, 5])) ** 2 / k4ar2 * dt * k1ra
            mq = 0.5 * (mb + ma)
            s = 1.0 - 0.5 * (na + nb) / nsyn
            q1h += m1 * r1w * (0.5 * (mnline.i1pu[j - 1] + mnline.i1pu[j]) * i1n) ** 2 * dt
            q2h += s * 2.0 * PI * nsyn * mq * dt
            iqt += mnline.i1pu[j - 1] ** 2 * dt
        else:
            # ---- intermediate interval ----
            nend = nb
            dt = fak * (nb - na) / (0.5 * (ma - mga + mb - mgb))
            if imerk == 1:
                dt = 0.0
            th += dt
            erw1 += (0.5 * (mnline.g[j - 1, 1] + mnline.g[j, 1])) ** 2 / k41 * dt
            erw2 += (0.5 * (mnline.g[j - 1, 2] + mnline.g[j, 2])) ** 2 / k42 * dt * k1
            erwa2 += (0.5 * (mnline.g[j - 1, 3] + mnline.g[j, 3])) ** 2 / k4a2 * dt
            erwr += (0.5 * (mnline.g[j - 1, 4] + mnline.g[j, 4])) ** 2 / k4r2 * dt * k1r
            erwra += (0.5 * (mnline.g[j - 1, 5] + mnline.g[j, 5])) ** 2 / k4ar2 * dt * k1ra
            mq = 0.5 * (mb + ma)
            s = 1.0 - 0.5 * (na + nb) / nsyn
            q1h += m1 * r1w * (0.5 * (mnline.i1pu[j - 1] + mnline.i1pu[j]) * i1n) ** 2 * dt
            q2h += s * 2.0 * PI * nsyn * mq * dt
            iqt += mnline.i1pu[j - 1] ** 2 * dt

    if verbose and idruck == 1:
        print(f"MITTLERES MOTORMOMENT: {mmit:5.2f} * MN   I**2 x t = {iqt:7.0f}")
        print(f"NEND = {nends * 60.0:7.1f} 1/MIN   TH = {th:7.1f} S")
        print(f"ERW1 = {erw1:7.1f} K   ERW2 = {erw2:7.1f} K   ERWA2 = {erwa2:7.1f} K")
        print(f"ERWR = {erwr:7.1f} K   ERWRA = {erwra:7.1f} K")
        print(f"Q1H = {q1h / 1000.0:7.2f} KJ   Q2H = {q2h / 1000.0:7.2f} KJ")

    return {
        "th": th, "erw1": erw1, "erw2": erw2, "erwa2": erwa2,
        "erwr": erwr, "erwra": erwra, "nends": nends,
        "q1h": q1h, "q2h": q2h, "iqt": iqt, "mmit": mmit,
    }
