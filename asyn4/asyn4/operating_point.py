"""
Seventh batch of ported FORTRAN subroutines from ASYN4.FOR - the
single-operating-point solver group. This is the computational core that
SLAST (slip-load calculation) and the starting-curve generation call
repeatedly, once per slip value, to get currents/powers/torque at that
point.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    MKWERT (7497-7548)  Interpolates the KREIS-generated magnetization
                         characteristic (st.kennli.mkpool) at a given
                         induced voltage UI. INTART=0 gives just
                         IMY/IFE (fast path used inside SPKT's iteration
                         loop); INTART=1 gives the full set of derived
                         quantities (used once the operating point has
                         converged, for reporting).
    SORT   (8479-8504)  Sorts the pre-generated (st.svorga) slip values
                         into descending order (selection sort, exactly
                         as the FORTRAN does it - not the fastest
                         algorithm, but matching it avoids any ordering
                         surprises later).
    SPKT   (8789-8978)  THE operating-point solver: for a given slip S,
                         iterates (simple successive substitution) on the
                         induced voltage UI, the stator slot-saturation
                         reactance correction DXS1, and (for cage rotors)
                         the slot-bridge saturation LAMSTG until self-
                         consistent, then computes current, power, and
                         torque contributions summed over every damped
                         harmonic pole-pair number. Updates st.dp and
                         st.esb when IS1PKT==1 (the synchronous/reference
                         point). Updates st.maxi (per-harmonic peak
                         torque tracking) as a side effect, exactly as
                         the FORTRAN COMMON /MAXI/ does.

Added `Dp` (COMMON /DP/) and `Svorga` (COMMON /SVORGA/) to state.py -
**replace state.py again**, see state_NEW.py.

Not yet ported: GENERA, DFMAG, PRAXIS, PRUEF, SK, SKGEN, SKNY, SPKTSP
(other members of the performance-calculation group), SLAST, HOLAUF
(load/starting curves that call SPKT/SORT), output/reporting, plotting
driver, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import math
import cmath

from .state import MachineState
from .winding_functions import dxs1, lastg, snue

PI = 3.14159
MY0 = 1.2566e-6


def mkwert(intart: int, ui: float, st: MachineState):
    """
    Interpolates the KREIS-generated magnetization characteristic at
    induced voltage UI. Raises ValueError if UI falls outside the
    generated curve (mirrors the FORTRAN `STOP` - a hard error, not
    something the caller is expected to recover from, since it means the
    machine design/voltage combination fell outside the range KREIS
    swept).

    Returns a dict. If intart != 1, only 'imy' and 'ife' are populated
    (fast path); if intart == 1, all fields are populated (matches the
    FORTRAN's early RETURN after the fast-path fields).
    """
    mkpool = st.kennli.mkpool
    eanz = st.kennli.eanz

    if ui < mkpool[1, 1] or ui > mkpool[eanz, 1]:
        raise ValueError("MKWERT: UI-WERT AUSSERHALB DER KENNLINIE")

    for j in range(2, eanz + 1):
        if ui <= mkpool[j, 1]:
            hg = (ui - mkpool[j - 1, 1]) / (mkpool[j, 1] - mkpool[j - 1, 1])

            def interp(col):
                return mkpool[j - 1, col] + (mkpool[j, col] - mkpool[j - 1, col]) * hg

            out = {"imy": interp(2), "ife": interp(8)}
            if intart != 1:
                return out
            out.update({
                "deltap": interp(3), "bp": interp(4), "b3p": interp(5),
                "pfej1": interp(6), "pfe1": interp(7),
                "bz1max": interp(9), "bz2max": interp(10),
                "bjoch1": interp(11), "bjoch2": interp(12),
                "vluft": interp(13), "vzahn1": interp(14), "vzahn2": interp(15),
                "vjoch1": interp(16), "vjoch2": interp(17), "alpha": interp(18),
            })
            return out
    # Only reached if the loop never finds UI <= mkpool[j,1], which the
    # range check above should already have excluded.
    raise ValueError("MKWERT: UI-WERT AUSSERHALB DER KENNLINIE")


def sort(st: MachineState):
    """
    Sorts st.svorga.svor (the pre-generated, unsorted slip values) into
    descending order via selection sort (matches the FORTRAN algorithm
    exactly, including its O(n^2) behaviour). A slip of exactly 0.0 is
    nudged to 1e-9 to avoid later division-by-zero, same as the FORTRAN.

    Returns (ssanz, ss) - ss is a 1-indexed np.ndarray of length 401.
    """
    import numpy as np
    sanz = st.svorga.sanz
    svor = st.svorga.svor.copy()
    ss = np.zeros(401)

    for j in range(1, sanz + 1):
        smax = -1.0
        jsave = None
        for jj in range(1, sanz + 1):
            if svor[jj] > smax:
                jsave = jj
                smax = svor[jj]
        ss[j] = svor[jsave]
        if ss[j] == 0.0:
            ss[j] = 1.0e-9
        svor[jsave] = -10.0

    return sanz, ss


def spkt(ilfr: int, s: float, p: float, fn: float, r1w: float,
         xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
         r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
         r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
         u1: float, m1: float, prbg0: float, rzus: float, is1pkt: int,
         st: MachineState):
    """
    Computes the full operating point (currents, powers, torque) at slip
    S. This is the workhorse called once per slip value by SLAST/HOLAUF's
    curve generation, and once more (IS1PKT=1) at synchronism to latch
    the equivalent-circuit reference values into st.esb/st.dp.

    IMPORTANT: st.dp.x1hp0 (the reference main-field reactance) is only
    ever written when is1pkt==1. Every other call divides by it, so you
    MUST call this once with is1pkt=1 (at the rated operating point,
    mirroring how PROGRAM ASYN4 does it before running SLAST's slip
    sweep) before calling it with is1pkt=0 anywhere else - otherwise
    x1hp0 stays at its state.py default of 0.0 and every harmonic
    impedance downstream divides by zero. This was discovered during
    integration testing of SLAST/HOLAUF (see load_curves.py) - not
    obvious from reading SPKT's FORTRAN source in isolation, since the
    FORTRAN COMMON /DP/ silently carries whatever the last IS1PKT=1 call
    left there.

    Returns a dict with i1, i2, i2o, i2u, imy, ui, p1, pfe1, pcu1, pmech,
    prbg, pcu2, pcu2ny, pzus, m, mny, n, cosph, bs2i, usd1, rsd1, isd1.
    """
    from .harmonics import z2nys  # local import: harmonics imports geometry_functions,
                                    # which does not import this module, so this is
                                    # not a true cycle, but kept local for clarity.

    nyinfo = st.nyinfo
    nyanz = nyinfo.nyanz

    if abs(s) < 1.0e-9:
        s = 1.0e-9
    jot = 1j
    om1 = 2.0 * PI * fn

    i1 = pn / (m1 * u1 * 0.9)
    i2 = i1 / nyinfo.uei[1]
    dxs1n = dxs1(i1, st)
    lamstg = None
    if ilfr != 1:
        lamstg = lastg(i2, st)
    ui = 0.97 * u1

    bs2i = (st.daten.bs2 if hs2 < 1.0e-9 else hs2 / lamstg) if (ilfr != 1) else None
    xssteg = 0.0

    for _iter in range(1, 1000):
        uiv = ui
        dxs1v = dxs1n

        if ilfr != 1:
            lastgv = lamstg
            # HS2=0 (no slot lip - confirmed as a real, valid input by
            # a real machine's data during validation of this port)
            # makes STEG's entire LAMDA table zero, so LAMSTG=0 and
            # BS2I=HS2/LAMSTG becomes a literal 0/0. The original
            # FORTRAN uses the exact same formula with no explicit
            # guard here, so this must have relied on that runtime's
            # non-strict floating-point exception handling rather than
            # any deliberate logic - Python's strict NaN propagation
            # doesn't replicate that. The real machine's own AUSASYN
            # output shows BS2I ending up exactly equal to BS2 in this
            # case, so that's applied here explicitly (mirroring the
            # ILFR==1 slip-ring bypass, which also just sets BS2I=BS2
            # directly) - and XSSTEG (the slot-bridge leakage
            # reactance) is 0 in this case since there's no separate
            # bridge to speak of when the slot has no lip.
            if hs2 < 1.0e-9:
                bs2i = st.daten.bs2
                xssteg = 0.0
            else:
                bs2i = hs2 / lamstg
                xssteg = om1 * MY0 * li * lamstg

        mk = mkwert(0, ui, st)
        imy = mk["imy"]
        x1hp = ui / imy
        if is1pkt == 1:
            st.dp.x1hp0 = x1hp
        x1hp0 = st.dp.x1hp0

        cz = ((r1w + rzus) +
              jot * (xssti1 + xsnut1 + nyinfo.sd1r * x1hp0 + dxs1n))
        czdvs = jot * nyinfo.sd1r * x1hp

        for j in range(nyanz, 0, -1):
            sny = snue(nyinfo.ny[j], p, s)
            x1hny = nyinfo.x1hbez[j] * x1hp
            x1hny0 = nyinfo.x1hbez[j] * x1hp0
            if j == 1:
                x1hny = x1hny0
            cz2nys, cdny, i2opu, i2upu = z2nys(
                ilfr, sny, nyinfo.ny[j], fn, x1hny, x1hny0, nyinfo.sf[j],
                nyinfo.uer[j], xssti2, xsnut2, nyinfo.sd2[j], r2w, xssteg,
                li, lamnu2, nyinfo.xsrny[j], nyinfo.rrny[j], r2s, r2kw,
                h42, br42, r2kwa, nyinfo.xsrnya[j], nyinfo.rrnya[j], 0, st)
            cz = cz + cz2nys
            if j != 1:
                czdvs = czdvs + cz2nys

        csd1 = czdvs / (jot * x1hp)
        rsd1 = csd1.real
        isd1 = csd1.imag
        usd1 = nyinfo.sd1
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

    # ---- operating point converged ----
    i1 = abs(ci1)
    pcu2 = 0.0
    mi = 0.0
    mny = 0.0

    for j in range(nyanz, 0, -1):
        sny = snue(nyinfo.ny[j], p, s)
        x1hny = nyinfo.x1hbez[j] * x1hp
        x1hny0 = nyinfo.x1hbez[j] * x1hp0
        if j == 1:
            x1hny = x1hny0
        cz2nys, cdny, i2opu, i2upu = z2nys(
            ilfr, sny, nyinfo.ny[j], fn, x1hny, x1hny0, nyinfo.sf[j],
            nyinfo.uer[j], xssti2, xsnut2, nyinfo.sd2[j], r2w, xssteg,
            li, lamnu2, nyinfo.xsrny[j], nyinfo.rrny[j], r2s, r2kw,
            h42, br42, r2kwa, nyinfo.xsrnya[j], nyinfo.rrnya[j], 0, st)
        mny = m1 * nyinfo.ny[j] * x1hny / om1 * i1 ** 2 * (1.0 - cdny).imag

        if abs(mny) > st.maxi.mmaxi[j] and s > 0.0:
            st.maxi.mmaxi[j] = abs(mny)
            st.maxi.smaxi[j] = s

        pdny = om1 / nyinfo.ny[j] * mny
        pcu2ny = sny * pdny
        mi += mny
        pcu2 += pcu2ny

    st.dp.cdp = cdny
    pzus = m1 * rzus * i1 ** 2
    prbg = abs(1.0 - s) ** 3 * prbg0
    pmech = om1 / p * (1.0 - s) * mi - prbg

    mk_full = mkwert(0, ui, st)  # IFE at converged UI (fast path is enough)
    ife = mk_full["ife"]
    ci1 = ci1 + ife
    i1 = abs(ci1)
    pcu1 = m1 * r1w * i1 ** 2
    n = (1.0 - s) * fn / p
    i2o = i2opu * i2
    i2u = i2upu * i2
    pfe1 = m1 * u1 * ife
    p1 = pcu1 + pfe1 + pzus + prbg + pmech + pcu2
    cosph = p1 / (m1 * u1 * i1)
    pcu2ny = pcu2 - pcu2ny
    mny = mi - mny
    if s == 1.0:
        m = mi
    else:
        m = mi - prbg / (om1 / p * abs(1.0 - s))

    if is1pkt == 1:
        cz2nys, cdny, i2opu, i2upu = z2nys(
            ilfr, s, p, fn, x1hp, x1hp0, nyinfo.sf[1], nyinfo.uer[1],
            xssti2, xsnut2, nyinfo.sd2[1], r2w, xssteg, li, lamnu2,
            nyinfo.xsrny[1], nyinfo.rrny[1], r2s, r2kw, h42, br42, r2kwa,
            nyinfo.xsrnya[1], nyinfo.rrnya[1], 1, st)
        st.esb.r1 = r1w
        st.esb.xsti1 = xssti1
        st.esb.xn1 = xsnut1
        st.esb.xsd1 = nyinfo.sd1 * x1hp0
        st.esb.xsd1r = nyinfo.sd1r * x1hp0

    return {
        "i1": i1, "i2": i2, "i2o": i2o, "i2u": i2u, "imy": imy, "ui": ui,
        "p1": p1, "pfe1": pfe1, "pcu1": pcu1, "pmech": pmech, "prbg": prbg,
        "pcu2": pcu2, "pcu2ny": pcu2ny, "pzus": pzus, "m": m, "mny": mny,
        "n": n, "cosph": cosph, "bs2i": bs2i, "usd1": usd1, "rsd1": rsd1,
        "isd1": isd1,
    }
