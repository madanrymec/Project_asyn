"""
Output/reporting functional group - the last major group before
PROGRAM ASYN4 itself. Most of this group (KDRUCK, PDRUCK, LISTE5,
KOPP1/4/6/7/8, KOK211, KOPASK, KOPKFE, XKEY) is FORTRAN `WRITE` statements
formatting already-computed COMMON-block values to line-printer reports
with HP PCL escape-code column positioning - little or no new
computation. This module ports the pieces with real logic in them first.

Ported so far (source line ranges refer to ASYN4.FOR):
    MWINFO (3023-3320)  "Vorabinformation an Materialwirtschaft" -
                         pre-production material-ordering info sheet.
                         Contains one genuine lookup table (V10 -> sheet-
                         steel grade name) and a solid-vs-segmented
                         lamination decision (DA1 >= 1286mm); the rest of
                         the subroutine is HP PCL report formatting to a
                         file 'MWDATEN', replaced here with a returned
                         dict (see the "output as data, not as printer
                         escape codes" pattern used throughout this port).
    SYNMOM (9353-9494)  Synchronous (cusp) torque analysis for cage
                         rotors: standstill synchronous-torque search
                         across load angle PHI0, then (if a synchronous
                         torque exists) the running-condition synchronous
                         torques at each qualifying harmonic pole-pair
                         number. This is real physics (not just
                         reporting) - determines whether the motor risks
                         getting "stuck" on a parasitic synchronous
                         torque during starting. Uses GGT, SPKT, MNYQMP,
                         NUE0, MGEGEN, and st.dp.cdp/x1hp0 (all already
                         ported/available).
    KDRUCK (1341-2054ish) "KST-Zusatzinformation" cooling/manufacturing
                         info sheet. Contains two genuine computation
                         blocks: cooling-air entry areas (LUNU1/LUSP/
                         LUR2/STDBOH/STDBOHK/LFRNG/LFRLEI and the IKLINK
                         "which cooling path dominates" decision) and
                         lamination-stack cross-section areas (AR1/AZ1/
                         AR2/AZ2, stator/rotor iron vs. copper split).
                         The rest (writing 'KSTDATEN') is report
                         formatting, replaced with a returned dict.
    KKWERT (6066-6157)  Starting-current correction factor KKA - a
                         lookup/formula selecting how much the locked-
                         rotor current/torque should be scaled down from
                         the idealized SPKT value, based on winding type
                         (Dahlander 2-zone vs. m-zone), skew, and bar
                         length vs. core length. TLAST's dependency.
    XKEY   (9614-9696)  Leakage-reactance breakdown (as % of total) and
                         the two time constants (THAUPT = rotor time
                         constant, TSTREU = leakage/transient time
                         constant) at standstill. Real computation
                         (differs materially between single-cage and
                         double-cage rotors via a complex parallel-
                         impedance combination), with the printer-report
                         write replaced by a returned dict.
    TLAST  (3906-5007)  THE motor-mode duty-cycle calculation - the
                         single largest subroutine in the program
                         (~1100 lines). Computes, in order: the no-load/
                         synchronism point, 5 partial-load points (1/4,
                         2/4, 3/4, rated, 5/4 x rated - FORTRAN J order,
                         matching GENERA's convention), the breakdown
                         (Kipp) point via SK, the theoretical locked-
                         rotor point (S=1), the KKWERT-corrected locked-
                         rotor point via SPKTSP, then XKEY's leakage
                         breakdown and time constants. Dispatches to
                         GENERA (if RKZ(2:2)=='G', generator duty) and
                         DFMAG (if RKZ(2:2)=='D', revolving-field-magnet
                         duty) exactly as the FORTRAN does. The call to
                         KOPASK (source line 4872, ported below) is
                         exposed as an optional `kopask_callback` hook,
                         since it's still useful to be able to bypass it
                         (e.g. when you don't need the compressor-duty
                         data file). All screen-position writes and the
                         many printer FORMAT statements were replaced
                         with a single returned dict.
    PDRUCK (3321-3658)  "Pruefeld-Information" (test-field data sheet):
                         torque/current at rated and reduced supply
                         voltage (with SPKTSP fallback via KKA when the
                         starting-current correction applies), then the
                         full short-circuit (locked-rotor) characteristic
                         curve across 15 fixed voltage ratios (10%-100%
                         of rated). Real, non-trivial SPKT-based
                         computation (a genuine test-certificate style
                         voltage sweep), not just formatting. The
                         FORTRAN writes two files (PRDATEN, VERTRIEB.DAT)
                         - both replaced with returned dict fields.
    KOPASK (6275-6306)  "Anbindung des Kolbenverdichterprogramms" -
                         builds a 50-point torque/current-vs-slip curve
                         for coupling to a separate piston-compressor
                         load-matching program. Real SPKT sweep (fine
                         resolution near synchronism, matching the
                         compressor program's needs), FORTRAN's
                         'ASYNASK.DAT' file write replaced with a
                         returned list of points.

*** SCOPE NOTE - external-tool data-transfer subroutines, not ported ***
KOK211, KOPKFE, KOPP1, KOPP4, KOPP6, KOPP7, KOPP8, and LISTE5 all write
input decks for SEPARATE FORTRAN programs not included in ASYN4.FOR
(a noise/vibration calculator "GERAEUSCHRECHNUNG", the "K211" mechanical
design tool, the "KFE" loss-calculation program, and similar). Their
only logic is unit conversion (m -> mm, etc.) and string formatting
before writing another program's expected input format - there is no
new physics or decision-making to port, and the receiving programs
aren't part of this codebase to validate against even if they were
ported. These are intentionally left unported; if you have those
downstream programs and need their input files generated, say so and
the specific ones you need can be ported on request.

THE ~1100-LINE TLAST SUBROUTINE - THE SINGLE LARGEST IN THE ENTIRE
PROGRAM - IS DONE, along with PDRUCK and KOPASK (this closes out every
subroutine in ASYN4.FOR with genuine physics/computation content - see
the scope note above for what remains and why it's out of scope). It's
ported without the CK211T, PINFO, CNLAST COMMON blocks (only ever
consumed by the external-tool bridges above) - their values are simply
included in TLAST's returned dict instead of being promoted to
state.py dataclasses prematurely.

Not yet ported: LESEN's screen-report tail (pure display, no new
computation), the plotting driver, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import cmath
import math

from .state import MachineState
from .winding_functions import ggt, mgegen, dxs1
from .geometry_functions import mnyqmp, nue0, k1k2
from .operating_point import spkt
from .special_points import sk, spktsp, genera, dfmag

PI = 3.14159
MY0 = 1.2566e-6


_V10_SHEET_GRADE = {
    1.15: "V290-50A", 1.25: "V310-50A", 1.35: "V330-50A", 1.50: "V350-50A",
    1.70: "V400-50A", 2.00: "V470-50A", 2.30: "V530-50A", 2.60: "V600-50A",
    3.00: "V700-50A", 3.60: "V800-50A",
}


def mwinfo(reke: str, v10: float, da1: float, di1: float,
           deltag: float) -> dict:
    """
    Pre-production material-ordering info: sheet-steel grade lookup and
    solid-vs-segmented lamination decision. Returns a dict instead of
    writing the FORTRAN's 'MWDATEN' HP-PCL report file: sheet_grade
    (str or None if V10 doesn't match a tabulated grade), stator_lam
    ('segmentiert' or 'Vollblech'), rotor_lam (same, based on the bore
    diameter net of air gap).
    """
    sheet_grade = _V10_SHEET_GRADE.get(round(v10, 2))
    stator_lam = "segmentiert" if da1 >= 1286.0 else "Vollblech"
    rotor_lam = "segmentiert" if (di1 - 2.0 * deltag) >= 1286.0 else "Vollblech"

    is_a5 = reke[0:2] == "A5"

    return {
        "sheet_grade": sheet_grade, "stator_lam": stator_lam,
        "rotor_lam": rotor_lam, "is_a5_frame": is_a5,
    }


def synmom(reke: str, ilfr: int, p: float, fn: float, r1w: float,
           xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
           r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
           r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
           u1: float, m1: float, mzone1: float, prbg0: float,
           mz1s: float, n1: float, n2: float, weite1: float, schr: float,
           mn: float, manzug: float, st: MachineState, verbose: bool = True):
    """
    Synchronous (cusp) torque analysis for cage rotors: how large a
    parasitic synchronous torque can appear at standstill (as a function
    of load angle), and, if significant, at which running speeds it
    recurs.

    Returns a dict: msymax (per-unit of MN), msydis (per-unit),
    phi0m (load angle at max), breakdown (list of (nyq, msy_per_mn) at
    standstill, only populated if msymax > 1e-6), starting_endangered
    (bool, mirrors "*** ANLAUF IST GEFAEHRDET ***"), ny0, running (list
    of dicts with nyq/s/n/msy_per_mn/m_per_mn - empty if NY0 has no
    solution, mirroring "*** KEINE SYNCHRONEN MOMENTE VORHANDEN ***").
    """
    jot = 1j
    ggteil = ggt(p * mz1s, n2)
    lam = p * mz1s * n2 / ggteil

    # ---- synchronous torque at standstill ----
    s = 1.0
    out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
              li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1, m1,
              prbg0, 0.0, 0, st)
    ui, imy, i1 = out["ui"], out["imy"], out["i1"]
    konst = m1 * ui / imy / (6.283 * fn) * i1 ** 2
    msymax = 0.0
    msydis = 0.0
    phi0m = 0.0

    for j in range(1, 73):
        phi0 = 0.08727 * j / lam
        msy = 0.0
        for g in range(1, 6):
            for gg in (1, 2):
                nyq = p + lam * g * (-1) ** gg
                if abs(nyq - p) > 5.001 * n2:
                    continue
                mnqp = mnyqmp(p, nyq, mzone1, n1, weite1, schr, n2, st)
                msyhg = konst * nyq * (mnqp * (1.0 - st.dp.cdp) *
                                        cmath.exp(jot * (nyq - p) * phi0)).imag
                msy += msyhg
                if abs(msyhg) > msydis:
                    msydis = abs(msyhg)
        if abs(msy) > msymax:
            msymax = abs(msy)
            phi0m = phi0

    if verbose:
        print(f"\n{reke}\nS Y N C H R O N E   M O M E N T E\nIM STILLSTAND")
        print(f"  MSYN/MN ={msymax / mn:7.3f}   M/MN     ={manzug / mn:7.3f}   "
              f"MSYDIS/MN ={msydis / mn:7.3f}")

    breakdown = []
    if msymax > 1.0e-6:
        for g in range(1, 6):
            for gg in (1, 2):
                nyq = p + lam * g * (-1) ** gg
                if abs(nyq - p) > 5.001 * n2:
                    continue
                mnqp = mnyqmp(p, nyq, mzone1, n1, weite1, schr, n2, st)
                msy = konst * nyq * (mnqp * (1.0 - st.dp.cdp) *
                                      cmath.exp(jot * (nyq - p) * phi0m)).imag
                breakdown.append((nyq, msy / mn))
                if verbose:
                    print(f"  NY      ={nyq:7.0f}   MSYNY/MN ={msy / mn:7.3f}")

    starting_endangered = bool(msymax >= manzug - mgegen(0.0, st))
    if verbose and starting_endangered:
        print("  *** ANLAUF IST GEFAEHRDET ***")

    ny0 = nue0(p, mz1s, n2)
    running = []
    if ny0 >= 1.0e6:
        if verbose:
            print("\nIM LAUF\n  *** KEINE SYNCHRONEN MOMENTE VORHANDEN ***")
        return {"msymax": msymax / mn, "msydis": msydis / mn,
                "phi0m": phi0m, "breakdown": breakdown,
                "starting_endangered": starting_endangered, "ny0": ny0,
                "running": running}

    for g in range(0, 6):
        for gg in (1, 2):
            if g == 0 and gg == 2:
                break
            nyq = ny0 + lam * g * (-1) ** gg
            if abs(nyq + p) > 5.001 * n2:
                continue
            s = (nyq - p) / (nyq + p)
            out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2,
                      xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42,
                      r2kwa, pn, u1, m1, prbg0, 0.0, 0, st)
            i1, ui, imy, n_val, m_val = (out["i1"], out["ui"], out["imy"],
                                          out["n"], out["m"])
            mnqp = mnyqmp(p, nyq, mzone1, n1, weite1, schr, n2, st)
            msy = (m1 * nyq * ui / imy / (6.283 * fn) * i1 ** 2 * mnqp *
                   abs(st.dp.cdp - 1.0))
            running.append({"nyq": nyq, "s": s, "n": n_val,
                             "msy_per_mn": abs(msy) / mn,
                             "m_per_mn": m_val / mn})
            if verbose:
                print(f"\nIM LAUF\n  NY ={nyq:7.0f}   S ={s:7.4f}   "
                      f"N ={n_val * 60.:7.1f} 1/MIN\n"
                      f"  MSYN/MN ={abs(msy) / mn:7.3f}   M/MN ={m_val / mn:7.3f}")

    return {"msymax": msymax / mn, "msydis": msydis / mn, "phi0m": phi0m,
            "breakdown": breakdown, "starting_endangered": starting_endangered,
            "ny0": ny0, "running": running}


def kdruck(ilfr: int, da1: float, di1: float, di2: float, l_: float,
           n1: float, n2: float, hs1: float, bs1: float, hk1: float,
           bn1: float, hs2: float, bn2s: float, br62: float,
           ho1: float, hcuo1: float, hzw1: float, hcuu1: float, hu1: float,
           ho2: float, hcuo2: float, hzw2: float, hcuu2: float, hu2: float,
           h22: float, h32: float, h42: float, h52: float, h62: float,
           h72: float, deltag: float, zk: float, bk: float, na: float,
           da: float, nk2: float, anut1: float, anut2: float,
           qstab: float, qstaba: float, drufi: float, k42: float,
           zn1: float, m1: float, ksys: float, eptauh: float,
           st: MachineState):
    """
    "KST-Zusatzinformation" cooling/manufacturing info: cooling-air
    entry areas (with the IKLINK "which cooling path is the bottleneck"
    decision) and lamination-stack cross-section areas. Returns a dict
    instead of writing the FORTRAN's 'KSTDATEN' report file.
    """
    rng1 = di1 / 2.0 + hs1 + hk1 + ho1 + hcuo1 + hzw1 + hcuu1 + hu1
    if ilfr == 1:
        # FIX: Pull hk2 from st.daten instead of trying to access it locally
        rng2 = di1 / 2.0 - deltag - hs2 - st.daten.hk2 - ho2 - hcuo2 - hzw2 - hcuu2 - hu2
    else:
        rng2 = di1 / 2.0 - deltag - hs2 - h22 - h32 - h42 - h52 - h62 - h72

    lunu1 = hs1 * bs1 * n1
    lusp = PI * di1 * deltag
    lur2 = na * PI / 4.0 * da ** 2
    azsteg = 23.0e-6 if round(bk * 1000.0) == 8 else 26.0e-6

    tn2ng = 2.0 * PI * rng2 / n2
    hg = bn2s if ilfr == 1 else br62
    lfrng = ((tn2ng - hg) * bk - azsteg) * nk2 * n2
    if round(na) == 0:
        lfrlei = di2 * PI * bk * nk2
    else:
        lfrlei = da * PI / 2.0 * na * nk2 * bk

    if hk1 < 0.0042:
        bnk1 = bn1 + 0.0034
    elif hk1 < 0.0052:
        bnk1 = bn1 + 0.0046
    elif hk1 < 0.0062:
        bnk1 = bn1 + 0.0054
    elif hk1 < 0.0082:
        bnk1 = bn1 + 0.0077
    elif hk1 < 0.0102:
        bnk1 = bn1 + 0.0100
    else:
        bnk1 = bn1 + 0.0123

    tn1k = PI * (di1 + 2.0 * hs1 + 2.0 * hk1) / n1
    stdboh = ((tn1k - bnk1) * bk - azsteg) * zk * n1
    tn1s = PI * (di1 + 2.0 * hs1) / n1
    stdbohk = ((tn1s - bn1) * bk - azsteg) * zk * n1

    iklink = 1 if (stdboh < lfrng and stdboh < lfrlei) else 0

    ar1 = PI / 4.0 * da1 ** 2 - PI * rng1 ** 2
    az1 = PI * rng1 ** 2 - PI / 4.0 * di1 ** 2 - n1 * anut1
    ar2 = PI * rng2 ** 2 - PI / 4.0 * di2 ** 2 - na * PI / 4.0 * da ** 2
    if ilfr == 1:
        az2 = PI / 4.0 * (di1 - 2.0 * deltag) ** 2 - PI * rng2 ** 2 - n2 * anut2
    else:
        az2 = PI / 4.0 * (di1 - 2.0 * deltag) ** 2 - PI * rng2 ** 2 - n2 * qstab
        if ilfr >= 3:
            az2 = az2 - n2 * qstaba

    result = {
        "lunu1": lunu1, "lusp": lusp, "lur2": lur2, "stdboh": stdboh,
        "stdbohk": stdbohk, "lfrng": lfrng, "lfrlei": lfrlei,
        "iklink": iklink, "ar1": ar1, "az1": az1, "ar2": ar2, "az2": az2,
        "rng1": rng1, "rng2": rng2,
    }

    if ilfr >= 2:
        k1, _k2, _k1r, _k1ra = k1k2(1.00, st)
        z1 = zn1 * n1 / m1 / ksys
        result.update({"k1_at_s1": k1, "cooling_time_const_inv": 1.0 / (k42 / 1.0e12),
                        "eptauh": eptauh, "z1": z1})

    return result


def kkwert(ilfr: int, m1: float, mzone1: float, n1: float, n2: float,
           lstab: float, l_: float, schr: float, p: float, weite1: float):
    """
    Starting-current correction factor KKA: how much the idealized
    locked-rotor current/torque from SPKT should be scaled down,
    based on winding type, skew, and bar length vs. core length.
    Returns kka (float).
    """
    if m1 == mzone1:
        pass  # Dahlander in high pole count -> m-zoned branch below
    else:
        if ilfr == 1:
            return 0.70 if n1 > n2 else 0.85
        if lstab - l_ - 0.001 < 0.0:
            return 0.50
        if schr == 0.0 and p > 1.001 and n2 > n1:
            return 1.00  # mirrors FORTRAN's RETURN 1 (alternate return, "no correction")
        if schr == 0.0 and p > 1.001 and n2 < n1:
            return 0.70
        if schr == 0.0 and p < 1.001:
            return 0.60
        if schr != 0.0 and n2 > n1:
            return 0.80
        if schr != 0.0 and n2 < n1:
            return 0.60
        return 1.00

    taup = n1 / (2.0 * p)
    if lstab - l_ - 0.001 < 0.0:
        return 0.50 - (weite1 - taup) / taup * 1.20
    if schr != 0.0:
        return 0.60 - (weite1 - taup) / taup * 0.80
    return 0.70 - (weite1 - taup) / taup * 0.80


def xkey(ilfr: int, i1th: float, bs2i0: float, bs2ia: float, om1: float,
         h42: float, br42: float, lamnu2: float, st: MachineState):
    """
    Leakage-reactance breakdown (% of total) and time constants at
    standstill (THAUPT = rotor time constant, TSTREU = leakage/
    transient time constant). Differs materially between single-cage
    (ILFR<3) and double-cage (ILFR>=3) rotors via a complex parallel-
    impedance combination for the latter.

    Returns a dict: stirn1, stirn2, dvk1, dvk2, schrg, steg2, nut1,
    nut2 (all as fractions, not percent - multiply by 100 for the
    FORTRAN report's convention), thaupt, tstreu.
    """
    esb = st.esb
    xs2sl = st.xs2sl
    jot = 1j

    k1, k2, k1r, k1ra = k1k2(1.000, st)
    xst2 = esb.xsteg * bs2i0 / bs2ia

    if ilfr < 3:
        if ilfr == 1:
            xnu2 = (xs2sl.xss2 / esb.xn2 * (esb.xn2 - esb.xsteg) +
                    xs2sl.xsn2 / esb.xn2 * (esb.xn2 - esb.xsteg) * k2)
        else:
            xnu2 = (esb.xn2 - esb.xsteg) * k2
        xk = (esb.xsti1 + esb.xsd1 + 2.0 * esb.xsch + esb.xn1 +
              dxs1(i1th, st) + esb.xring + esb.xsd2 + xst2 + xnu2)
        thaupt = esb.x1h / om1 / (esb.rk2 + esb.rr2 + esb.rs2)
        tstreu = xk / om1 / (esb.r1 + esb.rk2 + k1r * esb.rr2 + k1 * esb.rs2)
    else:
        hga = esb.xna - esb.xng
        hgb = ((esb.xn2 - esb.xng) * (lamnu2 * k2 + h42 / br42 + 0.79) /
               (lamnu2 * 1.0 + h42 / br42 + 0.79))
        hgz = esb.xng - esb.xsteg
        if ilfr == 3:
            cz = (1.0 / (1.0 / (esb.ra2 + jot * hga) +
                         1.0 / (esb.rk2 + esb.rs2 * k1 + jot * hgb)) +
                  k1r * esb.rr2)
        else:
            cz = (1.0 / (1.0 / (esb.ra2 + k1ra * esb.rra2 + jot * hga) +
                         1.0 / (esb.rk2 + esb.rs2 * k1 + k1r * esb.rr2 + jot * hgb)))

        xnu2 = cz.imag + hgz
        xk = (esb.xsti1 + esb.xsd1 + 2.0 * esb.xsch + esb.xn1 +
              dxs1(i1th, st) + esb.xring + esb.xsd2 + xst2 + xnu2)

        if ilfr == 3:
            r2g = 1.0 / (1.0 / esb.ra2 + 1.0 / (esb.rk2 + esb.rs2)) + esb.rr2
        else:
            r2g = 1.0 / (1.0 / (esb.ra2 + esb.rra2) +
                         1.0 / (esb.rk2 + esb.rs2 + esb.rr2))
        thaupt = esb.x1h / om1 / r2g
        tstreu = xk / om1 / (esb.r1 + cz.real)

    return {
        "stirn1": esb.xsti1 / xk, "stirn2": esb.xring / xk,
        "dvk1": esb.xsd1 / xk, "dvk2": esb.xsd2 / xk,
        "schrg": 2.0 * esb.xsch / xk, "steg2": xst2 / xk,
        "nut1": (esb.xn1 + dxs1(i1th, st)) / xk, "nut2": xnu2 / xk,
        "thaupt": thaupt, "tstreu": tstreu,
    }


def _locked_rotor_densities_and_heating(ilfr, i1, i2, i2o, qcu1, qcu2,
                                         qstab, qstaba, qring, qringa,
                                         n2, p, k41, k42, k4a2, k4r2,
                                         k4ar2, di1, zn1, n1, ksys, zn2,
                                         st):
    """
    Current densities and per-second adiabatic heating rates at a
    locked-rotor operating point (source ~4233-4298 and ~4650-4721,
    confirmed against the original FORTRAN during a cross-check against
    the Kirloskar/AEG output-parameter reference sheet - these were
    computed in TLAST but not previously surfaced by this port).

    Formulas (all confirmed directly from source):
        G1 = I1/QCU1, G2 = I2/QSTAB (or I2/QCU2 for slip-ring),
        G2A = I2O/QSTABA, GR = I2/QRING/(2*sin(P*PI/N2)) (cage) or
        I2/QCU2-based for slip-ring, GRA = I2O/QRINGA/(2*sin(P*PI/N2)).
        A1EFF = ZN1*N1*I1/(PI*DI1)/KSYS, A2EFF = ZN2*N2*I2/(PI*DI1)
        (slip-ring) or N2*I2/(PI*DI1) (cage).
        TA1 = G1**2/K41, TA2 = K1*G2**2/K42, TAA2 = G2A**2/K4A2,
        TAR = K1R*GR**2/K4R2, TARA = K1RA*GRA**2/K4AR2 (all K/sec).

    Returns a dict: g1, g2, g2a, gr, gra, a1eff, a2eff, ta1, ta2, taa2,
    tar, tara. Fields not applicable to the given ILFR (e.g. g2a/gra/
    taa2/tara for a single-cage rotor) are returned as 0.0, matching
    what the FORTRAN report would leave blank/zero in that case.
    """
    from .geometry_functions import k1k2
    k1, k2, k1r, k1ra = k1k2(1.0, st)

    g1 = i1 / qcu1 if qcu1 else 0.0
    a1eff = zn1 * n1 * i1 / (PI * di1) / ksys

    if ilfr == 1:
        g2 = i2 / qcu2 if qcu2 else 0.0
        a2eff = zn2 * n2 * i2 / (PI * di1)
        gr = g2a = gra = 0.0
    else:
        g2 = i2 / qstab if qstab else 0.0
        a2eff = n2 * i2 / (PI * di1)
        gr = (i2 / qring / (2.0 * math.sin(p * PI / n2))) if qring else 0.0
        if ilfr in (3, 4):
            g2a = i2o / qstaba if qstaba else 0.0
            gra = (i2o / qringa / (2.0 * math.sin(p * PI / n2))) if (ilfr == 4 and qringa) else 0.0
        else:
            g2a = gra = 0.0

    ta1 = g1 ** 2 / k41 if k41 else 0.0
    ta2 = k1 * g2 ** 2 / k42 if k42 else 0.0
    taa2 = g2a ** 2 / k4a2 if k4a2 else 0.0
    tar = k1r * gr ** 2 / k4r2 if k4r2 else 0.0
    tara = k1ra * gra ** 2 / k4ar2 if k4ar2 else 0.0

    # Display-unit conversion, confirmed directly from the FORTRAN
    # source's print statement (source ~4328: "GRA/1.E6,A1EFF/100."):
    # G1/G2/G2A/GR/GRA are computed here in SI (A/m^2, since qcu1 etc.
    # are already SI m^2 throughout this port) and are converted to
    # A/mm^2 (the AEG-documented unit) by /1e6; A1EFF/A2EFF are
    # computed in SI (A/m) and converted to A/cm (the documented unit)
    # by /100. TA1/TA2/TAA2/TAR/TARA (K/sec) are computed FROM the
    # already-SI g1/g2/etc (matching the FORTRAN's own internal use of
    # un-rescaled G1 for the heating-rate formulas) so those stay as-is.
    return {
        "g1": g1 / 1.0e6, "g2": g2 / 1.0e6, "g2a": g2a / 1.0e6,
        "gr": gr / 1.0e6, "gra": gra / 1.0e6,
        "a1eff": a1eff / 100.0, "a2eff": a2eff / 100.0,
        "ta1": ta1, "ta2": ta2, "taa2": taa2, "tar": tar, "tara": tara,
    }


def tlast(name: str, reke: str, ilfr: int, p: float, fn: float, r1w: float,
          xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
          r2w: float, li: float, hs2: float, bs2: float, lamnu2: float,
          r2s: float, r2kw: float, h42: float, br42: float, r2kwa: float,
          pn: float, u1: float, m1: float, prbg0: float, qcu1: float,
          qcu2: float, qstab: float, qstaba: float, qring: float,
          qringa: float, zn1: float, n1: float, zn2: float, n2: float,
          di1: float, i1n_unused: float, mn: float, kapst: float,
          k41: float, k42: float, k4a2: float, k4r2: float,
          k4ar2: float, zk: float, nwkz: float, taup: float,
          ku1_unused: float, kz: float, kt: float, ishalt: int,
          scha2: int, ksys: float, ski_unused: float, kma: float,
          kia: float, l_: float, lm1: float, lm2: float, lstab: float,
          mzone1: float, schr: float, weite1: float, st: MachineState,
          kopask_callback=None, verbose: bool = True):
    """
    Motor-mode duty-cycle calculation: no-load point, 5 partial-load
    points, breakdown (Kipp) point, theoretical + KKWERT-corrected
    locked-rotor points, and XKEY's leakage breakdown/time constants.
    Dispatches to GENERA (RKZ(2:2)=='G') and DFMAG (RKZ(2:2)=='D') as
    the FORTRAN does.

    Returns a dict with every quantity the FORTRAN passed back as a
    subroutine argument (sfeld, cfeld, efeld as 1-indexed length-6
    lists matching FORTRAN's (5) arrays, tk, eptauh, delta2, blm, alpha,
    pvsum, pvinn, ski, l/lm1/lm2/lstab/ksys pass-through, thaupt, mkipp,
    mak, i1ak, iabort), plus 'points' (list of 5 per-load-point dicts,
    ORDERED BY THE FORTRAN'S OWN ITERATION SEQUENCE - points[0] is the
    RATED (100%) point, then points[1..3] are 25%/50%/75%, then
    points[4] is 125% - NOT sorted by ascending load fraction. This
    tripped up this port's own testing once already; double-check
    which index you want rather than assuming points[3] is the rated
    point (it's actually the 75% point). Each dict has 'converged':
    True/False, and (only for the rated point, j==1) additional keys
    g1/g2/a1eff/a2eff (current densities and effective surface current
    densities, in A/mm^2 and A/cm respectively) alongside dt1/dt2/wkz/
    pvsum/pvinn), 'kipp' (the breakdown-point dict), 'theoretical' and
    'corrected' (the two locked-rotor point dicts, each also carrying a
    'densities_heating' sub-dict: g1/g2/g2a/gr/gra in A/mm^2, a1eff/
    a2eff in A/cm, ta1/ta2/taa2/tar/tara adiabatic heating rates in
    K/sec - newly added after cross-checking against the Kirloskar/AEG
    output-parameter reference sheet, confirmed against source lines
    ~4233-4298 and ~4650-4721), 'xkey_result', 'generator_result' and
    'dfmag_result' (None unless dispatched).
    """
    rkz = st.art.rkz
    esb = st.esb

    xyd1 = 0.0 if ishalt in (0, 2) else 1.0
    xyd2 = 0.0 if scha2 == 0 else 1.0

    esb.skipp = esb.r1 = esb.xsti1 = esb.xn1 = esb.xsd1 = esb.xsd1r = 0.0
    esb.xsch = esb.x1h = esb.xsd2 = esb.xring = esb.xn2 = 0.0
    esb.rk2 = esb.rr2 = esb.rs2 = esb.xsteg = esb.xng = esb.xna = 0.0
    esb.ra2 = esb.rra2 = esb.xringa = 0.0

    pmsv = [0.0, pn, 0.25 * pn, 0.50 * pn, 0.75 * pn, 1.25 * pn]

    bs2i_slfr = bs2 if ilfr == 1 else None

    if verbose:
        print(f"\n{reke}\n" + "=" * 70)

    # ---- Synchronism (no-load) ----
    s = 1.0e-9
    out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
              li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1, m1,
              prbg0, 0.0, 1, st)
    qkmax = 0.9 * m1 * u1 * out["i1"]
    bs2i0 = out["bs2i"] if ilfr != 1 else bs2i_slfr
    ui0 = out["ui"]
    dxk0 = dxs1(out["i1"], st)
    pfe0 = out["pfe1"]
    alpha0 = None  # populated from MKWERT below

    from .operating_point import mkwert
    mk = mkwert(1, out["ui"], st)
    delta2 = mk["deltap"]
    alpha0 = alpha = mk["alpha"]
    blm = mk["bp"] / (PI / 2.0)
    phi = blm * li * taup * u1 / out["ui"]
    imue = out["imy"]

    if verbose:
        print(f"SYNCHRONISMUS: UI={out['ui']:.1f} V  IMY={out['imy']:.1f} A  "
              f"PFE1={out['pfe1'] / 1000.:.3f} kW")

    # ---- Partial-load points ----
    i1n = pzusn = None
    pvsum = pvinn = None
    points = [None] * 6
    sfeld = [0.0] * 6
    cfeld = [0.0] * 6
    efeld = [0.0] * 6
    pcu2nf = [0.0] * 6
    j_to_k = {1: 4, 2: 1, 3: 2, 4: 3, 5: 5}
    iabort = 0
    nn = in_ = cosfin = bs2in = in2 = sn = None
    dt1 = dt2 = wkz = None

    for j in range(1, 6):
        pmmax = 0.0
        s = 0.0005
        ds = 0.012 if abs(pn) <= 33000.0 else 0.004
        dalt = 1.0
        pms = pmsv[j]

        converged = False
        out = None
        for _jj in range(1, 301):
            out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                      r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                      pn, u1, m1, prbg0, 0.0, 0, st)
            if out["pmech"] > pmmax:
                pmmax = out["pmech"]
            if abs((pms - out["pmech"]) / pms) > 0.0040:
                dneu = (pms - out["pmech"]) / pms
                if dneu * dalt < 0.0:
                    ds = -0.5 * ds
                s += ds
                dalt = dneu
                continue
            converged = True
            break

        if not converged:
            if verbose:
                print(f"\n{int(4.01 * pms / pn)}/4-LAST: KEINE KONVERGENZ "
                      f"(PSOLL={pms / 1000.:.1f} kW, PMAX={pmmax / 1000.:.1f} kW)")
            if j == 1:
                return {"iabort": 1}
            k = j_to_k[j]
            points[j] = {"converged": False, "j": j, "psoll": pms,
                         "pmax": pmmax}
            continue

        i1, pmech, pcu1, pfe1, m_val, s_val, cosph = (
            out["i1"], out["pmech"], out["pcu1"], out["pfe1"], out["m"],
            s, out["cosph"])
        if j == 1:
            i1n = i1
            pzusn = 0.005 * out["p1"]
        pdelta = out["p1"] - pcu1 - pfe0
        pku2 = s * pdelta
        pzus = pzusn * (i1 / i1n) ** 2
        eta = pmech / (pmech + pcu1 + pku2 + prbg0 + pfe0 + pzus)

        extra = {}
        if j == 1:
            nn, in_, cosfin, bs2in, in2, sn = (
                out["n"], i1, cosph, out["bs2i"], out["i2u"], s)
            g1 = i1 / qcu1
            a1eff = zn1 * n1 * i1 / (PI * di1) / ksys
            if kz > 0.0:
                dt1 = kz * g1 * a1eff * math.sqrt((pcu1 + pfe1) / pcu1)
            else:
                dt1 = kt * (2.0 * pcu1 + out["pcu2"] + pfe1)

            if ilfr == 1:
                g2 = out["i2u"] / qcu2
            elif ilfr == 2:
                g2 = out["i2u"] / qstab
            else:
                g2 = out["i2u"] / qstab
            a2eff = (zn2 * n2 * out["i2"] / (PI * di1) if ilfr == 1
                     else n2 * out["i2"] / (PI * di1))
            if kz > 0.0:
                dt2 = kz * g2 * a2eff if ilfr == 1 else kz * g2 * a2eff * (44.0e6 / kapst) ** 0.75
                ddt = dt2 - dt1
                if ddt < 0.0:
                    ddt = 0.0
                dt2 = dt2 - 0.67 * ddt
            else:
                hg = (1.0 + g2 * a2eff / g1 / a1eff) / 2.0
                hg = 1.3 * min(hg, 1.0)
                dt2 = hg * dt1 if ilfr == 1 else hg * dt1 * (44.0e6 / kapst) ** 0.75

            wkz = ((2.0 * pcu1 + out["pcu2"] + pfe1 + pzus) / nwkz
                   if round(zk) == 0 else (pcu1 + 0.2 * pfe1) / nwkz)
            pvsum = pfe0 + prbg0 + pcu1 + out["pcu2"] + pzus
            if round(zk) == 0:
                pvinn = pvsum - prbg0
            else:
                if (len(reke) > 2 and (reke[1] == "R" or reke[2] == "R") or
                        (len(reke) > 21 and (reke[16:18] == "LK" or reke[20:22] == "LK"))):
                    pvinn = pvsum - 0.15 * prbg0
                else:
                    pvinn = pvsum
            if ilfr == 1:
                pvinn = pvinn - 0.2 * prbg0
            extra = {"dt1": dt1, "dt2": dt2, "wkz": wkz, "pvsum": pvsum,
                     "pvinn": pvinn, "g1": g1 / 1.0e6, "g2": g2 / 1.0e6,
                     "a1eff": a1eff / 100.0, "a2eff": a2eff / 100.0}
            if verbose:
                print(f"NENNLAST: DT1={dt1:.1f} K  DT2={dt2:.1f} K  WKZ={wkz:.1f}")
            # pvsum/pvinn assigned to the outer-scope variables (declared
            # before this loop) so the final return statement can use
            # them without depending on points[]'s K-remapped indexing.

        result = dict(out)
        result.update({"converged": True, "j": j, "eta": eta, "pku2": pku2,
                        "pzus": pzus, **extra})
        k = j_to_k[j]
        points[j] = result
        sfeld[k] = s
        cfeld[k] = cosph
        efeld[k] = eta
        pcu2nf[k] = out["pcu2ny"] / 1000.0
        if verbose:
            print(f"  {int(4.01 * pms / pn)}/4-Last: P={pmech / 1000.:.1f} kW  "
                  f"I1={i1:.1f} A  ETA={eta:.4f}  COSPH={cosph:.3f}  S={s:.4f}")

    # ---- Breakdown (Kipp) point ----
    skipp = sk(bs2i0, bs2in, ilfr, p, fn, r1w, xssti1, xsnut1, xssti2,
              xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
              pn, u1, m1, prbg0, st)
    kipp_out = spkt(ilfr, skipp, p, fn, r1w, xssti1, xsnut1, xssti2,
                    xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42,
                    r2kwa, pn, u1, m1, prbg0, 0.0, 0, st)
    bs2ik = kipp_out["bs2i"]
    dxkk = dxs1(kipp_out["i1"], st)
    mkipp = kipp_out["m"]
    ski = skipp
    if verbose:
        print(f"\nKIPPUNKT: M/MN={mkipp / mn:.2f}  I1/I1N={kipp_out['i1'] / i1n:.2f}  "
              f"SKIPP={skipp * 100.:.2f}%")

    # ---- Theoretical locked-rotor point (S=1) ----
    s = 1.0
    th_out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                 r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                 u1, m1, prbg0, 0.0, 0, st)
    bs2ia = th_out["bs2i"]
    dxka = dxs1(th_out["i1"], st)
    mak = th_out["m"]
    i1ak = th_out["i1"]
    manzug = th_out["m"]
    i1a = th_out["i1"]
    i1th = th_out["i1"]

    theo_dh = _locked_rotor_densities_and_heating(
        ilfr, th_out["i1"], th_out["i2"], th_out["i2o"], qcu1, qcu2,
        qstab, qstaba, qring, qringa, n2, p, k41, k42, k4a2, k4r2,
        k4ar2, di1, zn1, n1, ksys, zn2, st)

    if verbose:
        print(f"STILLSTAND (THEORIE): M={manzug:.0f} Nm  I1={i1a:.1f} A")

    # ---- SKNY: populate slip list around each synchronous torque region ----
    from .special_points import skny
    skny(ilfr, skipp, th_out["ui"] / th_out["imy"], hs2, th_out["bs2i"],
         xssti2, xsnut2, r2w, li, lamnu2, r2s, r2kw, h42, br42, fn, p,
         st.dp.x1hp0, st)

    # ---- Corrected locked-rotor point (KKWERT + SPKTSP) ----
    kka = kkwert(ilfr, m1, mzone1, n1, n2, lstab, l_, schr, p, weite1)
    corr_out = spktsp(ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                      r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                      pn, u1, m1, kka, st)
    mak = corr_out["m"]
    i1ak = corr_out["i1"]
    if verbose:
        print(f"STILLSTAND (KORRIGIERT, KKA={kka:.2f}): M={mak:.0f} Nm  "
              f"I1={i1ak:.1f} A")

    corr_dh = _locked_rotor_densities_and_heating(
        ilfr, corr_out["i1"], corr_out["i2"], corr_out["i2o"], qcu1, qcu2,
        qstab, qstaba, qring, qringa, n2, p, k41, k42, k4a2, k4r2,
        k4ar2, di1, zn1, n1, ksys, zn2, st)

    kma_out = mak / manzug
    kia_out = i1ak / i1th

    # ---- XKEY: leakage breakdown / time constants ----
    xkey_result = xkey(ilfr, i1th, bs2i0, bs2ia, 2.0 * PI * fn, h42, br42,
                       lamnu2, st)
    thaupt = xkey_result["thaupt"]

    # ---- TK and 1+tau/2 for construction ----
    xsiga = u1 / th_out["i1"]
    xsign = (esb.xsti1 + esb.xn1 + dxs1(in_, st) + esb.xsd1 + 2.0 * esb.xsch +
             esb.xring + esb.xsd2 + esb.xn2 - esb.xsteg +
             esb.xsteg * bs2i0 / bs2in)
    r2strn = esb.rk2 + esb.rs2 + esb.rr2
    tk = (xsiga / esb.r1 + 1.414 * xsign / r2strn) / (2.0 * PI * fn * (1.0 + 1.414))
    xsig1a = esb.xsti1 + esb.xn1 + esb.xsd1
    eptauh = math.sqrt(esb.r1 ** 2 + (esb.x1h + xsig1a) ** 2) / esb.x1h

    # ---- Generator / Drehfeldmagnet dispatch ----
    generator_result = None
    dfmag_result = None
    if len(rkz) > 1 and rkz[1] == "G":
        generator_result = genera(
            reke, ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
            li, hs2, bs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, u1, m1,
            prbg0, pn, pmsv, skipp, pfe0, kz, kt, qcu1, zn1, n1, di1,
            ishalt, ksys, st, verbose=verbose)
    if len(rkz) > 1 and rkz[1] == "D":
        dfmag_result = dfmag(
            reke, ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
            li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1, m1,
            esb.x1h, st, verbose=verbose)

    if kopask_callback is not None:
        kopask_callback(name, reke, ilfr, p, fn, r1w, xssti1, xsnut1,
                        xssti2, xsnut2, r2w, li, hs2, bs2, lamnu2, r2s,
                        r2kw, h42, br42, r2kwa, u1, m1, prbg0, sn, mkipp,
                        skipp)

    return {
        "iabort": iabort, "points": points[1:], "sfeld": sfeld[1:6],
        "cfeld": cfeld[1:6], "efeld": efeld[1:6], "pcu2nf": pcu2nf[1:6],
        "kipp": {"skipp": skipp, "m": mkipp, "i1": kipp_out["i1"],
                 "bs2i": bs2ik, "out": kipp_out},
        "theoretical": {"m": manzug, "i1": i1a, "bs2i": bs2ia, "out": th_out,
                        "densities_heating": theo_dh},
        "corrected": {"m": mak, "i1": i1ak, "kka": kka, "out": corr_out,
                     "densities_heating": corr_dh},
        "kma": kma_out, "kia": kia_out, "i1a": i1a, "i1ak": i1ak,
        "taup": taup, "pfe0": pfe0, "ui0": ui0, "k1r1w": None,
        "kz": kz, "kt": kt, "sfeld_full": sfeld, "tk": tk,
        "eptauh": eptauh, "delta2": delta2, "blm": blm, "alpha": alpha,
        "pvsum": pvsum, "pvinn": pvinn, "ski": ski, "l": l_, "lm1": lm1,
        "lm2": lm2, "lstab": lstab, "ksys": ksys, "thaupt": thaupt,
        "mkipp": mkipp, "mak": mak, "i1ak_final": i1ak,
        "xkey_result": xkey_result, "generator_result": generator_result,
        "dfmag_result": dfmag_result, "dt1": dt1, "dt2": dt2, "wkz": wkz,
        "qkmax": qkmax,
    }


_KUFELD = [0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
           0.65, 0.70, 0.80, 0.90, 1.00]


def kopask(name: str, reke: str, ilfr: int, p: float, fn: float,
           r1w: float, xssti1: float, xsnut1: float, xssti2: float,
           xsnut2: float, r2w: float, li: float, hs2: float, bs2: float,
           lamnu2: float, r2s: float, r2kw: float, h42: float,
           br42: float, r2kwa: float, u1: float, m1: float, prbg0: float,
           sn: float, mkipp: float, skipp: float, pn: float,
           st: MachineState):
    """
    Builds a 50-point torque/current-vs-slip curve for coupling to a
    separate piston-compressor load-matching program ("Kolbenverdichter-
    programm"). Fine resolution near synchronism (slip steps of
    0.04*SN), matching what that program needs.

    Returns a list of 50 dicts: {'pmech', 'p1', 's', 'i1'}, in place of
    the FORTRAN's 'ASYNASK.DAT' file write.
    """
    bs2i = bs2 if ilfr == 1 else None
    points = []
    for j in range(1, 51):
        s = 0.04 * sn * j
        out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                  r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                  u1, m1, prbg0, 0.0, 0, st)
        points.append({"pmech": out["pmech"], "p1": out["p1"], "s": s,
                        "i1": out["i1"]})
    return points


def pdruck(irenr: int, ilfr: int, p: float, fn: float, r1w: float,
           xssti1: float, xsnut1: float, xssti2: float, xsnut2: float,
           r2w: float, li: float, hs2: float, lamnu2: float, r2s: float,
           r2kw: float, h42: float, br42: float, r2kwa: float, pn: float,
           u1: float, m1: float, prbg0: float, zn1: float, n1: float,
           xsip: float, n2: float, qstab: float, k42: float,
           eptauh: float, r1k: float, r2k: float, ku1: float,
           xnetz: float, ksys: float, kka_full: float, ishalt: int, st: MachineState):
    """
    "Pruefeld-Information" test-field data sheet: torque/current at
    rated and reduced supply voltage, then the full short-circuit
    characteristic.
    """
    from .special_points import spktsp
    rkz = st.art.rkz
    esb = st.esb

    is_gen = len(rkz) > 1 and rkz[1] == "G"

    # Exact mapping of Fortran XYD1 logic for Line/Phase scaling
    xyd1 = 0.0 if ishalt in (0, 2) else 1.0
    v_scale = math.sqrt(3.0) if xyd1 == 0.0 else 1.0
    i_scale = 1.0 if xyd1 == 0.0 else math.sqrt(3.0)

    s = -esb.skipp if is_gen else esb.skipp
    u1red = 1.00 * u1

    out = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
               r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
               u1red, m1, prbg0, 0.0, 0, st)

    s = 1.0
    if kka_full == 1.0:
        out_lr = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                      r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                      pn, u1red, m1, prbg0, 0.0, 0, st)
    else:
        out_lr = spktsp(ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                        r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42,
                        r2kwa, pn, u1red, m1, kka_full, st)

    result = {
        "kipp_m": out["m"], "kipp_i1": out["i1"] * i_scale,
        "locked_rotor_m": out_lr["m"], "locked_rotor_i1": out_lr["i1"] * i_scale,
        "reduced_voltage": None, "short_circuit_curve": [],
    }

    # ---- Reduced Voltage Starting Check ----
    if ku1 > 0.0:
        u1red2 = ku1 * u1
        s = 1.0e-8
        spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
             li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1red2,
             m1, prbg0, 0.0, 1, st)
        s = -esb.skipp if is_gen else esb.skipp
        out_r = spkt(ilfr, s, p, fn, r1w, xssti1 + xnetz, xsnut1, xssti2,
                     xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42,
                     r2kwa, pn, u1red2, m1, prbg0, 0.0, 0, st)

        s = 1.0
        dkka = 1.0 - kka_full
        kka_red = 1.0 - dkka * math.sqrt(ku1)

        if kka_full == 1.0:
            out_lrr = spkt(ilfr, s, p, fn, r1w, xssti1 + xnetz, xsnut1,
                           xssti2, xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw,
                           h42, br42, r2kwa, pn, u1red2, m1, prbg0, 0.0, 0, st)
        else:
            out_lrr = spktsp(ilfr, p, fn, r1w, xssti1 + xnetz, xsnut1,
                             xssti2, xsnut2, r2w, li, hs2, lamnu2, r2s, r2kw,
                             h42, br42, r2kwa, pn, u1red2, m1, kka_red, st)

        sini = -math.sin(math.acos(out_lrr["cosph"]))
        u1str = math.sqrt((u1red2 + xnetz * out_lrr["i1"] * sini) ** 2 +
                          (xnetz * out_lrr["i1"] * out_lrr["cosph"]) ** 2)
        result["reduced_voltage"] = {
            "u1str_over_u1": u1str / u1, "i1": out_lrr["i1"] * i_scale,
            "m": out_lrr["m"], "kipp_m": out_r["m"], "kipp_i1": out_r["i1"] * i_scale,
        }

    # ---- Short-circuit characteristic across 15 voltage ratios ----
    dkka = 1.0 - kka_full
    for ratio in _KUFELD:
        u1red3 = ratio * u1
        s = 1.0e-8
        spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
             li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn, u1red3,
             m1, prbg0, 0.0, 1, st)

        s = 1.0
        kka_step = 1.0 - dkka * math.sqrt(ratio)

        if kka_full == 1.0:
            out_sc = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                          r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                          pn, u1red3, m1, prbg0, 0.0, 0, st)
        else:
            out_sc = spktsp(ilfr, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                            r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa,
                            pn, u1red3, m1, kka_step, st)

        result["short_circuit_curve"].append({
            "u1": u1red3 * v_scale,
            "i1": out_sc["i1"] * i_scale,
            "m": out_sc["m"],
            "cosph": out_sc["cosph"],
            "p1_kw": out_sc["p1"] / 1000.0,
            "s_kva": out_sc["p1"] / 1000.0 / out_sc["cosph"] if out_sc["cosph"] else 0.0,
            "kka": kka_step
        })

    # ---- Restore nominal saturation, then no-load summary ----
    s = 1.0e-8
    out_nl = spkt(ilfr, s, p, fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                  r2w, li, hs2, lamnu2, r2s, r2kw, h42, br42, r2kwa, pn,
                  u1, m1, prbg0, 0.0, 1, st)
    p1 = out_nl["p1"] + prbg0
    cosph = p1 / (m1 * u1 * out_nl["i1"])
    result["no_load"] = {"i1": out_nl["i1"] * i_scale, "cosph": cosph,
                         "p1_kw": p1 / 1000.0, "pfe1_kw": out_nl["pfe1"] / 1000.0,
                         "prbg0_kw": prbg0 / 1000.0}

    if ilfr != 1:
        k1, k2, k1r, k1ra = k1k2(1.00, st)
        z1 = n1 * zn1 / m1 / ksys
        result["eexe"] = {"z1": z1, "xsip": xsip, "m1": m1, "n2": n2,
                          "qstab_mm2": qstab * 1.0e6, "k1": k1,
                          "inv_k42": 1.0 / (k42 / 1.0e12), "eptauh": eptauh}

    result["r1k"] = r1k
    result["r2k"] = r2k
    return result