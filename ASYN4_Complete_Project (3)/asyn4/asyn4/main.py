"""
PROGRAM ASYN4 (source 1-869 of ASYN4.FOR) - the main driver, ported as
`run_asyn4()`. Wires together, in the FORTRAN's own order:

    LESEN -> BRUVOR -> [geometry pre-calc] -> PAKET -> STD -> KIKR1 ->
    SLFR+KIKRSL (slip-ring) or KIKR (cage) -> KREIS -> STEG (cage) ->
    POPAZA -> NUT (stator, + rotor if slip-ring) -> [mass/inertia calc]
    -> [K1R1W form-coil skin-effect correction] -> TLAST -> MAGZUG ->
    KDRUCK -> SLAST -> HOLAUF (x1-3) -> PDRUCK -> SYNMOM (cage only)

Deliberately NOT called (see CONVERSION_STATUS.md / reports.py's scope
note): the plotting driver (CALL PLOTS/PLOT/PLOTT - not ported), and the
external-tool data-transfer bridges KOPP1/4/6/7/8, KOK211, KOPKFE,
LISTE5 (write input decks for separate programs not included in this
codebase).

This is a big, mostly-linear orchestration function rather than many
small ones, mirroring the FORTRAN main program's own structure - that
structure is what's being ported here, not reinvented. A few small
pieces of inline FORTRAN logic that live directly in PROGRAM ASYN4
(never their own subroutine) are extracted as private helpers so they
can be tested in isolation: `_geometry_precalc` (source 163-237, the
LI/LFE/TAUP/KSYS/W1/XSIP/UEU/HOJOC2/RNG2/double-cage-QSTABA/DRING
block) and `_stator_form_coil_correction` (source 599-633, K1R1W - the
skin-effect resistance increase for stator form coils, which turned out
to have no subroutine of its own in the original).
"""
from __future__ import annotations
import math

from .state import MachineState
from .lesen import lesen
from .geometry_functions import bruvor, paket, k1k2
from .geometry_functions import xsiz, xsis, nut
from .winding_functions import xsiz2
from .magnetic_circuit import std, steg
from .rotor_functions import slfr, kikr1, kikr, kikrsl
from .magnetics_solver import kreis
from .harmonics import popaza
from .load_curves import slast, holauf
from .reports import tlast, kdruck, pdruck
from .magzug import magzug

PI = 3.14159
MY0 = 1.2566e-6


def _geometry_precalc(d, rkz: str):
    """
    Source 163-237: pre-calculation of stack lengths, winding factor,
    and (depending on rotor type) the rotor tooth/back-iron reference
    radius, plus double-cage QSTABA/R2KWA defaults and DRING/DRINGA.
    Mutates `d` (st.daten) in place for the fields the rest of the
    program reads directly out of it, and returns the values only ever
    used locally in run_asyn4 (li, lfe, taup, ksys, w1, xsip, vzxsiz,
    fw1, z1, z2, ueu, hojoc2, rng2, dring, dringa, r2kwa, lz2).
    """
    if len(rkz) > 9 and rkz[9] == " ":
        li = d.l - d.zk * d.bk ** 2 / (d.bk + 2.5 * d.deltag)
    else:
        li = d.l - 2.0 * d.zk * d.bk ** 2 / (d.bk + 5.0 * d.deltag)
    lfe = d.l - d.zk * d.bk
    taup = PI * d.di1 / (2.0 * d.p)
    ksys = 2.0 if (len(rkz) > 10 and rkz[10] != " ") else 1.0

    w1 = d.zn1 * d.n1 / (2.0 * d.m1) / ksys
    q1 = d.n1 / (d.mzone1 * d.p)

    class _V:
        pass
    class _B:
        pass
    stub = MachineState()
    stub.vorzei.bs1di1 = 0.000001
    stub.vorzei.ksys = ksys

    xsip = xsiz(d.p, q1, d.n1, stub) * xsis(d.p, d.weite1, d.n1)
    fw1 = xsip
    z1 = d.n1 * d.zn1 / d.m1 / ksys
    z2 = d.n2 * d.zn2 / 3.0 if d.ilfr == 1 else 0.0

    vzxsiz = 1.0
    if xsip < 0.0:
        xsip = -xsip
        vzxsiz = -1.0

    if d.ilfr == 1:
        q2 = d.n2 / (6.0 * d.p)
        ueu = (xsip / (xsiz2(d.p, q2, d.n2) * xsis(d.p, d.weite2, d.n2)) *
               d.zn1 * d.n1 / (d.zn2 * d.n2))
        lz2 = d.hs2 + d.hk2 + d.ho2 + d.hcuo2 + d.hzw2 + d.hcuu2 + d.hu2
        hojoc2 = d.di1 / 2.0 - d.deltag - lz2 - d.di2 / 2.0
        rng2 = d.di2 / 2.0 + hojoc2
    else:
        ueu = xsip * 2.0 * w1
        da2 = d.di1 - 2.0 * d.deltag
        lz2 = d.hs2 + d.h22 + d.h32 + d.h42 + d.h52 + d.h62 + d.h72
        rng2 = (math.sqrt((da2 / 2.0) ** 2 + (lz2 - d.br62 / 2.0) ** 2 -
                (lz2 - d.br62 / 2.0) * da2 * math.cos(d.beta)) - d.br62 / 2.0)
        hojoc2 = rng2 - d.di2 / 2.0

    r2kwa = 0.0
    if d.ilfr in (3, 4):
        if d.qstaba <= 1.0e-6:
            if abs((d.h32 - d.br32) / d.br32) > 1.0e-6:
                d.qstaba = PI / 8.0 * d.br32 ** 2 + (d.h32 - 0.5 * d.br32) * d.br32
            else:
                d.qstaba = PI / 4.0 * d.br32 ** 2
        if d.ilfr == 3:
            r2kwa = d.lstab / (d.kappa3 * d.qstaba)
        else:
            r2kwa = d.lstaba / (d.kappa3 * d.qstaba)
    else:
        d.qstaba = 0.0

    dring = dringa = 0.0
    if d.ilfr in (2, 3):
        dring = d.di1 - 2.0 * (d.deltag + d.hs2) - (d.h22 + d.h32 + d.h42 + d.h52 + d.h62 + d.h72)
        dringa = 0.0
    if d.ilfr == 4:
        dring = d.di1 - 2.0 * (d.deltag + d.hs2 + d.h22 + d.h32 + d.h42) - (d.h52 + d.h62 + d.h72)
        dringa = d.di1 - 2.0 * (d.deltag + d.hs2 + d.h22) - d.h32

    return {
        "li": li, "lfe": lfe, "taup": taup, "ksys": ksys, "w1": w1,
        "xsip": xsip, "fw1": fw1, "z1": z1, "z2": z2, "vzxsiz": vzxsiz,
        "ueu": ueu, "hojoc2": hojoc2, "rng2": rng2, "lz2": lz2,
        "r2kwa": r2kwa, "dring": dring, "dringa": dringa,
    }


def _stator_form_coil_correction(d, lfe: float, fn: float) -> float:
    """
    Source 599-633: K1R1W, the skin-effect resistance-increase factor
    for stator form coils built from multiple parallel partial
    conductors - inline logic in the FORTRAN main program with no
    subroutine of its own.
    """
    if d.izgf != 1:
        return 1.0
    if d.mtl * d.ntl * d.utl * d.azweig < 0.001:
        return 1.0

    tpar = d.mtl * d.ntl / (d.zn1 * d.azweig)
    krp = 1.0
    if tpar / d.ntl >= 1.001:
        lambda_ = (d.lm1 - lfe) / lfe
        alpha = math.sqrt((d.ntl * d.btl * PI * fn * MY0 * 57.0e6 /
                           (1.0 + 0.0043 * d.theta1)) / (d.bn1 * (1.0 + lambda_)))
        ws = d.zn1 * d.azweig / 2.0
        hlage = d.mtl * d.htl * d.utl / (2.0 * ws)
        xsi = alpha * hlage
        k1 = xsi * (math.sinh(2.0 * xsi) + math.sin(2.0 * xsi)) / \
                   (math.cosh(2.0 * xsi) - math.cos(2.0 * xsi))
        k3 = 2.0 * xsi * (math.sinh(xsi) - math.sin(xsi)) / \
                         (math.cosh(xsi) + math.cos(xsi))
        krp = k1 + ((ws / d.utl) ** 2 - 1.0) * k3 / 4.0

    alpha1 = math.sqrt((d.ntl * d.btl * PI * fn * MY0 * 57.0e6 /
                        (1.0 + 0.0043 * d.theta1)) / d.bn1)
    xi1 = alpha1 * d.htl
    ws = d.zn1 * d.azweig / 2.0
    wzutau = d.weite1 / (d.n1 / 2.0 / d.p)
    hg = 0.67 + (wzutau - 0.37) / 0.63 * 0.33
    pt = 20.0 / 45.0 * ws ** 2 * hg * xi1 ** 4
    if round(tpar) == 1 and round(d.ntl) == 2:
        pt = pt / 4.0
    lambda_ = (d.lm1 - lfe) / lfe
    krstr = ((1.0 + pt) + lambda_) / (1.0 + lambda_)
    return krp + krstr - 1.0


def _finish_lesen_postprocessing(d, st, reke: str):
    """
    Source 2860-2931 of ASYN4.FOR: three genuine computational blocks
    that live in PROGRAM ASYN4's own body, past LESEN's documented
    card-reading boundary (source ~1586-2054) but BEFORE the pure-
    display report tail - these were missed in the original LESEN port
    and only found via this integration test, so they're applied here
    as a post-processing step rather than folded back into lesen.py
    (keeping lesen.py focused on the card reads themselves).

    1. K41-K4AR2: heating time constants (stator/rotor/ring), material-
       and rotor-type-dependent (source 2868-2897, "ERWAERNUNGSKONSTANTE").
    2. PRBG0 default: Nurnberg friction-loss formula, only applied if
       the input file left PRBG0 at 0 (source 2903-2914,
       "REIBVERLUSTE NACH NUERNBERG").
    3. RKZ(8)='X' override for fractional-slot slip-ring machines
       (source 2920-2925) - suppresses harmonic-field damping for a
       combination POPAZA's damping-selection logic isn't valid for.

    Also runs PRUEF (source 2931) - the input-data validation check
    that was ported long ago (special_points.pruef) but never actually
    wired into any pipeline until this integration test surfaced the
    gap. Returns (iabort, messages) from PRUEF; the caller should treat
    iabort==1 as "stop, don't trust the results".
    """
    from .special_points import pruef

    kapst = max(d.kappa5, d.kappa6, d.kappa7)
    d.k41 = 386.0 * 8900.0 * 57.0e6 / (1.0 + 0.00429 * d.theta1)
    if d.ilfr == 1:
        d.k42 = 386.0 * 8900.0 * 57.0e6 / (1.0 + 0.00429 * d.theta2)
        d.k4a2 = d.k4r2 = d.k4ar2 = 1.0
    else:
        k5 = 910.0 * 2700.0 if d.lue < 0.001 else 386.0 * 8900.0
        if d.ilfr == 2:
            d.k42 = k5 * kapst
            d.k4a2 = 1.0
            d.k4r2 = k5 * d.kapri
            d.k4ar2 = 1.0
        elif d.ilfr == 3:
            d.k42 = k5 * kapst
            d.k4a2 = k5 * d.kappa3
            d.k4r2 = k5 * d.kapri
            d.k4ar2 = 1.0
        elif d.ilfr == 4:
            d.k42 = k5 * kapst
            d.k4a2 = k5 * d.kappa3
            d.k4r2 = k5 * d.kapri
            d.k4ar2 = k5 * d.kapria

    prbgi = 8.0 * d.di1 * (d.l + 0.15) * (PI * d.di1 * d.fn / d.p) ** 2
    if (len(reke) > 21 and (reke[1] == "R" or reke[2] == "R" or
                             reke[16:18] == "LK" or reke[20:22] == "LK")):
        hg = 1.35
    else:
        hg = 1.00
    if d.ilfr == 1:
        hg = hg * 1.2
    prbgi = hg * prbgi
    if d.prbg0 == 0.0:
        d.prbg0 = prbgi

    rkz = st.art.rkz
    if len(rkz) > 0 and rkz[0] == "1":
        q1 = d.n1 / (2.0 * d.p * d.m1)
        if math.fmod(q1 + 0.001, 1.0) > 0.01:
            rkz = rkz[:7] + "X" + rkz[8:]
            st.art.rkz = rkz

    iabort, messages = pruef(st)
    return iabort, messages


def run_asyn4(einasyn_path: str, irenr: int = 1, verbose: bool = True,
              run_slip_curve: bool = True, run_reduced_voltage: bool = True,
              plot_path: str | None = None):
    """
    Runs the full ASYN4 pipeline against a real EINASYN-format input
    file: reads the machine data, runs the complete geometry/winding/
    magnetic-circuit/performance calculation chain, and returns a single
    consolidated results dict.

    run_slip_curve=False skips SLAST/HOLAUF/PDRUCK/SYNMOM (the slower,
    optional parts controlled by RKZ columns 3/6 in the original) and
    just returns the core TLAST/KDRUCK results - useful for a quick
    check that a data file loads and produces a sane operating point.
    run_reduced_voltage=False skips the KU1/XNETZ reduced-starting-
    voltage sweep inside SLAST/HOLAUF/PDRUCK even if the input file
    configures one.
    plot_path, if given, writes the torque-speed/current-speed chart
    (source: PLOTT, source line 779 - see plotting.py) to that path
    once SLAST has run, matching where the FORTRAN calls PLOTT in its
    own main-program flow. None (default) skips chart generation.

    Returns a dict: lesen_result, geometry (the _geometry_precalc dict),
    k1r1w, tlast_result, magzug (ce, fax), kdruck_result, gaktiv,
    jaktiv, jfe, jcu, anut1, anut2, slast_result (msat, nsat) or None,
    holauf_result or None, pdruck_result or None, synmom_result or None
    (only for cage rotors when RKZ column 6 is set).
    """
    st = MachineState()

    lesen_result = lesen(einasyn_path, irenr, st, verbose=verbose)
    d = st.daten
    rkz = st.art.rkz

    iabort, pruef_messages = _finish_lesen_postprocessing(d, st, lesen_result["reke"])
    rkz = st.art.rkz  # _finish_lesen_postprocessing may have amended RKZ(8)
    if iabort == 1:
        if verbose:
            for is_err, text in pruef_messages:
                print(("ERROR: " if is_err else "WARN: ") + text)
        return {"lesen_result": lesen_result, "iabort": 1,
                "pruef_messages": pruef_messages}

    st.vorzei.bs1di1 = 0.000001
    mz1s, weite1_corr, a1max = bruvor(d.n1, d.p, d.mzone1, d.m1, d.weite1, st)
    d.weite1 = weite1_corr

    geo = _geometry_precalc(d, rkz)
    li, lfe = geo["li"], geo["lfe"]
    nk2 = lesen_result.get("nk2", d.zk)
    lfe2 = d.l - nk2 * d.bk
    st.vorzei.vzxsiz = geo["vzxsiz"]

    paket(rkz, d.l, d.bk, d.zk, nk2, st)

    std_out = std(d.da1, d.di1, d.n1, d.p, d.weite1, d.mzone1, d.zn1,
                  d.lm1, d.qcu1, d.m1, d.fn, li, d.theta1, mz1s,
                  d.hs1, d.hk1, d.ho1, d.hcuo1, d.hzw1, d.hcuu1, d.hu1,
                  d.bs1, d.bn1, d.bn1s, d.lamst1, st)
    hojoc1, hozah1, brzah1 = std_out["hojoc1"], std_out["hozah1"], std_out["brzah1"]
    r1w, xssti1, xsnut1 = std_out["r1w"], std_out["xssti1"], std_out["xsnut1"]

    kikr1(d.izgf, d.mtl, d.ntl, d.azweig, d.btl, d.htl, d.fn, li, d.lm1,
          d.theta1, d.bn1, d.zn1)

    hojoc2 = geo["hojoc2"]
    if d.ilfr == 1:
        slfr_out = slfr(d.di2, hojoc2, d.n2, d.p, d.weite2, d.zn2, d.lm2,
                        d.qcu2, d.fn, li, d.theta2, d.hs2, d.hk2, d.ho2,
                        d.hcuo2, d.hzw2, d.hcuu2, d.hu2, d.bs2, d.bn2,
                        d.bn2s, st)
        hozah2 = slfr_out["hozah2"]
        brzah2 = slfr_out["brzah2"]
        lamnu2, r2k, r2w = slfr_out["lamnu2"], slfr_out["r2k"], slfr_out["r2w"]
        xssti2, xsnut2 = slfr_out["xssti2"], slfr_out["xsnut2"]
        r2s = r2kw = r2kwa = 0.0
        kikrsl(0.0, 0.0, 0.0, 0.0, 0.0, d.fn, d.n1, d.p, li, d.lm2,
               d.theta2, d.bn2, d.zn2, st)
    else:
        hozah2, brzah2, kappaf, lamnu2, r2s, r2kw, qstabb = kikr(
            d.ilfr, d.h22, d.br12, d.br22, d.kappa2, d.h32, d.br32,
            d.kappa3, d.h42, d.br42, d.kappa4, d.h52, d.bstr52, d.kappa5,
            d.h62, d.br52, d.br62, d.kappa6, d.h72, d.bstr62, d.kappa7,
            d.qring, d.kapri, d.qringa, d.kapria, d.fn, d.di2, hojoc2,
            d.n2, d.beta, d.n1, d.p, li, d.lstab, st)
        r2k = r2w = xssti2 = xsnut2 = 0.0
        r2kwa = geo["r2kwa"]
        if d.qstab > 1.0e-6:
            r2s = r2s * qstabb / d.qstab
            r2kw = r2kw * qstabb / d.qstab
        else:
            d.qstab = qstabb
        if d.ilfr in (3, 4):
            kikr(0, 0.0, d.br12, d.br22, d.kappa2, 0.0, d.br32,
                 d.kappa3, 0.0, d.br42, d.kappa4, d.h52, d.bstr52,
                 d.kappa5, d.h62, d.br52, d.br62, d.kappa6, d.h72,
                 d.bstr62, d.kappa7, d.qring, d.kapri, d.qringa,
                 d.kapria, d.fn, d.di2, hojoc2, d.n2, d.beta, d.n1, d.p,
                 li, d.lstab, st)

    kreis_out = kreis(d.n1, d.n2, d.p, d.kfe, d.di1 / 2.0, li, lfe, lfe2,
                      d.da1, d.di2, hojoc1, d.kj1, hojoc2, d.kj2,
                      d.deltag, d.bs1, d.bs2, d.hs1, hozah1, brzah1,
                      hozah2, brzah2, d.cfej, d.cfez, d.v10, d.u1,
                      d.fn, geo["w1"], geo["xsip"], d.m1, st)
    kc1, kc2 = kreis_out["kc1"], kreis_out["kc2"]
    gj1, gz1, gj2, gz2 = kreis_out["gj1"], kreis_out["gz1"], kreis_out["gj2"], kreis_out["gz2"]
    bs1s = kreis_out["bs1s"]
    d.cfej, d.cfez = kreis_out["cfej"], kreis_out["cfez"]

    if d.ilfr != 1:
        inut, lamda = steg(d.ilfr, li, lfe2, d.kfe, d.hs2, d.bs2,
                           d.br12, d.br22, d.br32, d.br42, d.br52,
                           d.bstr52, st)
        st.kennli.inut[1:21] = inut[1:21]
        st.kennli.lamda[1:21] = lamda[1:21]

    if lesen_result["iflag1"] == 0:
        st.vorzei.bs1di1 = bs1s / d.di1

    fw2 = popaza(d.ilfr, d.p, d.fn, d.zn1, d.zn2, d.n1, d.n2, d.m1,
                d.mzone1, mz1s, d.weite1, d.weite2, geo["dring"],
                d.kapri, d.qring, geo["dringa"], d.kapria, d.qringa,
                d.schr, st)

    anut1, kcuth1 = nut(d.bs1, d.bn1, d.bn1s, d.hk1, d.ho1, d.hcuo1,
                        d.hzw1, d.hcuu1, d.hu1, d.qcu1, d.zn1)
    anut2 = kcuth2 = 0.0
    if d.ilfr == 1:
        anut2, kcuth2 = nut(d.bs2, d.bn2, d.bn2s, d.hk2, d.ho2, d.hcuo2,
                            d.hzw2, d.hcuu2, d.hu2, d.qcu2, d.zn2)

    # ---- mass/inertia (source 480-531) ----
    gcu1 = d.n1 * d.zn1 * d.qcu1 * d.lm1 * 8900.0
    if d.ilfr == 1:
        gcu2 = d.n2 * d.zn2 * d.qcu2 * d.lm2 * 8900.0
    else:
        rho2 = 2700.0 if (d.lstab - d.l) < 0.001 else 8900.0
        gcust = d.qstab * d.lstab * d.n2 * rho2
        gcusta = d.qstaba * d.lstaba * d.n2 * rho2 if d.ilfr >= 3 else 0.0
        gcur = PI * geo["dring"] * d.qring * rho2 * 2.0
        gcura = PI * geo["dringa"] * d.qringa * rho2 * 2.0 if d.ilfr == 4 else 0.0
        gcu2 = gcust + gcusta + gcur + gcura

    na = lesen_result.get("na", 0.0)
    da_val = lesen_result.get("da", 0.0)
    gj2_adj = gj2 - na * PI / 4.0 * da_val ** 2 * lfe2 * d.kfe * 7750.0
    gfe1 = gj1 + gz1
    gfe2 = gj2_adj + gz2
    gaktiv = gcu1 + gcu2 + gfe1 + gfe2

    jaktiv = 0.5 * PI * lfe2 * d.kfe * 7750.0 * (geo["rng2"] ** 4 - (d.di2 / 2.0) ** 4)
    lz2 = geo["lz2"]
    if d.ilfr == 1:
        jaktiv += (geo["rng2"] + lz2 / 2.0) ** 2 * gz2
        jfe = jaktiv
        jaktiv += (geo["rng2"] + lz2 / 2.0) ** 2 * gcu2
    else:
        jaktiv += (geo["rng2"] + lz2 / 2.0) ** 2 * gz2
        jfe = jaktiv
        jaktiv += (geo["dring"] / 2.0) ** 2 * (gcust + gcur)
        if d.ilfr == 3:
            jaktiv += (geo["dring"] / 2.0) ** 2 * gcusta
        elif d.ilfr == 4:
            jaktiv += (geo["dringa"] / 2.0) ** 2 * (gcusta + gcura)
    jcu = jaktiv - jfe

    k1r1w = _stator_form_coil_correction(d, lfe, d.fn)

    # ---- TLAST: motor-mode duty cycle ----
    kapst = max(d.kappa5, d.kappa6, d.kappa7)
    nwkz = (PI * d.da1 * lfe * 80.0 if round(d.zk) == 0 else
            (d.lm1 - 0.5 * d.l) * (d.hcuo1 + d.hcuu1 + d.bn1 + d.bn1s) * 2.0 * d.n1 * 80.0)

    tlast_out = tlast(lesen_result["name"], lesen_result["reke"], d.ilfr,
                      d.p, d.fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                      r2w, li, d.hs2, d.bs2, lamnu2, r2s, r2kw, d.h42,
                      d.br42, r2kwa, d.pn, d.u1, d.m1, d.prbg0, d.qcu1,
                      d.qcu2, d.qstab, d.qstaba, d.qring, d.qringa,
                      d.zn1, d.n1, d.zn2, d.n2, d.di1, 0.0, None, kapst,
                      d.k41, d.k42, d.k4a2, d.k4r2, d.k4ar2, d.zk, nwkz,
                      geo["taup"], 0.0, d.kz, d.kt, lesen_result["ishalt"],
                      lesen_result.get("scha2", 0), geo["ksys"], 0.0,
                      1.0, 1.0, d.l, d.lm1, d.lm2, d.lstab, d.mzone1,
                      d.schr, d.weite1, st, verbose=verbose)

    if tlast_out.get("iabort") == 1:
        return {"lesen_result": lesen_result, "iabort": 1}

    ce, fax = magzug(d.di1, kc1 * kc2 * d.deltag, tlast_out["delta2"],
                     tlast_out["blm"], tlast_out["alpha"], nk2, d.bk, li,
                     d.p, tlast_out.get("kipp", {}).get("m", 0.0), d.schr)

    kdruck_out = kdruck(d.ilfr, d.da1, d.di1, d.di2, d.l, d.n1, d.n2,
                        d.hs1, d.bs1, d.hk1, d.bn1, d.hs2, d.bn2s,
                        d.br62, d.ho1, d.hcuo1, d.hzw1, d.hcuu1, d.hu1,
                        d.ho2, d.hcuo2, d.hzw2, d.hcuu2, d.hu2, d.h22,
                        d.h32, d.h42, d.h52, d.h62, d.h72, d.deltag,
                        d.zk, d.bk, na, da_val, nk2, anut1, anut2,
                        d.qstab, d.qstaba, lesen_result.get("drufi", 0.0),
                        d.k42, d.zn1, d.m1, geo["ksys"],
                        tlast_out["eptauh"], st)

    result = {
        "lesen_result": lesen_result, "geometry": geo, "k1r1w": k1r1w,
        "tlast_result": tlast_out, "magzug": {"ce": ce, "fax": fax},
        "kdruck_result": kdruck_out, "gaktiv": gaktiv, "jaktiv": jaktiv,
        "jfe": jfe, "jcu": jcu, "anut1": anut1, "anut2": anut2,
        "slast_result": None, "holauf_result": None,
        "pdruck_result": None, "synmom_result": None,
    }

    if not run_slip_curve or len(rkz) <= 2 or rkz[2] == " ":
        return result

    # ---- slip-load curve, run-up, test-field sheet ----
    i1n = tlast_out["theoretical"]["i1"]
    kma = tlast_out.get("kma", 1.0)
    kia = tlast_out.get("kia", 1.0)
    rated_pt = tlast_out["points"][3]
    mn_val = rated_pt["m"] if rated_pt.get("converged") else 1.0

    ku1_eff = d.ku1 if run_reduced_voltage else 0.0
    xnetz_eff = d.xnetz if run_reduced_voltage else 0.0

    msat, nsat = slast(lesen_result["reke"], d.ilfr, ku1_eff, xnetz_eff,
                       d.p, d.fn, r1w, xssti1, xsnut1, xssti2, xsnut2,
                       r2w, li, d.hs2, d.bs2, lamnu2, r2s, r2kw, d.h42,
                       d.br42, r2kwa, d.pn, d.u1, d.m1, d.prbg0, i1n,
                       mn_val, d.qcu1, d.qcu2, d.qstab, d.qstaba,
                       d.qring, d.qringa, d.n2, kma, kia, st,
                       verbose=verbose)
    result["slast_result"] = {"msat": msat, "nsat": nsat}

    if plot_path:
        from .plotting import plot_torque_speed_curve
        plot_torque_speed_curve(lesen_result["reke"], mn_val, ku1_eff,
                                xnetz_eff, st, plot_path)
        result["plot_path"] = plot_path

    jsum = 1.2 * jaktiv + d.jantr
    holauf_out = holauf(d.ilfr, d.fn, d.p, jsum, mn_val, d.k41, d.k42,
                        d.k4a2, d.k4r2, d.k4ar2, r1w, d.m1, i1n,
                        tlast_out["ski"], 1.00, 0.00, 1, st, verbose=verbose)
    result["holauf_result"] = holauf_out

    if run_reduced_voltage and (d.ku1 >= 0.000001 or d.xnetz >= 0.000001):
        slast(lesen_result["reke"], d.ilfr, d.ku1, d.xnetz, d.p, d.fn,
              r1w, xssti1, xsnut1, xssti2, xsnut2, r2w, li, d.hs2, d.bs2,
              lamnu2, r2s, r2kw, d.h42, d.br42, r2kwa, d.pn, d.u1, d.m1,
              d.prbg0, i1n, mn_val, d.qcu1, d.qcu2, d.qstab, d.qstaba,
              d.qring, d.qringa, d.n2, kma, kia, st, verbose=verbose)
        holauf(d.ilfr, d.fn, d.p, jsum, mn_val, d.k41, d.k42, d.k4a2,
               d.k4r2, d.k4ar2, r1w, d.m1, i1n, tlast_out["ski"], d.ku1,
               d.xnetz, 1, st, verbose=verbose)

    mgmerk = d.mg[0]
    d.mg[0] = 0.0
    th0_out = holauf(d.ilfr, d.fn, d.p, 1.2 * jaktiv, mn_val, d.k41,
                     d.k42, d.k4a2, d.k4r2, d.k4ar2, r1w, d.m1, i1n,
                     tlast_out["ski"], 1.00, 0.00, 0, st, verbose=False)
    d.mg[0] = mgmerk
    th0 = th0_out["th"]

    pdruck_out = pdruck(irenr, d.ilfr, d.p, d.fn, r1w, xssti1, xsnut1,
                        xssti2, xsnut2, r2w, li, d.hs2, lamnu2, r2s,
                        r2kw, d.h42, d.br42, r2kwa, d.pn, d.u1, d.m1,
                        d.prbg0, d.zn1, d.n1, geo["xsip"], d.n2, d.qstab,
                        d.k42, tlast_out["eptauh"], std_out["r1k"], r2k,
                        d.ku1, d.xnetz, geo["ksys"],
                        tlast_out["corrected"]["kka"], st)
    result["pdruck_result"] = pdruck_out

    if not (len(rkz) > 5 and rkz[5] == " ") and d.ilfr != 1:
        from .reports import synmom
        synmom_out = synmom(lesen_result["reke"], d.ilfr, d.p, d.fn,
                            r1w, xssti1, xsnut1, xssti2, xsnut2, r2w,
                            li, d.hs2, lamnu2, r2s, r2kw, d.h42, d.br42,
                            r2kwa, d.pn, d.u1, d.m1, d.mzone1, d.prbg0,
                            mz1s, d.n1, d.n2, d.weite1, d.schr, mn_val,
                            tlast_out["theoretical"]["m"], st, verbose=verbose)
        result["synmom_result"] = synmom_out

    return result
