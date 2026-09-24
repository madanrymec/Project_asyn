"""
Regression test suite for the ASYN4 Python port.

Runs with pytest (recommended: `pip install pytest`, then `pytest tests/`
from the project root) but every test is a plain function with `assert`
statements, so it will also run with the bundled `run_tests.py` script if
pytest isn't installed.

Expected values baked into these tests were captured from verified
smoke-test runs during development (see CONVERSION_STATUS.md). They are
regression checks (catch anything breaking as more of the program gets
ported) - they are NOT yet validated against a real machine's known-
correct FORTRAN output, since no sample EINASYN/AUSASYN file has been
available. See CONVERSION_STATUS.md "Open questions".
"""
import math
import os
import csv
import numpy as np
import pytest

from asyn4.state import MachineState
from asyn4 import winding_functions as wf
from asyn4 import geometry_functions as gf
from asyn4 import magnetic_circuit as mc
from asyn4 import rotor_functions as rf
from asyn4 import magnetics_solver as ms
from asyn4 import harmonics as hm
from asyn4 import operating_point as op
from asyn4 import load_curves as lc
from asyn4 import special_points as sp
from asyn4 import insulation as ins
from asyn4 import lesen as ls
from asyn4 import reports as rp
from asyn4 import report_export as rex
from asyn4 import main as asyn4_main
from asyn4 import plotting as pl


def make_state():
    st = MachineState()
    st.vorzei.bs1di1 = 0.000001
    st.vorzei.ksys = 1.0
    st.art.rkz = " " * 12
    return st


# ---------------------------------------------------------------- winding_functions

def test_xsis():
    assert wf.xsis(1, 0.9, 36) == pytest.approx(0.07845902959260323)


def test_xsisch_no_skew():
    assert wf.xsisch(1, 0) == pytest.approx(1.0)


def test_xsiz():
    st = make_state()
    assert wf.xsiz(1, 3, 36, st) == pytest.approx(0.9898718524076344)


def test_xsiz2():
    assert wf.xsiz2(1, 3, 36) == pytest.approx(0.9898718524077995)


def test_sgrund():
    assert wf.sgrund(0.03, 5, 2) == pytest.approx(0.612)


def test_snue():
    assert wf.snue(5, 2, 0.03) == pytest.approx(-1.4249999999999998)


def test_rmag_table_hit():
    # exact B=10000 Gauss table point -> exact H value from the table
    assert wf.rmag(10000) == pytest.approx(2.75)


def test_rmag_above_table_extrapolates():
    val = wf.rmag(40000)
    assert val > wf.rmag(31000)


def test_ggt():
    assert wf.ggt(12, 18) == 6.0
    assert wf.ggt(7, 13) == 1.0


def test_mgegen_constant_torque():
    st = make_state()
    st.daten.mg[0] = 1
    st.daten.mg[1] = 100
    st.daten.mg[3] = 0.2
    st.daten.mg[4] = 50
    st.daten.mgn = 1
    assert wf.mgegen(0.1, st) == pytest.approx(62.5)
    assert wf.mgegen(0.5, st) == pytest.approx(50.0)


def test_mgegen_disabled_returns_zero():
    st = make_state()
    assert wf.mgegen(0.5, st) == 0.0


# ---------------------------------------------------------------- geometry_functions

def test_bruvor_integer_slot_branch():
    st = make_state()
    mz1s, weite1, a1max = gf.bruvor(36, 2, 3, 3, 8, st)
    assert mz1s == pytest.approx(3)
    assert weite1 == pytest.approx(8)
    assert a1max == pytest.approx(1.0)
    assert st.brulo.ibruch == 0


def test_coil_known_frame():
    bcoil, da1max = gf.coil("A5XX500XXXXXXXXXXXXXXXXXXXX", 0, 0)
    assert (bcoil, da1max) == (800.0, 840.0)


def test_coil_unknown_frame_unchanged():
    bcoil, da1max = gf.coil("ZZXX999XXXXXXXXXXXXXXXXXXXX", 5, 6)
    assert (bcoil, da1max) == (5, 6)


def test_nut():
    anut, kcuth = gf.nut(bs=8, bn=8, bns=9, hk=1, ho=1, hcuo=5, hzw=1,
                          hcuu=10, hu=2, qcu=20, zn=40)
    assert anut == pytest.approx(161.5)
    assert kcuth == pytest.approx(5.223367697594502)


def test_vj():
    assert gf.vj(bj=1.4, lj=0.05, p=2) == pytest.approx(23.42, rel=1e-3)


def test_vlzfe():
    assert gf.vlzfe(v10=1.6, fn=60) == pytest.approx(1.9776)


def test_nue0():
    assert gf.nue0(2, 1.0, 36) == pytest.approx(34.0)


def test_paket_populates_cpaket():
    st = make_state()
    gf.paket("M-typ", l_=0.3, bk=0.01, nk1=5, nk2=5, st=st)
    assert st.cpaket.paket1[1] == pytest.approx(40.0)
    assert st.cpaket.paket1[3] == pytest.approx(45.0)


# ---------------------------------------------------------------- magnetic_circuit

def test_std_known_values():
    st = make_state()
    out = mc.std(da1=0.4, di1=0.25, n1=36, p=2, weite1=8, mzone1=3, zn1=40,
                 lm1=0.5, qcu1=1e-5, m1=3, fn=50, li=0.2, theta1=75, mz1s=3,
                 hs1=0.001, hk1=0.001, ho1=0.001, hcuo1=0.005, hzw1=0.001,
                 hcuu1=0.01, hu1=0.002, bs1=0.008, bn1=0.009, bn1s=0.01,
                 lamst1_in=0.0, st=st)
    assert out["r1k"] == pytest.approx(0.41379310344827586)
    assert out["lamst1"] == pytest.approx(0.15)
    assert out["xsnut1"] == pytest.approx(0.9608354115471317)


def test_steg_runs_and_returns_arrays():
    st = make_state()
    st.cbsort.h = np.linspace(0, 100, 46)
    st.cbsort.my0 = 1.2566e-6
    inut, lamda = mc.steg(ilfr=2, li=0.2, lfe=0.19, kfe=0.95, hs2=0.02,
                           bs2=0.004, br12=0.01, br22=0.01, br32=0.01,
                           br42=0.0, br52=0.0, bstr52=0.0, st=st)
    assert len(inut) == 21
    assert len(lamda) == 21
    assert lamda[1] == pytest.approx(5.0)


# ---------------------------------------------------------------- rotor_functions

def test_slfr():
    st = make_state()
    out = rf.slfr(di2=0.24, hojoc2=0.03, n2=24, p=2, weite2=6, zn2=20,
                   lm2=0.4, qcu2=1e-5, fn=50, li=0.2, theta2=75, hs2=0.001,
                   hk2=0.0005, ho2=0.001, hcuo2=0.004, hzw2=0.001,
                   hcuu2=0.008, hu2=0.001, bs2=0.006, bn2=0.007,
                   bn2s=0.008, st=st)
    assert out["r2k"] == pytest.approx(0.1103448275862069)
    assert st.xs2sl.xsn2 == pytest.approx(0.11420024165515455)


def test_kikr1_izgf_not_1_returns_all_ones():
    svp = rf.kikr1(izgf=0, mtl=2, ntl=10, azweig=1, btl=0.006, htl=0.002,
                    fn=50, li=0.2, lm1=0.5, theta1=75, bn1=0.009, zn1=40)
    assert np.all(svp[1:26, 2] == 1.0)


def test_kikr_cage_rotor():
    st = make_state()
    hozah2, brzah2, kappaf, lamnu2, r2s, r2kw, qstab = rf.kikr(
        ilfr=2, h22=0, br12=0, br22=0, kappa2=0,
        h32=0.02, br32=0.006, kappa3=30e6,
        h42=0.002, br42=0.005, kappa4=30e6,
        h52=0, bstr52=0, kappa5=0, h62=0, br52=0, br62=0, kappa6=0,
        h72=0, bstr62=0, kappa7=0,
        qring=0.0002, kapri=30e6, qringa=0, kapria=0,
        fn=50, di2=0.24, hojoc2=0.03, n2=24, beta=0,
        n1=36, p=2, li=0.2, lstab=0.25, st=st)
    assert kappaf == pytest.approx(29999999.999999985)
    assert lamnu2 == pytest.approx(1.4378824966385662)
    assert r2s == pytest.approx(5.2821198574701964e-05)
    # skin effect factors should start at 1.0 for S=0 and increase with slip
    assert st.kennli.svpool[1, 2] == pytest.approx(1.0)
    assert st.kennli.svpool[10, 2] > st.kennli.svpool[2, 2]


def test_kikrsx_converges():
    xsistr = rf.kikrsx(xsi=1.5, krtstr=1.2)
    assert xsistr >= 1.5


# ---------------------------------------------------------------- magnetics_solver

def test_sorte_default_curve():
    st = make_state()
    blech = ms.sorte(170, st)
    assert blech == " 1.7 W/KG "
    assert st.cbsort.h[1] == pytest.approx(24.0)


def test_sorte_unmatched_guete_falls_back_to_17():
    st = make_state()
    blech = ms.sorte(999, st)
    assert blech == " 1.7 W/KG "


def test_platt():
    assert ms.platt(1.0) == pytest.approx(1.3253262153484369)


def test_kzwert_known_lookups():
    assert ms.kzwert("ARX355" + "X" * 64, 2) == pytest.approx(28.0)
    assert ms.kzwert("AJX630" + "X" * 64, 4) == pytest.approx(18.0)
    assert ms.kzwert("AMX500" + "X" * 64, 2) == 0.0


def test_kreis_generates_magnetization_curve():
    st = make_state()
    brzah1 = np.zeros(6)
    brzah1[1:6] = [0.012, 0.0125, 0.013, 0.0135, 0.014]
    brzah2 = np.zeros(6)
    brzah2[1:6] = [0.010, 0.0105, 0.011, 0.0115, 0.012]

    out = ms.kreis(n1=36, n2=24, p=2, kfe=0.95, rb=0.13, li=0.2, lfe=0.19,
                    lfe2=0.19, da1=0.4, di2=0.25, hojoc1=0.03, kj1=1.0,
                    hojoc2=0.03, kj2=1.0, deltag=0.0006, bs1=0.006,
                    bs2=0.005, hs1=0.001, hozah1=0.02, brzah1=brzah1,
                    hozah2=0.02, brzah2=brzah2, cfej=0.0, cfez=0.0,
                    v10=1.65, u1=380.0, fn=50, w1=200, xsip=0.92, m1=3,
                    st=st)
    assert out["eanz"] == 14
    assert st.kennli.eanz == 14
    assert out["cfei"] == pytest.approx(4.127858315369495)
    # UI should be monotonically increasing across the generated curve
    ui_values = st.kennli.mkpool[1:15, 1]
    assert all(b > a for a, b in zip(ui_values, ui_values[1:]))


# ---------------------------------------------------------------- harmonics

def test_popaza_selects_harmonics():
    st = make_state()
    fw2 = hm.popaza(ilfr=2, p=2, fn=50, zn1=40, zn2=1, n1=36, n2=24, m1=3,
                     mzone1=3, mz1s=3, weite1=8, weite2=0, dring=0.3,
                     kapri=30e6, qring=0.0003, dringa=0, kapria=0,
                     qringa=0, schr=0, st=st)
    assert fw2 == pytest.approx(1.0)
    assert st.nyinfo.nyanz == 7
    assert list(st.nyinfo.ny[1:8]) == pytest.approx(
        [2., -4., 8., -10., -16., -34., 38.])


def test_z2nys_fundamental_updates_esb():
    st = make_state()
    rf.kikr(ilfr=2, h22=0, br12=0, br22=0, kappa2=0,
            h32=0.02, br32=0.006, kappa3=30e6,
            h42=0.002, br42=0.005, kappa4=30e6,
            h52=0, bstr52=0, kappa5=0, h62=0, br52=0, br62=0, kappa6=0,
            h72=0, bstr62=0, kappa7=0,
            qring=0.0003, kapri=30e6, qringa=0, kapria=0,
            fn=50, di2=0.24, hojoc2=0.03, n2=24, beta=0,
            n1=36, p=2, li=0.2, lstab=0.25, st=st)

    cz, cd, i2o, i2u = hm.z2nys(
        ilfr=2, sny=1.0, ny=2, fn=50, x1hny=5.0, x1hny0=5.0, sf=1.0,
        uer=600.0, xssti2=0.05, xsnut2=0.1, sd2=0.02, r2w=0.02,
        xssteg=0.01, li=0.2, lamnu2=1.4, xsrny=0.02, rrny=0.05,
        r2s=0.0001, r2kw=0.00005, h42=0.002, br42=0.005, r2kwa=0,
        xsrnya=0, rrnya=0, is1pkt=1, st=st)
    assert i2o == 0.0
    assert i2u == 1.0
    assert st.esb.x1h == pytest.approx(5.0)
    assert st.esb.xn2 > 0.0


# ---------------------------------------------------------------- operating_point

def _build_full_state():
    """Shared fixture: builds a MachineState with a full magnetization
    curve, cage-rotor skin-effect table, harmonic info, and slot-bridge
    saturation table populated - everything SPKT needs to run."""
    st = make_state()
    st.keil.xskeil = 0.0
    st.keil.teta0 = 2000.0 * 0.006 / 0.015 / 40

    brzah1 = np.zeros(6)
    brzah1[1:6] = [0.012, 0.0125, 0.013, 0.0135, 0.014]
    brzah2 = np.zeros(6)
    brzah2[1:6] = [0.010, 0.0105, 0.011, 0.0115, 0.012]
    ms.kreis(n1=36, n2=24, p=2, kfe=0.95, rb=0.13, li=0.2, lfe=0.19,
              lfe2=0.19, da1=0.4, di2=0.25, hojoc1=0.03, kj1=1.0,
              hojoc2=0.03, kj2=1.0, deltag=0.0006, bs1=0.006, bs2=0.005,
              hs1=0.001, hozah1=0.02, brzah1=brzah1, hozah2=0.02,
              brzah2=brzah2, cfej=0.0, cfez=0.0, v10=1.65, u1=380.0,
              fn=50, w1=200, xsip=0.92, m1=3, st=st)

    rf.kikr(ilfr=2, h22=0, br12=0, br22=0, kappa2=0,
            h32=0.02, br32=0.006, kappa3=30e6,
            h42=0.002, br42=0.005, kappa4=30e6,
            h52=0, bstr52=0, kappa5=0, h62=0, br52=0, br62=0, kappa6=0,
            h72=0, bstr62=0, kappa7=0,
            qring=0.0003, kapri=30e6, qringa=0, kapria=0,
            fn=50, di2=0.24, hojoc2=0.03, n2=24, beta=0,
            n1=36, p=2, li=0.2, lstab=0.25, st=st)
    hm.popaza(ilfr=2, p=2, fn=50, zn1=40, zn2=1, n1=36, n2=24, m1=3,
              mzone1=3, mz1s=3, weite1=8, weite2=0, dring=0.3,
              kapri=30e6, qring=0.0003, dringa=0, kapria=0, qringa=0,
              schr=0, st=st)

    st.cbsort.h = np.linspace(0, 100000, 46)
    st.cbsort.my0 = 1.2566e-6
    inut, lamda = mc.steg(ilfr=2, li=0.2, lfe=0.19, kfe=0.95, hs2=0.02,
                           bs2=0.004, br12=0.01, br22=0.01, br32=0.01,
                           br42=0.002, br52=0.0, bstr52=0.0, st=st)
    st.kennli.inut[1:21] = inut[1:21]
    st.kennli.lamda[1:21] = lamda[1:21]
    return st


def test_mkwert_interpolates_within_curve():
    st = _build_full_state()
    mk = op.mkwert(0, st.kennli.mkpool[5, 1], st)
    assert "imy" in mk and "ife" in mk
    assert mk["imy"] == pytest.approx(1.7815080740580804)


def test_mkwert_out_of_range_raises():
    st = _build_full_state()
    with pytest.raises(ValueError):
        op.mkwert(0, 1.0e9, st)


def test_sort_descending_and_zero_nudge():
    st = make_state()
    st.svorga.sanz = 5
    st.svorga.svor[1:6] = [0.03, 1.0, 0.1, 0.0, 0.5]
    ssanz, ss = op.sort(st)
    assert ssanz == 5
    assert list(ss[1:6]) == pytest.approx([1.0, 0.5, 0.1, 0.03, 1e-9])


def test_spkt_converges_and_updates_state():
    st = _build_full_state()
    out = op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05,
                   xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2,
                   hs2=0.02, lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5,
                   h42=0.002, br42=0.005, r2kwa=0, pn=200000, u1=380.0,
                   m1=3, prbg0=500, rzus=0, is1pkt=1, st=st)
    assert out["i1"] > 0
    assert out["n"] == pytest.approx(24.25)
    assert 0.0 < out["cosph"] < 1.0
    # IS1PKT=1 branch should have written the fundamental ESB values
    assert st.esb.r1 == pytest.approx(0.55)
    assert st.esb.xn1 == pytest.approx(0.96)


# ---------------------------------------------------------------- load_curves

def test_slast_and_holauf_full_pipeline():
    """
    End-to-end integration test spanning nearly the whole ported program:
    KREIS -> KIKR -> POPAZA -> STEG -> SPKT (rated point) -> SLAST ->
    HOLAUF. Not checked against a real machine's known-correct values
    (still no sample EINASYN available) - this confirms the pipeline
    runs to completion and produces internally-consistent results, which
    is what a synthetic test can promise.
    """
    st = _build_full_state()

    # SPKT must be called once with is1pkt=1 (the rated point) before any
    # is1pkt=0 sweep - see the note in operating_point.spkt's docstring.
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)

    slips = list(np.arange(0.05, 1.05, 0.05))
    st.svorga.sanz = len(slips)
    for i, s in enumerate(slips, start=1):
        st.svorga.svor[i] = s

    st.daten.mg[0] = 1
    st.daten.mg[1] = 0.3
    st.daten.mg[3] = 0.2
    st.daten.mg[4] = 0.5
    st.daten.mgn = 1.0

    msat, nsat = lc.slast(
        reke="TEST-MACHINE", ilfr=2, ku1=0, xnetz=0, p=2, fn=50, r1w=0.55,
        xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146,
        li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44, r2s=5.28e-5,
        r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0, pn=200000, u1=380.0,
        m1=3, prbg0=500, i1n=350, mn=800, qcu1=1e-5, qcu2=1e-4, qstab=1e-4,
        qstaba=0, qring=1e-4, qringa=0, n2=24, kma=1.0, kia=1.0, st=st,
        verbose=False)

    assert st.mnline.nanz == 20
    assert msat > 0

    holauf_out = lc.holauf(ilfr=2, fn=50, p=2, jsum=5.0, mn=800, k41=1e6,
                            k42=1e6, k4a2=1e6, k4r2=1e6, k4ar2=1e6,
                            r1w=0.55, m1=3, i1n=350, ski=0.03, ku1=0,
                            xnetz=0, idruck=1, st=st, verbose=False)

    assert holauf_out["th"] > 0
    assert holauf_out["erw1"] >= 0


# ---------------------------------------------------------------- special_points

def test_sk_finds_breakdown_slip():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    skipp = sp.sk(bs2i0=0.004, bs2in=0.004, ilfr=2, p=2, fn=50, r1w=0.55,
                   xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                   r2w=0.146, li=0.2, hs2=0.02, lamnu2=1.44, r2s=5.28e-5,
                   r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0,
                   pn=200000, u1=380.0, m1=3, prbg0=500, st=st)
    assert skipp == pytest.approx(0.09762319717364473)
    assert st.esb.skipp == pytest.approx(skipp)


def test_skgen_finds_negative_slip():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    skippg = sp.skgen(skipp=0.0976, ilfr=2, p=2, fn=50, r1w=0.55,
                       xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                       r2w=0.146, li=0.2, hs2=0.02, lamnu2=1.44,
                       r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                       r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, st=st)
    assert skippg < 0


def test_skny_populates_svorga():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    x1hp0 = st.dp.x1hp0
    sanz = sp.skny(ilfr=2, skipp=0.0976, x1hp=x1hp0, hs2=0.02, bs2i=0.004,
                    xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, lamnu2=1.44,
                    r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                    fn=50, p=2, x1hp0=x1hp0, st=st)
    assert sanz == st.svorga.sanz
    assert sanz > 60  # at least the fixed 0.025-step grid (61 points)


def test_praxis_converges():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    out = sp.praxis(ilfr=2, s=1.0, p=2, fn=50, r1w=0.55, xssti1=2.05,
                     xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146,
                     li=0.2, hs2=0.02, lamnu2=1.44, r2s=5.28e-5,
                     r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0,
                     pn=200000, u1=380.0, m1=3, x1hp0=st.dp.x1hp0,
                     mnyth=1.0, i1th=300.0, pc2nth=10.0, st=st)
    assert out["i1"] > 0
    assert out["n"] == 0.0  # S=1 -> synchronous-slip speed is 0


def test_dfmag_runs_five_points():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    results = sp.dfmag(reke="TEST-MACHINE", ilfr=2, p=2, fn=50, r1w=0.55,
                        xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                        r2w=0.146, li=0.2, hs2=0.02, lamnu2=1.44,
                        r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                        r2kwa=0, pn=200000, u1=380.0, m1=3, x1h=5.0,
                        st=st, verbose=False)
    assert len(results) == 5
    assert [r["s"] for r in results] == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_pruef_clean_data_no_error():
    st = make_state()
    d = st.daten
    d.ilfr = 2
    d.m1 = 3
    d.mzone1 = 3
    d.lm1 = 0.5
    d.l = 0.2
    d.n2 = 24
    d.p = 2
    d.lstab = 0.25
    d.zk = 0.1
    d.bk = 0.01
    d.kfe = 0.95
    d.di1 = 0.25
    d.da1 = 0.4
    d.di2 = 0.24
    d.mtl = 2
    d.izgf = 0
    iabort, messages = sp.pruef(st)
    assert iabort == 0
    assert not any(is_err for is_err, _ in messages)


def test_pruef_catches_di2_greater_than_di1():
    st = make_state()
    d = st.daten
    d.ilfr = 2
    d.m1 = 3
    d.mzone1 = 3
    d.lm1 = 0.5
    d.l = 0.2
    d.n2 = 24
    d.p = 2
    d.lstab = 0.25
    d.kfe = 0.95
    d.di1 = 0.25
    d.da1 = 0.4
    d.di2 = 0.5  # invalid: DI2 > DI1
    d.mtl = 2
    d.izgf = 0
    iabort, messages = sp.pruef(st)
    assert iabort == 1
    assert any("DI2 GROESSER DI1" in text for _, text in messages)


# ---------------------------------------------------------------- insulation

def test_ald_in_range():
    out = ins.ald(z=0.008, r=0.05, rd=0.02, l1=0.05, hstreu=0.001,
                   hzw=0.001, n1=36, weite1=8, di1=0.25, l=0.2,
                   bkopf=0.01, hkopf=0.02)
    assert out["in_range"] is True
    assert out["lges1"] == pytest.approx(0.6456220212180396)


def test_ald_out_of_range_flagged():
    # a much smaller Z pushes ALPHA (and hence HG1) well outside range
    out = ins.ald(z=0.0001, r=0.05, rd=0.02, l1=0.05, hstreu=0.001,
                   hzw=0.001, n1=36, weite1=8, di1=0.25, l=0.2,
                   bkopf=0.01, hkopf=0.02)
    assert out["in_range"] is False


def test_nariss_random_wound_low_voltage():
    st = make_state()
    bnmind, hnmind, hzw, iokay, dik = ins.nariss(
        iss=3, i2l1gs=0, unx=400.0, btl=1.5, htl=1.5, mtl=4, ntl=2, tpar=1,
        hs=0.5, hk=4.0, bn=15.0, hn=40.0, st=st)
    assert iokay == 1
    assert bnmind == pytest.approx(4.999999999999999)
    assert hnmind == pytest.approx(18.58)
    assert st.spkopf.bkopf == pytest.approx(4.5)


def test_fspule_symmetric_coil():
    out = ins.fspule(iart=0, l=200.0, di1=250.0, n1=36, hn1min=40.0,
                      bn1=15.0, hstreu=1.0, hzw1=3.5, diso1=3.0, weite=8,
                      a=0.0, b=6.0, sh=20.0, dstr=0.0, l1=65.0, l2=65.0,
                      lstr=5.0, r=35.0, rd=15.0)
    assert out["beta1"] == pytest.approx(40.0)
    assert out["auslad"] == pytest.approx(511.0278663434302)


def test_zgf_symmetric_and_asymmetric():
    lw, ald1s, ald1a, sym, asym = ins.zgf(
        reke="TEST", l=200.0, di1=250.0, en1=36, hn1min=40.0, bn1=15.0,
        hstreu=1.0, hzw1=3.5, diso1=3.0, eweite=8, a=0.0, b=6.0, sh=20.0,
        dstr=0.0, l1=65.0, l2=65.0, lstr=5.0, r=35.0, rd=15.0)
    assert ald1s == pytest.approx(511.0278663434302)
    assert ald1a < ald1s  # asymmetric coil has a shorter overhang here


def test_wvar_finds_conductor_configuration():
    st = make_state()
    hzw, mtl, ntl, tpar, qcu, btl, htl, dik = ins.wvar(
        iss=3, i2l1gs=0, unx=400.0, zn=40, azweig=1, hs=0.5, hk=4.0,
        bn=15.0, hn=40.0, st=st)
    assert qcu > 0
    assert btl > 0 and htl > 0


def test_isostd_conductor_given_designs_slot():
    st = make_state()
    reke = "A5XX500" + "X" * 63
    out = ins.isostd(
        reke=reke, l=200.0, di1=250.0, n1=36, weite1=8, isoart=0,
        unx=400.0, qcu1=0.0, bs1=0.0, bn1=0.0, bn1s=0.0, hs1=0.0,
        hk1=0.0, ho1=0.0, hcuo1=0.0, hzw1=3.5, hcuu1=0.0, hu1=0.0,
        zn1=40, azweig=1, lw1=0.0, mtl=4, ntl=2, utl=0.0, btl=1.5,
        htl=1.5, diso1=0.0, l1=0.0, l2=0.0, lstr=5.0, r=0.0, rd=0.0,
        b=0.0, dstr=0.0, m1=3, p=2, da1=400.0, zk=0.1, bk=0.01,
        schr=0.0, st=st, verbose=False)
    assert out["bn1"] == pytest.approx(7.0)
    assert out["drufi"] == pytest.approx(25.0)  # A5/500 frame override
    assert out["iss"] == 3  # random-wound, UNX <= 1500


def test_isostd_slot_given_designs_conductor():
    st = make_state()
    out = ins.isostd(
        reke="M-TYPE-FRAME", l=200.0, di1=250.0, n1=36, weite1=8,
        isoart=0, unx=400.0, qcu1=0.0, bs1=15.5, bn1=15.0, bn1s=15.0,
        hs1=0.5, hk1=4.0, ho1=0.0, hcuo1=15.0, hzw1=3.5, hcuu1=15.0,
        hu1=0.0, zn1=40, azweig=1, lw1=0.0, mtl=0.0, ntl=0.0, utl=0.0,
        btl=0.0, htl=0.0, diso1=0.0, l1=0.0, l2=0.0, lstr=5.0, r=0.0,
        rd=0.0, b=0.0, dstr=0.0, m1=3, p=2, da1=400.0, zk=0.1, bk=0.01,
        schr=0.0, st=st, verbose=False)
    assert out["qcu1"] > 0
    assert out["mtl"] > 0 and out["ntl"] > 0


def test_isostd_rejects_zero_azweig():
    st = make_state()
    with pytest.raises(ValueError):
        ins.isostd(
            reke="X", l=1, di1=1, n1=36, weite1=8, isoart=0, unx=400,
            qcu1=0, bs1=0, bn1=0, bn1s=0, hs1=0, hk1=0, ho1=0, hcuo1=0,
            hzw1=0, hcuu1=0, hu1=0, zn1=40, azweig=0.0, lw1=0, mtl=1,
            ntl=1, utl=0, btl=1, htl=1, diso1=0, l1=0, l2=0, lstr=0,
            r=0, rd=0, b=0, dstr=0, m1=3, p=2, da1=1, zk=0, bk=0,
            schr=0, st=st)


def test_cufl():
    acu = ins.cufl(bs=8, bn=8, bns=9, hk=1, ho=1, hcuo=5, hzw=1, hcuu=10, hu=2)
    assert acu == pytest.approx(153.1578947368421)


def test_zplan_no_overlap():
    a = [0, 1, 0, 0]  # q=3, one turn at position 1
    assert ins.zplan(3, a, svk=1) == 1


def test_teilen_uneven_split_fallback():
    a, iokay = ins.teilen(q1=3.0, wdg=3.5, gwdg=3.0, bwdg=0.5, weite1=8)
    assert a == [3, 4, 4]
    assert iokay == 0


def test_teilen_finds_symmetric_pattern():
    a, iokay = ins.teilen(q1=4.0, wdg=4.25, gwdg=4.0, bwdg=0.25, weite1=10)
    assert iokay == 1
    assert sum(a) == pytest.approx(4 * 4 + 1)  # 3 slots @4 turns + 1 @5


def test_isotr_returns_suggestions_without_tpar():
    out = ins.isotr(rkz=" " * 12, bs1=8, bn1=8, bn1s=9, hk1=1, ho1=1,
                     hcuo1=5, hzw1=1, hcuu1=10, hu1=2, zn1=40, azweig=1,
                     n1=36, p=2, mzone1=3, btl=1.5, htl=1.5, weite1=8,
                     qcu1=0.0, tpar=None, verbose=False)
    assert out["mtl"] is None
    assert len(out["suggestions"]) > 0


def test_isotr_finalizes_with_chosen_tpar():
    out = ins.isotr(rkz=" " * 12, bs1=8, bn1=8, bn1s=9, hk1=1, ho1=1,
                     hcuo1=5, hzw1=1, hcuu1=10, hu1=2, zn1=40, azweig=1,
                     n1=36, p=2, mzone1=3, btl=1.5, htl=1.5, weite1=8,
                     qcu1=0.0, tpar=1.0, verbose=False)
    assert out["mtl"] == pytest.approx(1.0)
    assert out["qcu1"] == pytest.approx(1.767144375)


def test_isotr_rejects_zero_azweig():
    with pytest.raises(ValueError):
        ins.isotr(rkz=" " * 12, bs1=8, bn1=8, bn1s=9, hk1=1, ho1=1,
                  hcuo1=5, hzw1=1, hcuu1=10, hu1=2, zn1=40, azweig=0.0,
                  n1=36, p=2, mzone1=3, btl=1.5, htl=1.5, weite1=8)


def test_nar2_search_and_print_modes():
    st = make_state()
    bnmind, hnmind, hzw, iokay, unx2 = ins.nar2(
        iprint=0, imerk=0, u20=400.0, btl=2.0, htl=8.0, mtl=4, ntl=2,
        tpar=1, hs=0.5, hk=4.0, bs=15.5, bn=15.0, hn=40.0, hzw=0.0,
        qcue=0.0, st=st)
    assert unx2 == pytest.approx(999.99)

    bnmind2, hnmind2, hzw2, iokay2, unx2b = ins.nar2(
        iprint=1, imerk=0, u20=400.0, btl=2.0, htl=8.0, mtl=4, ntl=2,
        tpar=1, hs=0.5, hk=4.0, bs=15.5, bn=15.0, hn=40.0, hzw=hzw,
        qcue=1.5, st=st)
    assert st.ck211l.delb == pytest.approx(3.5)
    assert st.ck211l.del1l == pytest.approx(0.25)


def test_isolfr_finds_conductor():
    st = make_state()
    out = ins.isolfr(u20=400.0, zn=4, n=1500, lm=500.0, azweig=1,
                      weite=6, bs=15.5, bn=15.0, hs=0.5, hk=4.0, ho=0.0,
                      hcuo=15.0, hzw=3.5, hcuu=15.0, hu=0.0, st=st)
    assert out["iokay"] == 1
    assert out["mtl"] == pytest.approx(4.0)
    assert out["wdg"] == pytest.approx(2.0)


def test_isolfr_reversible_duty_doubles_voltage():
    st = make_state()
    out = ins.isolfr(u20=400.0, zn=4, n=1500, lm=500.0, azweig=1,
                      weite=6, bs=15.5, bn=15.0, hs=0.5, hk=4.0, ho=0.0,
                      hcuo=15.0, hzw=3.5, hcuu=15.0, hu=0.0, st=st,
                      reversible_duty=True)
    assert out["u20"] == pytest.approx(800.0)


def test_isolfr_raises_when_no_conductor_fits():
    st = make_state()
    with pytest.raises(ValueError):
        ins.isolfr(u20=400.0, zn=40, n=1500, lm=500.0, azweig=1,
                   weite=6, bs=15.5, bn=15.0, hs=0.5, hk=4.0, ho=0.0,
                   hcuo=15.0, hzw=3.5, hcuu=15.0, hu=0.0, st=st)


# ---------------------------------------------------------------- SPKTSP / GENERA

def test_spktsp_converges():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    out = sp.spktsp(ilfr=2, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
                     xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
                     lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002,
                     br42=0.005, r2kwa=0, pn=200000, u1=380.0, m1=3,
                     kka=1.0, st=st)
    assert out["i1"] > 0
    assert 0.0 < abs(out["cosph"]) <= 1.0


def test_genera_full_load_sweep_converges():
    """
    Uses smaller, self-consistent PN/PMSV targets than the rest of the
    test suite's synthetic machine, since GENERA's slip search needs the
    target power to actually be reachable within +/-0.4% in 300
    iterations - see the module docstring's note on GENERA's RuntimeError
    for the (deliberately loud) alternative when it isn't.
    """
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=2000, u1=380.0, m1=3, prbg0=5, rzus=0,
             is1pkt=1, st=st)
    skipp = sp.sk(bs2i0=0.004, bs2in=0.004, ilfr=2, p=2, fn=50, r1w=0.55,
                   xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                   r2w=0.146, li=0.2, hs2=0.02, lamnu2=1.44, r2s=5.28e-5,
                   r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0,
                   pn=2000, u1=380.0, m1=3, prbg0=5, st=st)

    pn = 2000.0
    pmsv = [0.0, 0.25 * pn, 0.5 * pn, 0.75 * pn, pn, 1.25 * pn]
    out = sp.genera(reke="TEST-GEN", ilfr=2, p=2, fn=50, r1w=0.55,
                     xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                     r2w=0.146, li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44,
                     r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                     r2kwa=0, u1=380.0, m1=3, prbg0=5, pn=pn, pmsv=pmsv,
                     skipp=skipp, pfe0=15.0, kz=1.0, kt=1.0, qcu1=1e-5,
                     zn1=40, n1=36, di1=0.25, ishalt=0, ksys=1.0, st=st,
                     verbose=False)
    assert all(p["converged"] for p in out["points"])
    assert out["dtkz"] is not None and out["dtkt"] is not None
    assert out["breakdown"]["m_over_mngen"] > 0


# ---------------------------------------------------------------- lesen

def _make_synthetic_einasyn(path):
    """
    Builds a synthetic EINASYN-format file using lesen.py's own
    _CardWriter (which shares its skip-N-lines logic with _CardReader,
    so this can't drift out of sync with the reader the way a hand-
    rolled version did earlier in this port - see CONVERSION_STATUS.md
    round 21 for that story). This is a PARSER round-trip test, not a
    real machine data file - no sample EINASYN has been available
    during this port, so this only proves the fixed-column reads land
    on the columns lesen.py's docstring claims they do, and that the
    UMWANDLUNG IN SI-EINHEITEN (mm/kW -> m/W) conversion block produces
    the expected SI values, not that these specific columns are
    correct for a real file.

    All geometry fields here are in MILLIMETRES and PN is in KILOWATTS,
    matching EINASYN's raw convention (LESEN converts both to SI/metres
    and watts internally - see the UNITS note in lesen.py's docstring).
    """
    from asyn4.lesen import _CardWriter
    w = _CardWriter()
    w.write_mixed(2, [("a", 16, "TEST-MACHINE"), ("a", 8, "01072026"),
                       ("a", 25, "KENNWORT-TEXT")])
    w.write_mixed(2, [("a", 27, "A5-315M04"), ("a", 17, "VANR-1"),
                       ("a", 15, "FBNR-1"), ("i", 5, 12345), ("a", 12, "MONR-1")])
    w.write_mixed(2, [("a", 11, "2 X     E")])
    w.write_f(3, [50, 0.95, 1.0, 1.0, 2, 400, 200, 0, 0, 0.05, 0], decimals=2)
    w.write_mixed(2, [("i", 7, 0)] + [("f", 7, v) for v in
                  [380, 50, 1.65, 0, 0, 0, 0.6, 0, 0, 0]])
    w.write_mixed(2, [("f", 7, v) for v in [3, 3, 250, 36, 40, 8, 0, 0, 75]] +
                  [("i", 7, 0), ("f", 7, 400)])
    w.write_mixed(2, [("f", 7, 6.0), ("f", 7, 1.0), ("i", 7, 0), ("a", 7, "N1"),
                       ("f", 7, 0), ("f", 7, 0), ("f", 7, 0), ("f", 7, 0)])
    w.write_f(2, [9.0, 10.0, 1.0, 1.0, 5.0, 1.0, 10.0, 2.0, 0, 0, 0], decimals=1)
    w.write_mixed(2, [("f", 7, v) for v in
                  [150, 24, 5.0, 20.0, 4, 2, 0, 1, 1.2, 1.2]] + [("a", 5, "N2")])
    w.write_f(3, [1, 0, 0, 0, 75, 0, 0, 0, 0, 0])
    w.write_mixed(2, [("f", 7, v) for v in
                  [7.0, 8.0, 1.0, 1.0, 5.0, 1.0, 5.0, 0]] + [("i", 7, 1)])
    w.write_f(3, [20.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_f(2, [6.0, 5.0, 0, 0, 0, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_f(2, [30, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_mixed(2, [("f", 7, 0, 4), ("f", 7, 300, 1), ("f", 7, 30, 1)] +
                  [("f", 7, 0, 4) for _ in range(8)])
    w.write_f(3, [0, 0, 0, 0, 0])
    w.write_mixed(3, [("f", 14, 0, 4), ("f", 14, 1, 4)] +
                  [("f", 7, v) for v in [1, 100, 0, 0.2, 50]])
    w.write_f(3, [0] * 11)
    w.write_mixed(3, [("a", 10, "WP1TEXT"), ("a", 4, ""),
                       ("a", 10, "WP2TEXT"), ("a", 4, "")])
    w.write(path)


def test_lesen_parses_fixed_columns_correctly(tmp_path):
    path = tmp_path / "test_einasyn.txt"
    _make_synthetic_einasyn(str(path))

    st = MachineState()
    out = ls.lesen(str(path), irenr=1, st=st, verbose=False)

    assert out["name"].strip() == "TEST-MACHINE"
    assert out["reke"].strip() == "A5-315M04"
    assert out["wanr"] == 12345
    assert out["monr"].strip() == "MONR-1"
    assert out["rkz"][8] == "E"

    d = st.daten
    assert d.ilfr == 2
    assert d.fn == 50.0
    assert d.p == 2.0
    assert d.da1 == pytest.approx(0.400)
    assert d.l == pytest.approx(0.200)
    assert d.pn == pytest.approx(50.0)  # 0.05 kW raw -> 50 W after *K
    assert d.m1 == 3.0
    assert d.di1 == pytest.approx(0.250)
    assert d.n1 == 36.0
    assert d.bs1 == pytest.approx(0.006)
    assert d.bn1 == pytest.approx(0.009)
    assert d.bn1s == pytest.approx(0.010)
    assert d.di2 == pytest.approx(0.150)
    assert d.n2 == 24.0
    # This test geometry (H32=20mm != BR32=6mm) triggers the hammer-
    # head-slot transform (source 2780-2799), which moves KAPPA3's
    # value into KAPPA2 and zeros KAPPA3 - see lesen.py's docstring.
    assert d.kappa2 == pytest.approx(30.0e6)  # 30 raw -> 30e6 S/m after *M
    assert d.kappa3 == 0.0


def test_lesen_routes_izgf_to_isotr_not_isostd(tmp_path):
    path = tmp_path / "test_einasyn.txt"
    _make_synthetic_einasyn(str(path))
    st = MachineState()
    out = ls.lesen(str(path), irenr=1, st=st, verbose=False)
    # IZGF=0 in the synthetic file -> random-wound path (ISOTR), not
    # form-coil (ISOSTD); ILFR=2 (not slip-ring) -> ISOLFR skipped.
    assert out["isostd_result"] is None
    assert out["isotr_result"] is not None
    assert out["isolfr_result"] is None


# ---------------------------------------------------------------- reports

def test_mwinfo_sheet_grade_lookup_and_lamination():
    out = rp.mwinfo(reke="A5XX500XX", v10=1.70, da1=1300.0, di1=1290.0,
                     deltag=1.0)
    assert out["sheet_grade"] == "V400-50A"
    assert out["stator_lam"] == "segmentiert"
    assert out["is_a5_frame"] is True


def test_mwinfo_unmatched_grade_and_solid_lamination():
    out = rp.mwinfo(reke="M-TYPE", v10=1.65, da1=400.0, di1=250.0,
                     deltag=0.6)
    assert out["sheet_grade"] is None
    assert out["stator_lam"] == "Vollblech"
    assert out["is_a5_frame"] is False


def test_synmom_runs_and_returns_breakdown():
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=200000, u1=380.0, m1=3, prbg0=500, rzus=0,
             is1pkt=1, st=st)
    out = rp.synmom(reke="TEST-MACHINE", ilfr=2, p=2, fn=50, r1w=0.55,
                     xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                     r2w=0.146, li=0.2, hs2=0.02, lamnu2=1.44,
                     r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                     r2kwa=0, pn=200000, u1=380.0, m1=3, mzone1=3,
                     prbg0=500, mz1s=3, n1=36, n2=24, weite1=8, schr=0,
                     mn=800, manzug=500, st=st, verbose=False)
    assert out["msymax"] >= 0
    assert isinstance(out["breakdown"], list)
    assert len(out["breakdown"]) > 0
    assert out["starting_endangered"] in (True, False)


def test_kdruck_cooling_and_lamination_areas():
    st = _build_full_state()
    out = rp.kdruck(ilfr=2, da1=0.4, di1=0.25, di2=0.24, l_=0.2, n1=36,
                     n2=24, hs1=0.001, bs1=0.006, hk1=0.001, bn1=0.009,
                     hs2=0.001, bn2s=0.008, br62=0.0, ho1=0.001,
                     hcuo1=0.005, hzw1=0.001, hcuu1=0.01, hu1=0.002,
                     ho2=0.0, hcuo2=0.005, hzw2=0.001, hcuu2=0.005,
                     hu2=0.0, h22=0, h32=0.02, h42=0.002, h52=0, h62=0,
                     h72=0, deltag=0.0006, zk=4, bk=0.01, na=1, da=0.02,
                     nk2=4, anut1=0.0002, anut2=0.0002, qstab=0.0002,
                     qstaba=0.0, drufi=0.025, k42=1e12, zn1=40, m1=3,
                     ksys=1.0, eptauh=1.5, st=st)
    assert out["iklink"] in (0, 1)
    assert out["k1_at_s1"] > 0
    assert out["ar1"] > 0


def test_kkwert_lookup():
    kka = rp.kkwert(ilfr=2, m1=3, mzone1=3, n1=36, n2=24, lstab=0.25,
                     l_=0.2, schr=0, p=2, weite1=8)
    assert 0.0 < kka <= 1.0


def _run_tlast(st, reke=" " * 12, pn=2000.0, prbg0=5.0):
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=pn, u1=380.0, m1=3, prbg0=prbg0, rzus=0,
             is1pkt=1, st=st)
    return rp.tlast(name="TEST", reke=reke, ilfr=2, p=2, fn=50, r1w=0.55,
                     xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
                     r2w=0.146, li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44,
                     r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
                     r2kwa=0, pn=pn, u1=380.0, m1=3, prbg0=prbg0,
                     qcu1=1e-5, qcu2=1e-4, qstab=1e-4, qstaba=0,
                     qring=1e-4, qringa=0, zn1=40, n1=36, zn2=1, n2=24,
                     di1=0.25, i1n_unused=0, mn=800, kapst=30e6,
                     k41=1e6, k42=1e6, k4a2=1e6, k4r2=1e6, k4ar2=1e6,
                     zk=4, nwkz=1.0, taup=0.09, ku1_unused=0, kz=1.0,
                     kt=1.0, ishalt=0, scha2=1, ksys=1.0, ski_unused=0,
                     kma=1.0, kia=1.0, l_=0.2, lm1=0.5, lm2=0.4,
                     lstab=0.25, mzone1=3, schr=0, weite1=8, st=st,
                     verbose=False)


def test_tlast_full_pipeline():
    """
    End-to-end test of the largest subroutine in the program: no-load,
    5 partial-load points, breakdown, theoretical + corrected locked-
    rotor points, XKEY time constants. Uses the same small self-
    consistent PN as GENERA's test, for the same reason (see that
    test's docstring).
    """
    st = _build_full_state()
    out = _run_tlast(st)
    assert out["iabort"] == 0
    assert all(p["converged"] for p in out["points"])
    assert out["kipp"]["m"] > 0
    assert out["theoretical"]["m"] > 0
    assert out["corrected"]["m"] > 0
    assert out["tk"] > 0
    assert out["xkey_result"]["thaupt"] > 0
    assert out["pvsum"] is not None and out["pvinn"] is not None
    assert out["generator_result"] is None
    assert out["dfmag_result"] is None


def test_tlast_dispatches_generator_mode():
    st = _build_full_state()
    st.art.rkz = " G" + " " * 10  # TLAST checks st.art.rkz[1], not reke
    out = _run_tlast(st)
    assert out["generator_result"] is not None
    assert all(p["converged"] for p in out["generator_result"]["points"])


def test_tlast_kopask_callback_fires():
    st = _build_full_state()
    calls = []
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=2000.0, u1=380.0, m1=3, prbg0=5.0, rzus=0,
             is1pkt=1, st=st)
    rp.tlast(name="TEST", reke=" " * 12, ilfr=2, p=2, fn=50, r1w=0.55,
             xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146,
             li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44, r2s=5.28e-5,
             r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0, pn=2000.0,
             u1=380.0, m1=3, prbg0=5.0, qcu1=1e-5, qcu2=1e-4, qstab=1e-4,
             qstaba=0, qring=1e-4, qringa=0, zn1=40, n1=36, zn2=1, n2=24,
             di1=0.25, i1n_unused=0, mn=800, kapst=30e6, k41=1e6,
             k42=1e6, k4a2=1e6, k4r2=1e6, k4ar2=1e6, zk=4, nwkz=1.0,
             taup=0.09, ku1_unused=0, kz=1.0, kt=1.0, ishalt=0, scha2=1,
             ksys=1.0, ski_unused=0, kma=1.0, kia=1.0, l_=0.2, lm1=0.5,
             lm2=0.4, lstab=0.25, mzone1=3, schr=0, weite1=8, st=st,
             kopask_callback=lambda *a: calls.append(a), verbose=False)
    assert len(calls) == 1


def test_kopask_builds_50_point_curve():
    st = _build_full_state()
    tlast_out = _run_tlast(st)
    points = rp.kopask(
        name="TEST", reke="M-TYPE-FRAME", ilfr=2, p=2, fn=50, r1w=0.55,
        xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146,
        li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44, r2s=5.28e-5,
        r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0, u1=380.0, m1=3,
        prbg0=5, sn=tlast_out["kipp"]["skipp"], mkipp=tlast_out["kipp"]["m"],
        skipp=tlast_out["kipp"]["skipp"], pn=2000, st=st)
    assert len(points) == 50
    assert all("pmech" in pt and "i1" in pt for pt in points)


def test_pdruck_produces_short_circuit_curve_and_reduced_voltage():
    st = _build_full_state()
    tlast_out = _run_tlast(st)
    out = rp.pdruck(
        irenr=1, ilfr=2, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
        xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
        lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
        r2kwa=0, pn=2000, u1=380.0, m1=3, prbg0=5, zn1=40, n1=36,
        xsip=0.92, n2=24, qstab=1e-4, k42=1e12, eptauh=1.5, r1k=0.5,
        r2k=0.1, ku1=0.7, xnetz=0.1, ksys=1.0,
        kka=tlast_out["corrected"]["kka"], st=st)
    assert len(out["short_circuit_curve"]) == 15
    assert out["reduced_voltage"] is not None
    assert out["reduced_voltage"]["u1str_over_u1"] > 0
    assert out["no_load"]["i1"] > 0
    assert out["eexe"] is not None  # ILFR=2 -> cage rotor -> Ex-e data present


# ---------------------------------------------------------------- report_export

def test_write_csv_report_produces_expected_files(tmp_path):
    st = _build_full_state()
    tlast_out = _run_tlast(st)
    synmom_out = rp.synmom(
        reke="TEST-MACHINE", ilfr=2, p=2, fn=50, r1w=0.55, xssti1=2.05,
        xsnut1=0.96, xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
        lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
        r2kwa=0, pn=2000, u1=380.0, m1=3, mzone1=3, prbg0=5, mz1s=3,
        n1=36, n2=24, weite1=8, schr=0, mn=800, manzug=500, st=st,
        verbose=False)
    kdruck_out = rp.kdruck(
        ilfr=2, da1=0.4, di1=0.25, di2=0.24, l_=0.2, n1=36, n2=24,
        hs1=0.001, bs1=0.006, hk1=0.001, bn1=0.009, hs2=0.001,
        bn2s=0.008, br62=0.0, ho1=0.001, hcuo1=0.005, hzw1=0.001,
        hcuu1=0.01, hu1=0.002, ho2=0.0, hcuo2=0.005, hzw2=0.001,
        hcuu2=0.005, hu2=0.0, h22=0, h32=0.02, h42=0.002, h52=0, h62=0,
        h72=0, deltag=0.0006, zk=4, bk=0.01, na=1, da=0.02, nk2=4,
        anut1=0.0002, anut2=0.0002, qstab=0.0002, qstaba=0.0,
        drufi=0.025, k42=1e12, zn1=40, m1=3, ksys=1.0, eptauh=1.5, st=st)
    mwinfo_out = rp.mwinfo(reke="M-TYPE-FRAME", v10=1.65, da1=0.4,
                            di1=0.25, deltag=0.0006)

    prefix = str(tmp_path / "report")
    paths = rex.write_csv_report(tlast_out, st.kennli.mkpool, st.kennli.eanz,
                                  synmom_out, kdruck_out, mwinfo_out,
                                  path_prefix=prefix)
    assert len(paths) == 6
    for p in paths:
        assert os.path.exists(p)
        with open(p) as f:
            rows = list(csv.reader(f))
        assert len(rows) >= 2  # header + at least one data row


def test_write_csv_report_skips_missing_sections(tmp_path):
    prefix = str(tmp_path / "report")
    paths = rex.write_csv_report(path_prefix=prefix)
    assert paths == []


def test_write_pdf_report_produces_valid_pdf(tmp_path):
    st = _build_full_state()
    tlast_out = _run_tlast(st)
    header = {"name": "TEST-MACHINE", "reke": "M-TYPE-FRAME",
              "datum": "10072026"}
    out_path = str(tmp_path / "report.pdf")
    result_path = rex.write_pdf_report(
        header_info=header, tlast_result=tlast_out,
        mkpool_rows=st.kennli.mkpool, eanz=st.kennli.eanz,
        output_path=out_path)
    assert result_path == out_path
    assert os.path.exists(out_path)
    assert os.path.getsize(out_path) > 1000
    with open(out_path, "rb") as f:
        assert f.read(5) == b"%PDF-"


# ---------------------------------------------------------------- main (run_asyn4)

def test_run_asyn4_full_pipeline_from_file(tmp_path):
    """
    The top-level integration test: runs the ENTIRE pipeline (LESEN
    through TLAST/KDRUCK/MAGZUG) starting from an EINASYN-format file
    on disk, exactly as a user would invoke it. This is what actually
    exercises the mm/kW -> SI conversion block, the heating-constant
    and friction-loss post-processing, PRUEF validation, and the
    hammerhead-slot/die-cast-rotor special cases - none of which are
    reachable by calling individual functions directly with pre-built
    MachineState objects the way every other test in this suite does.

    Not checked against a real machine's known-correct values (still no
    sample EINASYN available - see CONVERSION_STATUS.md) - this
    confirms the full pipeline runs to completion from a real input
    file and produces physically sane (positive mass/inertia, all
    partial-load points converged) results.
    """
    path = tmp_path / "test_einasyn.txt"
    _make_synthetic_einasyn(str(path))

    result = asyn4_main.run_asyn4(str(path), irenr=1, verbose=False,
                                   run_slip_curve=False)

    assert result.get("iabort") in (None, 0)
    tl = result["tlast_result"]
    assert all(p["converged"] for p in tl["points"])
    assert result["jaktiv"] > 0  # rotor inertia must be physically positive
    assert result["gaktiv"] > 0  # total active mass must be positive
    assert result["kdruck_result"]["iklink"] in (0, 1)
    assert result["magzug"]["ce"] > 0


def test_run_asyn4_rejects_invalid_data(tmp_path):
    """PRUEF's validation is now actually wired into the pipeline (it
    wasn't, until this integration test found the gap - see
    CONVERSION_STATUS.md round 21). A DI2 > DI1 data error should stop
    the pipeline with iabort=1 rather than silently continuing."""
    from asyn4.lesen import _CardWriter
    w = _CardWriter()
    w.write_mixed(2, [("a", 16, "BAD"), ("a", 8, ""), ("a", 25, "")])
    w.write_mixed(2, [("a", 27, "X"), ("a", 17, ""), ("a", 15, ""),
                       ("i", 5, 1), ("a", 12, "")])
    w.write_mixed(2, [("a", 11, "2")])
    w.write_f(3, [50, 0.95, 1.0, 1.0, 2, 400, 200, 0, 0, 0.05, 0], decimals=2)
    w.write_mixed(2, [("i", 7, 0)] + [("f", 7, v) for v in
                  [380, 50, 1.65, 0, 0, 0, 0.6, 0, 0, 0]])
    w.write_mixed(2, [("f", 7, v) for v in [3, 3, 250, 36, 40, 8, 0, 0, 75]] +
                  [("i", 7, 0), ("f", 7, 400)])
    w.write_mixed(2, [("f", 7, 6.0), ("f", 7, 1.0), ("i", 7, 0), ("a", 7, "N1"),
                       ("f", 7, 0), ("f", 7, 0), ("f", 7, 0), ("f", 7, 0)])
    w.write_f(2, [9.0, 10.0, 1.0, 1.0, 5.0, 1.0, 10.0, 2.0, 0, 0, 0], decimals=1)
    # DI2=500mm > DI1=250mm -> invalid, PRUEF should catch this
    w.write_mixed(2, [("f", 7, v) for v in
                  [500, 24, 5.0, 20.0, 4, 2, 0, 1, 1.2, 1.2]] + [("a", 5, "N2")])
    w.write_f(3, [1, 0, 0, 0, 75, 0, 0, 0, 0, 0])
    w.write_mixed(2, [("f", 7, v) for v in
                  [7.0, 8.0, 1.0, 1.0, 5.0, 1.0, 5.0, 0]] + [("i", 7, 1)])
    w.write_f(3, [20.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_f(2, [6.0, 5.0, 0, 0, 0, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_f(2, [30, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0], decimals=1)
    w.write_mixed(2, [("f", 7, 0, 4), ("f", 7, 300, 1), ("f", 7, 30, 1)] +
                  [("f", 7, 0, 4) for _ in range(8)])
    w.write_f(3, [0, 0, 0, 0, 0])
    w.write_mixed(3, [("f", 14, 0, 4), ("f", 14, 1, 4)] +
                  [("f", 7, v) for v in [1, 100, 0, 0.2, 50]])
    w.write_f(3, [0] * 11)
    w.write_mixed(3, [("a", 10, ""), ("a", 4, ""), ("a", 10, ""), ("a", 4, "")])
    path = tmp_path / "bad_einasyn.txt"
    w.write(str(path))

    result = asyn4_main.run_asyn4(str(path), irenr=1, verbose=False,
                                   run_slip_curve=False)
    assert result["iabort"] == 1
    assert any("DI2 GROESSER DI1" in text for _, text in result["pruef_messages"])


# ---------------------------------------------------------------- plotting

def test_plot_torque_speed_curve_requires_slast_first():
    st = make_state()
    with pytest.raises(ValueError):
        pl.plot_torque_speed_curve("X", mn=1, ku1=0, xnetz=0, st=st,
                                    output_path="/tmp/should_not_exist.png")


def test_plot_torque_speed_curve_produces_valid_png(tmp_path):
    st = _build_full_state()
    op.spkt(ilfr=2, s=0.03, p=2, fn=50, r1w=0.55, xssti1=2.05, xsnut1=0.96,
             xssti2=0.1, xsnut2=0.21, r2w=0.146, li=0.2, hs2=0.02,
             lamnu2=1.44, r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005,
             r2kwa=0, pn=2000, u1=380.0, m1=3, prbg0=5, rzus=0,
             is1pkt=1, st=st)
    slips = list(np.arange(0.05, 1.05, 0.05))
    st.svorga.sanz = len(slips)
    for i, s in enumerate(slips, start=1):
        st.svorga.svor[i] = s
    st.daten.mg[0] = 1
    st.daten.mg[1] = 0.3
    st.daten.mg[3] = 0.2
    st.daten.mg[4] = 0.5
    st.daten.mgn = 1.0

    lc.slast(reke="TEST-MACHINE", ilfr=2, ku1=0, xnetz=0, p=2, fn=50,
             r1w=0.55, xssti1=2.05, xsnut1=0.96, xssti2=0.1, xsnut2=0.21,
             r2w=0.146, li=0.2, hs2=0.02, bs2=0.004, lamnu2=1.44,
             r2s=5.28e-5, r2kw=1.32e-5, h42=0.002, br42=0.005, r2kwa=0,
             pn=2000, u1=380.0, m1=3, prbg0=5, i1n=15, mn=8, qcu1=1e-5,
             qcu2=1e-4, qstab=1e-4, qstaba=0, qring=1e-4, qringa=0, n2=24,
             kma=1.0, kia=1.0, st=st, verbose=False)

    out_path = str(tmp_path / "chart.png")
    result_path = pl.plot_torque_speed_curve(
        "TEST-MACHINE", mn=8.0, ku1=0.0, xnetz=0.0, st=st,
        output_path=out_path)
    assert result_path == out_path
    assert os.path.exists(out_path)
    with open(out_path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"


def test_run_asyn4_generates_plot_when_requested(tmp_path):
    """
    Confirms run_asyn4's plot_path parameter actually wires SLAST's
    output into plot_torque_speed_curve(), and that RKZ column 3 (the
    flag controlling whether the slip curve/plot runs at all, source
    "IF (RKZ(3:3) .EQ. ' ') GOTO 500") is respected - reusing the same
    synthetic file as the other run_asyn4 tests, which sets that column
    to a non-blank value.
    """
    path = tmp_path / "test_einasyn.txt"
    _make_synthetic_einasyn(str(path))
    plot_path = str(tmp_path / "chart.png")

    result = asyn4_main.run_asyn4(str(path), irenr=1, verbose=False,
                                   run_slip_curve=True, plot_path=plot_path)

    assert result.get("iabort") in (None, 0)
    if result.get("slast_result") is not None:
        assert result.get("plot_path") == plot_path
        assert os.path.exists(plot_path)
