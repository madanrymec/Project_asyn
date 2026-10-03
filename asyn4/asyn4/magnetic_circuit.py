"""
Third batch of ported FORTRAN subroutines from ASYN4.FOR - stator geometry
and cage-rotor slot-bridge (Steg) leakage curve.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    STD  (9158-9285)  Berechnung einiger Staendergroessen
                       (stator slot leakage, tooth widths, resistance, yoke
                       height, end-winding/slot leakage reactances).
                       Updates st.keil (XSKEIL/TETA0).
    STEG (9286-9352)  Generierung der Streuschlitzaufweitungskurve beim
                       Kaefiglaeufer (cage-rotor slot-bridge saturation
                       curve). Uses HEISEN from winding_functions, needs
                       st.cbsort.h / st.cbsort.my0 populated.
    SORTE(8505-8697)  Laedt eine der 13 eingebauten Magnetisierungskurven
                       in st.cbsort.h / st.cbsort.my0 (this is what
                       populates the table STEG/VZ/HEISEN read from).
    VZ   (9575-9599)  Magnetische Zahnspannungen und Zahneisenverluste an
                       5 Stuetzstellen (uses HEISEN).
    PLATT(7724-7734)  Abplattungsfaktor (flattening factor) fuer die
                       Luftspaltfeldkurve.
    KREIS(6947-7146)  Durchrechnung des gesamten magnetischen Kreises:
                       generiert die Magnetisierungskennlinie UI=f(IMY)
                       ueber 50 Luftspaltinduktionswerte (fills MKPOOL /
                       st.kennli.mkpool). The central magnetic-circuit
                       solver; calls SORTE, VZ, PLATT, VJ, VLZFE.
    KZWERT(7147-7232) Erwaermungsziffer (heating coefficient) lookup table
                       for internally-cooled machines, frame sizes 355-900,
                       2-12 poles.

Not yet ported from this functional area: MKWERT.
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState
from .winding_functions import xsiz, heisen
from .geometry_functions import vj, vlzfe

PI = 3.14159
MY0 = 1.2566e-6


def _rkz_char(rkz: str, pos: int) -> str:
    if pos - 1 < len(rkz):
        return rkz[pos - 1]
    return " "


def std(da1: float, di1: float, n1: float, p: float, weite1: float,
        mzone1: float, zn1: float, lm1: float, qcu1: float, m1: float,
        fn: float, li: float, theta1: float, mz1s: float,
        hs1: float, hk1: float, ho1: float, hcuo1: float, hzw1: float,
        hcuu1: float, hu1: float,
        bs1: float, bn1: float, bn1s: float, lamst1_in: float,
        st: MachineState):
    """
    Stator geometric/electric quantities: yoke height, magnetically
    effective tooth length, tooth widths at 5 points, slot leakage
    permeance, ohmic resistance, end-winding leakage LAMST1 (computed
    from empirical tables if not already supplied > 0.001), and the
    stator end-winding (XSSTI1) and slot (XSNUT1) leakage reactances.

    lamst1_in mirrors the FORTRAN in/out argument LAMST1: pass 0 (or any
    value <= 0.001) to have it computed from the built-in empirical table;
    pass a value > 0.001 to keep it fixed (as the FORTRAN GOTO 800 does).

    Returns a dict of all computed outputs (brzah1 is a 1-indexed
    np.ndarray of length 6, slots 1..5 used, matching FORTRAN BRZAH1(5)).
    """
    rkz = st.art.rkz
    ksys = st.vorzei.ksys
    bs1di1 = st.vorzei.bs1di1

    hojoc1 = 0.5 * (da1 - di1) - hs1 - hk1 - ho1 - hcuo1 - hzw1 - hcuu1 - hu1

    if bs1 > 0.9 * bn1:
        hozah1 = hs1 + hk1 + ho1 + hcuo1 + hzw1 + hcuu1 + hu1
    else:
        hozah1 = 0.5 * hk1 + ho1 + hcuo1 + hzw1 + hcuu1 + hu1

    brzah1 = np.zeros(6)
    if bs1 > 0.9 * bn1:
        brzah1[1] = PI * (di1 + hs1 + hk1) / n1 - (bs1 + bn1) / 2.0
    else:
        brzah1[1] = PI * (di1 + 2.0 * hs1 + 2.0 * hk1) / n1 - bn1
    hdiff = ho1 + hcuo1 + hzw1 + hcuu1 + hu1
    for j in range(2, 6):
        brzah1[j] = (PI * (di1 + 2.0 * hs1 + 2.0 * hk1 + 2.0 * hdiff *
                            (j - 2.0 + 0.5) / 4.0) / n1 -
                      (bn1 + (bn1s - bn1) * (j - 2.0 + 0.5) / 4.0))

    lamstg = hs1 / bs1 + hk1 / (bs1 + (bn1 - bs1) / 3.0)
    lamo = (0.33 * hcuo1 / (bn1 + (bn1s - bn1) * ho1 / hdiff) +
            lamstg + ho1 / bn1)
    lamu = (0.33 * hcuu1 / (bn1 + (bn1s - bn1) * (ho1 + hcuo1 + hzw1) / hdiff) +
            (ho1 + hcuo1 + hzw1) / (bn1 + 0.33 * (bn1s - bn1) *
                                     (ho1 + hcuo1 + hzw1) / hdiff) + lamstg)
    lamou = (0.5 * hcuo1 / (bn1 + (bn1s - bn1) * ho1 / hdiff) +
             lamstg + ho1 / bn1)

    q1 = n1 / (p * mzone1)
    hgsum = (xsiz(p, q1, n1, st) ** 2 * math.cos(2.0 * PI * p * weite1 / n1) *
              (p * bs1di1 / math.sin(p * bs1di1)) ** 2)
    for g in range(1, 101):
        for gg in (1, 2):
            ny = p * (1.0 + mz1s * g * (-1) ** gg)
            if ny > n1 / 2.0 or ny <= -n1 / 2.0:
                continue
            hgsum += (xsiz(ny, q1, n1, st) ** 2 * math.cos(2.0 * PI * ny * weite1 / n1) *
                      (ny * bs1di1 / math.sin(ny * bs1di1)) ** 2)
    if _rkz_char(rkz, 9) == "E":
        hgsum = -1.0
    lamnu1 = (lamo + lamu - 2.0 * lamou * hgsum) / 4.0

    r1k = zn1 * n1 * lm1 / (qcu1 * 58.0e6 * m1) / ksys
    r1w = r1k * (1.0 + 0.00429 * theta1)

    lamst1 = lamst1_in
    if lamst1 <= 0.001:
        if abs(bn1 - bn1s) < 0.0001:
            p_i = round(p)
            if p_i == 1:
                lamst1 = 0.150
            elif p_i == 2:
                lamst1 = 0.175
            else:
                lamst1 = 0.200
        else:
            p_i = round(p)
            if p_i == 1:
                lamst1 = 0.200
            else:
                lamst1 = 0.150
        m1_i = round(m1)
        if m1_i == 5:
            lamst1 *= 1.70
        if m1_i == 6:
            lamst1 *= 2.00
        if _rkz_char(rkz, 9) == "E":
            lamst1 *= 1.75

    xssti1 = (2.0 * PI * fn * MY0 * (lm1 - li) * q1 * zn1 ** 2 * n1 / m1 *
              lamst1 / ksys ** 2)

    xsnut1 = 2.0 * PI * fn * MY0 * li * zn1 ** 2 * n1 / m1 * lamnu1 / ksys
    keil = st.keil
    if _rkz_char(rkz, 4) == "M":
        lakeil = hk1 / bs1
        keil.xskeil = xsnut1 * lakeil / lamnu1
    else:
        keil.xskeil = 0.0
    keil.teta0 = 2000.0 * bs1 / 0.015 / zn1

    return {
        "hojoc1": hojoc1, "hozah1": hozah1, "brzah1": brzah1,
        "lamu": lamu, "lamo": lamo, "lamou": lamou, "lamnu1": lamnu1,
        "r1k": r1k, "r1w": r1w, "lamst1": lamst1,
        "xssti1": xssti1, "xsnut1": xsnut1,
    }


_BSTEG = np.array([0.0, 0.000001, 0.8, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5,
                    2.7, 3.0, 3.3, 3.6, 4.0, 4.5, 5.0, 5.7, 6.5, 7.5, 8.5])


def steg(ilfr: int, li: float, lfe: float, kfe: float, hs2: float, bs2: float,
         br12: float, br22: float, br32: float, br42: float, br52: float,
         bstr52: float, st: MachineState):
    """
    Generates the slot-bridge (Steg) widening curve for a cage rotor:
    20 points of (INUT, LAMDA) - slot induction/leakage vs. bridge current
    - used later for interpolation by LASTG. Requires st.cbsort.h and
    st.cbsort.my0 to be populated (HEISEN dependency).

    Returns (inut, lamda), 1-indexed np.ndarrays of length 21 (slots 1..20).
    """
    if ilfr == 3 or ilfr == 4:
        brauf = 1.2 * br32
    else:
        if br32 > 0.001 or br22 > 0.001:
            brauf = 1.2 * max(br12, br22, br32)
        elif br42 > 0.001:
            brauf = 1.2 * br42
        else:
            brauf = 1.2 * max(br52, bstr52)

    bs2s = 0.05 * hs2 if bs2 < 0.05 * hs2 else bs2

    f = MY0 * (1.0 - kfe) / kfe

    inut = np.zeros(21)
    lamda = np.zeros(21)
    for j in range(1, 21):
        vlu = _BSTEG[j] / MY0 * bs2s
        bfe = _BSTEG[j] * li / (lfe * kfe)
        vfe = heisen(bfe, f, st) * (brauf - bs2s)
        if j == 1:
            vfe = 0.0
        bs2i = bs2s * (vlu + vfe) / vlu
        inut[j] = vlu + vfe
        lamda[j] = hs2 / bs2i
    return inut, lamda


# ---------------------------------------------------------------------------
# SORTE - the 13 built-in magnetization curves, H in A/m at 0.05 T steps.
# ---------------------------------------------------------------------------
_H00 = [4.8, 9.7, 15.1, 20.9, 27.5, 35., 43.4, 53., 63.6, 75.4, 88.3,
        102., 117., 132., 148., 165., 181., 199., 218., 240., 266., 300.,
        348., 414., 510., 647., 843., 1119., 1502., 2029., 2744., 3703.,
        4972., 6637., 8797., 11570., 15100., 19600., 25200., 32200.,
        40900., 51600., 64700., 80600., 99800.]
_H01 = [10., 25., 38., 48., 58., 78., 93., 100., 118., 132., 155., 180.,
        208., 225., 270., 300., 335., 385., 450., 500., 620., 700.,
        850., 1000., 1200., 1430., 1750., 2050., 2500., 3000., 4000.,
        5000., 6850., 8600., 10600., 13100., 15000., 20000., 25500.,
        31500., 44000., 56000., 78000., 100000., 136000.]
_H11 = [12., 18., 24., 28., 32., 36., 40., 44., 48., 52., 55., 58., 62.,
        68., 72., 78., 84., 90., 98., 108., 118., 132., 151., 168., 214.,
        360., 560., 860., 1280., 2000., 3300., 4600., 6200., 8000.,
        10600., 13400., 16800., 21600., 28800., 42000., 73000.,
        113000., 153000., 193000., 233000.]
_H125 = [14., 24., 30., 36., 40., 44., 47., 49., 52., 56., 60., 64., 70.,
         75., 80., 86., 94., 102., 114., 125., 140., 158., 178., 208.,
         246., 320., 440., 640., 1280., 2100., 3000., 4800., 6200.,
         8200., 10800., 13700., 17800., 22800., 29200., 39000.,
         65000., 105000., 145000., 185000., 225000.]
_H15 = [22., 32., 40., 46., 50., 54., 58., 62., 64., 68., 72., 77., 82.,
        86., 92., 98., 106., 116., 127., 140., 156., 176., 204., 240.,
        300., 380., 490., 680., 1120., 2000., 3000., 4000., 5300., 7200.,
        9600., 12600., 16000., 20000., 27600., 34000., 42000., 52000.,
        65000., 80000., 97000.]
_H17 = [24., 36., 45., 52., 57., 62., 66., 71., 74., 77., 81., 85., 89.,
        94., 99., 106., 114., 121., 132., 144., 158., 176., 198., 228.,
        276., 340., 430., 600., 880., 1400., 2100., 3200., 4600.,
        6400., 8800., 11600., 14800., 18800., 23200., 28400., 35000.,
        41500., 49800., 62300., 87300.]
_H20 = [28., 43., 54., 61., 66., 71., 74., 78., 82., 86., 89., 93., 97.,
        102., 108., 114., 122., 131., 142., 153., 168., 184., 210., 245.,
        290., 360., 440., 580., 790., 1200., 1840., 2800., 4000.,
        5600., 8000., 10600., 13800., 17600., 22200., 27000., 34000.,
        39400., 47800., 60300., 85300.]
_H23 = [40., 53., 64., 72., 80., 86., 93., 98., 104., 110., 116., 122.,
        128., 134., 140., 146., 154., 162., 170., 183., 197., 214., 238.,
        274., 340., 400., 480., 560., 640., 960., 1460., 2130., 3200.,
        4800., 6700., 9500., 12500., 16000., 20600., 25400., 31000.,
        37100., 45500., 57900., 82900.]
_H26 = [38., 58., 72., 82., 90., 97., 102., 108., 113., 118., 124., 130.,
        136., 142., 148., 156., 163., 172., 182., 195., 210., 226.,
        246., 265., 290., 340., 410., 500., 700., 1040., 1560., 2320.,
        3400., 4800., 6700., 9300., 12400., 15600., 20000., 25000.,
        34000., 45100., 59400., 79400., 104400.]
_H30 = [48., 68., 84., 94., 102., 108., 114., 119., 124., 130., 136.,
        142., 148., 154., 162., 170., 178., 188., 198., 212., 226., 244.,
        260., 284., 330., 380., 460., 560., 700., 940., 1320., 1860.,
        2620., 3700., 5300., 7400., 10000., 13200., 17000., 22000.,
        30000., 41000., 57600., 80300., 108000.]
_H36 = [66., 90., 110., 122., 132., 140., 146., 152., 157., 162., 167.,
        173., 178., 184., 190., 198., 206., 214., 223., 233., 244., 256.,
        270., 284., 300., 340., 400., 500., 620., 840., 1180., 1650., 2400.,
        3400., 4800., 6600., 9200., 12200., 16000., 20800., 28000.,
        41300., 59200., 80000., 106400.]
_H42 = [40., 70., 90., 120., 140., 160., 170., 185., 195., 200., 204.,
        208., 212., 216., 220., 224., 228., 232., 236., 240., 246., 255.,
        265., 290., 320., 380., 450., 530., 660., 840., 1120., 1600.,
        2360., 3600., 5000., 6400., 8600., 11600., 16000., 22000.,
        29000., 41000., 57600., 80300., 108000.]
_H60 = [40., 70., 90., 120., 140., 160., 180., 210., 230., 245., 260.,
        280., 300., 310., 325., 340., 360., 380., 400., 420., 445., 480.,
        520., 560., 630., 700., 780., 900., 1040., 1280., 1620., 2240.,
        3400., 4500., 6000., 7800., 10000., 13000., 18200., 24000.,
        32000., 43000., 59600., 82300., 110000.]

_GUETE_TABLE = {
    110: (_H11, " 1.1 W/KG "), 115: (_H11, " 1.1 W/KG "), 120: (_H11, " 1.1 W/KG "),
    125: (_H125, " 1.25 W/KG"), 135: (_H125, " 1.25 W/KG"),
    150: (_H15, " 1.5 W/KG "),
    200: (_H20, " 2.0 W/KG "),
    230: (_H23, " 2.3 W/KG "),
    260: (_H26, " 2.6 W/KG "),
    300: (_H30, " 3.0 W/KG "),
    360: (_H36, " 3.6 W/KG "),
    420: (_H42, " 4.2 W/KG "),
    600: (_H60, " 6.0 W/KG "),
    0: (_H00, " SCHORCH  "),
    10: (_H01, "NUERNBERG "),
}


def sorte(guete: int, st: MachineState) -> str:
    """
    Loads one of the 13 built-in magnetization curves into
    st.cbsort.h / st.cbsort.my0. Falls back to the 1.7 W/kg curve
    (the FORTRAN default ELSE branch) for any unrecognized GUETE.
    Returns the sheet-steel grade label (BLECH).
    """
    st.cbsort.my0 = math.asin(1.0) * 8.0e-7
    table, blech = _GUETE_TABLE.get(guete, (_H17, " 1.7 W/KG "))
    st.cbsort.h[0] = 0.0
    st.cbsort.h[1:46] = table
    return blech


def platt(vzvl: float):
    """Flattening factor (Abplattungsfaktor) for the air-gap field curve."""
    hg = 0.5 * vzvl - PI / 4.0
    return math.sqrt(hg * hg + vzvl) - hg


def vz(bluft: float, li: float, lfe: float, kfe: float, tn: float,
       brzah, brnut, hozah: float, n: float, cfez: float, st: MachineState):
    """
    Magnetic tooth voltages and tooth iron losses at 5 support points.
    brzah/brnut are 1-indexed arrays/sequences of length >= 6 (slots 1..5).
    Returns (vzahn, bzm, pfez, gz).
    """
    hg_weights = {1: 0.125, 2: 0.25, 3: 0.25, 4: 0.25, 5: 0.125}
    fak1 = bluft * tn * li / (lfe * kfe)
    fak2 = 7750.0 * lfe * kfe * hozah / 5.0 * n
    gz = 0.0
    pfez = 0.0
    vzahn = 0.0
    bzm = 0.0
    for j in range(1, 6):
        bz = fak1 / brzah[j]
        if bz > bzm:
            bzm = bz
        f = MY0 * (1.0 - kfe + brnut[j] / brzah[j]) / kfe
        vzahn += heisen(bz, f, st) * hozah * hg_weights[j]
        pfez += fak2 * brzah[j] * bz ** 2 * cfez
        gz += fak2 * brzah[j]
    return vzahn, bzm, pfez, gz


def kreis(n1: float, n2: float, p: float, kfe: float, rb: float, li: float,
          lfe: float, lfe2: float, da1: float, di2: float,
          hojoc1: float, kj1: float, hojoc2: float, kj2: float,
          deltag: float, bs1: float, bs2: float, hs1: float,
          hozah1: float, brzah1, hozah2: float, brzah2,
          cfej: float, cfez: float, v10: float,
          u1: float, fn: float, w1: float, xsip: float, m1: float,
          st: MachineState):
    """
    Runs the whole magnetic circuit: generates the magnetization
    characteristic UI=f(IMY) over up to 50 air-gap induction points and
    fills st.kennli.mkpool / st.kennli.eanz.

    brzah1/brzah2 are 1-indexed arrays/sequences of length >= 6.

    Returns a dict with kc1, kc2, gj1, gz1, gj2, gz2, bs1s, cfei, plus the
    (possibly-updated) cfej/cfez (the FORTRAN in/out CFEJ/CFEZ), matching
    the FORTRAN OUT arguments. st.kennli.mkpool / eanz hold the generated
    curve, one row per air-gap induction step (18 columns, see MKPOOL
    layout comments in the FORTRAN source).
    """
    rkz = st.art.rkz

    guete = int(v10 * 100.0 + 0.1)
    sorte(guete, st)

    tn1 = 2.0 * PI * rb / n1
    tn2 = 2.0 * PI * rb / n2
    kc1 = tn1 / (tn1 - bs1 * bs1 / (5.0 * deltag + bs1))
    kc1fe = kc1
    bs1s = bs1
    if _rkz_char(rkz, 4) == "M":
        bs1s = 0.67 * bs1
        dts = deltag + max(hs1 - 0.00025, 0.0)
        kcn = bs1 / (bs1 - bs1s * bs1s / (bs1s + 5.0 * dts))
        dtss = dts * kcn
        dt1 = deltag + bs1 / tn1 * (dtss - deltag)
        kc1o = kc1
        kc1m = dt1 / deltag
        kc1 = min(kc1o, kc1m)
        kc1fe = tn1 / (tn1 - 0.49 * bs1 * bs1 / (5.0 * deltag + 0.7 * bs1))

    kc2 = tn2 / (tn2 - bs2 * bs2 / (5.0 * deltag + bs2))
    kcfe = kc1fe * kc2

    if _rkz_char(rkz, 4) == "M":
        for j in range(0, 1001):
            bs1s_try = 0.0001 * j + 1e-8
            kc1s = tn1 / (tn1 - bs1s_try * bs1s_try / (5.0 * deltag + bs1s_try))
            if kc1s >= kc1:
                bs1s = bs1s_try
                break
    else:
        bs1s = 1.00 * bs1

    lamda0 = MY0 / (kc1 * kc2 * deltag)

    v10fn = vlzfe(v10, fn)
    if v10 <= 1.35:
        cfei = 1.65 * kcfe ** 2 * v10fn
    elif v10 <= 1.70:
        cfei = 1.40 * kcfe ** 2 * v10fn
    elif v10 <= 2.00:
        cfei = 1.25 * kcfe ** 2 * v10fn
    elif v10 <= 2.30:
        cfei = 1.20 * kcfe ** 2 * v10fn
    else:
        cfei = 1.15 * kcfe ** 2 * v10fn
    if cfej + cfez == 0.0:
        cfej = cfei
        cfez = cfei

    lj1 = PI * (da1 - hojoc1) / (4.0 * p)
    lj2 = PI * (di2 + hojoc2) / (4.0 * p)
    gj1 = 7750.0 * 4.0 * p * lj1 * hojoc1 * lfe * kfe
    gj2 = 7750.0 * 4.0 * p * lj2 * hojoc2 * lfe2 * kfe

    brnut1 = np.zeros(6)
    brnut2 = np.zeros(6)
    for j in range(1, 6):
        brnut1[j] = (PI * (da1 - 2.0 * hojoc1 - 2.0 * hozah1 + 2.0 * (j - 0.5) / 5.0 * hozah1) / n1
                     - brzah1[j])
        brnut2[j] = (PI * (di2 + 2.0 * hojoc2 + 2.0 * hozah2 - 2.0 * (j - 0.5) / 5.0 * hozah2) / n2
                     - brzah2[j])

    imyfak = PI * p / (1.414 * xsip * m1 * w1)
    uifak = 4.443 * fn * w1 * xsip * PI * rb * li / p
    fakj1 = li * PI * rb / p / (2.0 * lfe * kfe * hojoc1 * kj1)
    fakj2 = li * PI * rb / p / (2.0 * lfe2 * kfe * hojoc2 * kj2)

    bl = -0.0749
    ui = 0.0
    eanz = 0
    mkpool = st.kennli.mkpool
    gz1 = gz2 = 0.0

    for _j in range(1, 51):
        if ui > 1.50 * u1:
            break
        elif ui > 1.02 * u1:
            bl += 0.05
        elif ui > 0.90 * u1:
            bl += 0.008
        elif ui > 0.80 * u1:
            bl += 0.025
        elif ui > 0.70 * u1:
            bl += 0.040
        else:
            bl += 0.075

        vluft = bl / lamda0

        vzahn1, bz1max, pfez1, gz1 = vz(bl, li, lfe, kfe, tn1, brzah1, brnut1, hozah1, n1, cfez, st)
        vzahn2, bz2max, pfez2, gz2 = vz(bl, li, lfe2, kfe, tn2, brzah2, brnut2, hozah2, n2, cfez, st)
        alpha = platt((vzahn1 + vzahn2) / vluft)

        blm = bl / alpha
        bjoch1 = fakj1 * blm
        bjoch2 = fakj2 * blm
        vjoch1 = vj(bjoch1, lj1, p)
        vjoch2 = vj(bjoch2, lj2, p)

        pfe1 = pfez1 + cfej * gj1 * bjoch1 ** 2
        vsum = vluft + vzahn1 + vzahn2 + vjoch1 + vjoch2
        imy = vsum * imyfak
        ui = blm * uifak
        ife = pfe1 / (m1 * u1)
        hg1 = vluft / vsum
        hg2 = vluft / (vluft + vjoch1 + vjoch2)
        bp = blm * PI / 2.0
        b3p = (hg1 - hg2) / (hg2 + 3.0 * hg1) * bp
        pfej1 = pfe1 - pfez1
        deltap = deltag * kc1 * kc2 * vsum / vluft * alpha / 1.571

        eanz += 1
        mkpool[eanz, 1] = ui
        mkpool[eanz, 2] = imy
        mkpool[eanz, 3] = deltap
        mkpool[eanz, 4] = bp
        mkpool[eanz, 5] = b3p
        mkpool[eanz, 6] = pfej1
        mkpool[eanz, 7] = pfe1
        mkpool[eanz, 8] = ife
        mkpool[eanz, 9] = bz1max
        mkpool[eanz, 10] = bz2max
        mkpool[eanz, 11] = bjoch1
        mkpool[eanz, 12] = bjoch2
        mkpool[eanz, 13] = vluft
        mkpool[eanz, 14] = vzahn1
        mkpool[eanz, 15] = vzahn2
        mkpool[eanz, 16] = vjoch1
        mkpool[eanz, 17] = vjoch2
        mkpool[eanz, 18] = alpha

    st.kennli.eanz = eanz

    return {
        "kc1": kc1, "kc2": kc2, "gj1": gj1, "gz1": gz1, "gj2": gj2, "gz2": gz2,
        "bs1s": bs1s, "cfei": cfei, "cfej": cfej, "cfez": cfez, "eanz": eanz,
    }


_ARKZ = np.array([
    [0, 0, 0, 0, 0, 0, 0],  # row 0 unused (1-indexed pole columns 1..6)
    [0, 25., 28., 30., 30., 30., 30.],
    [0, 25., 28., 30., 30., 30., 30.],
    [0, 25., 28., 30., 30., 30., 30.],
    [0, 25., 28., 30., 28., 30., 30.],
    [0, 25., 28., 30., 28., 30., 30.],
    [0, 23., 28., 30., 30., 25., 30.],
    [0, 23., 23., 30., 25., 30., 30.],
    [0, 25., 23., 30., 25., 25., 30.],
    [0, 25., 28., 30., 28., 30., 30.],
])


_AJKZ = np.array([
    [0, 0, 0, 0, 0, 0, 0],
    [0, 15., 20., 23., 20., 25., 25.],
    [0, 14., 18., 23., 20., 25., 25.],
    [0, 15., 20., 23., 20., 25., 25.],
    [0, 13., 14., 20., 18., 25., 25.],
    [0, 15., 18., 20., 19., 25., 25.],
    [0, 15., 15., 15., 18., 25., 25.],
    [0, 15., 13., 15., 20., 25., 25.],
    [0, 15., 18., 20., 20., 25., 25.],
    [0, 15., 18., 20., 20., 25., 25.],
])


def kzwert(reke: str, p: float) -> float:
    """
    Heating coefficient (Erwaermungsziffer) for internally-cooled
    machines, frame sizes 355-900, 2-12 poles (June 1990 report).
    reke should be the frame designation string (e.g. 'AJ  500...').
    """
    if round(p) > 6:
        return 0.0
    if reke[0:2] == "AM":
        return 0.0

    frame_code = reke[3:6]
    frame_to_row = {"355": 1, "400": 2, "450": 3, "500": 4, "560": 5,
                    "630": 6, "710": 7, "800": 8, "900": 9}
    if frame_code not in frame_to_row:
        return 0.0
    iah = frame_to_row[frame_code]
    ip = round(p)
    if ip < 1 or ip > 6:
        return 0.0

    prefix = reke[0:2]
    if prefix in ("AJ", "AW", "AO", "AN"):
        return float(_AJKZ[iah, ip])
    elif prefix == "AR":
        return float(_ARKZ[iah, ip])
    return 0.0
