"""
Fifth batch of ported FORTRAN subroutines from ASYN4.FOR - the magnetic
circuit solver and its direct dependencies.

Ported in this batch (source line ranges refer to ASYN4.FOR):
    SORTE  (8505-8731)  Lookup of one of 13 built-in magnetization curves
                         (B in 0.05 T steps, H in A/m) into st.cbsort.h,
                         selected by GUETE = round(V10*100). This is what
                         makes HEISEN/STEG/VZ physically meaningful - until
                         SORTE has been called, st.cbsort.h is all zeros.
    VZ     (9575-9614)  Magnetic tooth-voltage calculation at 5 points
                         along the tooth height (uses HEISEN).
    PLATT  (7724-7734)  "Abplattungsfaktor" (flattening factor) for the
                         combined tooth+air-gap MMF wave.
    KZWERT (7147-7232)  Heating-coefficient lookup table (frame
                         size/pole-count dependent), report from June 1990,
                         internal-air-circuit machines only, frame sizes
                         355-900mm, 2-12 poles.
    KREIS  (6947-7146)  THE magnetic circuit solver: generates the full
                         magnetization characteristic (MKPOOL table) over
                         50 points of increasing air-gap induction, calling
                         SORTE, VZ, PLATT, VJ, and VLZFE. This is the core
                         computation the whole rest of the program is built
                         around.

Not yet ported: LESEN (input reader), SLAST/HOLAUF (load/starting curves),
output/report subroutines, plotting driver, PROGRAM ASYN4 itself.
"""
from __future__ import annotations
import math
import numpy as np

from .state import MachineState
from .winding_functions import heisen
from .geometry_functions import vj, vlzfe

PI = 3.14159
MY0 = 1.2566e-6


def _rkz_char(rkz: str, pos: int) -> str:
    if pos - 1 < len(rkz):
        return rkz[pos - 1]
    return " "


# --- 13 built-in magnetization curves, B = 0.05*I Tesla (I=1..45), H in A/m ---
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

_CURVES = {
    110: (" 1.1 W/KG ", _H11), 115: (" 1.1 W/KG ", _H11), 120: (" 1.1 W/KG ", _H11),
    125: (" 1.25 W/KG", _H125), 135: (" 1.25 W/KG", _H125),
    150: (" 1.5 W/KG ", _H15),
    200: (" 2.0 W/KG ", _H20),
    230: (" 2.3 W/KG ", _H23),
    260: (" 2.6 W/KG ", _H26),
    300: (" 3.0 W/KG ", _H30),
    360: (" 3.6 W/KG ", _H36),
    420: (" 4.2 W/KG ", _H42),
    600: (" 6.0 W/KG ", _H60),
    0: (" SCHORCH  ", _H00),
    10: ("NUERNBERG ", _H01),
}


def sorte(guete: int, st: MachineState) -> str:
    """
    Loads the magnetization curve selected by GUETE into st.cbsort.h
    (H(0:45), index 0 = 0.0, indices 1..45 = H at B = 0.05*i Tesla).
    Falls back to the 1.7 W/kg curve (GUETE not matched) like the FORTRAN
    ELSE branch. Returns the BLECH label string.
    """
    st.cbsort.my0 = math.asin(1.0) * 8.0e-7
    st.cbsort.h[0] = 0.0
    blech, table = _CURVES.get(guete, (" 1.7 W/KG ", _H17))
    for i in range(1, 46):
        st.cbsort.h[i] = table[i - 1]
    return blech


def platt(vzvl: float):
    """Flattening ("Abplattung") factor for the combined tooth+gap MMF wave."""
    hg = 0.5 * vzvl - PI / 4.0
    return math.sqrt(hg * hg + vzvl) - hg


def vz(bluft: float, li: float, lfe: float, kfe: float, tn: float,
       brzah, brnut, hozah: float, n: float, cfez: float, st: MachineState):
    """
    Magnetic tooth-voltage at 5 support points along the tooth height.
    brzah/brnut are 1-indexed arrays/sequences (index 1..5 used, matching
    FORTRAN BRZAH(5)/BRNUT(5); index 0 ignored if present).
    Returns (vzahn, bzm, pfez, gz).
    """
    hg_weights = [0.0, 0.125, 0.25, 0.25, 0.25, 0.125]  # index 0 unused
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


_ARKZ_FLAT = [25., 25., 25., 25., 25., 23., 23., 25., 25.,
              28., 28., 28., 28., 28., 28., 23., 23., 28.,
              30., 30., 30., 30., 30., 30., 30., 30., 30.,
              30., 30., 30., 28., 28., 30., 25., 25., 28.,
              30., 30., 30., 30., 30., 25., 30., 25., 30.,
              30., 30., 30., 30., 30., 30., 30., 30., 30.]
_AJKZ_FLAT = [15., 14., 15., 13., 15., 15., 15., 15., 15.,
              20., 18., 20., 14., 18., 15., 13., 18., 18.,
              23., 23., 23., 20., 20., 15., 15., 20., 20.,
              20., 20., 20., 18., 19., 18., 20., 20., 20.,
              25., 25., 25., 25., 25., 25., 25., 25., 25.,
              25., 25., 25., 25., 25., 25., 25., 25., 25.]
# FORTRAN DATA fills ARKZ(9,6) column-major (first index fastest) - use
# order='F' so _ARKZ[iah, ip] matches ARKZ(IAH,IP) with iah,ip 0-based;
# pad with a dummy row/column at index 0 so 1-based FORTRAN indices work
# directly (_ARKZ[iah][ip] for iah=1..9, ip=1..6).
_ARKZ = np.zeros((10, 7))
_ARKZ[1:10, 1:7] = np.array(_ARKZ_FLAT).reshape(9, 6, order="F")
_AJKZ = np.zeros((10, 7))
_AJKZ[1:10, 1:7] = np.array(_AJKZ_FLAT).reshape(9, 6, order="F")

_AH_CODES = {"355": 1, "400": 2, "450": 3, "500": 4, "560": 5,
             "630": 6, "710": 7, "800": 8, "900": 9}


def kzwert(reke: str, p: float) -> float:
    """
    Heating-coefficient lookup (June 1990 report), internal-air-circuit
    machines only, frame sizes 355-900mm, 2-12 poles (P=1..6 pole pairs).
    reke is the 70-char frame-type code string (FORTRAN 1-based
    substrings REKE(1:2) and REKE(4:6)).
    """
    ip = round(p)
    if ip > 6:
        return 0.0
    if reke[0:2] == "AM":
        return 0.0
    iah = _AH_CODES.get(reke[3:6])
    if iah is None:
        return 0.0
    if reke[0:2] in ("AJ", "AW", "AO", "AN"):
        return float(_AJKZ[iah][ip])
    elif reke[0:2] == "AR":
        return float(_ARKZ[iah][ip])
    return 0.0


def kreis(n1: float, n2: float, p: float, kfe: float, rb: float, li: float,
          lfe: float, lfe2: float, da1: float, di2: float,
          hojoc1: float, kj1: float, hojoc2: float, kj2: float,
          deltag: float, bs1: float, bs2: float, hs1: float,
          hozah1: float, brzah1, hozah2: float, brzah2,
          cfej: float, cfez: float, v10: float,
          u1: float, fn: float, w1: float, xsip: float, m1: float,
          st: MachineState):
    """
    Magnetic circuit solver: generates the magnetization characteristic
    (populates st.kennli.mkpool, up to 50 rows, 18 columns) as air-gap
    induction BL is swept up from -0.0749 T until the induced voltage UI
    exceeds 1.5x the rated voltage U1.

    Returns a dict of the outputs that used to be plain arguments:
    kc1, kc2, gj1, gz1, gj2, gz2, eanz, bs1s, cfei (cfej/cfez are updated
    in place via the returned dict since the FORTRAN mutates them too).

    brzah1/brzah2 are 1-indexed sequences (index 1..5 used).
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
            bs1s = 0.0001 * j + 1.0e-8
            kc1s = tn1 / (tn1 - bs1s * bs1s / (5.0 * deltag + bs1s))
            if kc1s >= kc1:
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
        brnut1[j] = (PI * (da1 - 2.0 * hojoc1 - 2.0 * hozah1 + 2.0 * (j - 0.5) / 5.0
                           * hozah1) / n1 - brzah1[j])
        brnut2[j] = (PI * (di2 + 2.0 * hojoc2 + 2.0 * hozah2 - 2.0 * (j - 0.5) / 5.0
                           * hozah2) / n2 - brzah2[j])

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

        vzahn1, bz1max, pfez1, gz1 = vz(bl, li, lfe, kfe, tn1, brzah1, brnut1,
                                         hozah1, n1, cfez, st)
        vzahn2, bz2max, pfez2, gz2 = vz(bl, li, lfe2, kfe, tn2, brzah2, brnut2,
                                         hozah2, n2, cfez, st)
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
        "eanz": eanz, "bs1s": bs1s, "cfei": cfei, "cfej": cfej, "cfez": cfez,
    }
