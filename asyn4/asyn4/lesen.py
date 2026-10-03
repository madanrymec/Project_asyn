"""
LESEN (source 1586-2054 of ASYN4.FOR) - reads the EINASYN input file and
populates the machine-data record (st.daten / st.art.rkz / st.materi),
then runs the stator/rotor winding design (ISOSTD or ISOTR, and ISOLFR
for slip-ring rotors) exactly as the FORTRAN does at the end of reading.

*** VALIDATION STATUS - READ THIS FIRST ***
This is a best-effort reconstruction of LESEN's fixed-column FORTRAN
READ/FORMAT statements (documented field-by-field below, each traceable
to a specific ASYN4.FOR line range). It has NOT been validated against a
real EINASYN file, because none has been available during this port -
see CONVERSION_STATUS.md. Fixed-format FORTRAN reads are notoriously
easy to get subtly wrong (off-by-one column, wrong field width) in ways
that silently produce wrong numbers instead of an error. Treat this as a
strong starting point, not a verified reader, until you can run it
against a real EINASYN and check the resulting st.daten fields against
known-correct values (or against a print of what the original ASYN4.FOR
read for the same file, if you still have access to a working copy).

FILE FORMAT (19 "cards", labelled Z1-Z19 below, matching the FORTRAN's
own C-comments) - each FORMAT string is quoted from the source:

  Z1  '(//,A16,a8,a25)'        NAME(16), DATUM(8), KENNWORT(25)
  Z2  '(//,A27,a17,a15,i5,a12)' REKE(27), VANR(17), FBNR(15), WANR(i5), MONR(12)
  Z3  '(//,A11)'               RKZ(11)
  Z4  '(///,11f7.0)'           FN,KFE,KJ1,KJ2,P,DA1,L,ZK,BK,PN,KZ
  Z5  '(//,i7,10f7.0)'         ISHALT(i7),UN,PRBG0,V10,CFEJ,CFEZ,SCHR,DELTAG,KU1,XNETZ,KT
  Z6  '(//,9f7.0,i7,f7.0)'     M1,MZONE1,DI1,N1,ZN1,WEITE1,QCU1,LW1,THETA1,IKOPP(i7),UNX
  Z7  '(//,2f7.0,i7,a7,4f7.0)' BS1,HS1,IZGF(i7),NUTNR1(a7),NUTF1,NA,DA,DTA
  Z8  '(//,11f7.0)'            BN1,BN1S,HK1,HO1,HCUO1,HZW1,HCUU1,HU1,LAMST1,NUTF2 (11th field unused)
  Z9  '(//,10f7.0,a5)'         DI2,N2,BS2,HS2,MTL,NTL,UTL,AZWEIG,BTL,HTL,NUTNR2(a5)
  Z10 '(///,10f7.0)'           ZN2,WEITE2,QCU2,LW2,THETA2,KP2,LNEB2,LUEB2,HLTR2,BLTR2
  Z11 '(//,8f7.0,i7)'          BN2,BN2S,HK2,HO2,HCUO2,HZW2,HCUU2,HU2,SCHA2(i7)
  Z12 '(///,10f7.0)'           H32,H42,H52,H62,H72,MATSO,MATRO,MATSU,MATRU (10th field unused)
  Z13 '(//,11f7.0)'            BR32,BR42,BR52,BR62,HUSTAB,BUSTAB,AUSTAB (4 trailing fields unused)
  Z14 '(//,11f7.0)'            KAPPA3,KAPPA4,KAPPA5,KAPPA6,KAPPA7,HOSTAB,BOSTAB,AOSTAB,AN2,VN2 (1 unused)
  Z15 '(//,11f7.0)'            LUE,QRING,KAPRI,BETA,QSTAB,QSTABA (5 trailing fields unused)
  Z16 '(///,5f7.0)'            LUEA,QRINGA,KAPRIA,NK2,BK2
  Z17 '(///,2F14.0,5f7.0)'     JANTR,MGN,MG(0),MG(1),MG(2),MG(3),MG(4)
  Z18 '(///,11F7.0)'           MG(5)..MG(15)
  Z19 '(///,A10,4X,A10,4X)'    WP1(10),WP2(10)  (each followed by 4 skipped columns)

Each Zn's leading '/' characters mean "skip N-1 blank/header lines before
this one" - i.e. Z1's '(//,...)' skips 2 lines then reads the 3rd.

NOT PORTED: the ~800-line tail of the original LESEN (source 2055-3023)
that follows the data reads. That tail is entirely screen/printer report
generation (writes to units 3/4/50/51 with HP terminal escape codes -
the "cover sheet" summary display) with no new persisted computation -
it belongs conceptually to the not-yet-ported "output/reporting"
functional group, not to input parsing. 'lesen' here stops once the
machine-data record is fully populated and the winding design
(ISOSTD/ISOTR/ISOLFR) has run, which is everything downstream
calculation code actually depends on.

*** UNITS - CORRECTED (see CONVERSION_STATUS.md for the full story) ***
An earlier round of this port flagged a "millimetres vs metres" concern
after spotting FORTRAN print statements labelled ' MM' (e.g. source line
816, 'HOJOC1 =',F7.2,' MM'). That concern was WRONG, and is corrected
here: checking those same WRITE statements shows every one of them
multiplies the variable by 1000 immediately before printing it (e.g.
source line 388: `HOJOC1*1000.,HOJOC2*1000.,LI*1000.` - and the same
pattern for DA1*1000., DI1*1000. at every external-file-transfer call
site). That `*1000.` only happens at the display/output boundary, which
proves the INTERNAL representation throughout ASYN4.FOR - and therefore
throughout every field LESEN reads from EINASYN - is METRES, exactly
matching the metre-scale values every function in this port has been
tested with since round 1. No unit conversion is needed anywhere in
this codebase. If you obtain a real EINASYN file, its numeric fields
(DA1, DI1, L, geometry dimensions generally) will already be in metres
(e.g. "0.400" for a 400mm bore, not "400").
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field

from .state import MachineState
from .magnetics_solver import kzwert
from .insulation import isostd, isotr, isolfr

PI = 3.14159


# --------------------------------------------------------------------- #
# Fixed-width "card" parsing helpers, mirroring FORTRAN F/I/A edit
# descriptors: blank/short fields read as 0 (numeric) or space-padded
# (character), matching FORTRAN's list-directed blank-field behaviour.
# --------------------------------------------------------------------- #

def _field(line: str, start: int, width: int) -> str:
    return line[start:start + width] if len(line) >= start else ""


def _f(line: str, start: int, width: int) -> float:
    s = _field(line, start, width).strip()
    return float(s) if s else 0.0


def _i(line: str, start: int, width: int) -> int:
    s = _field(line, start, width).strip()
    return int(s) if s else 0


def _a(line: str, start: int, width: int) -> str:
    return _field(line, start, width)


class _CardReader:
    """Sequential reader over a EINASYN file's lines, with FORTRAN-style
    '/' skip-N-lines-then-read semantics built into each read_* call."""

    def __init__(self, lines: list):
        self.lines = lines
        self.pos = 0

    def _next_line(self, skip: int) -> str:
        self.pos += skip
        line = self.lines[self.pos] if self.pos < len(self.lines) else ""
        self.pos += 1
        return line

    def read_f(self, skip: int, n: int, width: int = 7) -> list:
        line = self._next_line(skip)
        return [_f(line, j * width, width) for j in range(n)]

    def read_mixed(self, skip: int, spec: list) -> list:
        """spec: list of ('f'|'i'|'a', width) tuples, read in order from
        one line."""
        line = self._next_line(skip)
        out = []
        pos = 0
        for kind, width in spec:
            if kind == "f":
                out.append(_f(line, pos, width))
            elif kind == "i":
                out.append(_i(line, pos, width))
            else:
                out.append(_a(line, pos, width))
            pos += width
        return out


def _fmt_field(kind: str, width: int, value, decimals: int = 4) -> str:
    if kind == "f":
        s = f"{float(value):{width}.{decimals}f}"
        if len(s) > width:
            # too many decimals to fit - drop decimals until it fits,
            # matching how a real EINASYN field would have been punched
            for d in range(decimals - 1, -1, -1):
                s = f"{float(value):{width}.{d}f}"
                if len(s) <= width:
                    break
        return s.rjust(width)[:width]
    if kind == "i":
        return f"{int(value)}".rjust(width)[:width]
    return str(value).ljust(width)[:width]


class _CardWriter:
    """
    Companion to _CardReader: writes an EINASYN-format file using the
    exact same skip-N-lines-then-write semantics, so a synthetic test
    file can never drift out of sync with the reader (which is exactly
    what caused a debugging detour during this port - hand-counting
    blank lines between cards is error-prone; this makes the skip
    counts shared code instead of duplicated-by-hand).
    """

    def __init__(self):
        self.lines = []

    def _pad_to(self, skip: int):
        for _ in range(skip):
            self.lines.append("")

    def write_f(self, skip: int, values, width: int = 7, decimals: int = 4):
        self._pad_to(skip)
        self.lines.append("".join(_fmt_field("f", width, v, decimals) for v in values))

    def write_mixed(self, skip: int, fields):
        """fields: list of (kind, width, value[, decimals]) - decimals
        only used for kind='f', defaults to a sensible value for the
        given width if omitted."""
        self._pad_to(skip)
        parts = []
        for field in fields:
            kind, width, value = field[0], field[1], field[2]
            decimals = field[3] if len(field) > 3 else max(width - 3, 0)
            parts.append(_fmt_field(kind, width, value, decimals))
        self.lines.append("".join(parts))

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"

    def write(self, path: str):
        with open(path, "w") as f:
            f.write(self.text())


def lesen(filepath: str, irenr: int, st: MachineState,
          run_isostd: bool = True, run_isolfr: bool = True,
          notbetrieb: bool = False, verbose: bool = True):
    """
    Reads an EINASYN-format input file and populates st.daten, st.art.rkz,
    and st.materi, then runs the stator winding design (ISOSTD for
    form-coil, ISOTR for random-wound) and, for slip-ring rotors (ILFR==1),
    the rotor winding design (ISOLFR).

    run_isostd/run_isolfr mirror the FORTRAN's "... ISO durchlaufen?
    (J/n)" prompts (only asked when QCU1/QCU2 already != 0 in the file) -
    default True matches not answering 'n'. notbetrieb mirrors the
    "SOLL NOTBETRIEB GERECHNET WERDEN?" prompt for 6-phase machines
    (default False = 'N').

    Returns a dict with every value the FORTRAN passed back as a
    subroutine argument: name, reke, kennwort, datum, vanr, fbnr, wanr,
    monr, rkz, iflag1, ikopp, ald1s, ald1a, lges1, drufi, nk2, bk2, na,
    da, dta, wp1, wp2, iokay, isoart, unx, unx2, lueb2, lneb2, kp2,
    bltr2, hltr2, ihk, plus the ISOSTD/ISOTR/ISOLFR result dicts under
    'isostd_result'/'isotr_result'/'isolfr_result' (whichever ran).
    """
    with open(filepath, "r") as fh:
        lines = fh.read().splitlines()

    d = st.daten
    r = _CardReader(lines)

    # Z1
    name, datum, kennwort = r.read_mixed(2, [("a", 16), ("a", 8), ("a", 25)])
    # Z2
    reke, vanr, fbnr, wanr, monr = r.read_mixed(
        2, [("a", 27), ("a", 17), ("a", 15), ("i", 5), ("a", 12)])
    # Z3
    (rkz,) = r.read_mixed(2, [("a", 11)])
    rkz = (rkz + " " * 12)[:12]  # pad to 12 like the FORTRAN CHARACTER*12

    # ILFR from RKZ(1:1)
    ilfr = int(rkz[0]) if rkz[0].strip() else 0
    if ilfr > 4:
        ilfr -= 4
        iflag1 = 1
    else:
        iflag1 = 0
    d.ilfr = ilfr

    # Z4
    (d.fn, d.kfe, d.kj1, d.kj2, d.p, d.da1, d.l, d.zk, d.bk, d.pn,
     d.kz) = r.read_f(3, 11)

    # Z5
    line5 = r.read_mixed(2, [("i", 7)] + [("f", 7)] * 10)
    (ishalt, un, d.prbg0, d.v10, d.cfej, d.cfez, d.schr, d.deltag, d.ku1,
     d.xnetz, d.kt) = line5

    # Z6
    line6 = r.read_mixed(2, [("f", 7)] * 9 + [("i", 7), ("f", 7)])
    (d.m1, d.mzone1, d.di1, d.n1, zn1_raw, d.weite1, d.qcu1, lw1, d.theta1,
     ikopp, unx) = line6
    if d.m1 < 0.0001:
        d.m1 = 3.0
    d.zn1 = round(zn1_raw * d.n1 / d.m1) * d.m1 / d.n1

    # Z7
    line7 = r.read_mixed(2, [("f", 7), ("f", 7), ("i", 7), ("a", 7),
                              ("f", 7), ("f", 7), ("f", 7), ("f", 7)])
    d.bs1, d.hs1, izgf_i, nutnr1, nutf1, na, da, dta = line7
    d.izgf = izgf_i
    d.diso1 = d.l1 = d.l2 = d.lstr = d.r = d.rd = d.b = d.dstr = 0.0

    # Z8 (format 2231, 11f7.0)
    (d.bn1, d.bn1s, d.hk1, d.ho1, d.hcuo1, d.hzw1, d.hcuu1, d.hu1,
     d.lamst1, nutf2, _spare8) = r.read_f(2, 11)
    hn1 = d.hs1 + d.hk1 + d.ho1 + d.hcuo1 + d.hzw1 + d.hcuu1 + d.hu1
    if rkz[8] == "E":
        d.hcuo1 = (d.hcuo1 + d.hcuu1) / 2.0
        d.hcuu1 = d.hcuo1

        # Z9
    line9 = r.read_mixed(2, [("f", 7)] * 10 + [("a", 5)])
    (d.di2, d.n2, d.bs2, d.hs2, d.mtl, d.ntl, d.utl, d.azweig, d.btl,
    d.htl, nutnr2) = line9

    if d.hs2 < 1.0e-8:
        d.hs2 = 1.0e-8

    # Z10
    (zn2_raw, d.weite2, d.qcu2, lw2, d.theta2, kp2, lneb2, lueb2, hltr2,
     bltr2) = r.read_f(3, 10)

    # Safely calculate zn2 to prevent ZeroDivisionError if n2 or m1 is missing/zero
    if d.n2 > 0.0 and d.m1 > 0.0:
        d.zn2 = round(zn2_raw * d.n2 / d.m1) * d.m1 / d.n2
    else:
        d.zn2 = zn2_raw

    # Z11
    line11 = r.read_mixed(2, [("f", 7)] * 8 + [("i", 7)])
    (d.bn2, d.bn2s, d.hk2, d.ho2, d.hcuo2, d.hzw2, d.hcuu2, d.hu2,
     scha2) = line11

    # Z12
    (d.h32, d.h42, d.h52, d.h62, d.h72, st.materi.matso, st.materi.matro,
     st.materi.matsu, st.materi.matru, _spare12) = r.read_f(3, 10)
    d.h22 = 0.0

    # Z13
    (d.br32, d.br42, d.br52, d.br62, hustab, bustab, austab, *_spare13
     ) = r.read_f(2, 11)
    d.br12 = 0.0
    d.br22 = 0.0
    d.bstr52 = d.br52
    d.bstr62 = d.br62

    # Z14 (format 2231)
    (d.kappa3, d.kappa4, d.kappa5, d.kappa6, d.kappa7, hostab, bostab,
     aostab, an2, vn2, _spare14) = r.read_f(2, 11)
    d.kappa2 = 0.0

    # Z15 (format 2231)
    (d.lue, d.qring, d.kapri, d.beta, d.qstab, d.qstaba, *_spare15
     ) = r.read_f(2, 11)
    d.k42 = 0.0
    d.k4a2 = 0.0

    # Z16
    luea, d.qringa, d.kapria, nk2, bk2 = r.read_f(3, 5)

    # Z17
    line17 = r.read_mixed(3, [("f", 14), ("f", 14)] + [("f", 7)] * 5)
    d.jantr, d.mgn = line17[0], line17[1]
    for j in range(0, 5):
        d.mg[j] = line17[2 + j]

    # Z18
    mg_rest = r.read_f(3, 11)
    for j in range(5, 16):
        d.mg[j] = mg_rest[j - 5]

    # Z19
    wp1, wp2 = r.read_mixed(3, [("a", 10), ("a", 4), ("a", 10), ("a", 4)])[0::2]

    if ilfr != 1:
        wp2 = " " * 10

    ihk = 0
    if ilfr != 1:
        d.qstab = (austab + bustab) / 2.0 * hustab
        if d.h52 + d.h62 > 0.0:
            qstabo = PI / 2.0 * d.h52 ** 2
            qstabu = PI / 2.0 * d.h72 ** 2
            d.qstab = d.qstab - d.h52 * austab + qstabo
            d.qstab = d.qstab - d.h72 * bustab + qstabu
        if austab * bustab * hustab == 0.0:
            d.qstab = 0.0
        if d.qstab == 0.0:
            austab = d.br52
            bustab = d.br62
            hustab = d.h52 + d.h62 + d.h72

        if aostab == hostab:
            d.qstaba = PI / 4.0 * aostab ** 2
        else:
            d.qstaba = (aostab + bostab) / 2.0 * hostab
        if aostab * bostab * hostab == 0.0:
            d.qstaba = 0.0
        if d.qstaba == 0.0:
            aostab = d.br32
            bostab = d.br32
            hostab = d.h32

        ihk = 0 if d.h32 == d.br32 else 1
    if ilfr == 1:
        ihk = 0

    # ---- default values ----
    if d.fn < 0.0001:
        d.fn = 50.0
    if d.kfe < 0.0001:
        d.kfe = 0.92
    if d.kj1 < 0.0001:
        d.kj1 = 1.0
    if d.bk < 0.0001:
        d.bk = 10.0
    if d.zk < 0.0001:
        d.bk = 0.0
    if nk2 < 0.0001:
        nk2 = d.zk
    if bk2 < 0.0001:
        bk2 = d.bk
    if nk2 < 0.0001:
        bk2 = 0.0
    if d.v10 < 0.0001:
        d.v10 = 2.00
    if d.mzone1 < 0.0001:
        d.mzone1 = 6.0
    if d.bn1s < 0.0001:
        brz1 = PI * (d.di1 + 2.0 * d.hs1 + 2.0 * d.hk1) / d.n1 - d.bn1
        d.bn1s = (PI * (d.di1 + 2.0 * d.hs1 + 2.0 * d.hk1 + 2.0 * d.ho1 +
                         2.0 * d.hcuo1 + 2.0 * d.hzw1 + 2.0 * d.hcuu1 +
                         2.0 * d.hu1) / d.n1 - brz1)
    if d.utl < 0.0001:
        d.utl = 1.0
    if ilfr == 1 and d.bn2s < 0.00001:
        brz2 = (PI * (d.di1 - 2.0 * d.deltag - 2.0 * d.hs2 - 2.0 * d.hk2) /
                d.n2 - d.bn2)
        d.bn2s = (PI * (d.di1 - 2.0 * d.deltag - 2.0 * d.hs2 - 2.0 * d.hk2 -
                         2.0 * d.ho2 - 2.0 * d.hcuo2 - 2.0 * d.hzw2 -
                         2.0 * d.hcuu2 - 2.0 * d.hu2) / d.n2 - brz2)
    if unx == 0.0:
        unx = un
    isoart = 0
    if unx < 0.0:
        isoart = 1
        unx = -unx
    i2l1gs = 1

    ihg = round(d.zk - nk2)
    if ihg % 2 != 0:
        raise ValueError("LESEN: Differenz der Kuehlkanaele ist ungerade")

    if d.kz == 0.0:
        d.kz = kzwert(reke, d.p)

    if ishalt in (1, 3):
        d.u1 = un
    else:
        if round(d.m1) == 6:
            d.u1 = un / 1.73205
        else:
            d.u1 = un / (2.0 * math.sin(PI / d.m1))

    st.art.rkz = rkz

    # ---- stator winding design ----
    isostd_result = isotr_result = isolfr_result = None
    chtyp = reke[0:22]
    chauf = fbnr[0:12]

    if d.izgf == 1:
        if run_isostd:
            isostd_result = isostd(
                reke, d.l, d.di1, d.n1, d.weite1, isoart, unx, d.qcu1,
                d.bs1, d.bn1, d.bn1s, d.hs1, d.hk1, d.ho1, d.hcuo1,
                d.hzw1, d.hcuu1, d.hu1, d.zn1, d.azweig, lw1, d.mtl,
                d.ntl, d.utl, d.btl, d.htl, d.diso1, d.l1, d.l2, d.lstr,
                d.r, d.rd, d.b, d.dstr, d.m1, d.p, d.da1, d.zk, d.bk,
                d.schr, st, verbose=verbose)
            for key in ("bs1", "bn1", "bn1s", "hs1", "hk1", "ho1",
                        "hcuo1", "hzw1", "hcuu1", "hu1", "qcu1", "mtl",
                        "ntl", "btl", "htl", "l1", "l2", "r", "rd", "b",
                        "diso1"):
                setattr(d, key, isostd_result[key])
            lges1 = isostd_result["lges1"]
            ald1s = isostd_result["ald1s"]
            ald1a = isostd_result["ald1a"]
            drufi = isostd_result["drufi"]
        else:
            lges1 = ald1s = ald1a = drufi = 0.0
    else:
        if d.qcu1 != 0.0:
            # Caller already specified the copper area directly - ISOTR's
            # QCU1!=0 branch just back-derives MTL from it, no interactive
            # decision needed.
            isotr_result = isotr(
                rkz, d.bs1, d.bn1, d.bn1s, d.hk1, d.ho1, d.hcuo1, d.hzw1,
                d.hcuu1, d.hu1, d.zn1, d.azweig, d.n1, d.p, d.mzone1,
                d.btl, d.htl, d.weite1, d.qcu1, tpar=None, verbose=verbose)
        else:
            # Mirrors the FORTRAN's "ANZAHL PARALLELER TEILLEITER = ?"
            # prompt, which run_asyn4() (and any other non-interactive
            # caller) cannot answer interactively. First get the
            # suggestion table, then auto-pick the first entry whose
            # theoretical (blank-wire) fill factor is <= 0.68 - i.e. the
            # smallest TPAR that still physically fits, mirroring what an
            # operator would sensibly choose absent other guidance. If
            # ISOTR finds no suggestion at all (BTL/HTL don't fit the
            # slot at any TPAR), d.qcu1 is left at 0 and a warning is
            # printed - the caller should treat that as a data problem to
            # fix in the input file, not something this port can silently
            # paper over.
            probe = isotr(
                rkz, d.bs1, d.bn1, d.bn1s, d.hk1, d.ho1, d.hcuo1, d.hzw1,
                d.hcuu1, d.hu1, d.zn1, d.azweig, d.n1, d.p, d.mzone1,
                d.btl, d.htl, d.weite1, 0.0, tpar=None, verbose=False)
            chosen_tpar = None
            for tp, kcuth, kcup, tp2, kcuth2, kcup2 in probe["suggestions"]:
                if kcup <= 0.68:
                    chosen_tpar = tp
                    break
                if kcup2 <= 0.68:
                    chosen_tpar = tp2
                    break
            if chosen_tpar is None:
                if verbose:
                    print("\nISOTR: KEIN TPAR GEFUNDEN, DER DIE NUT AUSFUELLT "
                          "- QCU1 BLEIBT 0, DATENSATZ PRUEFEN")
                isotr_result = probe
            else:
                isotr_result = isotr(
                    rkz, d.bs1, d.bn1, d.bn1s, d.hk1, d.ho1, d.hcuo1,
                    d.hzw1, d.hcuu1, d.hu1, d.zn1, d.azweig, d.n1, d.p,
                    d.mzone1, d.btl, d.htl, d.weite1, 0.0,
                    tpar=chosen_tpar, verbose=verbose)
        if isotr_result["mtl"] is not None:
            d.mtl = isotr_result["mtl"]
            d.qcu1 = isotr_result["qcu1"]
        lges1 = ald1s = ald1a = drufi = 0.0

    # ---- slip-ring rotor winding design ----
    u20 = 0.0
    unx2 = 0.0
    if ilfr == 1:
        proceed = run_isolfr or d.qcu2 <= 0.0
        if proceed:
            u20 = d.u1 * d.zn2 * d.n2 / (d.zn1 * d.n1)
            if scha2 == 0:
                u20 = u20 * 1.73
            if lw2 == 0.0:
                dm2 = (d.di1 - 2.0 * (d.deltag + d.hs2 + d.hk2) - d.ho2 -
                       d.hcuo2 - d.hzw2 - d.hcuu2 - d.hu2)
                tnm2 = PI * dm2 / d.n2
                lw2 = (150.0 + d.n2 / (2.0 * d.p) * tnm2 /
                       math.sqrt(1.0 - ((d.bn2 + 3.0) / tnm2) ** 2))
            isolfr_result = isolfr(
                u20, d.zn2, d.n2, d.l + lw2, kp2, d.weite2, d.bs2, d.bn2,
                d.hs2, d.hk2, d.ho2, d.hcuo2, d.hzw2, d.hcuu2, d.hu2, st)
            lueb2, lneb2 = isolfr_result["tpar"], isolfr_result["ntl"]
            d.qcu2 = isolfr_result["qcu"]
            bltr2, hltr2 = isolfr_result["btl"], isolfr_result["htl"]
            unx2 = isolfr_result["unx2"]

    if round(d.m1) == 6 and notbetrieb:
        rkz = rkz[:10] + "X" + rkz[11:]
        st.art.rkz = rkz
        d.m1 = 3.0
        d.mzone1 = 6.0

    if ilfr == 1:
        d.h22 = d.h32 = d.h42 = d.h52 = d.h62 = d.h72 = 0.0
        d.br12 = d.br22 = d.br32 = d.br42 = d.br52 = d.bstr52 = 0.0
        d.br62 = d.bstr62 = 0.0

    if ilfr >= 2 and d.qring == 0.0:
        qsthg = (d.h52 + d.h62 + d.h72) * (d.br52 + d.br62) / 2.0 + d.h32 * d.br32
        qrhg = qsthg * d.n2 / (2.0 * d.p * PI)
        d.qring = math.trunc(qrhg / 100.0) * 100.0 + 100.0

    # ---- MAGN. WIRKS. STREUSCHLITZ BEI SONDERSTREUNUTFORMEN (source
    # ~2676-2683): special slot-shape codes override HS2 with AN2.
    if an2 > 0.0001 and round(nutf2) in (2156, 2157, 2158, 2159):
        d.hs2 = an2

    # ---- UMWANDLUNG IN SI-EINHEITEN (source 2689-2774) ----
    # *** CRITICAL: this is where the mm-vs-metres question this port
    # earlier got wrong (twice - see CONVERSION_STATUS.md) is settled
    # for good. K=1000, M=1e6, defined at the very top of the original
    # LESEN. Every length field read from EINASYN is in MILLIMETRES and
    # is divided by K here; every area field (QCU1/QCU2/QSTAB/QSTABA/
    # QRING/QRINGA) is in mm^2 and divided by M; every conductivity
    # field (KAPPA2-7, KAPRI, KAPRIA) is in S/mm^2-equivalent and
    # multiplied by M. This conversion is what makes the internal
    # representation - and therefore every ported function in this
    # package, all tested with metre-scale values - correct: EINASYN's
    # raw fields are millimetres, and this block converts them to
    # metres before anything else in the program touches them. ***
    k = 1.0e3
    m = 1.0e6
    d.da1 /= k
    d.l /= k
    d.bk /= k
    d.pn *= k
    d.prbg0 *= k
    d.deltag /= k
    d.di1 /= k
    d.qcu1 /= m
    lw1 /= k
    d.lm1 = lw1 + d.l
    d.bs1 /= k
    d.hs1 /= k
    d.diso1 /= k
    d.l1 /= k
    d.l2 /= k
    d.lstr /= k
    d.r /= k
    d.rd /= k
    d.b /= k
    d.bn1 /= k
    d.bn1s /= k
    d.hk1 /= k
    d.ho1 /= k
    d.hcuo1 /= k
    d.hzw1 /= k
    d.hcuu1 /= k
    d.hu1 /= k
    d.dstr /= k
    d.di2 /= k
    d.bs2 /= k
    d.hs2 /= k
    d.btl /= k
    d.htl /= k
    d.qcu2 /= m
    lw2 /= k
    d.lm2 = lw2 + d.l
    d.bn2 /= k
    d.bn2s /= k
    d.hk2 /= k
    d.ho2 /= k
    d.hcuo2 /= k
    d.hzw2 /= k
    d.hcuu2 /= k
    d.hu2 /= k
    d.h22 /= k
    d.h32 /= k
    d.h42 /= k
    d.h52 /= k
    d.h62 /= k
    d.h72 /= k
    d.br12 /= k
    d.br22 /= k
    d.br32 /= k
    d.br42 /= k
    d.br52 /= k
    d.bstr52 /= k
    d.br62 /= k
    d.bstr62 /= k
    d.kappa2 *= m
    d.kappa3 *= m
    d.kappa4 *= m
    d.kappa5 *= m
    d.kappa6 *= m
    d.kappa7 *= m
    d.lue /= k
    d.lstab = 2.0 * d.lue + d.l
    d.qring /= m
    d.kapri *= m
    d.beta *= PI / 180.0
    d.qstab /= m
    d.qstaba /= m
    luea /= k
    d.lstaba = 2.0 * luea + d.l
    d.qringa /= m
    d.kapria *= m
    d.kz /= 1.0e11
    d.kt /= 1.0e3
    ald1s /= k
    ald1a /= k
    lges1 /= k
    drufi /= k
    bk2 /= k
    da /= k
    dta /= k
    bltr2 /= k
    hltr2 /= k

    # ---- UMWANDLUNG VON DOPPELNUT IN HAMMERKOPF (source 2780-2799) ----
    if ihk == 1:
        d.kappa4 = d.kappa6 if d.kappa3 == d.kappa6 else 0.0
        rkz = "2" + rkz[1:]
        st.art.rkz = rkz
        ilfr = 2
        d.ilfr = 2
        d.h22 = d.h32
        d.h32 = 0.0
        d.br12 = d.br32
        d.br22 = d.br32
        d.br32 = 0.0
        d.kappa2 = d.kappa3
        d.kappa3 = 0.0
        if d.qstab + d.qstaba > 0.0:
            d.qstab = d.qstab + d.qstaba + d.h42 * d.br42
            d.qstaba = 0.0

    # ---- DRUCKGUSS-LAEUFER (source 2805-2810): die-cast rotor with no
    # ring-overhang length configured downgrades from double-cage
    # common-ring (ILFR=3) to single-cage (ILFR=2). ----
    if ilfr == 3 and d.lue == 0.0:
        ilfr = 2
        d.ilfr = 2
        rkz = "2" + rkz[1:]
        st.art.rkz = rkz
        d.qstab = 0.0
        d.qstaba = 0.0

    return {
        "name": name, "reke": reke, "kennwort": kennwort, "datum": datum,
        "vanr": vanr, "fbnr": fbnr, "wanr": wanr, "monr": monr,
        "rkz": st.art.rkz, "iflag1": iflag1, "ikopp": ikopp,
        "ald1s": ald1s, "ald1a": ald1a, "lges1": lges1, "drufi": drufi,
        "nk2": nk2, "bk2": bk2, "na": na, "da": da, "dta": dta,
        "wp1": wp1, "wp2": wp2, "isoart": isoart, "unx": unx,
        "unx2": unx2, "lueb2": lueb2 if ilfr == 1 else 0.0,
        "lneb2": lneb2 if ilfr == 1 else 0.0, "kp2": kp2 if ilfr == 1 else 0.0,
        "bltr2": bltr2 if ilfr == 1 else 0.0,
        "hltr2": hltr2 if ilfr == 1 else 0.0, "ihk": ihk,
        "ishalt": ishalt, "i2l1gs": i2l1gs,
        "isostd_result": isostd_result, "isotr_result": isotr_result,
        "isolfr_result": isolfr_result,
    }
