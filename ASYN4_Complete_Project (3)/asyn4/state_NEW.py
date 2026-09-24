"""
Replacement for the FORTRAN COMMON blocks used throughout ASYN4.FOR.

The original program shares state between ~100 subroutines purely through
COMMON blocks (no argument passing). Rather than trying to force that into
Python globals, every COMMON block becomes one dataclass here, and all of
them are bundled into a single `MachineState` object that gets threaded
through the ported functions explicitly. This keeps the translation
traceable back to the FORTRAN (block names and variable names are kept
identical, just lower-cased) while making the data flow visible instead of
implicit.

INDEXING CONVENTION
--------------------
FORTRAN arrays in this program are 1-indexed (e.g. NY(25), MG(0:15)).
To avoid introducing off-by-one bugs while porting ~11,000 more lines that
constantly use these arrays with their original FORTRAN indices, arrays are
allocated ONE ELEMENT LARGER than their FORTRAN declaration and index 0 is
simply left unused (except where FORTRAN itself declares a 0 lower bound,
e.g. MG(0:15), which is a natural 0-based array of 16 elements).

Example: FORTRAN `NY(25)` -> `np.zeros(26)`, valid Python indices 1..25.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


def _f(n: int) -> np.ndarray:
    """Helper: float64 zeros array with FORTRAN 1-based indexing (n usable slots)."""
    return np.zeros(n + 1, dtype=float)


@dataclass
class Brulo:
    """COMMON /BRULO/ - fractional-slot ("Bruchlochwicklung") winding data."""
    bruhi1: float = 0.0
    bruhi2: float = 0.0
    bruhi3: float = 0.0
    bruhi4: float = 0.0
    bruhi5: float = 0.0
    ibruch: int = 0
    zv: float = 0.0


@dataclass
class Vorzei:
    """COMMON /VORZEI/ - sign/system conventions."""
    vzxsiz: float = 1.0
    ksys: float = 1.0
    bs1di1: float = 0.0


@dataclass
class Keil:
    """COMMON /KEIL/ - slot wedge saturation data."""
    xskeil: float = 0.0
    teta0: float = 0.0


@dataclass
class CBSort:
    """COMMON /CBSORT/ - magnetization curve table used by HEISEN."""
    h: np.ndarray = field(default_factory=lambda: _f(45))  # H(0:45) -> 46 slots, 0 used
    my0: float = 1.2566e-6


@dataclass
class NyInfo:
    """COMMON /NYINFO/ - harmonic order (Oberfeld) bookkeeping, 25 slots each."""
    sd1: float = 0.0
    sd1r: float = 0.0
    nyanz: int = 0
    ny: np.ndarray = field(default_factory=lambda: _f(25))
    sd2: np.ndarray = field(default_factory=lambda: _f(25))
    uer: np.ndarray = field(default_factory=lambda: _f(25))
    uei: np.ndarray = field(default_factory=lambda: _f(25))
    sf: np.ndarray = field(default_factory=lambda: _f(25))
    x1hbez: np.ndarray = field(default_factory=lambda: _f(25))
    rrny: np.ndarray = field(default_factory=lambda: _f(25))
    rrnya: np.ndarray = field(default_factory=lambda: _f(25))
    xsrny: np.ndarray = field(default_factory=lambda: _f(25))
    xsrnya: np.ndarray = field(default_factory=lambda: _f(25))


@dataclass
class Maxi:
    smaxi: np.ndarray = field(default_factory=lambda: _f(25))
    mmaxi: np.ndarray = field(default_factory=lambda: _f(25))


@dataclass
class Kennli:
    """COMMON /KENNLI/ - magnetization/leakage characteristic curve pools."""
    eanz: int = 0
    mkpool: np.ndarray = field(default_factory=lambda: np.zeros((51, 19)))   # (50,18)
    svpool: np.ndarray = field(default_factory=lambda: np.zeros((26, 6)))    # (25,5)
    inut: np.ndarray = field(default_factory=lambda: _f(20))
    lamda: np.ndarray = field(default_factory=lambda: _f(20))


@dataclass
class Esb:
    """COMMON /ESB/ - equivalent-circuit (Ersatzschaltbild) geometry."""
    skipp: float = 0.0
    r1: float = 0.0
    xsti1: float = 0.0
    xn1: float = 0.0
    xsd1: float = 0.0
    xsd1r: float = 0.0
    xsch: float = 0.0
    x1h: float = 0.0
    xsd2: float = 0.0
    xring: float = 0.0
    xn2: float = 0.0
    rk2: float = 0.0
    rr2: float = 0.0
    rs2: float = 0.0
    xsteg: float = 0.0
    xng: float = 0.0
    xna: float = 0.0
    ra2: float = 0.0
    rra2: float = 0.0
    xringa: float = 0.0


@dataclass
class Art:
    rkz: str = " " * 12


@dataclass
class Laerm:
    """COMMON /LAERM/ - noise/loss-related geometry."""
    nn: float = 0.0
    in_: float = 0.0
    cosfin: float = 0.0
    brz1mi: float = 0.0
    brz1ma: float = 0.0
    dring: float = 0.0
    bs2in: float = 0.0
    lamn2: float = 0.0
    in2: float = 0.0
    imue: float = 0.0
    bp: float = 0.0
    alpha: float = 0.0
    vluft: float = 0.0
    vzahn1: float = 0.0
    vzahn2: float = 0.0
    vjoch1: float = 0.0
    vjoch2: float = 0.0
    deltap: float = 0.0


@dataclass
class Mnline:
    """COMMON /MNLINE/ - torque-speed (n,M) load curve, 400 points, 5 columns."""
    nanz: int = 0
    npu: np.ndarray = field(default_factory=lambda: _f(400))
    mpu: np.ndarray = field(default_factory=lambda: _f(400))
    i1pu: np.ndarray = field(default_factory=lambda: _f(400))
    cosph: np.ndarray = field(default_factory=lambda: _f(400))
    g: np.ndarray = field(default_factory=lambda: np.zeros((401, 6)))
    upu: np.ndarray = field(default_factory=lambda: _f(400))


@dataclass
class Mnlred(Mnline):
    """COMMON /MNLRED/ - same layout as MNLINE, for the reduced-voltage run."""
    pass


@dataclass
class Daten:
    """
    COMMON /DATEN/ - the big one: main machine design data, ~100 scalars
    plus the MG(0:15) counter-torque table. Populated by LESEN.
    """
    ilfr: int = 0
    fn: float = 0.0
    kfe: float = 0.0
    kj1: float = 0.0
    kj2: float = 0.0
    p: float = 0.0
    da1: float = 0.0
    l: float = 0.0
    zk: float = 0.0
    bk: float = 0.0
    pn: float = 0.0
    un: float = 0.0
    u1: float = 0.0
    prbg0: float = 0.0
    v10: float = 0.0
    cfej: float = 0.0
    cfez: float = 0.0
    schr: float = 0.0
    deltag: float = 0.0
    ku1: float = 0.0
    xnetz: float = 0.0
    m1: float = 0.0
    mzone1: float = 0.0
    di1: float = 0.0
    n1: float = 0.0
    zn1: float = 0.0
    weite1: float = 0.0
    qcu1: float = 0.0
    lm1: float = 0.0
    theta1: float = 0.0
    bs1: float = 0.0
    hs1: float = 0.0
    izgf: int = 0
    diso1: float = 0.0
    l1: float = 0.0
    l2: float = 0.0
    lstr: float = 0.0
    r: float = 0.0
    rd: float = 0.0
    b: float = 0.0
    bn1: float = 0.0
    bn1s: float = 0.0
    hk1: float = 0.0
    ho1: float = 0.0
    hcuo1: float = 0.0
    hzw1: float = 0.0
    hcuu1: float = 0.0
    hu1: float = 0.0
    lamst1: float = 0.0
    dstr: float = 0.0
    di2: float = 0.0
    n2: float = 0.0
    bs2: float = 0.0
    hs2: float = 0.0
    mtl: float = 0.0
    ntl: float = 0.0
    utl: float = 0.0
    azweig: float = 0.0
    btl: float = 0.0
    htl: float = 0.0
    zn2: float = 0.0
    weite2: float = 0.0
    qcu2: float = 0.0
    lm2: float = 0.0
    theta2: float = 0.0
    bn2: float = 0.0
    bn2s: float = 0.0
    hk2: float = 0.0
    ho2: float = 0.0
    hcuo2: float = 0.0
    hzw2: float = 0.0
    hcuu2: float = 0.0
    hu2: float = 0.0
    h22: float = 0.0
    h32: float = 0.0
    h42: float = 0.0
    h52: float = 0.0
    h62: float = 0.0
    h72: float = 0.0
    br12: float = 0.0
    br22: float = 0.0
    br32: float = 0.0
    br42: float = 0.0
    br52: float = 0.0
    bstr52: float = 0.0
    br62: float = 0.0
    bstr62: float = 0.0
    kappa2: float = 0.0
    kappa3: float = 0.0
    kappa4: float = 0.0
    kappa5: float = 0.0
    kappa6: float = 0.0
    kappa7: float = 0.0
    lstab: float = 0.0
    qring: float = 0.0
    kapri: float = 0.0
    beta: float = 0.0
    qstab: float = 0.0
    qstaba: float = 0.0
    k41: float = 0.0
    k42: float = 0.0
    k4a2: float = 0.0
    k4r2: float = 0.0
    k4ar2: float = 0.0
    lstaba: float = 0.0
    qringa: float = 0.0
    kapria: float = 0.0
    jantr: float = 0.0
    mgn: float = 0.0
    mg: np.ndarray = field(default_factory=lambda: np.zeros(16))  # MG(0:15), 0-based
    kz: float = 0.0
    kt: float = 0.0
    lue: float = 0.0
    k41: float = 0.0
    k42: float = 0.0
    k4a2: float = 0.0
    k4r2: float = 0.0
    k4ar2: float = 0.0


@dataclass
class CPaket:
    """COMMON /CPAKET/ - lamination-stack (Blechpaket) split, 23 slots each."""
    paket1: np.ndarray = field(default_factory=lambda: _f(23))
    paket2: np.ndarray = field(default_factory=lambda: _f(23))


@dataclass
class Ctest:
    """COMMON /CTEST/ - flag for the reduced ("abgespeckt") insulation system."""
    itest: int = 0


@dataclass
class Spkopf:
    """COMMON /SPKOPF/ - coil-head clearance dimensions, set by NARISS."""
    bkopf: float = 0.0
    hkopf: float = 0.0


@dataclass
class Ck211s:
    """COMMON /CK211S/ - insulation clearance breakdown, set by NARISS,
    consumed by the (not yet ported) K211-family drawing/report routines."""
    hwk: float = 0.0
    bwk: float = 0.0
    del1: float = 0.0
    del2: float = 0.0
    del3: float = 0.0
    del4: float = 0.0
    del5: float = 0.0
    del6: float = 0.0
    del7: float = 0.0
    del8: float = 0.0
    del9: float = 0.0
    del10: float = 0.0
    del11: float = 0.0
    swk: float = 0.0
    ausldg: float = 0.0


@dataclass
class Ck211l:
    """COMMON /CK211L/ - rotor-side insulation clearance breakdown, set
    by NAR2, consumed by the (not yet ported) K211-family drawing/report
    routines - the slip-ring-rotor counterpart to Ck211s."""
    del1l: float = 0.0
    del2l: float = 0.0
    del3l: float = 0.0
    del4l: float = 0.0
    del5l: float = 0.0
    del6l: float = 0.0
    del7l: float = 0.0
    del8l: float = 0.0
    del9l: float = 0.0
    del10l: float = 0.0
    del11l: float = 0.0
    delb: float = 0.0


@dataclass
class Materi:
    """COMMON /MATERI/ - rotor bar/ring material codes, set by LESEN."""
    matro: float = 0.0
    matru: float = 0.0
    matso: float = 0.0
    matsu: float = 0.0


@dataclass
class Dp:
    """COMMON /DP/ - synchronous-point damping factor and reference X1H."""
    cdp: complex = 0j
    x1hp0: float = 0.0


@dataclass
class Svorga:
    """COMMON /SVORGA/ - pre-generated (unsorted) slip value list."""
    sanz: int = 0
    svor: np.ndarray = field(default_factory=lambda: _f(400))


@dataclass
class Xs2sl:
    """COMMON /XS2SL/ - slip-ring rotor slot-leakage reactance split."""
    xsn2: float = 0.0
    xss2: float = 0.0


@dataclass
class MachineState:
    """Bundle of every COMMON block, passed around instead of module globals."""
    brulo: Brulo = field(default_factory=Brulo)
    vorzei: Vorzei = field(default_factory=Vorzei)
    keil: Keil = field(default_factory=Keil)
    cbsort: CBSort = field(default_factory=CBSort)
    nyinfo: NyInfo = field(default_factory=NyInfo)
    maxi: Maxi = field(default_factory=Maxi)
    kennli: Kennli = field(default_factory=Kennli)
    esb: Esb = field(default_factory=Esb)
    art: Art = field(default_factory=Art)
    laerm: Laerm = field(default_factory=Laerm)
    mnline: Mnline = field(default_factory=Mnline)
    mnlred: Mnlred = field(default_factory=Mnlred)
    daten: Daten = field(default_factory=Daten)
    cpaket: CPaket = field(default_factory=CPaket)
    xs2sl: Xs2sl = field(default_factory=Xs2sl)
    dp: Dp = field(default_factory=Dp)
    svorga: Svorga = field(default_factory=Svorga)
    ctest: Ctest = field(default_factory=Ctest)
    spkopf: Spkopf = field(default_factory=Spkopf)
    ck211s: Ck211s = field(default_factory=Ck211s)
    ck211l: Ck211l = field(default_factory=Ck211l)
    materi: Materi = field(default_factory=Materi)
