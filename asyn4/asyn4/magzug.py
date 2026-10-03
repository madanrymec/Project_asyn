"""
MAGZUG / MAGZU1 (source 7331-7394 of ASYN4.FOR) - axial magnetic pull
calculation ("Axialer magnetischer Zug"), source: report by Dr. Meister.
Small, self-contained, genuine physics - kept in its own module since
it doesn't naturally belong to any of the existing functional-group
files and is only needed by the main driver's call into KDRUCK.
"""
from __future__ import annotations
import math

PI = 3.14159
MY0 = 1.2566e-6


def magzu1(gam: float) -> float:
    """Helper for MAGZUG: finds A such that a monotonic function of A
    first reaches GAM, scanning A = 0.01, 0.02, ..., 1.99."""
    a = 0.0
    for j in range(1, 200):
        a = 0.01 * j
        wert = (2.0 / PI * (3.0 * math.sqrt(1.0 + a) / (2.0 - a) +
                math.log((math.sqrt(1.0 + a) - 1.0) / math.sqrt(a))) - gam)
        if wert >= 0.0:
            break
    return a


def magzug(di1: float, delta1: float, delta2: float, blm: float,
           alpha: float, nk2: float, bk: float, li: float, p: float,
           mn: float, schr: float):
    """
    Axial magnetic pull as a function of axial rotor displacement X
    (0.001m steps, 8 points), plus the radial unbalanced-pull
    coefficient CE.

    Returns (ce, fax) where fax is a list of 8 (x, force) tuples,
    matching the FORTRAN's FAX(8,2) array.
    """
    af = PI * di1 / (4.0 * MY0) * delta1 * (alpha * blm) ** 2
    fax = []
    for j in range(1, 9):
        x = 0.001 * j
        y = x / delta1
        g = 1.0 / PI * (2.0 * y * math.atan(y) - math.log(1.0 + y ** 2))
        lam = 1.0 - delta1 / li * g
        dlam = -1.0 / (PI * li) * 2.0 * math.atan(x / delta1)
        fs = af / lam ** 2 * li * dlam
        gam1 = (bk - x) / delta1
        gam2 = (bk + x) / delta1
        a1 = magzu1(gam1)
        a2 = magzu1(gam2)
        fk = (nk2 * af * (-2.0 / PI * math.atan(x / delta1) +
              (2.0 - a1) / (2.0 * math.sqrt(1.0 + a1)) -
              (2.0 - a2) / (2.0 * math.sqrt(1.0 + a2))))
        if schr == 0.0:
            fn = 0.0
        else:
            gammas = 2.0 * PI / abs(schr)
            fn = 2.0 * abs(mn) / di1 * math.tan(gammas)
        fax.append((x, -fs - fk + fn))

    ce = (PI * di1 * li / (4.0 * MY0 * delta2) * (1.57 * blm) ** 2)
    if round(p) == 1:
        ce = ce / 2.0

    return ce, fax
