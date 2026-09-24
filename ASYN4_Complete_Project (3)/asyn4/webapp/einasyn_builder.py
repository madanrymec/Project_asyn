"""
Converts a dict of web-form field values (keyed by field_schema.py's
field names) into an EINASYN-format file, using asyn4.lesen's own
_CardWriter - the same utility that generates this port's own test
fixtures - so the file this produces can never drift out of sync with
what asyn4.lesen.lesen() actually reads. See field_schema.py for what
each field means and its unit.

Updated to build the full 11-column RKZ control code (previously only
4 of 11 columns were wired up - see field_schema.py's module docstring
for the corrections made after cross-checking against the original
Kirloskar/AEG "INPUT DATA FOR ASYN4" specification sheet), and to
write the MG00-MG10 11-point load-torque table into Z18 (previously
hardcoded to all zeros, silently ignoring the eleven mg00..mg10 fields
even if a user set IMG=4 and filled them in).
"""
from asyn4.lesen import _CardWriter
from .field_schema import default_values


def _v(values: dict, name: str):
    """Look up a form value, falling back to the schema default if the
    form omitted it or left it blank."""
    val = values.get(name)
    if val is None or val == "":
        return default_values()[name]
    return val


def build_rkz(values: dict) -> str:
    """
    Builds the 11-character RKZ control-code string from the individual
    per-column form fields - see field_schema.py's Z3 entries and
    asyn4/lesen.py's docstring for what each column controls.
    """
    def col(name):
        s = str(_v(values, name))
        return s[0] if s else " "

    cols = [" "] * 11
    cols[0] = col("rkz_ilfr") or "2"
    cols[1] = col("rkz_dispatch")
    cols[2] = col("rkz_slipcurve")
    cols[3] = col("rkz_wedge")
    cols[4] = col("rkz_ring_dist")
    cols[5] = col("rkz_harmonic")
    cols[6] = col("rkz_plot")
    cols[7] = col("rkz_plot_no_harmonic")
    cols[8] = col("rkz_formcoil")
    cols[9] = col("rkz_duct_align")
    cols[10] = col("rkz_6phase")
    return "".join(cols)


def build_einasyn(values: dict, output_path: str) -> str:
    """
    Writes an EINASYN-format file at output_path from the given form
    values dict. Returns output_path.
    """
    w = _CardWriter()

    # Z1
    w.write_mixed(2, [
        ("a", 16, _v(values, "name")),
        ("a", 8, _v(values, "datum")),
        ("a", 25, _v(values, "kennwort")),
    ])
    # Z2
    w.write_mixed(2, [
        ("a", 27, _v(values, "reke")),
        ("a", 17, _v(values, "vanr")),
        ("a", 15, _v(values, "fbnr")),
        ("i", 5, int(_v(values, "wanr"))),
        ("a", 12, _v(values, "monr")),
    ])
    # Z3
    w.write_mixed(2, [("a", 11, build_rkz(values))])
    # Z4
    w.write_f(3, [
        _v(values, "fn"), _v(values, "kfe"), _v(values, "kj1"),
        _v(values, "kj2"), _v(values, "p"), _v(values, "da1"),
        _v(values, "l"), _v(values, "zk"), _v(values, "bk"),
        _v(values, "pn"), _v(values, "kz"),
    ], decimals=2)
    # Z5
    w.write_mixed(2, [("i", 7, int(_v(values, "ishalt")))] + [
        ("f", 7, _v(values, "un"), 1), ("f", 7, _v(values, "prbg0"), 3),
        ("f", 7, _v(values, "v10"), 2), ("f", 7, _v(values, "cfej"), 3),
        ("f", 7, _v(values, "cfez"), 3), ("f", 7, _v(values, "schr"), 2),
        ("f", 7, _v(values, "deltag"), 2), ("f", 7, _v(values, "ku1"), 3),
        ("f", 7, _v(values, "xnetz"), 4), ("f", 7, _v(values, "kt"), 3),
    ])
    # Z6
    w.write_mixed(2, [
        ("f", 7, _v(values, "m1"), 1), ("f", 7, _v(values, "mzone1"), 1),
        ("f", 7, _v(values, "di1"), 1), ("f", 7, _v(values, "n1"), 0),
        ("f", 7, _v(values, "zn1"), 1), ("f", 7, _v(values, "weite1"), 1),
        ("f", 7, _v(values, "qcu1"), 2), ("f", 7, _v(values, "lw1"), 1),
        ("f", 7, _v(values, "theta1"), 1),
        ("i", 7, int(_v(values, "ikopp"))),
        ("f", 7, _v(values, "unx"), 1),
    ])
    # Z7
    nutnr1 = str(_v(values, "nutnr1"))[:7] or "N1"
    w.write_mixed(2, [
        ("f", 7, _v(values, "bs1"), 1), ("f", 7, _v(values, "hs1"), 1),
        ("i", 7, int(_v(values, "izgf"))), ("a", 7, nutnr1),
        ("f", 7, _v(values, "nutf1"), 1), ("f", 7, _v(values, "na"), 1),
        ("f", 7, _v(values, "da"), 1), ("f", 7, _v(values, "dta"), 1),
    ])
    # Z8
    w.write_f(2, [
        _v(values, "bn1"), _v(values, "bn1s"), _v(values, "hk1"),
        _v(values, "ho1"), _v(values, "hcuo1"), _v(values, "hzw1"),
        _v(values, "hcuu1"), _v(values, "hu1"), _v(values, "lamst1"),
        0, 0,
    ], decimals=2)
    # Z9
    nutnr2 = str(_v(values, "nutnr2"))[:5] or "N2"
    w.write_mixed(2, [
        ("f", 7, _v(values, "di2"), 1), ("f", 7, _v(values, "n2"), 0),
        ("f", 7, _v(values, "bs2"), 1), ("f", 7, _v(values, "hs2"), 1),
        ("f", 7, _v(values, "mtl"), 1), ("f", 7, _v(values, "ntl"), 1),
        ("f", 7, _v(values, "utl"), 1), ("f", 7, _v(values, "azweig"), 1),
        ("f", 7, _v(values, "btl"), 2), ("f", 7, _v(values, "htl"), 2),
        ("a", 5, nutnr2),
    ])
    # Z10
    w.write_f(3, [
        _v(values, "zn2"), _v(values, "weite2"), _v(values, "qcu2"),
        _v(values, "lw2"), _v(values, "theta2"), _v(values, "kp2"),
        _v(values, "lneb2"), _v(values, "lueb2"), _v(values, "hltr2"),
        _v(values, "bltr2"),
    ], decimals=2)
    # Z11
    w.write_mixed(2, [
        ("f", 7, _v(values, "bn2"), 1), ("f", 7, _v(values, "bn2s"), 1),
        ("f", 7, _v(values, "hk2"), 1), ("f", 7, _v(values, "ho2"), 1),
        ("f", 7, _v(values, "hcuo2"), 1), ("f", 7, _v(values, "hzw2"), 1),
        ("f", 7, _v(values, "hcuu2"), 1), ("f", 7, _v(values, "hu2"), 1),
        ("i", 7, int(_v(values, "scha2"))),
    ])
    # Z12
    w.write_f(3, [
        _v(values, "h32"), _v(values, "h42"), _v(values, "h52"),
        _v(values, "h62"), _v(values, "h72"), _v(values, "matso"),
        _v(values, "matro"), _v(values, "matsu"), _v(values, "matru"), 0,
    ], decimals=2)
    # Z13
    w.write_f(2, [
        _v(values, "br32"), _v(values, "br42"), _v(values, "br52"),
        _v(values, "br62"), _v(values, "hustab"), _v(values, "bustab"),
        _v(values, "austab"), 0, 0, 0, 0,
    ], decimals=2)
    # Z14 (KAPPA values need enough decimal precision to preserve
    # temperature-derated conductivities like 41.665 m/(ohm.mm^2))
    w.write_f(2, [
        _v(values, "kappa3"), _v(values, "kappa4"), _v(values, "kappa5"),
        _v(values, "kappa6"), _v(values, "kappa7"), _v(values, "hostab"),
        _v(values, "bostab"), _v(values, "aostab"), _v(values, "an2"),
        _v(values, "vn2"), 0,
    ], decimals=1)
    # Z15
    w.write_mixed(2, [
        ("f", 7, 0, 4), ("f", 7, _v(values, "qring"), 1),
        ("f", 7, _v(values, "kapri"), 1), ("f", 7, _v(values, "beta"), 2),
        ("f", 7, _v(values, "qstab"), 1), ("f", 7, _v(values, "qstaba"), 1),
        ("f", 7, 0, 4), ("f", 7, 0, 4), ("f", 7, 0, 4), ("f", 7, 0, 4),
        ("f", 7, 0, 4),
    ])
    # Z16
    w.write_f(3, [
        _v(values, "luea"), _v(values, "qringa"), _v(values, "kapria"),
        _v(values, "nk2"), _v(values, "bk2"),
    ], decimals=1)
    # Z17
    mg_shape = int(_v(values, "mg_shape"))
    w.write_mixed(3, [
        ("f", 14, _v(values, "jantr"), 4), ("f", 14, _v(values, "mgn"), 4),
        ("f", 7, mg_shape, 0), ("f", 7, _v(values, "mg1"), 2),
        ("f", 7, _v(values, "mg2"), 2), ("f", 7, _v(values, "mg3"), 2),
        ("f", 7, _v(values, "mg4"), 2),
    ])
    # Z18: MG00-MG10, the 11-point load-torque table (only meaningful
    # when mg_shape/IMG == 4, but always written so the field is
    # available if the user switches IMG later without resubmitting)
    w.write_f(3, [
        _v(values, "mg00"), _v(values, "mg01"), _v(values, "mg02"),
        _v(values, "mg03"), _v(values, "mg04"), _v(values, "mg05"),
        _v(values, "mg06"), _v(values, "mg07"), _v(values, "mg08"),
        _v(values, "mg09"), _v(values, "mg10"),
    ], decimals=3)
    # Z19
    w.write_mixed(3, [
        ("a", 10, _v(values, "wp1")), ("a", 4, ""),
        ("a", 10, _v(values, "wp2")), ("a", 4, ""),
    ])

    w.write(output_path)
    return output_path
