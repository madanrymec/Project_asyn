"""
Field schema for the ASYN4 web application.

This is the single source of truth mapping every field the LESEN input
reader (asyn4/lesen.py) expects from an EINASYN file to a web-form
field: name, label, unit, data type, default value, and which section
of the form it belongs in. Both the HTML form (rendered from this list)
and the EINASYN file builder (einasyn_builder.py) are driven from here,
so the form can never drift out of sync with what LESEN actually reads.

Units shown to the user match EINASYN's own raw-file convention
(millimetres, kilowatts, kilowatts; conductivity in m/(ohm.mm^2) before
the internal *1e6-to-S/m scaling) - see asyn4/lesen.py's docstring for
the full explanation of the unit-conversion block this mirrors. A
blank/0 value for most fields lets LESEN's own default-filling logic
take over, exactly as leaving a field blank in a real EINASYN file
would.

Cross-checked (this revision) against the original Kirloskar Electric /
AEG "INPUT DATA FOR ASYN4" specification sheet (AEG-A-01-004 R1,
1992) and the "INDUCTION MOTOR DESIGN OUTPUT PARAMETERS (AEG GERMANY)"
reference sheet. Corrections made from that cross-check versus this
port's earlier (inferred-from-source-code-only) labels are noted
inline where they matter - e.g. BN1/BN1S were swapped, ISHALT's real
meaning is star/delta/double-star/double-delta (not "direct line
voltage"), KAPPA fields need temperature-derated conductivity in
m/(ohm.mm^2) rather than a generic "value*1e6" placeholder, MATSO/
MATRO/MATSU/MATRU are a coded 1-6 material list, and the MG00-MG10
11-point load-torque table (used when IMG=4) is now exposed as
individual fields, which it wasn't before.
"""

FIELD_TYPE_FLOAT = "float"
FIELD_TYPE_INT = "int"
FIELD_TYPE_STR = "str"

# Each entry: (name, card, label, unit, type, default, help_text)
# `card` groups fields by their EINASYN card (Z1-Z19) for the builder;
# `section` (added below) groups them for the form UI.
FIELDS = [
    # ---- Z1: identification ----
    dict(name="name", card="Z1", section="Identification", label="Machine name",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=16),
    dict(name="datum", card="Z1", section="Identification", label="Date",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=8),
    dict(name="kennwort", card="Z1", section="Identification", label="Keyword / project code",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=25),

    # ---- Z2: identification continued ----
    dict(name="reke", card="Z2", section="Identification", label="Frame type code (REKE)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=27,
         help="e.g. 'A5-315M04'. First 2 chars matter for frame-size lookups."),
    dict(name="vanr", card="Z2", section="Identification", label="Order number (VANR)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=17),
    dict(name="fbnr", card="Z2", section="Identification", label="Design number (FBNR)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=15),
    dict(name="wanr", card="Z2", section="Identification", label="Drawing number (WANR)",
         unit="", type=FIELD_TYPE_INT, default=0),
    dict(name="monr", card="Z2", section="Identification", label="Motor number (MONR)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=12),

    # ---- Z3: calculation control code (RK1 / RKZ, 11 characters) ----
    # Column-by-column meaning per the original AEG-A-01-004 spec sheet
    # (Note 1 confirms the ILFR values below):
    #   1: ILFR (rotor type)              6: X = calculate harmonic torque
    #   2: blank=motor, G=generator       7: X = plot current and torque
    #   3: X = calculate torque-speed     8: X = plot WITHOUT harmonic content
    #   4: M = magnetic wedge used        9: blank = double-layer winding
    #   5: ring current distribution     10: blank = rotor/stator radial ducts opposite
    #                                     11: X = 6-phase winding (converter appln.)
    dict(name="rkz_ilfr", card="Z3", section="Rotor Type", label="Rotor type (ILFR) - col.1",
         unit="", type=FIELD_TYPE_STR, default="2",
         choices=[("1", "Wound rotor (slip-ring)"), ("2", "Single cage"),
                   ("3", "Double cage, brazed/diecast/sash bar, ONE ring"),
                   ("4", "Double cage, brazed, TWO rings")]),
    dict(name="rkz_dispatch", card="Z3", section="Rotor Type", label="Motor/generator - col.2",
         unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Motor duty"), ("G", "Generator duty")],
         help="RKZ column 2. Note: this port's revolving-field-magnet (braking) duty dispatch "
              "also uses this column with value 'D' - not on the original spec sheet, kept for "
              "backward compatibility with earlier calculations run through this port."),
    dict(name="rkz_slipcurve", card="Z3", section="Rotor Type",
         label="Calculate torque-speed characteristic - col.3",
         unit="", type=FIELD_TYPE_STR, default="X",
         choices=[("", "No (quick check only)"), ("X", "Yes")]),
    dict(name="rkz_wedge", card="Z3", section="Stator Slot", label="Magnetic wedge - col.4",
         unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "No magnetic wedge"), ("M", "Magnetic wedge used")]),
    dict(name="rkz_ring_dist", card="Z3", section="Rotor Ring (Cage)",
         label="Ring current distribution - col.5", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Uniformly distributed"), ("O", "Non-uniformly distributed"),
                   ("H", "Non-uniform, blended below S=0.7"),
                   ("J", "Non-uniform, blended above S=0.7")],
         help="Exact column-5 code semantics are the one item on the original spec sheet this "
              "port could not fully reconcile with the FORTRAN source (the sheet's OCR scan was "
              "ambiguous here) - the O/H/J values are confirmed directly from source code."),
    dict(name="rkz_harmonic", card="Z3", section="Rotor Type",
         label="Calculate harmonic torque - col.6", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "No"), ("X", "Yes - also runs synchronous-torque (SYNMOM) analysis")]),
    dict(name="rkz_plot", card="Z3", section="Rotor Type",
         label="Plot current and torque - col.7", unit="", type=FIELD_TYPE_STR, default="X",
         choices=[("", "No plot"), ("X", "Generate torque/current-speed chart")]),
    dict(name="rkz_plot_no_harmonic", card="Z3", section="Rotor Type",
         label="Plot without harmonic content - col.8", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Include harmonic content in plot"), ("X", "Plot fundamental only")]),
    dict(name="rkz_formcoil", card="Z3", section="Stator Winding",
         label="Winding layer - col.9", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Double-layer winding"),
                   ("E", "Fractional-slot form-coil special case")],
         help="Spec sheet says 'blank for double-layer winding' without enumerating the "
              "non-blank options; 'E' (fractional-slot/form-coil) is confirmed from source code."),
    dict(name="rkz_duct_align", card="Z3", section="Cooling",
         label="Radial cooling ducts alignment - col.10", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Rotor/stator ducts opposite (aligned)"),
                   ("X", "Ducts not opposite")]),
    dict(name="rkz_6phase", card="Z3", section="Stator Winding",
         label="6-phase converter winding - col.11", unit="", type=FIELD_TYPE_STR, default="",
         choices=[("", "Normal 3-phase winding"),
                   ("X", "6-phase winding for converter application (3ph in/3ph out)")]),

    # ---- Z4: main ratings ----
    dict(name="fn", card="Z4", section="Rating", label="Frequency (FN)", unit="Hz",
         type=FIELD_TYPE_FLOAT, default=50.0),
    dict(name="kfe", card="Z4", section="Rating", label="Iron stacking factor (KFE)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0,
         help="0 = default 0.92. Guidance: 0.92 with radial vent ducts, 0.96 without"),
    dict(name="kj1", card="Z4", section="Rating", label="Core-height factor for magnetising current (KJ1)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = default 1.0, otherwise 0"),
    dict(name="kj2", card="Z4", section="Rating", label="Rotor yoke-height factor (KJ2)",
         unit="", type=FIELD_TYPE_FLOAT, default=1.0,
         help="Default 1.0; >1 for 2-pole rotor or rotor with ribbed shaft; 0 otherwise"),
    dict(name="p", card="Z4", section="Rating", label="Number of pole pairs (P)", unit="",
         type=FIELD_TYPE_FLOAT, default=2.0),
    dict(name="da1", card="Z4", section="Stator Geometry", label="Outer diameter of stator stamping (DA1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=400.0),
    dict(name="l", card="Z4", section="Stator Geometry", label="Gross core length (L)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=200.0),
    dict(name="zk", card="Z4", section="Cooling", label="Number of radial ducts (ZK)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="bk", card="Z4", section="Cooling", label="Width of radial ducts (BK)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = default 10mm"),
    dict(name="pn", card="Z4", section="Rating", label="Output of motor (PN)", unit="kW",
         type=FIELD_TYPE_FLOAT, default=2.0),
    dict(name="kz", card="Z4", section="Cooling", label="Thermal coefficient, internally cooled motor (KZ)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0,
         help="0 = auto (frame-size lookup); see design writeup. Can be back-derived from test "
              "data: KZ = (dT x 1000) / (g x A x sqrt((Pcu_stator+Pfe)/Pcu_stator)), where dT is "
              "measured temperature rise (K), g is current density (A/mm^2), A is slot loading "
              "(A/cm). Example worked calc: dT=54K, g=3.88 A/mm^2, A=680 A/cm, Pcu=30kW, "
              "Pfe=28kW -> KZ = 54000 / (3.88 x 680 x sqrt(58/30)) = 14.7, rounded to 15."),

    # ---- Z5: supply / operating conditions ----
    dict(name="ishalt", card="Z5", section="Rating", label="Winding connection (ISHALT)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "Star"), (1, "Delta"),
                   (2, "Double star (variable speed)"), (3, "Double delta (variable speed)")]),
    dict(name="un", card="Z5", section="Rating", label="Rated voltage (UN)",
         unit="V", type=FIELD_TYPE_FLOAT, default=400.0),
    dict(name="prbg0", card="Z5", section="Rating", label="Friction and windage loss (PRBG0)",
         unit="kW", type=FIELD_TYPE_FLOAT, default=0.0,
         help="0 = auto (Nurnberg formula). Independent of bearing type (brush + fans)"),
    dict(name="v10", card="Z5", section="Materials", label="Theoretical iron loss @1T/50Hz (V10)",
         unit="W/kg", type=FIELD_TYPE_FLOAT, default=1.65,
         help="Must be one of: 1.25, 1.35, 1.5, 1.7, 2.0, 2.3, 2.6, 3.0, 3.6"),
    dict(name="cfej", card="Z5", section="Materials", label="Actual iron loss in yoke (CFEJ)",
         unit="W/kg", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = to be computed"),
    dict(name="cfez", card="Z5", section="Materials", label="Actual iron loss in teeth (CFEZ)",
         unit="W/kg", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = to be computed"),
    dict(name="schr", card="Z5", section="Rotor Geometry", label="Skew, oblong of circumference (SCHR)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0,
         help="High value = no skew (typical for AW machines and size >400 frame); 0 otherwise. "
              "No. of slots skewed by 1 slot, or half that value if skew = 1/2 slot"),
    dict(name="deltag", card="Z5", section="Stator Geometry", label="Radial air gap (DELTAG)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.6),
    dict(name="ku1", card="Z5", section="Starting", label="Fall in voltage during start (KU1)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0,
         help="0 = no reduced-voltage start"),
    dict(name="xnetz", card="Z5", section="Starting", label="External reactance to motor terminals (XNETZ)",
         unit="Ohm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="kt", card="Z5", section="Cooling", label="Thermal coefficient, surface cooled motor (KT)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto"),

    # ---- Z6: stator winding ----
    dict(name="m1", card="Z6", section="Stator Winding", label="Number of phases (M1)",
         unit="", type=FIELD_TYPE_FLOAT, default=3.0, help="Default 3"),
    dict(name="mzone1", card="Z6", section="Stator Winding", label="Number of leads (MZONE1)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = default 6"),
    dict(name="di1", card="Z6", section="Stator Geometry", label="Inner diameter of stator stamping (DI1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=250.0),
    dict(name="n1", card="Z6", section="Stator Geometry", label="Number of stator slots (N1)",
         unit="", type=FIELD_TYPE_FLOAT, default=36.0),
    dict(name="zn1", card="Z6", section="Stator Winding",
         label="Conductors/slot x parallel paths (ZN1)",
         unit="", type=FIELD_TYPE_FLOAT, default=40.0,
         help="This is conductors/slot MULTIPLIED BY number of parallel paths, not simply "
              "turns per slot"),
    dict(name="weite1", card="Z6", section="Stator Winding", label="Pitch of stator winding (WEITE1)",
         unit="slots", type=FIELD_TYPE_FLOAT, default=8.0),
    dict(name="qcu1", card="Z6", section="Stator Winding",
         label="Area of each copper conductor x conductors/slot x parallel paths (QCU1)",
         unit="mm^2", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Total copper cross-section per slot circuit, not a single strand's area. "
              "0 = auto-select (random-wound) or must be set for form-coil"),
    dict(name="lw1", card="Z6", section="Stator Winding", label="Overhang length (LW1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto"),
    dict(name="theta1", card="Z6", section="Materials", label="Temperature rise in stator (THETA1)",
         unit="degC", type=FIELD_TYPE_FLOAT, default=75.0),
    dict(name="ikopp", card="Z6", section="Stator Winding", label="IKOPP",
         unit="", type=FIELD_TYPE_INT, default=0, help="0 or 1, not used"),
    dict(name="unx", card="Z6", section="Materials", label="Insulation voltage (UNX)",
         unit="V", type=FIELD_TYPE_FLOAT, default=0.0,
         help="0 = use UN. A 5% safety margin is built into the program automatically"),

    # ---- Z7: stator slot (part 1) ----
    dict(name="bs1", card="Z7", section="Stator Slot", label="Slot opening width of stator (BS1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=6.0, help="See enclosed slot-geometry drawing"),
    dict(name="hs1", card="Z7", section="Stator Slot", label="Lip height of stator slot (HS1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0, help="See enclosed slot-geometry drawing"),
    dict(name="izgf", card="Z7", section="Stator Winding", label="Coil type (IZGF)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "Mush windings (random-wound)"), (1, "Formed coils")]),
    dict(name="dta", card="Z7", section="Cooling", label="Vent-hole pitch diameter (DTA)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),

    # ---- Z8: stator slot (part 2) ----
    dict(name="bn1", card="Z8", section="Stator Slot", label="Stator slot width at opening (BN1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=9.0),
    dict(name="bn1s", card="Z8", section="Stator Slot", label="Stator slot width at bottom (BN1S)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto"),
    dict(name="hk1", card="Z8", section="Stator Slot", label="Wedge groove height (HK1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="ho1", card="Z8", section="Stator Slot", label="Space for 2nd winding, upper (HO1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="hcuo1", card="Z8", section="Stator Slot", label="Outer bar height (HCUO1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=5.0),
    dict(name="hzw1", card="Z8", section="Stator Slot", label="Thickness of separator (HZW1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="hcuu1", card="Z8", section="Stator Slot", label="Inner bar height (HCUU1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=10.0),
    dict(name="hu1", card="Z8", section="Stator Slot", label="Space for 2nd winding, lower (HU1)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=2.0),
    dict(name="lamst1", card="Z8", section="Stator Winding", label="Overhang leakage permeance (LAMST1)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = default"),
    dict(name="nutnr1", card="Z7", section="Stator Slot", label="Shape of stator slots (NUTNR1)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=7,
         help="For printing/drawing labels only - does not affect the calculation"),
    dict(name="nutf1", card="Z7", section="Stator Slot", label="Type of stator slots (NUTF1)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0,
         help="For printing/drawing labels only - does not affect the calculation"),

    # ---- Z9: rotor geometry + bar (part 1) ----
    dict(name="na", card="Z7", section="Cooling", label="Number of axial vent holes (NA)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="da", card="Z7", section="Cooling", label="Diameter of axial vent holes in rotor (DA)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="di2", card="Z9", section="Rotor Geometry", label="Rotor inner diameter (DI2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=150.0),
    dict(name="n2", card="Z9", section="Rotor Geometry", label="Number of rotor slots (N2)",
         unit="", type=FIELD_TYPE_FLOAT, default=24.0),
    dict(name="bs2", card="Z9", section="Rotor Slot", label="Width of rotor slot at opening (BS2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=5.0),
    dict(name="hs2", card="Z9", section="Rotor Slot", label="Lip height of rotor slot (HS2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=20.0),
    dict(name="mtl", card="Z9", section="Rotor Winding (Slip-Ring)",
         label="No. of parallel conductors in slot, depthwise (MTL)",
         unit="", type=FIELD_TYPE_FLOAT, default=4.0),
    dict(name="ntl", card="Z9", section="Rotor Winding (Slip-Ring)",
         label="No. of sections widthwise (NTL)", unit="", type=FIELD_TYPE_FLOAT, default=2.0),
    dict(name="utl", card="Z9", section="Rotor Winding (Slip-Ring)",
         label="No. of electrical sections (UTL)", unit="", type=FIELD_TYPE_FLOAT, default=1.0,
         help="Always 1"),
    dict(name="azweig", card="Z9", section="Stator Winding", label="No. of parallel paths in stator (AZWEIG)",
         unit="", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="btl", card="Z9", section="Stator Winding", label="Width of stator conductor (BTL)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.2),
    dict(name="htl", card="Z9", section="Stator Winding", label="Depth of stator conductor (HTL)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.2),
    dict(name="nutnr2", card="Z9", section="Rotor Slot", label="Shape type of rotor slot (NUTNR2)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=5,
         help="For printing/drawing labels only - does not affect the calculation"),

    # ---- Z10: slip-ring rotor winding ----
    dict(name="zn2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="No. of bars/turns per slot in rotor (ZN2)", unit="", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="weite2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="Pitch of rotor winding (WEITE2)", unit="slots", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="qcu2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="Area of copper in rotor, each slot (QCU2)", unit="mm^2",
         type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto-select"),
    dict(name="lw2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="Overhang length in rotor (LW2)", unit="mm", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Similar convention to LW1"),
    dict(name="theta2", card="Z10", section="Materials", label="Rotor temperature rise (THETA2)",
         unit="degC", type=FIELD_TYPE_FLOAT, default=75.0),
    dict(name="kp2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="No. of parallel paths in rotor (KP2)", unit="", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="lneb2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="No. of sections of rotor cond. widthwise (LNEB2)", unit="",
         type=FIELD_TYPE_FLOAT, default=0.0, help="Rotor equivalent of NTL"),
    dict(name="lueb2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="No. of cond. in rotor slot depthwise (LUEB2)", unit="",
         type=FIELD_TYPE_FLOAT, default=0.0, help="Rotor equivalent of MTL"),
    dict(name="hltr2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="Depth of rotor conductor (HLTR2)", unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="bltr2", card="Z10", section="Rotor Winding (Slip-Ring)",
         label="Width of rotor conductor (BLTR2)", unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),

    # ---- Z11: rotor slot (part 2) / slip-ring only ----
    dict(name="bn2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Rotor slot width at bottom (BN2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=7.0),
    dict(name="bn2s", card="Z11", section="Rotor Slot (Slip-Ring)", label="Rotor slot width at centre (BN2S)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto"),
    dict(name="hk2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Rotor slot height (HK2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="ho2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Space for 2nd winding, upper (HO2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="hcuo2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Outer bar height (HCUO2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=5.0),
    dict(name="hzw2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Thickness of separator (HZW2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="hcuu2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Inner bar position (HCUU2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=5.0),
    dict(name="hu2", card="Z11", section="Rotor Slot (Slip-Ring)", label="Space for 2nd winding, lower (HU2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="scha2", card="Z11", section="Rotor Winding (Slip-Ring)",
         label="Rotor connection (SCHA2)", unit="", type=FIELD_TYPE_INT, default=1,
         choices=[(0, "Star"), (1, "Delta")]),

    # ---- Z12: cage bar geometry ----
    dict(name="h32", card="Z12", section="Rotor Bar (Cage)", label="Dimension of rotor bar, H32",
         unit="mm", type=FIELD_TYPE_FLOAT, default=20.0, help="See enclosed slot-geometry drawing"),
    dict(name="h42", card="Z12", section="Rotor Bar (Cage)", label="Dimension of rotor bar, H42",
         unit="mm", type=FIELD_TYPE_FLOAT, default=2.0, help="See enclosed slot-geometry drawing"),
    dict(name="h52", card="Z12", section="Rotor Bar (Cage, Double)", label="Double-cage rotor dimension, H52",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="See enclosed slot-geometry drawing"),
    dict(name="h62", card="Z12", section="Rotor Bar (Cage, Double)", label="Double-cage rotor dimension, H62",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="See enclosed slot-geometry drawing"),
    dict(name="h72", card="Z12", section="Rotor Bar (Cage, Double)", label="Double-cage rotor dimension, H72",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="See enclosed slot-geometry drawing"),
    dict(name="matso", card="Z12", section="Materials", label="Material for outer bar (MATSO)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "(unset)"), (1, "Copper"), (2, "Brass 95 (CuZn5)"), (3, "Brass 63 (CuZn37)"),
                   (4, "Bronze"), (5, "Aluminium"), (6, "GDAlSi - 2x aluminium resistivity")]),
    dict(name="matro", card="Z12", section="Materials", label="Material for outer ring (MATRO)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "(unset)"), (1, "Copper"), (2, "Brass 95 (CuZn5)"), (3, "Brass 63 (CuZn37)"),
                   (4, "Bronze"), (5, "Aluminium"), (6, "GDAlSi - 2x aluminium resistivity")]),
    dict(name="matsu", card="Z12", section="Materials", label="Material for inner bar (MATSU)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "(unset)"), (1, "Copper"), (2, "Brass 95 (CuZn5)"), (3, "Brass 63 (CuZn37)"),
                   (4, "Bronze"), (5, "Aluminium"), (6, "GDAlSi - 2x aluminium resistivity")]),
    dict(name="matru", card="Z12", section="Materials", label="Material for inner ring (MATRU)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "(unset)"), (1, "Copper"), (2, "Brass 95 (CuZn5)"), (3, "Brass 63 (CuZn37)"),
                   (4, "Bronze"), (5, "Aluminium"), (6, "GDAlSi - 2x aluminium resistivity")]),

    # ---- Z13: cage bar widths ----
    dict(name="br32", card="Z13", section="Rotor Bar (Cage)", label="Dimension of double-cage rotor, BR32",
         unit="mm", type=FIELD_TYPE_FLOAT, default=6.0, help="See enclosed slot-geometry drawing"),
    dict(name="br42", card="Z13", section="Rotor Bar (Cage)", label="Dimension of double-cage rotor, BR42",
         unit="mm", type=FIELD_TYPE_FLOAT, default=5.0, help="See enclosed slot-geometry drawing"),
    dict(name="br52", card="Z13", section="Rotor Bar (Cage, Double)", label="Dimension of double-cage rotor, BR52",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="See enclosed slot-geometry drawing"),
    dict(name="br62", card="Z13", section="Rotor Bar (Cage, Double)", label="Dimension of double-cage rotor, BR62",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="See enclosed slot-geometry drawing"),
    dict(name="hustab", card="Z13", section="Rotor Bar (Cage)", label="Depth of inner bar (HUSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="bustab", card="Z13", section="Rotor Bar (Cage)", label="Width of inner bar at bottom (BUSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="austab", card="Z13", section="Rotor Bar (Cage)", label="Width of inner bar at top (AUSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),

    # ---- Z14: conductivities ----
    dict(name="kappa3", card="Z14", section="Materials", label="Conductivity in area 3 (KAPPA3)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Enter conductivity AT OPERATING TEMPERATURE (already temperature-derated), in "
              "m/(ohm.mm^2) - the program multiplies by 1e6 internally to get S/m. Base (20C) "
              "reference values: Copper=58, Brass CuZn5=32, Brass CuZn37=15 m/(ohm.mm^2). "
              "Derate to operating temperature via: kappa = 1 / (1/kappa_20C + tempco x rise), "
              "with tempco (ohm.mm^2/m per K) approx: Cu=0.000068, CuZn5=0.000071, "
              "CuZn37=0.000112. Worked example (copper bar, Class B rise = 80K over 40C "
              "ambient): resistivity(20C)=0.01724 ohm.mm^2/m -> resistivity(operating) = "
              "(235+120)/(235+20) x 0.01724 = 0.024 -> kappa = 41.665 m/(ohm.mm^2)."),
    dict(name="kappa4", card="Z14", section="Materials", label="Conductivity in area 4 (KAPPA4)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Same convention as KAPPA3 - enter value at operating temperature"),
    dict(name="kappa5", card="Z14", section="Materials", label="Conductivity in area 5 (KAPPA5)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Same convention as KAPPA3 - enter value at operating temperature"),
    dict(name="kappa6", card="Z14", section="Materials", label="Conductivity in area 6 (KAPPA6)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Same convention as KAPPA3 - enter value at operating temperature"),
    dict(name="kappa7", card="Z14", section="Materials", label="Conductivity in area 7 (KAPPA7)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Same convention as KAPPA3 - enter value at operating temperature"),
    dict(name="hostab", card="Z14", section="Rotor Bar (Cage)", label="Depth of outer bar (HOSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="bostab", card="Z14", section="Rotor Bar (Cage)", label="Width of outer bar at bottom (BOSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="aostab", card="Z14", section="Rotor Bar (Cage)", label="Width of outer bar at top (AOSTAB)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="an2", card="Z14", section="Rotor Slot", label="Dimension of slot type, AN2",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Cross-referenced to a specific slot-type code (e.g. type 159) on the original drawing"),
    dict(name="vn2", card="Z14", section="Rotor Slot", label="Dimension of slot type, VN2",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="Paired with AN2"),

    # ---- Z15: ring geometry ----
    dict(name="lue", card="Z15", section="Rotor Ring (Cage)", label="Length of bar outside the core (LUE)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=15.0),
    dict(name="qring", card="Z15", section="Rotor Ring (Cage)", label="Area of ring (QRING)",
         unit="mm^2", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = auto-estimate"),
    dict(name="kapri", card="Z15", section="Materials", label="Ring conductivity (KAPRI)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=30.0,
         help="Same derating convention as KAPPA3 (operating temperature). Typical: ring IACS "
              "is 85% (min) of the derated bar value per drawing, i.e. KAPRI approx 0.85 x "
              "KAPPA(bar)"),
    dict(name="beta", card="Z15", section="Rotor Ring (Cage)", label="Skew angle (BETA)",
         unit="deg", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="qstab", card="Z15", section="Rotor Bar (Cage)", label="Operating area of bar (QSTAB)",
         unit="mm^2", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = computed from bar geometry"),
    dict(name="qstaba", card="Z15", section="Rotor Bar (Cage, Double)",
         label="Area of cross bar (QSTABA)", unit="mm^2", type=FIELD_TYPE_FLOAT, default=0.0),

    # ---- Z16: outer ring (double cage) + lamination stack ----
    dict(name="luea", card="Z16", section="Rotor Ring (Cage, Double)",
         label="Length of cross bar (LUEA)", unit="mm", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Suffix A stands for outer cage"),
    dict(name="qringa", card="Z16", section="Rotor Ring (Cage, Double)",
         label="Area of ring during starting (QRINGA)", unit="mm^2", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Suffix A stands for outer cage"),
    dict(name="kapria", card="Z16", section="Materials", label="Conductivity of ring, outer cage (KAPRIA)",
         unit="m/(ohm.mm^2)", type=FIELD_TYPE_FLOAT, default=0.0,
         help="Same derating convention as KAPPA3 (operating temperature). Suffix A stands for outer cage"),
    dict(name="nk2", card="Z16", section="Cooling", label="Number of radial ducts in rotor (NK2)",
         unit="", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = same as ZK"),
    dict(name="bk2", card="Z16", section="Cooling", label="Width of radial ducts in rotor (BK2)",
         unit="mm", type=FIELD_TYPE_FLOAT, default=0.0, help="0 = same as BK"),

    # ---- Z17/Z18: counter-torque (load) curve ----
    dict(name="jantr", card="Z17", section="Load Curve", label="Load inertia (JANTR)",
         unit="kg*m^2", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mgn", card="Z17", section="Load Curve", label="Nominal load torque (MGN)",
         unit="", type=FIELD_TYPE_FLOAT, default=1.0),
    dict(name="mg_shape", card="Z17", section="Load Curve", label="Type of load during starting (IMG)",
         unit="", type=FIELD_TYPE_INT, default=0,
         choices=[(0, "No load"), (1, "Constant load"),
                   (2, "Linear load"), (3, "Quadrilateral load"),
                   (4, "Any other load - specify 11 points below")]),
    dict(name="mg1", card="Z17", section="Load Curve", label="Load torque parameter MG(1)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0,
         help="MGLOS - per-unit starting torque, when IMG is 1/2/3"),
    dict(name="mg2", card="Z17", section="Load Curve", label="Load torque parameter MG(2)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0,
         help="MGSATT - saddle point during starting, when IMG is 1/2/3"),
    dict(name="mg3", card="Z17", section="Load Curve", label="Load torque parameter MG(3)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg4", card="Z17", section="Load Curve", label="Load torque parameter MG(4)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0,
         help="MGEND - load torque at synchronous speed, when IMG is 1/2/3"),

    # ---- Z18: 11-point load-torque table (only used when IMG=4 above) ----
    dict(name="mg00", card="Z18", section="Load Curve", label="Load torque at N=1.0 p.u. (MG00)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0, help="Synchronous speed"),
    dict(name="mg01", card="Z18", section="Load Curve", label="Load torque at N=0.9 p.u. (MG01)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg02", card="Z18", section="Load Curve", label="Load torque at N=0.8 p.u. (MG02)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg03", card="Z18", section="Load Curve", label="Load torque at N=0.7 p.u. (MG03)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg04", card="Z18", section="Load Curve", label="Load torque at N=0.6 p.u. (MG04)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg05", card="Z18", section="Load Curve", label="Load torque at N=0.5 p.u. (MG05)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg06", card="Z18", section="Load Curve", label="Load torque at N=0.4 p.u. (MG06)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg07", card="Z18", section="Load Curve", label="Load torque at N=0.3 p.u. (MG07)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg08", card="Z18", section="Load Curve", label="Load torque at N=0.2 p.u. (MG08)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg09", card="Z18", section="Load Curve", label="Load torque at N=0.1 p.u. (MG09)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0),
    dict(name="mg10", card="Z18", section="Load Curve", label="Load torque at N=0 (locked rotor) (MG10)",
         unit="p.u.", type=FIELD_TYPE_FLOAT, default=0.0, help="Standstill"),

    # ---- Z19: winding-diagram cross-reference (informational only) ----
    dict(name="wp1", card="Z19", section="Identification", label="Winding reference number, stator (WP1)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=10,
         help="Printed only if entered"),
    dict(name="wp2", card="Z19", section="Identification", label="Winding reference number, rotor (WP2)",
         unit="", type=FIELD_TYPE_STR, default="", maxlen=10,
         help="Printed only if entered"),
]

# Ordered list of section names, controlling the order tabs/groups appear
# in the form. Any field whose section isn't listed here falls back to
# appearing at the end, so this list doesn't need to be kept in lockstep
# with FIELDS - it's just display ordering.
SECTION_ORDER = [
    "Identification", "Rotor Type", "Rating", "Stator Geometry",
    "Stator Winding", "Stator Slot", "Rotor Geometry", "Rotor Slot",
    "Rotor Slot (Slip-Ring)", "Rotor Winding (Slip-Ring)",
    "Rotor Bar (Cage)", "Rotor Bar (Cage, Double)", "Rotor Ring (Cage)",
    "Rotor Ring (Cage, Double)", "Materials", "Cooling", "Starting",
    "Load Curve",
]

# ---------------------------------------------------------------------
# Dynamic visibility rules for the frontend: which sections/fields only
# make sense for certain rotor types or load-curve shapes. Rather than
# showing all ~150 fields at once, the form only shows what's relevant
# to the current ROTOR TYPE (RKZ column 1) and LOAD-CURVE SHAPE (IMG)
# selections, updating live as the person changes those two controls.
#
# SECTION_VISIBILITY: section name -> {"field": controlling field name,
# "values": [...allowed values of that field, as strings...]}. A
# section with no entry here is always shown.
#
# FIELD_VISIBILITY: same shape, but for individual fields within an
# always-shown section (used for mg1-4 vs mg00-10, which share the
# "Load Curve" section but depend on which load shape is selected).
# ---------------------------------------------------------------------
SECTION_VISIBILITY = {
    "Rotor Slot (Slip-Ring)": {"field": "rkz_ilfr", "values": ["1"]},
    "Rotor Winding (Slip-Ring)": {"field": "rkz_ilfr", "values": ["1"]},
    "Rotor Bar (Cage)": {"field": "rkz_ilfr", "values": ["2", "3", "4"]},
    "Rotor Bar (Cage, Double)": {"field": "rkz_ilfr", "values": ["3", "4"]},
    "Rotor Ring (Cage)": {"field": "rkz_ilfr", "values": ["2", "3", "4"]},
    "Rotor Ring (Cage, Double)": {"field": "rkz_ilfr", "values": ["4"]},
}

FIELD_VISIBILITY = {
    "mg1": {"field": "mg_shape", "values": ["1", "2", "3"]},
    "mg2": {"field": "mg_shape", "values": ["1", "2", "3"]},
    "mg3": {"field": "mg_shape", "values": ["1", "2", "3"]},
    "mg4": {"field": "mg_shape", "values": ["1", "2", "3"]},
    "mg00": {"field": "mg_shape", "values": ["4"]},
    "mg01": {"field": "mg_shape", "values": ["4"]},
    "mg02": {"field": "mg_shape", "values": ["4"]},
    "mg03": {"field": "mg_shape", "values": ["4"]},
    "mg04": {"field": "mg_shape", "values": ["4"]},
    "mg05": {"field": "mg_shape", "values": ["4"]},
    "mg06": {"field": "mg_shape", "values": ["4"]},
    "mg07": {"field": "mg_shape", "values": ["4"]},
    "mg08": {"field": "mg_shape", "values": ["4"]},
    "mg09": {"field": "mg_shape", "values": ["4"]},
    "mg10": {"field": "mg_shape", "values": ["4"]},
}


def visibility_rules():
    """Returns a plain dict (JSON-serializable) the frontend can use to
    show/hide sections and fields as the person changes rkz_ilfr or
    mg_shape - see SECTION_VISIBILITY/FIELD_VISIBILITY above."""
    return {"sections": SECTION_VISIBILITY, "fields": FIELD_VISIBILITY}


def fields_by_section():
    """Returns an ordered dict: section name -> list of field dicts."""
    sections = {name: [] for name in SECTION_ORDER}
    for f in FIELDS:
        sections.setdefault(f["section"], []).append(f)
    # drop empty sections while preserving order
    return {k: v for k, v in sections.items() if v}


def default_values():
    """Returns a dict of field name -> default value, for pre-filling the form."""
    return {f["name"]: f["default"] for f in FIELDS}


def field_names():
    return [f["name"] for f in FIELDS]
