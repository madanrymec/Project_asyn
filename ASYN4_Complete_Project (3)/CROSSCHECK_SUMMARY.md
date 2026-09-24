# ASYN4 — Cross-Check Against Original Kirloskar/AEG Documents

This note summarizes a review of two real documents you provided:
`Asyn_Design_input_sheet.pdf` (the original "INPUT DATA FOR ASYN4" spec
sheet, Kirloskar Electric / AEG-A-01-004 R1, 1992) and
`ASYN_DESIGN_input_kec.pdf` (the "INDUCTION MOTOR DESIGN OUTPUT
PARAMETERS (AEG GERMANY)" reference sheet, plus worked calculation
examples for KAPPA and KZ). These are the actual original design-
engineering documentation for ASYN4 — the most authoritative source
available for this port, and genuinely valuable: they caught real
errors and filled real gaps that no amount of FORTRAN-source-only
reading would have caught.

**Files changed:** `asyn4/reports.py`, `webapp/field_schema.py`,
`webapp/einasyn_builder.py`. All changes tested — 86/86 tests still
pass, full pipeline runs end-to-end from a synthetic file.

---

## 1. Input side (`webapp/field_schema.py`, `webapp/einasyn_builder.py`)

### Real errors fixed (would have produced wrong results)
- **`BN1`/`BN1S` and `BN2`/`BN2S` were reversed.** `BN1` is slot width
  *at opening*, not bottom; `BN1S` is at the bottom. Confirmed against
  the enclosed slot-geometry drawing. `BN2S` (rotor) is at *centre*,
  not top — a different convention from the stator side.
- **`ISHALT`** — real meaning is Star / Delta / Double-Star / Double-
  Delta (for variable-speed motors). Earlier version had guessed
  "direct line voltage" framing, which was wrong.
- **`UTL`** — must always be 1. Default was 0.

### A genuine functional gap, now fixed
- The **11-point load-torque table** (`MG00`–`MG10`, used when
  `IMG=4`) was never exposed as form fields — `Z18` was silently
  hardcoded to all zeros. Now wired through as 11 real fields and
  tested end-to-end (`MG00`→`mg[5]` … `MG10`→`mg[15]`, confirmed
  against the FORTRAN's own indexing convention).

### The most valuable correction — `KAPPA` conductivity fields
The worked "Kappa value by calculation" and "Rotor bar conductivity"
sheets revealed that `KAPPA3`–`KAPPA7`, `KAPRI`, `KAPRIA` must be
entered **already temperature-derated** for the winding's expected
operating rise — not a generic room-temperature material constant.
The correct formula and reference values are now built into the field
help text:

```
kappa = 1 / (1/kappa_20C + tempco x temperature_rise)
```

Base (20°C) values, m/(Ω·mm²): Copper=58, Brass CuZn5=32, Brass
CuZn37=15. Temperature coefficients (Ω·mm²/m per K): Cu=0.000068,
CuZn5=0.000071, CuZn37=0.000112.

Worked example from the sheet (copper bar, Class B rise = 80K over
40°C ambient): resistivity at 20°C = 0.01724 Ω·mm²/m → at operating
temperature = (235+120)/(235+20) × 0.01724 = 0.024 Ω·mm²/m → KAPPA =
41.665 m/(Ω·mm²).

This also **independently confirmed** an earlier finding from this
port's FORTRAN-source review: the `×1,000,000` conversion applied to
KAPPA fields internally (mm-based units → S/m) is dimensionally
correct — the sheet's own numbers work out to the same conversion
factor when checked by hand.

### Also corrected
- `MATSO`/`MATRO`/`MATSU`/`MATRU` now have real dropdown choices
  (1=Copper, 2=Brass95, 3=Brass63, 4=Bronze, 5=Aluminium,
  6=GDAlSi/2× aluminium resistivity) instead of unexplained numbers.
- `HO1`/`HU1`/`HO2`/`HU2` correctly relabeled "space reserved for 2nd
  winding" — not insulation, as previously guessed.
- `THETA1`/`THETA2` correctly labeled **temperature rise**, not
  absolute operating temperature.
- Added previously-missing fields the reader (`LESEN`) already
  consumes but the form never exposed: `DTA` (vent-hole pitch
  diameter), `NA`/`DA` (axial vent holes), `NUTNR1`/`NUTF1`/`NUTNR2`
  (cosmetic slot-shape codes for drawing labels only).
- The 11-character `RKZ` control code now has **all 11 columns**
  properly separated and labeled (previously columns 3 and 7 — "compute
  torque-speed curve" vs. "plot current and torque" — were incorrectly
  treated as one flag).
- `KZ` field help text now includes the back-derivation formula from
  the "K2 calculation" worked example (`KZ = (ΔT×1000)/(g×A×√((Pcu+Pfe)/Pcu))`),
  useful if you ever need to calibrate KZ from test data for a new
  frame/cooling type.

---

## 2. Output side (`asyn4/reports.py`)

The AEG "Output Parameters" reference sheet is a comprehensive glossary
of what the original program reports. Cross-checking it against my
port confirmed the overwhelming majority of output values already
match correctly (KC1/KC2, HOJOC1/HOJOC2, all the ESB equivalent-circuit
values, weights, moment of inertia, lamination areas, XKEY's leakage
breakdown, HOLAUF's heating results — all consistent). Two real gaps
were found and fixed:

### New: current densities and heating rates at locked rotor
The sheet documents `G1`/`G2`/`G2A`/`GR`/`GRA` (current densities,
A/mm²), `A1EFF`/`A2EFF` (effective surface current density, A/cm), and
`TA1`/`TA2`/`TAA2`/`TAR`/`TARA` (adiabatic heating rates, K/sec) at
the locked-rotor operating point. I confirmed the exact formulas
directly in the FORTRAN source (previously not ported):

```
G1 = I1/QCU1, G2 = I2/QSTAB, G2A = I2O/QSTABA
GR = I2/QRING/(2·sin(P·π/N2)), GRA = I2O/QRINGA/(2·sin(P·π/N2))
TA1 = G1²/K41, TA2 = K1·G2²/K42, TAA2 = G2A²/K4A2
TAR = K1R·GR²/K4R2, TARA = K1RA·GRA²/K4AR2
```

Added as a new `_locked_rotor_densities_and_heating()` helper in
`reports.py`, called at both the theoretical and KKA-corrected
locked-rotor points inside `tlast()`. Results land in
`tlast_result['theoretical']['densities_heating']` and
`['corrected']['densities_heating']`.

**A units subtlety worth knowing if you use these values directly:**
the FORTRAN's raw internal values are SI (A/m², A/m); the AEG-
documented display units are A/mm² and A/cm. I confirmed the exact
scaling directly from the source's own print statement
(`GRA/1.E6, A1EFF/100.`) and applied the same `/1,000,000` and `/100`
conversions before returning these values, so what you get back
already matches the documented units — no further conversion needed
on your end.

### Also newly surfaced (already computed internally, just not returned)
- `G1`/`G2`/`A1EFF`/`A2EFF` at the **rated** operating point (25%–125%
  load sweep) — these were already computed for the `DT1`/`DT2`
  heating-rise calculation but never included in the returned dict.
  Now available on the rated-load point (`tlast_result['points'][0]`
  — see the indexing note below).

### One thing I want to flag clearly: a real indexing gotcha
While testing this, I initially checked the wrong list index myself.
**`tlast_result['points']` is ordered by the FORTRAN's own loop
sequence, not by ascending load fraction:**

| Index | Load |
|---|---|
| `points[0]` | **Rated (100%)** |
| `points[1]` | 25% |
| `points[2]` | 50% |
| `points[3]` | 75% |
| `points[4]` | 125% |

I added this as an explicit, prominent note in `tlast()`'s docstring
so it can't be missed again.

---

## 3. What's still open

- **`RKZ` column 5** (ring current-distribution flag) — the scanned
  spec sheet's wording here was ambiguous (partly garbled OCR). The
  `O`/`H`/`J` values used in the code are confirmed directly from
  source; the sheet's own description of this column couldn't be
  fully reconciled and is noted as such in the field's help text.
- **No real `EINASYN` file or known-correct `AUSASYN` output has
  still been available** to validate actual numbers against — this
  session's documents were specification sheets (what each field
  means), not a filled-in example file with verified results. That
  remains the single biggest opportunity to fully validate this port.
