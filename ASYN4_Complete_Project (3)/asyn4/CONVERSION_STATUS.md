# ASYN4 FORTRAN -> Python conversion: status

Source: ASYN4.FOR, 12,835 lines, 110 subroutines/functions
(Dr.-Ing. W. Janssen, IEM Hannover / AEG-LDW Bremen, 1990/1996).

## Done (this round)
- `asyn4/state.py` — every COMMON block turned into a dataclass
  (BRULO, VORZEI, KEIL, CBSORT, NYINFO, MAXI, KENNLI, ESB, ART, LAERM,
  MNLINE, MNLRED, DATEN), bundled into one `MachineState` object that gets
  passed explicitly instead of relying on FORTRAN's implicit global COMMON.
  FORTRAN's 1-based array indices are preserved (arrays allocated one slot
  larger, index 0 unused) to avoid off-by-one bugs while porting the rest.
- `asyn4/winding_functions.py` — 13 standalone FUNCTIONs/SUBROUTINEs ported
  and smoke-tested: XSIS, XSISCH, XSIZ, XSIZB, XSIZ2, ZFAK, SGRUND, SNUE,
  RMAG, HEISEN, GGT, DXS1, LASTG, MGEGEN.

## Done (round 2)
- `asyn4/geometry_functions.py` — 11 more subroutines ported and
  smoke-tested: BRUVOR, COIL, NUT, PAKET, SPAK, VJ, UMRI, VLZFE, K1K2,
  NUE0, MNYQMP. Added `CPaket` (COMMON /CPAKET/) to `state.py`.
  **24 of 110 subroutines/functions ported so far.**

## Done (round 3)
- `asyn4/magnetic_circuit.py` — 2 more subroutines ported and smoke-tested:
  STD (stator geometry: yoke height, tooth widths, slot/end-winding
  leakage, resistance, updates st.keil), STEG (cage-rotor slot-bridge
  saturation curve, uses HEISEN). **26 of 110 subroutines/functions ported
  so far.**

## Done (round 4)
- `asyn4/rotor_functions.py` — 5 more subroutines ported and smoke-tested:
  SLFR (slip-ring rotor geometry/leakage/resistance, updates st.xs2sl),
  KIKR1 (stator form-coil skin-effect vs. drive frequency), KIKR (cage
  rotor bar skin-effect via the partial-conductor method - the biggest
  subroutine ported so far, ~300 lines, updates st.kennli.svpool), KIKRSL
  (slip-ring rotor skin-effect), KIKRSX (helper for KIKRSL). Added `Xs2sl`
  (COMMON /XS2SL/) to `state.py` - **you need to replace state.py again**,
  see state_NEW.py. **31 of 110 subroutines/functions ported so far.**

## Done (round 5)
- `asyn4/magnetics_solver.py` — 5 more subroutines ported and smoke-tested:
  SORTE (13 built-in magnetization curves - this is what makes HEISEN/STEG/
  VZ physically meaningful, previously they only worked with a caller-
  supplied placeholder table), VZ (tooth-voltage calc), PLATT (flattening
  factor), KZWERT (heating-coefficient lookup table), and **KREIS** - the
  core magnetic-circuit solver (~200 lines) that generates the full
  magnetization characteristic table (st.kennli.mkpool) the rest of the
  program is built around. **This completes the entire "Magnetic circuit"
  functional group.** No new COMMON blocks needed - state.py unchanged.
  **36 of 110 subroutines/functions ported so far.**

## Done (round 6)
- `asyn4/harmonics.py` — 2 more subroutines ported and smoke-tested:
  POPAZA (selects which harmonic pole-pair numbers need active damping,
  computes stator double-coupled leakage SD1/SD1R and per-harmonic
  auxiliary quantities, populates st.nyinfo) and Z2NYS (damped harmonic
  field impedance + complex damping factor, for all 4 rotor types -
  slip-ring, single-cage, double-cage common-ring, double-cage separate
  rings - updates st.esb when computing the fundamental branch). **This
  completes the "Harmonic-field & winding-factor engine" functional
  group.** No new COMMON blocks needed. **38 of 110 subroutines/functions
  ported so far.**

## Done (round 7)
- `asyn4/operating_point.py` — 3 more subroutines ported and tested:
  MKWERT (interpolates the KREIS magnetization curve at a given induced
  voltage), SORT (descending selection sort of pre-generated slip
  values, distinct from SORTE), and **SPKT** - the single-operating-point
  solver (~190 lines): given a slip value, iterates to a self-consistent
  induced voltage/saturation state and computes current, power, and
  torque summed over every damped harmonic. This is what SLAST/HOLAUF
  call once per point to build the full torque-speed curve. Updates
  st.esb/st.dp (synchronous reference point) and st.maxi (per-harmonic
  peak torque tracking) as side effects, matching the FORTRAN COMMON
  behaviour. Added `Dp` (COMMON /DP/) and `Svorga` (COMMON /SVORGA/) to
  `state.py` - **replace state.py again**, see state_NEW.py.
  **41 of 110 subroutines/functions ported so far. Test suite: 36/36
  passing** (tests/test_asyn4.py, run with `python run_tests.py`).

## Done (round 8)
- `asyn4/load_curves.py` — 2 more subroutines ported and tested: **SLAST**
  (slip-load calculation: sorts slip values, sweeps SPKT across all of
  them to build the full M-N torque-speed characteristic in st.mnline,
  finds the saddle/pull-up torque, and optionally repeats at reduced
  starting voltage into st.mnlred) and **HOLAUF** (run-up/starting
  calculation: integrates run-up time and adiabatic heating - stator,
  cage-bar bottom/top, ring - from the M-N curve and a counter-torque
  curve via MGEGEN). Printing was rewritten as plain console output
  instead of line-printer FORMATs / HP terminal escape codes (per the
  earlier "Open questions" decision); pass `verbose=False` to silence it.
  **This is the first true end-to-end integration test**, exercising
  nearly the whole pipeline in one run: KREIS -> KIKR -> POPAZA -> STEG ->
  SPKT -> SLAST -> HOLAUF - it runs to completion and produces
  internally-consistent output. One important discovery from this
  integration test, now documented in SPKT's docstring: **SPKT must be
  called once with `is1pkt=1` (the rated operating point) before any
  `is1pkt=0` sweep call**, or its reference reactance `st.dp.x1hp0` stays
  at 0 and every harmonic impedance downstream divides by zero. This
  wasn't obvious from the FORTRAN source in isolation, since COMMON
  silently carries state between calls.
  **43 of 110 subroutines/functions ported so far. Test suite: 37/37
  passing.**

## Done (round 9)
- `asyn4/special_points.py` — 6 more subroutines ported and tested: SK
  (finds the motor breakdown/pull-out slip, updates st.esb.skipp), SKGEN
  (generator-mode breakdown slip), **SKNY** (generates the slip-value list
  in st.svorga that SORT/SLAST sweep over - this was a missing link:
  earlier SLAST/HOLAUF tests had to populate st.svorga by hand; now the
  real generator exists), PRAXIS (practical starting/plugging-region
  operating point for revolving-field-magnet duty), DFMAG
  (revolving-field-magnet braking behaviour across S=1..5, combines SPKT
  + PRAXIS), and PRUEF (validates the read-in machine data for internal
  consistency - returns error/warning messages instead of FORTRAN WRITEs;
  first subroutine that meaningfully exercises st.daten/st.art.rkz as a
  whole record, ahead of LESEN existing to populate it from a real file).
  No new COMMON blocks needed. **49 of 110 subroutines/functions ported
  so far. Test suite: 44/44 passing.**
  Not ported this round: GENERA (generator-mode performance report across
  5 load points - large, mostly screen/report formatting, lower priority
  since SK/SKGEN/SPKT/PRAXIS - its actual computational dependencies -
  are now all done).

## Done (round 10)
- `asyn4/insulation.py` — first piece of the insulation-design group:
  **ALD** (stator coil-end overhang length, Jungbluth formulas). The
  FORTRAN original is an interactive REPL (computes, warns if angles are
  out of range, reads a terminal Y/N, loops on a new Z if requested) -
  ported as a single-pass function returning `in_range` instead of
  blocking on terminal input; a plain Python `while` loop around it
  reproduces the original REPL exactly if ever wanted. **50 of 110
  subroutines/functions ported so far. Test suite: 46/46 passing.**

  **Size flag on what's left in this group:** ISOSTD, ISOTR, ISOLFR,
  NAR2, NARISS, WVAR, FSPULE, CUFL, TEILEN, SUCH1-6, ZPLAN (source
  9961-12836, ~2900 lines) is the single largest remaining functional
  group - bigger than everything ported so far combined. It's mostly
  form-coil/random-wound insulation buildup and wire-gauge selection
  logic (manufacturing detail) rather than electromagnetic performance
  calculation, and several of its subroutines are interactive terminal
  REPLs needing the same treatment as ALD above. Continuing through it
  is straightforward in the same style as everything so far - noting the
  size here so it's visible, not as a reason to stop.

## Done (round 10)
- `asyn4/insulation.py` — first 4 pieces of the insulation-design group:
  ALD (stator coil-end overhang length), **NARISS** (minimum slot width/
  height from empirical insulation tables - pure computation, updates
  st.spkopf/st.ck211s), **FSPULE** (two-layer full-form-coil end-winding
  geometry per Schorch W01481), **ZGF** (combines FSPULE symmetric +
  asymmetric to get both overhang lengths), and **WVAR** (searches the
  standard wire-gauge table for the best-fit conductor configuration).
  NARISS/FSPULE are pure computation (no interactive I/O in the
  original); ALD and WVAR were interactive REPLs in the FORTRAN and were
  ported the same way as before - single deterministic pass, with the
  yes/no/override decisions exposed as optional parameters instead of
  blocking reads. Added `Ctest`, `Spkopf`, `Ck211s` (COMMON /CTEST/,
  /SPKOPF/, /CK211S/) to `state.py` - **replace state.py again**, see
  state_NEW.py. **54 of 110 subroutines/functions ported so far. Test
  suite: 50/50 passing.**

  Still not ported in this group: ISOSTD (the orchestrator that calls
  ALD/NARISS/WVAR/ZGF together and writes the MAISO/KASIMIRS transfer
  file - depends on all 4 of the above, now portable), ISOTR, ISOLFR,
  NAR2, CUFL, TEILEN, SUCH1-6, ZPLAN.

## Done (round 11)
- `asyn4/insulation.py` — **ISOSTD** added, the orchestrator that ties
  ALD/NARISS/WVAR/ZGF together: given either a target conductor size, a
  target slot size, or neither, works out the other, sizes the
  press-finger, computes coil-head clearance defaults, and converges on
  the stator overhang length. Tested across both its "conductor given"
  and "slot given" branches, plus its AZWEIG=0 validation error. The
  MAISO/KASIMIRS transfer-file write was kept as data returned in the
  result dict rather than writing 'FILE.ISO' directly (same optional-
  file-write pattern as ZGF's 'WKFISO'). No new COMMON blocks needed.
  **55 of 110 subroutines/functions ported so far. Test suite: 53/53
  passing.**

  Still not ported in this group: ISOTR, ISOLFR, NAR2, CUFL, TEILEN,
  SUCH1-6, ZPLAN.

## Done (round 12)
- `asyn4/insulation.py` — 6 more pieces added: **ISOTR** (random-wound
  stator winding design - copper fill-factor suggestions across
  candidate parallel-conductor counts, finalized once a count is chosen),
  **CUFL** (copper area per slot, standalone), **TEILEN** (splits uneven
  turns/slot counts across a zone's slots), **SUCH1-6** (collapsed into
  one generic `such(q, svk, loch)` using itertools.product for the exact
  same brute-force enumeration order as the original 6 near-identical
  subroutines), and **ZPLAN** (symmetric-pattern checker, implemented as
  an O(Q) periodicity check instead of literally building the FORTRAN's
  6*Q-length working arrays - same result, cheaper). No new COMMON blocks
  needed. **62 of 110 subroutines/functions ported so far** (SUCH1-6
  counted as the 6 originals they replace). **Test suite: 60/60 passing.**

  Still not ported in this group: ISOLFR (cage/slip-ring rotor bar
  insulation), NAR2 (a second slot-drawing variant). Once those two are
  done, the entire ~2900-line insulation-design group is complete.

## Done (round 13)
- `asyn4/insulation.py` — final 2 pieces added: **NAR2** (slot layout for
  a slip-ring rotor bar, the rotor-side counterpart to NARISS - pure
  computation, updates st.ck211l) and **ISOLFR** (slip-ring rotor winding
  design - searches the standard conductor tables via NAR2, with the
  FORTRAN's 3 interactive prompts exposed as optional parameters:
  reversible_duty, ntl_override, tpar_override). Added `Ck211l` (COMMON
  /CK211L/) to `state.py` - **replace state.py again**, see state_NEW.py.
  Also added a bounds guard in ISOLFR's table search that raises a clear
  `ValueError` instead of an `IndexError` when no standard conductor size
  fits the given slot (the FORTRAN would have read off the end of its
  DATA arrays in that case - undefined behaviour there, a clean error
  here). **64 of 110 subroutines/functions ported so far** (66 counting
  SUCH1-6 as the 6 originals). **Test suite: 64/64 passing.**

  **THE ENTIRE ~2900-LINE "INSULATION DESIGN" FUNCTIONAL GROUP IS NOW
  COMPLETE**: ALD, NARISS, FSPULE, ZGF, WVAR, ISOSTD, ISOTR, CUFL,
  TEILEN, SUCH1-6, ZPLAN, NAR2, ISOLFR - 15 subroutines across 6 rounds,
  all tested. This was the single largest remaining functional group in
  the program.

## Done (round 14)
- `asyn4/special_points.py` — final 2 pieces added: **SPKTSP** (locked-
  rotor/S=1 "special point" solver, like SPKT but with a KKA damping
  factor scaling harmonic contributions differently near vs. far from
  the fundamental pole-pair number) and **GENERA** (generator-mode
  performance: 5 partial-load points via SPKT + the generator breakdown
  point via SKGEN+SPKT; all screen-position/printer output replaced with
  a returned dict). Caught and hardened one real fragility during
  testing: if GENERA's rated-load point never converges, the original
  FORTRAN would silently use uninitialized memory for a later
  calculation - this port raises a clear `RuntimeError` instead. No new
  COMMON blocks needed (COMMON /PINFO/ and /CNLAST/, only read by other
  not-yet-ported report routines, are returned as plain dict fields for
  now rather than promoted to state.py dataclasses prematurely).
  **68 of 110 subroutines/functions ported (70 counting SUCH1-6
  individually). Test suite: 66/66 passing.**

  **THE "PERFORMANCE CALCULATIONS" FUNCTIONAL GROUP IS NOW COMPLETE.**

## Done (round 15)
- `asyn4/lesen.py` — **LESEN**, the EINASYN input-file reader, ported as
  a best-effort reconstruction of its 19 fixed-column FORTRAN READ/
  FORMAT statements (documented card-by-card in the module docstring,
  each traceable to a source line range). Populates st.daten/st.art.rkz/
  st.materi, then runs ISOSTD/ISOTR (stator) and ISOLFR (slip-ring rotor)
  exactly as the FORTRAN does at the end of reading. Added `Materi`
  (COMMON /MATERI/) to `state.py` - **replace state.py again**, see
  state_NEW.py. The ~800-line screen/printer-report tail that follows
  the data reads in the original (source 2055-3023) was NOT ported here
  - it's pure display formatting with no new computation, and belongs to
  the not-yet-ported "output/reporting" group instead.

  **Tested via a synthetic file** (tests/test_asyn4.py generates one
  matching the documented column layout) confirming the parser lands on
  the right columns - this is a **self-consistency test, not validation
  against a real EINASYN file**, since none has been available during
  this port. `run_tests.py` gained a `tmp_path` fixture shim to support
  this.

  ⚠️ **UNITS DISCOVERY (round 15)** - later corrected, then corrected
  again with actual proof; see round 21 above for the final, source-
  backed answer: EINASYN's raw fields are in millimetres/kW, and LESEN
  converts them to metres/W via an explicit `UMWANDLUNG IN SI-EINHEITEN`
  block (K=1000, M=1e6) that this round initially missed. If you're
  reading this file top-to-bottom, skip ahead to round 21's write-up
  rather than trusting this round's version of events.

  **73 of 110 subroutines/functions ported (75 counting SUCH1-6
  individually). Test suite: 68/68 passing.**

## Done (round 16)
- `asyn4/reports.py` — new module starting the "output/reporting"
  functional group, with the 2 pieces that have real logic in them:
  **MWINFO** (sheet-steel grade lookup + solid-vs-segmented lamination
  decision - the report-formatting parts replaced with a returned dict,
  per the established pattern) and **SYNMOM** (synchronous/cusp torque
  analysis for cage rotors - genuine physics, not just reporting:
  determines whether a parasitic synchronous torque could stall the
  motor during starting, and at what running speeds it would recur).
  No new COMMON blocks needed. **75 of 110 subroutines/functions ported
  (77 counting SUCH1-6 individually). Test suite: 71/71 passing.**

  **Size flag: TLAST** (source 3906-5007, ~1100 lines) is now clearly
  the single largest remaining subroutine in the entire program - bigger
  than KIKR. It's the motor-mode counterpart to GENERA (no-load + 5
  partial-load points + breakdown + locked-rotor point), interleaved
  with report writes and 3 more COMMON blocks not yet modelled
  (CK211T, PINFO, CNLAST - only consumed by other unported report
  routines). All of its computational dependencies are already ported,
  so it's ready to tackle, just sizeable - flagged rather than rushed.

  Still not ported: TLAST, KDRUCK, PDRUCK, LISTE5, KOPP1/4/6/7/8,
  KOK211, KOPASK, KOPKFE, XKEY (mostly pure report formatting drawing on
  already-computed values), LESEN's deferred screen-report tail.

## Done (round 17)
- `asyn4/reports.py` — **KDRUCK** added: "KST-Zusatzinformation"
  cooling/manufacturing info sheet. Two genuine computation blocks:
  cooling-air entry areas (with the IKLINK "which cooling path is the
  bottleneck" decision) and lamination-stack cross-section areas
  (stator/rotor iron vs. copper split). Report-file writing ('KSTDATEN')
  replaced with a returned dict, same pattern as everywhere else. No new
  COMMON blocks needed. **76 of 110 subroutines/functions ported (78
  counting SUCH1-6 individually). Test suite: 72/72 passing.**

  Still not ported: TLAST (~1100 lines, largest remaining subroutine -
  see round 16's size flag), PDRUCK, LISTE5, KOPP1/4/6/7/8, KOK211,
  KOPASK, KOPKFE, XKEY, LESEN's deferred screen-report tail.

## Done (round 18) — TLAST, the largest subroutine in the program, is done
- `asyn4/reports.py` — 3 more pieces added: **KKWERT** (starting-current
  correction factor lookup, TLAST's dependency), **XKEY** (leakage-
  reactance breakdown + time constants at standstill, differs materially
  for single vs. double-cage rotors via a complex parallel-impedance
  combination), and **TLAST** itself - the ~1100-line motor-mode duty-
  cycle calculation, the single largest subroutine in the entire
  program. Computes the no-load point, 5 partial-load points, the
  breakdown (Kipp) point, the theoretical locked-rotor point, the
  KKWERT-corrected locked-rotor point (via SPKTSP), and XKEY's time
  constants, then dispatches to GENERA (generator duty) or DFMAG
  (revolving-field-magnet duty) exactly as the FORTRAN does based on
  st.art.rkz's 2nd character. The not-yet-ported KOPASK report call was
  exposed as an optional `kopask_callback` hook rather than silently
  dropped. **Tested end-to-end**: the full pipeline (no-load through
  locked-rotor plus both dispatch branches) runs clean and produces
  internally-consistent results. No new COMMON blocks needed (CK211T/
  PINFO/CNLAST - only consumed by other not-yet-ported report routines -
  are returned as plain dict fields instead of being promoted to
  state.py dataclasses prematurely).

  **80 of 110 subroutines/functions ported (82 counting SUCH1-6
  individually). Test suite: 76/76 passing.**

  Still not ported: PDRUCK, LISTE5, KOPP1/4/6/7/8, KOK211, KOPASK,
  KOPKFE, LESEN's deferred screen-report tail — all of these are pure
  report/file-writing on already-computed values, no remaining
  substantial computation anywhere in the program.

## Done (round 19) — CSV and PDF report export
- `asyn4/report_export.py` — new module, NOT a FORTRAN port (this
  replaces what PDRUCK/LISTE5/KOPP1-8/KOK211/KOPASK/KOPKFE would have
  written as HP-PCL printer streams). Two entry points:
  - `write_csv_report(...)` — writes one CSV file per available section
    (partial-load performance, breakdown/locked-rotor summary,
    magnetization curve, synchronous-torque breakdown, cooling areas,
    material info). Every argument is optional; pass only whichever
    results you've computed.
  - `write_pdf_report(...)` — a single formatted multi-page PDF with a
    title page (machine ID from LESEN's header fields), styled tables
    for every section, page breaks between major sections. Built with
    reportlab (already available in this environment).
  Both take the result dicts already returned by TLAST/SYNMOM/KDRUCK/
  MWINFO plus the raw KREIS magnetization-curve array (st.kennli.mkpool)
  directly - no new data model needed. **Tested end-to-end**, including
  a real example report (6 CSVs + 1 PDF) generated from the same
  synthetic test machine used throughout this test suite - see
  `example_report/` for a live sample of what the output looks like.

  **80 of 110 FORTRAN subroutines/functions ported (82 counting SUCH1-6
  individually), plus this new non-FORTRAN reporting module. Test suite:
  79/79 passing.**

## Done (round 20) — PDRUCK and KOPASK; every subroutine with real
## computational content in ASYN4.FOR is now ported
- `asyn4/reports.py` — 2 more pieces added: **PDRUCK** (Pruefeld-
  Information test-field data sheet: torque/current at rated and
  reduced supply voltage with SPKTSP/KKA fallback, plus the full
  15-point short-circuit characteristic curve from 10%-100% of rated
  voltage - a genuine test-certificate style voltage sweep) and
  **KOPASK** (builds a 50-point torque/current-vs-slip curve for
  coupling to a separate piston-compressor load-matching program). Both
  tested end-to-end against the same synthetic machine fixture used
  throughout, chained after a full TLAST run (which is what latches the
  st.esb reference values these two depend on, matching the FORTRAN's
  COMMON /PINFO/ dependency chain).

  **Scope decision, documented in the module docstring**: KOK211,
  KOPKFE, KOPP1, KOPP4, KOPP6, KOPP7, KOPP8, and LISTE5 were reviewed
  and found to be pure data-transfer bridges to SEPARATE FORTRAN
  programs not included in ASYN4.FOR (a noise/vibration calculator, the
  "K211" mechanical design tool, the "KFE" loss-calculation program).
  Their only logic is unit conversion before writing another program's
  input format - no physics or decisions to port, and the receiving
  programs aren't available to validate against. Left unported by
  design; happy to port specific ones on request if those downstream
  programs are ever needed.

  **82 of 110 FORTRAN subroutines/functions ported (84 counting SUCH1-6
  individually). Test suite: 81/81 passing.**

  **This closes out every subroutine in ASYN4.FOR with genuine physics,
  engineering logic, or decision-making.** What's left in the whole
  project is: LESEN's screen-report tail (pure display, no computation),
  the plotting driver (Calcomp/HP-GL -> matplotlib), and PROGRAM ASYN4
  itself (the main driver wiring everything together).

## Done (round 21) — 🎉 PROGRAM ASYN4 itself is ported: the pipeline runs
## end-to-end from a real EINASYN-format file for the first time
- `asyn4/main.py` — **`run_asyn4()`**, the main driver (source 1-869).
  Wires together, in the FORTRAN's own order: LESEN -> BRUVOR ->
  geometry pre-calc -> PAKET -> STD -> KIKR1 -> SLFR+KIKRSL (slip-ring)
  or KIKR (cage) -> KREIS -> STEG -> POPAZA -> NUT -> mass/inertia calc
  -> TLAST -> MAGZUG -> KDRUCK -> SLAST -> HOLAUF -> PDRUCK -> SYNMOM.
  Two inline FORTRAN blocks that were never their own subroutine are
  extracted as testable helpers: `_geometry_precalc` (stack lengths,
  winding factor, rotor back-iron radius) and
  `_stator_form_coil_correction` (K1R1W skin-effect factor).
- `asyn4/magzug.py` — **MAGZUG/MAGZU1** (axial magnetic pull), a small
  self-contained piece needed by the main driver's KDRUCK call.

- **This is the first time any real EINASYN-format file has been run
  through the pipeline**, and it was extremely productive: integration
  testing at this level found and fixed 8 real bugs that no amount of
  function-level testing with hand-built MachineState objects could
  have caught, because they only manifest when data actually flows
  through LESEN from a file:

  1. **ISOTR never finalized for random-wound stators** - LESEN called
     ISOTR with tpar=None (get suggestions) but never made the follow-
     up call to actually finalize MTL/QCU1, silently leaving QCU1 at 0.
     Fixed: LESEN now auto-selects the first suggestion with a
     practical fill factor <=0.68, mirroring a sensible non-interactive
     answer to the FORTRAN's "ANZAHL PARALLELER TEILLEITER?" prompt.
  2. **Fabricated H32/BR32 zeroing block** - an `elif ilfr == 2: h32 =
     h42 = 0.0; br32 = br42 = 0.0` block in lesen.py did not correspond
     to anything in the real FORTRAN source (confirmed by the very next
     line using h32*br32, which would be nonsensically zero if that
     block were real) and was removed.
  3. **Unbound `r2w`** in run_asyn4's cage-rotor branch (legitimately
     unused for cage rotors, but needed a default value).
  4. **Unbound `g1`/`a1eff`** in TLAST - these were nested inside `if kz
     > 0.0:` but used unconditionally in the `else` branch below;
     FORTRAN source confirms they're computed unconditionally.
  5. **K41-K4AR2 (heating time constants) never computed** - this
     entire block (source 2868-2897, "ERWAERNUNGSKONSTANTE") lives in
     PROGRAM ASYN4's own body, past where this port's LESEN boundary
     was drawn, and was missed entirely. Now computed by a new
     `_finish_lesen_postprocessing` helper in main.py.
  6. **PRBG0 default (Nurnberg friction-loss formula) never applied** -
     same location, same fix.
  7. **PRUEF (data validation) was ported back in round 9 but never
     actually called from anywhere** - now wired into run_asyn4()
     immediately after LESEN, with a proper `iabort` short-circuit.
  8. **LUE was read but discarded** (not stored on st.daten) - added
     `Daten.lue` and fixed the assignment. Also added `Daten.k41`
     through `Daten.k4ar2`.

- **The units question, resolved for real this time.** This port's
  understanding of EINASYN's units has now been wrong twice and is
  corrected a third time, with actual source evidence this time:
  searching for `UMWANDLUNG IN SI-EINHEITEN` in ASYN4.FOR (source
  2689-2774) reveals LESEN defines `K = 1.E3` and `M = 1.E6` and
  systematically divides every length field by K, area field by M, and
  multiplies every conductivity field by M, immediately after the raw
  card reads and before anything else touches the data. This is
  conclusive: **EINASYN's raw fields are in millimetres (mm), mm^2, and
  a conductivity unit that becomes S/m after *1e6, and PN/PRBG0 are in
  kilowatts** - LESEN converts everything to metres/W/S-per-metre
  before the rest of the program (everything ported in rounds 1-20,
  all along) ever sees it. This is now implemented in lesen.py, along
  with the two special-case transforms that immediately follow it in
  the source: **UMWANDLUNG VON DOPPELNUT IN HAMMERKOPF** (double-slot
  to hammerhead-slot geometric model conversion) and **DRUCKGUSS-
  LAEUFER** (die-cast rotor downgrade from double-cage to single-cage
  when no ring-overhang length is configured). Also added: LM1/LM2
  (source LW+L) and LSTAB/LSTABA (source 2*LUE+L / 2*LUEA+L)
  computation, and the AN2/NUTF2 special-slot-shape HS2 override.

- **New `_CardWriter` utility in lesen.py**, a companion to
  `_CardReader` sharing the same skip-N-lines logic, so synthetic test
  files can never drift out of sync with the reader the way an early
  hand-rolled version did during this round's debugging (a blank-line
  miscount between two cards silently shifted every field after it by
  one column - exactly the class of error fixed-format FORTRAN parsing
  is notorious for, and exactly why this integration test was worth
  doing even without a real reference file).

- **Honest open item found, not swept under the rug**: with valid
  (if synthetic) geometry, `SK`'s breakdown-torque search converges to
  a large, wrong-signed value (~-230 Nm against a ~6-10 Nm working
  range) for this test machine. This may be a real bug in SK's slip-
  scanning logic, or - more likely given everything else in this round
  - another artifact of synthetic geometry landing in a region SK's
  simple local-search algorithm wasn't designed for. Flagged here
  rather than either hand-tuning test data until it disappears
  (masking a possible real bug) or claiming it's fixed without knowing
  why. **This needs a real EINASYN file to resolve properly.**

- Test suite: **83/83 passing**, including two new top-level
  integration tests: `test_run_asyn4_full_pipeline_from_file` (the
  first test in this whole port to go file-in to results-out) and
  `test_run_asyn4_rejects_invalid_data` (confirms PRUEF's validation is
  actually reachable now). The existing LESEN parser test was also
  fixed to expect properly SI-converted values instead of raw mm/kW
  numbers.

- **87 of 110 FORTRAN subroutines/functions ported (89 counting
  SUCH1-6 individually).** Every subroutine with real computational
  content in ASYN4.FOR - and now PROGRAM ASYN4 itself - is ported and
  integration-tested.

## Done (round 22) — THE PLOTTING DRIVER IS DONE. EVERY FUNCTIONAL GROUP
## IN THE ENTIRE PROGRAM IS NOW PORTED.
- `asyn4/plotting.py` — replaces the ~1000-line Calcomp/HP-GL pen-
  plotter subsystem (PLOTT, FRAME, HEAD, GITTER, and the low-level
  primitives AXIS/SYMBOL/NUMBER/SCALE/DELTI/LINE/DASHL/DASHP/CIRCL/
  PLOT/PLOTS/NEWPEN/NEWPLOT/BMASKE/WHERE/FACTOR) with a single
  matplotlib chart function, NOT a line-by-line port. Tracing every
  caller of the plotting system shows PLOTT is called exactly ONCE in
  the whole program (right after SLAST), and despite the machinery's
  size it exists to draw exactly one chart: the combined torque-speed /
  current-speed characteristic (FRAME's IACHS=3 mode). Reimplementing
  ~20 generic Calcomp primitives that would only ever be exercised by
  this one caller would have been pure busywork - matplotlib already
  does axis scaling and text rendering better than 1990s HP-GL
  emulation. plot_torque_speed_curve() reproduces exactly what the
  original chart showed (read directly from FRAME's own axis/curve
  setup, source 5306-5334): X=N/N1, left Y=M/MN (motor + load-torque
  curves), right Y=I1/I1N, with a reduced-voltage overlay when KU1/
  XNETZ are configured - see example_report/torque_speed_chart.png for
  a live sample generated from the same synthetic test machine used
  throughout this suite.
- Wired into run_asyn4() via a new plot_path parameter, positioned
  exactly where PLOTT is called in the original main-program flow
  (immediately after SLAST). Also surfaced (while wiring the test)
  that RKZ(3:3) is the flag controlling whether the whole slip-curve/
  plot block runs at all (source "IF (RKZ(3:3).EQ.' ') GOTO 500") -
  worth knowing if a real EINASYN file's slip curve/chart isn't
  appearing when expected.
- The chart also visually confirms the already-flagged SK/breakdown-
  torque anomaly from round 21 (a sharp, unphysical spike near
  synchronism) - not a new bug, just the existing open item becoming
  visible now that there's a chart to look at. Still needs a real
  EINASYN file to resolve.
- Test suite: 86/86 passing, including
  test_run_asyn4_generates_plot_when_requested (confirms the plot
  wiring end-to-end from a file) and 2 direct plotting.py tests.
- 89 of 110 FORTRAN subroutines/functions ported (91 counting SUCH1-6
  individually) - and every one of the 10 functional groups in the
  original roadmap is now complete, including the main program and the
  plotting driver. What remains unported by design (not by gap): the
  external-tool data-transfer bridges (KOK211, KOPKFE, KOPP1/4/6/7/8,
  LISTE5 - see round 20's scope note) and LESEN's pure-display
  cover-sheet print.

## Not started yet (roadmap, in the order the FORTRAN calls them)
1. **I/O layer** — LESEN: fully done as of round 21, including the
   SI-unit conversion, hammerhead-slot and die-cast-rotor special
   cases, and PRUEF validation wiring. The only thing not ported is the
   pure-display "AUSDRUCK DER EINGELESENEN DATEN" cover-sheet print
   (source ~2816 onward) - no remaining computation there.
2. **Core geometry & winding pre-calc** — all done (BRUVOR, PAKET/SPAK,
   NUT, COIL, STD, SLFR, KIKR/KIKR1/KIKRSL/KIKRSX).
3. **Magnetic circuit** — all done (SORTE, VZ, PLATT, KZWERT, KREIS, STEG,
   RMAG/HEISEN/VJ/UMRI/VLZFE).
4. **Harmonic-field & winding-factor engine** — all done (NUE0, MNYQMP,
   K1K2, ZFAK/XSIZ family, POPAZA, Z2NYS).
5. **Insulation design** — all done (ISOSTD, ISOTR, ISOLFR, NAR2, NARISS,
   WVAR, ALD, FSPULE, ZGF, CUFL, TEILEN, SUCH1-6, ZPLAN).
6. **Performance calculations** — all done (SPKT, MKWERT, SK, SKGEN,
   SKNY, PRAXIS, PRUEF, DFMAG, SPKTSP, GENERA).
7. **Load / starting / heating** — all done (SLAST, HOLAUF, SORT).
8. **Output/reporting** — MWINFO, SYNMOM, KDRUCK, KKWERT, XKEY, TLAST,
   PDRUCK, KOPASK, MAGZUG all done - this closes out every subroutine
   with real computational content anywhere in ASYN4.FOR. Intentionally
   not ported (external-tool data-transfer bridges to separate programs
   not included in this codebase, see round 20): KOK211, KOPKFE,
   KOPP1/4/6/7/8, LISTE5.
9. **Plotting driver** — PLOTS/PLOT/AXIS/SYMBOL/NUMBER/SCALE/DELTI/LINE/
   DASHL/DASHP/CIRCL/FRAME/HEAD/GITTER/BMASKE/WHERE/FACTOR/NEWPEN/NEWPLOT/
   PLOTT — this whole family is a Calcomp/HP-GL pen-plotter emulation.
   Not yet ported; would replace it with a small `PlotContext` class
   backed by matplotlib. **This is now the only unstarted functional
   group in the entire program.**
10. **`PROGRAM ASYN4` itself** — **DONE** as of round 21
    (`asyn4/main.py`'s `run_asyn4()`). Runs end-to-end from a real
    EINASYN-format file, minus the plotting calls (not ported, see #9)
    and the external-tool bridges (out of scope, see #8).

## Open questions I'll need from you as we go
- **A sample real `EINASYN` file is now more valuable than ever.**
  Round 21's integration testing (with synthetic-but-format-correct
  data) found and fixed 8 real bugs unreachable by any amount of
  function-level testing - but it also surfaced one genuine open
  question a real file could resolve outright: `SK`'s breakdown-torque
  search produces a large, wrong-signed result for the current
  synthetic test machine (~-230 Nm against a ~6-10 Nm working torque
  range) - this may be a real bug in SK's slip-scanning algorithm, or
  an artifact of synthetic geometry outside its intended operating
  envelope. A real file (or a known-correct AUSASYN output to compare
  against) would settle this and validate every numeric result in the
  package for the first time, not just its internal consistency.
- Whether the terminal-control-code screens (`ESC[...H` cursor positioning
  written straight to stdout) should be kept, or replaced with normal
  console/log output — this port defaults to normal output throughout.
- Whether you want the printer report layouts preserved character-for-
  character (T-column FORMATs) or just the same data in a normal
  Python/CSV/PDF report (the CSV/PDF export in `report_export.py`
  already covers the latter).

## How to run the full pipeline
```python
from asyn4.main import run_asyn4

result = run_asyn4("path/to/your/EINASYN", verbose=True)
if result.get("iabort") == 1:
    for is_error, text in result["pruef_messages"]:
        print(("ERROR: " if is_error else "WARN: ") + text)
else:
    tl = result["tlast_result"]
    print("Rated point:", tl["points"][3])
    print("Breakdown:", tl["kipp"])
```

## How to run individual functions (as in earlier rounds)
```
cd asyn4
python3 -c "
from asyn4.state import MachineState
from asyn4 import winding_functions as wf
st = MachineState()
st.vorzei.bs1di1 = 0.000001
print(wf.xsiz(1, 3, 36, st))
"
```

