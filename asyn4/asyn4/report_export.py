"""
Report export - turns the result dicts already returned by the ported
calculation functions (TLAST, GENERA, SYNMOM, KDRUCK, MWINFO, and the
KREIS magnetization curve in st.kennli.mkpool) into two concrete
deliverables: a CSV file (one per logical table, or a combined multi-
section CSV) and a single formatted PDF report.

This is NOT a port of any FORTRAN subroutine - none of PDRUCK, LISTE5,
KOPP1/4/6/7/8, KOK211, KOPASK, KOPKFE (the FORTRAN's own line-printer
report writers, still unported - see CONVERSION_STATUS.md) work this
way; they format directly to HP PCL escape-code printer streams. This
module is a new, from-scratch replacement built on the same underlying
data those FORTRAN routines would have printed, targeting normal
CSV/PDF output instead.

Usage
-----
    from asyn4 import report_export as rex

    rex.write_csv_report(tlast_result, mkpool_rows, path_prefix="machine1")
    rex.write_pdf_report(header_info, tlast_result, synmom_result,
                          kdruck_result, mwinfo_result, mkpool_rows,
                          "machine1_report.pdf")

Every argument other than the output path is optional (pass None to
skip that section) so this can be used with whatever subset of results
you've actually computed.
"""
from __future__ import annotations
import csv
import os

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, PageBreak)


# --------------------------------------------------------------------- #
# Shared table-building helpers
# --------------------------------------------------------------------- #

_LOAD_LABELS = ["1/4", "2/4", "3/4", "4/4 (rated)", "5/4"]


def _tlast_points_table(tlast_result: dict):
    """
    Builds a list-of-lists table (header + one row per converged partial-
    load point) from a TLAST result dict's 'points' list, reordering them
    from computation order (Rated, 25%, 50%, 75%, 125%) to ascending order.
    """
    header = ["Load", "P_mech (kW)", "Eta", "cos(phi)", "I1 (A)",
              "Slip (%)", "Pcu1 (kW)", "Pcu2 (kW)", "M (Nm)"]
    rows = [header]

    # Reorder from (Rated, 1/4, 2/4, 3/4, 5/4) -> (1/4, 2/4, 3/4, Rated, 5/4)
    raw_points = tlast_result.get("points", [])
    if len(raw_points) == 5:
        ordered_points = [raw_points[1], raw_points[2], raw_points[3], raw_points[0], raw_points[4]]
    else:
        ordered_points = raw_points

    for j, pt in enumerate(ordered_points, start=1):
        label = _LOAD_LABELS[j - 1]
        if not pt.get("converged", False):
            rows.append([label, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a"])
            continue
        rows.append([
            label,
            f"{pt['pmech'] / 1000.0:.2f}",
            f"{pt['eta']:.4f}",
            f"{pt['cosph']:.3f}",
            f"{pt['i1']:.2f}",
            f"{pt.get('s', 0.0) * 100.0:.3f}" if "s" in pt else "",
            f"{pt['pcu1'] / 1000.0:.3f}",
            f"{pt['pcu2'] / 1000.0:.3f}",
            f"{pt['m']:.1f}",
        ])
    return rows


def _tlast_summary_table(tlast_result: dict):
    rows = [["Quantity", "Value"]]
    kipp = tlast_result.get("kipp")
    theo = tlast_result.get("theoretical")
    corr = tlast_result.get("corrected")
    if kipp:
        rows.append(["Breakdown slip (Kipp)", f"{kipp['skipp'] * 100.0:.3f} %"])
        rows.append(["Breakdown torque", f"{kipp['m']:.1f} Nm"])
    if theo:
        rows.append(["Locked-rotor torque (theoretical)", f"{theo['m']:.1f} Nm"])
        rows.append(["Locked-rotor current (theoretical)", f"{theo['i1']:.1f} A"])
    if corr:
        rows.append(["Locked-rotor torque (corrected, KKA="
                      f"{corr['kka']:.2f})", f"{corr['m']:.1f} Nm"])
        rows.append(["Locked-rotor current (corrected)", f"{corr['i1']:.1f} A"])
    if tlast_result.get("tk") is not None:
        rows.append(["Time constant TK", f"{tlast_result['tk']:.4f} s"])
    if tlast_result.get("thaupt") is not None:
        rows.append(["Time constant THAUPT (rotor)", f"{tlast_result['thaupt']:.4f} s"])
    xk = tlast_result.get("xkey_result")
    if xk:
        rows.append(["Time constant TSTREU (leakage)", f"{xk['tstreu']:.4f} s"])
    return rows


def _mkpool_table(mkpool_rows, eanz: int):
    """
    mkpool_rows: st.kennli.mkpool, a 1-indexed np.ndarray(51,19) as
    populated by KREIS. eanz: st.kennli.eanz (number of valid rows).
    """
    header = ["UI (V)", "Imy (A)", "Bp (T)", "Pfe1 (kW)", "Bjoch1 (T)",
              "Bjoch2 (T)", "Bz1max (T)", "Bz2max (T)"]
    rows = [header]
    for j in range(1, eanz + 1):
        r = mkpool_rows[j]
        rows.append([
            f"{r[1]:.1f}", f"{r[2]:.2f}", f"{r[4]:.3f}", f"{r[7] / 1000.0:.3f}",
            f"{r[11]:.3f}", f"{r[12]:.3f}", f"{r[9]:.3f}", f"{r[10]:.3f}",
        ])
    return rows


def _synmom_table(synmom_result: dict):
    header = ["Pole-pair number NY", "Msy/Mn (standstill)"]
    rows = [header]
    for nyq, msy in synmom_result.get("breakdown", []):
        rows.append([f"{nyq:.0f}", f"{msy:.4f}"])
    return rows


def _kv(d: dict, keys_labels: list):
    """Builds a [label, value] table from a dict, skipping missing keys."""
    rows = [["Quantity", "Value"]]
    for key, label, fmt in keys_labels:
        if key in d and d[key] is not None:
            val = d[key]
            rows.append([label, fmt(val) if callable(fmt) else f"{val}"])
    return rows


# --------------------------------------------------------------------- #
# CSV export
# --------------------------------------------------------------------- #

def write_csv_report(tlast_result: dict | None = None,
                      mkpool_rows=None, eanz: int = 0,
                      synmom_result: dict | None = None,
                      kdruck_result: dict | None = None,
                      mwinfo_result: dict | None = None,
                      path_prefix: str = "asyn4_report"):
    """
    Writes one CSV file per available section, named
    '{path_prefix}_<section>.csv'. Returns the list of paths written.
    Any argument left as None (or eanz=0) skips that section - so this
    can be called with only whichever results you've actually computed.
    """
    written = []

    if tlast_result is not None:
        path = f"{path_prefix}_partial_load.csv"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(_tlast_points_table(tlast_result))
        written.append(path)

        path = f"{path_prefix}_summary.csv"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(_tlast_summary_table(tlast_result))
        written.append(path)

    if mkpool_rows is not None and eanz > 0:
        path = f"{path_prefix}_magnetization_curve.csv"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(_mkpool_table(mkpool_rows, eanz))
        written.append(path)

    if synmom_result is not None:
        path = f"{path_prefix}_synchronous_torque.csv"
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(_synmom_table(synmom_result))
        written.append(path)

    if kdruck_result is not None:
        path = f"{path_prefix}_cooling.csv"
        rows = _kv(kdruck_result, [
            ("iklink", "Dominant cooling path (IKLINK)", None),
            ("ar1", "Stator yoke iron area AR1 (m2)", None),
            ("az1", "Stator tooth iron area AZ1 (m2)", None),
            ("ar2", "Rotor yoke iron area AR2 (m2)", None),
            ("az2", "Rotor tooth iron area AZ2 (m2)", None),
            ("k1_at_s1", "Skin-effect factor K1 at S=1", None),
        ])
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(rows)
        written.append(path)

    if mwinfo_result is not None:
        path = f"{path_prefix}_material.csv"
        rows = _kv(mwinfo_result, [
            ("sheet_grade", "Sheet steel grade", None),
            ("stator_lam", "Stator lamination", None),
            ("rotor_lam", "Rotor lamination", None),
        ])
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(rows)
        written.append(path)

    return written


# --------------------------------------------------------------------- #
# PDF export
# --------------------------------------------------------------------- #

def _styled_table(data, col_widths=None):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2b3a55")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f0f2f6")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#c0c4cc")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


def write_pdf_report(header_info: dict | None = None,
                      tlast_result: dict | None = None,
                      synmom_result: dict | None = None,
                      kdruck_result: dict | None = None,
                      mwinfo_result: dict | None = None,
                      mkpool_rows=None, eanz: int = 0,
                      output_path: str = "asyn4_report.pdf"):
    """
    Builds a single formatted PDF report covering whichever sections are
    supplied (pass None to skip a section). header_info is a plain dict
    of machine identification (e.g. from lesen()'s returned 'name',
    'reke', 'datum', ...) - shown as a title-page key/value table.

    Returns output_path.
    """
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("ReportTitle", parent=styles["Title"],
                                  fontSize=20, spaceAfter=6)
    section_style = ParagraphStyle("Section", parent=styles["Heading2"],
                                    spaceBefore=14, spaceAfter=6,
                                    textColor=colors.HexColor("#2b3a55"))
    body_style = styles["Normal"]

    doc = SimpleDocTemplate(output_path, pagesize=A4,
                             leftMargin=18 * mm, rightMargin=18 * mm,
                             topMargin=16 * mm, bottomMargin=16 * mm)
    story = []

    story.append(Paragraph("ASYN4 Machine Performance Report", title_style))
    if header_info:
        rows = [["Field", "Value"]]
        for key in ("name", "reke", "datum", "vanr", "fbnr", "monr"):
            if key in header_info and str(header_info[key]).strip():
                rows.append([key.upper(), str(header_info[key]).strip()])
        if len(rows) > 1:
            story.append(_styled_table(rows, col_widths=[40 * mm, 120 * mm]))
    story.append(Spacer(1, 8))

    if tlast_result is not None:
        story.append(Paragraph("Partial-Load Performance", section_style))
        story.append(_styled_table(_tlast_points_table(tlast_result)))
        story.append(Spacer(1, 6))

        story.append(Paragraph("Breakdown / Locked-Rotor Summary", section_style))
        story.append(_styled_table(_tlast_summary_table(tlast_result),
                                    col_widths=[90 * mm, 60 * mm]))

    if mkpool_rows is not None and eanz > 0:
        story.append(PageBreak())
        story.append(Paragraph("Magnetization Characteristic", section_style))
        story.append(Paragraph(
            "Generated by the magnetic-circuit solver (KREIS) as air-gap "
            "induction is swept from below rated voltage up to 150% of "
            "rated voltage.", body_style))
        story.append(Spacer(1, 4))
        story.append(_styled_table(_mkpool_table(mkpool_rows, eanz)))

    if synmom_result is not None and synmom_result.get("breakdown"):
        story.append(PageBreak())
        story.append(Paragraph("Synchronous (Cusp) Torque Analysis", section_style))
        story.append(Paragraph(
            f"Maximum standstill synchronous torque: "
            f"{synmom_result['msymax']:.4f} x Mn. "
            f"Starting endangered: {synmom_result['starting_endangered']}.",
            body_style))
        story.append(Spacer(1, 4))
        story.append(_styled_table(_synmom_table(synmom_result),
                                    col_widths=[70 * mm, 60 * mm]))

    if kdruck_result is not None or mwinfo_result is not None:
        story.append(PageBreak())
        story.append(Paragraph("Cooling & Material Information", section_style))
        if kdruck_result is not None:
            rows = _kv(kdruck_result, [
                ("iklink", "Dominant cooling path (IKLINK)", None),
                ("ar1", "Stator yoke iron area AR1", lambda v: f"{v:.5f} m2"),
                ("az1", "Stator tooth iron area AZ1", lambda v: f"{v:.5f} m2"),
                ("ar2", "Rotor yoke iron area AR2", lambda v: f"{v:.5f} m2"),
                ("az2", "Rotor tooth iron area AZ2", lambda v: f"{v:.5f} m2"),
            ])
            story.append(_styled_table(rows, col_widths=[90 * mm, 60 * mm]))
            story.append(Spacer(1, 6))
        if mwinfo_result is not None:
            rows = _kv(mwinfo_result, [
                ("sheet_grade", "Sheet steel grade", None),
                ("stator_lam", "Stator lamination", None),
                ("rotor_lam", "Rotor lamination", None),
            ])
            story.append(_styled_table(rows, col_widths=[90 * mm, 60 * mm]))

    doc.build(story)
    return output_path
