"""
Ties together einasyn_builder -> asyn4.main.run_asyn4 -> report_export
-> plotting for one web request: takes the submitted form values,
produces every output artifact (EINASYN file, CSV report(s), PDF
report, torque-speed chart), and returns a plain-dict summary suitable
both for rendering the results page and for storing in MySQL.
"""
import os
import traceback
import uuid

from asyn4.main import run_asyn4
from asyn4 import report_export as rex
from asyn4.prdaten_generator import generate_prdaten_report
from .einasyn_builder import build_einasyn


def run_calculation(values: dict, output_dir: str):
    """
    Runs one full calculation from submitted form values.

    Returns a dict:
        status: "ok" | "validation_error" | "error"
        iabort: 0 or 1
        error_message: str or None
        pruef_messages: list of (is_error, text) or []
        summary: dict of headline numbers for the results page (only
                 present if status == "ok")
        csv_prefix, pdf_path, plot_path: file paths (only if generated)
        run_dir: the per-run output directory these files live in
    """
    run_id_slug = uuid.uuid4().hex[:12]
    run_dir = os.path.join(output_dir, run_id_slug)
    os.makedirs(run_dir, exist_ok=True)

    einasyn_path = os.path.join(run_dir, "EINASYN.txt")
    try:
        build_einasyn(values, einasyn_path)
    except Exception as exc:
        return {
            "status": "error", "iabort": 1,
            "error_message": f"Could not build input file: {exc}",
            "pruef_messages": [], "run_dir": run_dir,
        }

    try:
        want_plot = str(values.get("rkz_slipcurve", "")).strip() != ""
        plot_path = os.path.join(run_dir, "torque_speed_chart.png") if want_plot else None
        result = run_asyn4(einasyn_path, irenr=1, verbose=False,
                            run_slip_curve=True, plot_path=plot_path)
    except Exception as exc:
        return {
            "status": "error", "iabort": 1,
            "error_message": f"Calculation failed: {exc}\n\n{traceback.format_exc()}",
            "pruef_messages": [], "run_dir": run_dir,
        }

    if result.get("iabort") == 1:
        return {
            "status": "validation_error", "iabort": 1,
            "error_message": "Input data failed validation (see messages below).",
            "pruef_messages": result.get("pruef_messages", []),
            "run_dir": run_dir,
        }

    tlast_out = result.get("tlast_result")
    synmom_out = result.get("synmom_result")
    kdruck_out = result.get("kdruck_result")
    pdruck_out = result.get("pdruck_result")

    from asyn4.state import MachineState  # noqa: F401 (documents dependency, no instantiation here)

    csv_prefix = os.path.join(run_dir, "report")
    csv_paths = []
    pdf_path = None
    try:
        csv_paths = rex.write_csv_report(
            tlast_result=tlast_out, synmom_result=synmom_out,
            kdruck_result=kdruck_out, path_prefix=csv_prefix)
        pdf_path = os.path.join(run_dir, "report.pdf")
        rex.write_pdf_report(
            header_info={"name": values.get("name", ""),
                         "reke": values.get("reke", ""),
                         "datum": values.get("datum", "")},
            tlast_result=tlast_out, synmom_result=synmom_out,
            kdruck_result=kdruck_out, output_path=pdf_path)
    except Exception:
        # Reports are a nice-to-have on top of a successful calculation;
        # don't fail the whole request if report generation stumbles -
        # surface it, but keep the numeric results available.
        pdf_path = None

    summary = _build_summary(tlast_out, kdruck_out, synmom_out, result, values)

    if tlast_out and pdruck_out:
        prdaten_text = generate_prdaten_report(result, values)
        prdaten_path = os.path.join(run_dir, "PRDATEN_EN.TXT")
        with open(prdaten_path, "w", encoding="utf-8") as f:
            f.write(prdaten_text)
        summary["prdaten_text"] = prdaten_text
        summary["prdaten_path"] = prdaten_path

    return {
        "status": "ok", "iabort": 0, "error_message": None,
        "pruef_messages": [], "run_dir": run_dir,
        "summary": summary,
        "csv_paths": csv_paths,
        "csv_prefix": csv_prefix if csv_paths else None,
        "pdf_path": pdf_path,
        "plot_path": result.get("plot_path"),
    }


def _build_summary(tlast_out, kdruck_out, synmom_out, result, values=None):
    """Extracts the headline numbers shown on the results page."""
    import math
    summary = {}
    if tlast_out:
        rated = tlast_out["points"][0] if tlast_out["points"] else None

        # Scale phase current to line current for display (Delta = 1 or 3)
        ishalt = int(values.get("ishalt", 1)) if values else 1
        i_scale = 1.0 if ishalt in (0, 2) else math.sqrt(3.0)

        mn_rated = rated.get("m") if rated and rated.get("m") else 1.0
        i_rated = (rated.get("i1") * i_scale) if rated and rated.get("i1") else 1.0

        summary["rated_point"] = (
            {
                "converged": rated.get("converged", False),
                "p_mech_kw": rated.get("pmech", 0.0) / 1000.0 if rated.get("converged") else None,
                "current_a": i_rated if rated.get("converged") else None,
                "cos_phi": rated.get("cosph") if rated.get("converged") else None,
                "efficiency": rated.get("eta") if rated.get("converged") else None,
                "torque_nm": mn_rated if rated.get("converged") else None,
            } if rated else None
        )
        if values is not None:
            summary["rated_voltage_v"] = values.get("un")

        # Breakdown Torque
        mk_m = tlast_out.get("kipp", {}).get("m", 0.0)
        summary["breakdown_torque_nm"] = mk_m
        summary["breakdown_torque_pu"] = mk_m / mn_rated if mn_rated else 0.0
        summary["breakdown_slip_pct"] = (tlast_out.get("kipp", {}).get("skipp", 0) or 0) * 100.0

        # Locked-Rotor Torque and Current (MUST USE 'corrected', NOT 'theoretical')
        lr_m = tlast_out.get("corrected", {}).get("m", 0.0)
        lr_i1 = tlast_out.get("corrected", {}).get("i1", 0.0)

        summary["locked_rotor_torque_nm"] = lr_m
        summary["locked_rotor_torque_pu"] = lr_m / mn_rated if mn_rated else 0.0
        summary["locked_rotor_current_a"] = (lr_i1 * i_scale) if lr_i1 else None
        summary["locked_rotor_current_pu"] = (lr_i1 * i_scale) / i_rated if i_rated else 0.0

        summary["thermal_time_constant_s"] = tlast_out.get("thaupt")

    if kdruck_out:
        summary["dominant_cooling_path"] = kdruck_out.get("iklink")
    if synmom_out:
        summary["max_synchronous_torque_pu"] = synmom_out.get("msymax")
        summary["starting_endangered"] = synmom_out.get("starting_endangered")
    if result.get("slast_result"):
        summary["saddle_torque_nm"] = result["slast_result"].get("msat")
    summary["active_mass_kg"] = result.get("gaktiv")
    summary["rotor_inertia_kg_m2"] = result.get("jaktiv")
    return summary
