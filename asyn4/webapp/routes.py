"""
Routes for the ASYN4 web application.
"""
import os
import json
from flask import (Blueprint, render_template, request, redirect,
                    url_for, current_app, send_file, flash, abort)

from .field_schema import fields_by_section, default_values, FIELDS, visibility_rules
from .calc_runner import run_calculation
from . import db as dbmod

bp = Blueprint("asyn4", __name__)

_FIELD_BY_NAME = {f["name"]: f for f in FIELDS}


def _parse_form_values(form) -> dict:
    """Reads every field_schema field out of the submitted form,
    coercing to the declared type. Blank/invalid numeric fields fall
    back to 0 (which then triggers LESEN's own default-filling logic -
    see field_schema.py's module docstring), never raise."""
    values = {}
    for f in FIELDS:
        raw = form.get(f["name"], "")
        if f["type"] == "str":
            values[f["name"]] = raw.strip()[:f.get("maxlen", 255)]
        elif f["type"] == "int":
            try:
                values[f["name"]] = int(float(raw)) if raw.strip() != "" else f["default"]
            except ValueError:
                values[f["name"]] = f["default"]
        else:  # float
            try:
                values[f["name"]] = float(raw) if raw.strip() != "" else f["default"]
            except ValueError:
                values[f["name"]] = f["default"]
    return values


def _render_form(values: dict, edit_of_run_id=None):
    """Shared form renderer for both a fresh (blank/default) calculation
    and editing a past run's inputs. edit_of_run_id, when set, is shown
    in the form so the person knows they're starting from a previous
    run's values rather than a blank sheet - it does NOT overwrite that
    run when submitted; submitting always creates a new run, keeping
    the original in history for audit/comparison."""
    return render_template("form.html", sections=fields_by_section(),
                            values=values, field_by_name=_FIELD_BY_NAME,
                            visibility_rules_json=json.dumps(visibility_rules()),
                            edit_of_run_id=edit_of_run_id)


@bp.route("/", methods=["GET"])
def index():
    return _render_form(default_values())


@bp.route("/runs/<int:run_id>/edit", methods=["GET"])
def edit_run(run_id):
    """Reopens the input form pre-filled with a past run's exact
    submitted values, so the person can tweak a design and recalculate
    instead of re-entering ~150 fields from scratch. Submitting the
    form from here creates a NEW run (via the normal /calculate route)
    - the original run this was opened from is left untouched in
    history, so both the "before" and "after" calculations stay
    available to compare."""
    row = dbmod.get_run(run_id)
    if row is None:
        abort(404)
    # Merge onto the schema's defaults first, so any field added to
    # field_schema.py *after* this run was saved still gets a sensible
    # default instead of appearing blank/missing in the reopened form.
    values = default_values()
    values.update(row["input_values"])
    return _render_form(values, edit_of_run_id=run_id)


@bp.route("/calculate", methods=["POST"])
def calculate():
    values = _parse_form_values(request.form)
    output_dir = current_app.config["OUTPUT_DIR"]

    result = run_calculation(values, output_dir)

    csv_prefix = result.get("csv_prefix")
    pdf_path = result.get("pdf_path")
    plot_path = result.get("plot_path")

    # PRUEF validation messages have no dedicated DB column - fold them
    # into error_message as formatted text so view_run() can still show
    # them after a redirect/reload, instead of only being available on
    # this first response.
    error_message = result.get("error_message")
    if result.get("pruef_messages"):
        lines = [("ERROR: " if is_err else "WARN: ") + text
                 for is_err, text in result["pruef_messages"]]
        error_message = (error_message or "") + "\n" + "\n".join(lines)

    try:
        run_id = dbmod.save_run(
            values, status=result["status"], iabort=result["iabort"],
            error_message=error_message,
            result_summary=result.get("summary"),
            csv_prefix=csv_prefix, pdf_path=pdf_path, plot_path=plot_path,
        )
    except Exception as exc:
        # The calculation itself may have succeeded even if the DB
        # write failed (e.g. MySQL not reachable) - show the person
        # their results anyway, with a clear warning, rather than
        # discarding a successful calculation because of a storage
        # problem.
        flash(f"Warning: results were computed but could not be saved "
              f"to the database ({exc}). They are shown below but will "
              f"not appear in history.", "warning")
        return render_template("result.html", result=result, values=values,
                                run_id=None)

    return redirect(url_for("asyn4.view_run", run_id=run_id))


@bp.route("/runs/<int:run_id>", methods=["GET"])
def view_run(run_id):
    row = dbmod.get_run(run_id)
    if row is None:
        abort(404)

    # Recover the individual (is_error, text) messages that were folded
    # into error_message at save time (see calculate() above) so the
    # results page can render them the same way whether viewed right
    # after submission or reloaded later from history.
    pruef_messages = []
    error_message = row["error_message"]
    if row["status"] == "validation_error" and error_message:
        lines = error_message.split("\n")[1:]  # skip the summary line
        for line in lines:
            if line.startswith("ERROR: "):
                pruef_messages.append((True, line[len("ERROR: "):]))
            elif line.startswith("WARN: "):
                pruef_messages.append((False, line[len("WARN: "):]))

    result = {
        "status": row["status"], "iabort": row["iabort"],
        "error_message": error_message,
        "pruef_messages": pruef_messages,
        "summary": row["result_summary"],
        "csv_prefix": row["csv_prefix"],
        "pdf_path": row["pdf_path"],
        "plot_path": row["plot_path"],
    }
    return render_template("result.html", result=result,
                            values=row["input_values"], run_id=run_id,
                            created_at=row["created_at"])


@bp.route("/history", methods=["GET"])
def history():
    page = max(int(request.args.get("page", 1)), 1)
    per_page = 25
    runs = dbmod.list_runs(limit=per_page, offset=(page - 1) * per_page)
    return render_template("history.html", runs=runs, page=page)


@bp.route("/runs/<int:run_id>/delete", methods=["POST"])
def delete_run(run_id):
    dbmod.delete_run(run_id)
    flash("Run deleted.", "info")
    return redirect(url_for("asyn4.history"))


@bp.route("/runs/<int:run_id>/download/<kind>")
def download(run_id, kind):
    row = dbmod.get_run(run_id)
    if row is None:
        abort(404)

    if kind == "pdf" and row["pdf_path"]:
        return send_file(row["pdf_path"], as_attachment=True)
    if kind == "plot" and row["plot_path"]:
        return send_file(row["plot_path"], as_attachment=True)
    if kind == "prdaten" and row.get("result_summary", {}).get("prdaten_path"):
        path = row["result_summary"]["prdaten_path"]
        if os.path.exists(path):
            return send_file(path, as_attachment=True, download_name="PRDATEN_EN.TXT")
    if kind.startswith("csv:") and row["csv_prefix"]:
        section = kind.split(":", 1)[1]
        path = f"{row['csv_prefix']}_{section}.csv"
        if os.path.exists(path):
            return send_file(path, as_attachment=True)
    abort(404)


@bp.route("/runs/<int:run_id>/plot")
def view_plot(run_id):
    """Serves the chart inline (not as a download) for embedding in
    the results page with an <img> tag."""
    row = dbmod.get_run(run_id)
    if row is None or not row["plot_path"] or not os.path.exists(row["plot_path"]):
        abort(404)
    return send_file(row["plot_path"], mimetype="image/png")
