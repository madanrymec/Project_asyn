"""
MySQL persistence layer for the ASYN4 web application, using pymysql
directly (no ORM, per the project's requirements).

Two tables:
  runs        - one row per calculation: every input field (as
                submitted), plus top-level status/summary output
                fields and file paths for the generated CSV/PDF/PNG.
  run_inputs  - a normalized key/value table mirroring every field in
                field_schema.FIELDS, so individual input values can be
                queried/filtered without parsing the `runs.input_json`
                blob. (Kept alongside the JSON blob rather than
                instead of it: the blob is what repopulates the form
                exactly as submitted; the normalized table is what
                makes "find all runs with PN between X and Y" a plain
                SQL query instead of a JSON scan.)

Connections are opened per-request and closed immediately (Flask's
request lifecycle already bounds this correctly) rather than pooled -
pymysql doesn't include a pool itself, and for a form-driven app like
this one, connection-per-request is simple and correct. If this needs
to scale beyond that, swap in DBUtils.PooledDB or PyMySQL's own
`Connection` reuse pattern without changing any of the call sites
below, since they all go through `get_connection()`.
"""
import json
import pymysql
import pymysql.cursors

from .config import get_db_config


def _json_safe(obj):
    """
    Fallback encoder for json.dumps(default=...): converts numpy scalar
    types (bool_, integer, floating - and anything else numpy-ish with
    an .item() method) to native Python types.

    Added after a real bug: numpy.bool_ does NOT subclass Python's
    bool (unlike numpy.float64, which does subclass float and so
    serializes fine on its own) - a numpy.bool_ produced by a numpy-
    array comparison (synmom's starting_endangered field) reached
    json.dumps() here unconverted and raised "Object of type bool is
    not JSON serializable". That specific case is now fixed at the
    source (reports.py casts it to a native bool before returning),
    but this fallback exists so any other stray numpy scalar - now or
    added later - degrades to a clear, correct JSON value instead of
    crashing the whole save.
    """
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _dumps(obj):
    return json.dumps(obj, default=_json_safe)


def get_connection():
    """Opens a new pymysql connection using the app's configured
    credentials. Caller is responsible for closing it (use as a
    context manager: `with get_connection() as conn:`)."""
    cfg = get_db_config()
    return pymysql.connect(
        host=cfg["host"],
        port=cfg["port"],
        user=cfg["user"],
        password=cfg["password"],
        database=cfg["database"],
        charset="utf8mb4",
        cursorclass=pymysql.cursors.DictCursor,
        autocommit=False,
    )


SCHEMA_STATEMENTS = [
    """
    CREATE TABLE IF NOT EXISTS runs (
        id INT AUTO_INCREMENT PRIMARY KEY,
        created_at DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
        machine_name VARCHAR(64),
        reke VARCHAR(64),
        status VARCHAR(16) NOT NULL,
        iabort TINYINT NOT NULL DEFAULT 0,
        error_message TEXT,
        input_json LONGTEXT NOT NULL,
        result_summary_json LONGTEXT,
        csv_prefix VARCHAR(255),
        pdf_path VARCHAR(255),
        plot_path VARCHAR(255),
        INDEX idx_created_at (created_at),
        INDEX idx_machine_name (machine_name)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """,
    """
    CREATE TABLE IF NOT EXISTS run_inputs (
        id INT AUTO_INCREMENT PRIMARY KEY,
        run_id INT NOT NULL,
        field_name VARCHAR(64) NOT NULL,
        field_value VARCHAR(255),
        FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE,
        INDEX idx_run_id (run_id),
        INDEX idx_field_name (field_name)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """,
]


def init_schema():
    """Creates the runs/run_inputs tables if they don't already exist.
    Safe to call on every app startup."""
    with get_connection() as conn:
        with conn.cursor() as cur:
            for stmt in SCHEMA_STATEMENTS:
                cur.execute(stmt)
        conn.commit()


def save_run(values: dict, status: str, iabort: int = 0,
             error_message: str = None, result_summary: dict = None,
             csv_prefix: str = None, pdf_path: str = None,
             plot_path: str = None) -> int:
    """
    Inserts one run: the full submitted input (as JSON, for exact form
    repopulation) plus a normalized copy in run_inputs (for querying),
    plus the result summary and generated-file paths. Returns the new
    run's id.
    """
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO runs
                    (machine_name, reke, status, iabort, error_message,
                     input_json, result_summary_json, csv_prefix,
                     pdf_path, plot_path)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                """,
                (
                    values.get("name", ""), values.get("reke", ""),
                    status, iabort, error_message,
                    _dumps(values),
                    _dumps(result_summary) if result_summary else None,
                    csv_prefix, pdf_path, plot_path,
                ),
            )
            run_id = cur.lastrowid

            rows = [(run_id, k, "" if v is None else str(v))
                    for k, v in values.items()]
            if rows:
                cur.executemany(
                    "INSERT INTO run_inputs (run_id, field_name, field_value) "
                    "VALUES (%s, %s, %s)",
                    rows,
                )
        conn.commit()
    return run_id


def get_run(run_id: int):
    """Returns the run row (dict) for run_id, with input_json and
    result_summary_json already parsed back into dicts, or None if not
    found."""
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT * FROM runs WHERE id = %s", (run_id,))
            row = cur.fetchone()
    if row is None:
        return None
    row["input_values"] = json.loads(row["input_json"]) if row["input_json"] else {}
    row["result_summary"] = (json.loads(row["result_summary_json"])
                              if row["result_summary_json"] else None)
    return row


def list_runs(limit: int = 50, offset: int = 0):
    """Returns the most recent runs (without the full input/result JSON
    blobs, to keep the history listing light), most recent first."""
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT id, created_at, machine_name, reke, status, iabort,
                       csv_prefix, pdf_path, plot_path
                FROM runs
                ORDER BY created_at DESC, id DESC
                LIMIT %s OFFSET %s
                """,
                (limit, offset),
            )
            return cur.fetchall()


def delete_run(run_id: int):
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("DELETE FROM runs WHERE id = %s", (run_id,))
        conn.commit()
