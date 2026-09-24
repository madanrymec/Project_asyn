# ASYN4 Web Application

A full-stack Flask web app wrapping the ASYN4 induction-machine
calculation engine: fill in a machine-data form in the browser, run
the full calculation, and get back a results page with headline
numbers, a torque-speed chart, and downloadable PDF/CSV reports.
Every submission is stored in MySQL (via `pymysql`, no ORM) so you can
revisit any past calculation from the History page.

## Architecture

```
asyn4/              the calculation engine (already built, unchanged)
webapp/
  __init__.py        Flask app factory
  config.py           settings, all overridable via environment variables
  db.py                pymysql schema + CRUD (runs, run_inputs tables)
  field_schema.py     the ~140 EINASYN input fields: name/label/unit/default
  einasyn_builder.py  form values -> EINASYN file (via asyn4.lesen's own
                       _CardWriter, so it can't drift out of sync with the reader)
  calc_runner.py       orchestrates: build EINASYN -> run_asyn4() ->
                        CSV/PDF/chart generation -> summary dict
  routes.py            HTTP routes
  templates/           Jinja2 templates
  static/style.css     styling
run_webapp.py          dev-server entry point
schema.sql              MySQL DDL (for manual setup; the app also
                         creates these tables automatically on startup)
requirements_webapp.txt Python dependencies
```

## 1. Install dependencies

```bash
pip install -r requirements_webapp.txt --break-system-packages
```
(drop `--break-system-packages` if you're using a virtualenv, which is recommended)

## 2. Set up MySQL

```bash
mysql -u root -p -e "CREATE DATABASE asyn4 CHARACTER SET utf8mb4;"
```

The application creates its two tables (`runs`, `run_inputs`)
automatically on first startup - you don't need to run `schema.sql` by
hand unless you want to review or customize the schema first.

## 3. Configure

Set these environment variables (defaults shown are for a typical
local MySQL install; change as needed):

```bash
export ASYN4_DB_HOST=localhost
export ASYN4_DB_PORT=3306
export ASYN4_DB_USER=root
export ASYN4_DB_PASSWORD=yourpassword
export ASYN4_DB_NAME=asyn4
export ASYN4_SECRET_KEY=change-this-to-something-random
```

(`ASYN4_OUTPUT_DIR` is optional - defaults to `webapp/run_outputs/`,
where generated EINASYN files, CSVs, PDFs, and charts are written, one
subfolder per run.)

## 4. Run it

```bash
python run_webapp.py
```

Then open **http://localhost:5000** in a browser.

For a production deployment, don't use the Flask dev server - point a
real WSGI server at the app factory instead, e.g.:

```bash
pip install gunicorn --break-system-packages
gunicorn -w 4 -b 0.0.0.0:8000 "webapp:create_app()"
```

## How it works

1. **Form** (`/`) — organized into sections matching the original
   EINASYN input file's own card layout (Z1 identification, Z4 rating,
   Z8 stator slot, Z12 rotor bar, etc.) — each section header shows its
   card tag. Leave a field at 0/blank to use the same automatic default
   the original design program would compute.
2. **Submit** (`POST /calculate`) — the form values are written into a
   real EINASYN-format file (using the exact same `_CardWriter` utility
   the calculation engine's own test suite uses), then run through the
   full pipeline (`asyn4.main.run_asyn4`): geometry, magnetic circuit,
   windings, performance, load curve, and (if enabled) the torque-speed
   chart.
3. **Results** — headline numbers (rated output, current, efficiency,
   breakdown torque, locked-rotor values), the chart inline, and
   download buttons for the PDF report and CSV exports. Everything is
   saved to MySQL — the exact input you submitted (for exact
   re-editing/audit later) and a summary of the results.
4. **History** (`/history`) — every past run, with status (success /
   validation error / calculation error), and a link back to its full
   results page.

## Data validation

Invalid input (e.g. rotor bore larger than stator bore) is caught by
the calculation engine's own `PRUEF` data-validation step before any
calculation runs — you'll see the specific error message(s) on the
results page instead of a wrong or crashed calculation.

## Notes on the calculation engine itself

This web app is a thin, additive layer — it does not modify
`asyn4/` in any way. See `asyn4/CONVERSION_STATUS.md` for the full
history of the underlying ASYN4 FORTRAN-to-Python port, including one
important open item: the breakdown-torque search (`SK`) has shown an
anomalous result on synthetic test data and would benefit from
validation against a real motor's known-correct data. This doesn't
affect the web app's functionality, but is worth knowing if a
breakdown-torque figure looks implausible for a particular design.
