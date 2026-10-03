"""
Entry point for running the ASYN4 web application's development server.

Usage:
    python run_webapp.py

Configure the MySQL connection and app secret via environment
variables first - see WEBAPP_README.md for the full list
(ASYN4_DB_HOST, ASYN4_DB_PORT, ASYN4_DB_USER, ASYN4_DB_PASSWORD,
ASYN4_DB_NAME, ASYN4_SECRET_KEY, ASYN4_OUTPUT_DIR).

For a production deployment, use a real WSGI server (gunicorn,
waitress, etc.) pointed at `webapp:create_app()` instead of this
script's built-in Flask dev server.
"""
from webapp import create_app

app = create_app()

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
