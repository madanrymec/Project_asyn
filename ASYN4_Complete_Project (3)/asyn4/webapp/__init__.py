"""
Flask application factory for the ASYN4 web application.
"""
import os
from flask import Flask

from .config import get_app_config


def create_app():
    app = Flask(__name__)
    app_config = get_app_config()
    app.config.update(app_config)

    os.makedirs(app_config["OUTPUT_DIR"], exist_ok=True)

    from . import routes
    app.register_blueprint(routes.bp)

    # Initialize the MySQL schema on startup. If MySQL isn't reachable
    # yet (e.g. first-ever run before the DB exists), don't crash app
    # startup - the error will surface clearly the first time a route
    # actually needs the database, which is easier to diagnose than a
    # failed `flask run`.
    try:
        from .db import init_schema
        init_schema()
    except Exception as exc:
        app.logger.warning(
            "Could not initialize MySQL schema at startup (%s). "
            "Check ASYN4_DB_HOST/PORT/USER/PASSWORD/NAME and that the "
            "database exists - see WEBAPP_README.md.", exc)

    return app
