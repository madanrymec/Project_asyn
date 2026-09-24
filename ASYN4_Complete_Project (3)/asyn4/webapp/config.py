"""
Configuration for the ASYN4 web application. All settings are
overridable via environment variables so the same code runs in dev
and in a real deployment without editing source files - set these in
your shell, a `.env` file (if you use python-dotenv), or your hosting
platform's environment settings before starting the app.
"""
import os


def get_db_config():
    return {
        "host": os.environ.get("ASYN4_DB_HOST", "localhost"),
        "port": int(os.environ.get("ASYN4_DB_PORT", "3306")),
        "user": os.environ.get("ASYN4_DB_USER", "root"),
        "password": os.environ.get("ASYN4_DB_PASSWORD", ""),
        "database": os.environ.get("ASYN4_DB_NAME", "asyn4"),
    }


def get_app_config():
    return {
        "SECRET_KEY": os.environ.get("ASYN4_SECRET_KEY", "dev-secret-change-me"),
        "OUTPUT_DIR": os.environ.get("ASYN4_OUTPUT_DIR",
                                      os.path.join(os.path.dirname(__file__), "run_outputs")),
        "MAX_CONTENT_LENGTH": 2 * 1024 * 1024,  # 2 MB form-post safety cap
    }
