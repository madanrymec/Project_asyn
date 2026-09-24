-- ASYN4 web application - MySQL schema.
--
-- This file is provided for manual setup / review; the application
-- also creates these tables automatically on startup (see
-- webapp/db.py's init_schema()), so running this by hand is optional.
--
-- Setup:
--   mysql -u root -p -e "CREATE DATABASE asyn4 CHARACTER SET utf8mb4;"
--   mysql -u root -p asyn4 < schema.sql

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

CREATE TABLE IF NOT EXISTS run_inputs (
    id INT AUTO_INCREMENT PRIMARY KEY,
    run_id INT NOT NULL,
    field_name VARCHAR(64) NOT NULL,
    field_value VARCHAR(255),
    FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE,
    INDEX idx_run_id (run_id),
    INDEX idx_field_name (field_name)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
