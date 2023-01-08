#!/bin/bash
APP_PORT=${PORT:-8080}
cd /app
/opt/venv/bin/gunicorn --worker-tmp-dir /dev/shm hic_server.wsgi:application --bind "0.0.0.0:${APP_PORT}"
