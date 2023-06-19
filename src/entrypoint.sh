#!/bin/bash
PORT=${PORT:-8020}
cd /app
/opt/venv/bin/gunicorn -w 8 --worker-tmp-dir /dev/shm hic_server.wsgi:application --bind "0.0.0.0:${PORT}" --worker-class gevent
