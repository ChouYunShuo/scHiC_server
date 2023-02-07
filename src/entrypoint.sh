#!/bin/bash
cd /app
/opt/venv/bin/gunicorn --worker-tmp-dir /dev/shm hic_server.wsgi:application
