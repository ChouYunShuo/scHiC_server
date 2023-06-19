#!/bin/bash

SUPERUSER_EMAIL=${DJANGO_SUPERUSER_EMAIL:-"yunshuoc@andrew.cmu.edu"}
cd /app/

# Collect static files
/opt/venv/bin/python manage.py collectstatic --noinput

/opt/venv/bin/python manage.py migrate --noinput 
/opt/venv/bin/python manage.py createsuperuser --email $SUPERUSER_EMAIL --noinput || true