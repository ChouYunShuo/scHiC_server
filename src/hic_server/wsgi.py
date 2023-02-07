"""
WSGI config for hic_server project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.1/howto/deployment/wsgi/
"""

import os
import pathlib
from django.core.wsgi import get_wsgi_application

CURRENT_DIR = pathlib.Path(__file__).resolve().parent


os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'hic_server.settings')

application = get_wsgi_application()
