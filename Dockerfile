FROM python:3.10.0-slim 

COPY . /app
WORKDIR /app

RUN python3 -m venv /opt/venv

RUN pip install pip --upgrade
RUN /opt/venv/bin/pip install -r requirements.txt

RUN chmod +x entrypoint.sh

CMD ["/app/entrypoint.sh"]