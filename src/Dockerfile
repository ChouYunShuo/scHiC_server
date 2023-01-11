FROM python:3.10.0-slim 

# copy your local files to your
# docker container
COPY . /app

# update your environment to work
# within the folder you copied your 
# files above into
WORKDIR /app

# /opt: reserved for the installation of add-on application software packages.
# We'll use this to create & store our virtual environment

# Create a virtual environment in /opt
RUN python3 -m venv /opt/venv


#install h5py wheel dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends --yes build-essential
RUN apt-get install -y pkg-config libhdf5-dev libcurl4-gnutls-dev default-libmysqlclient-dev

# install psycopg2 mysql dependencies
RUN apt-get install -y postgresql-server-dev-all gcc python3-dev musl-dev


# Install requirments to new virtual environment
# requirements.txt must have gunicorn & django
RUN /opt/venv/bin/pip install pip --upgrade && \
    /opt/venv/bin/pip install -r requirements.txt --use-pep517 && \
    chmod +x entrypoint.sh

# CMD ["/app/entrypoint.sh"]