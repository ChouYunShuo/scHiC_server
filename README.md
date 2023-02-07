# A Django server to host single cell HiC data.

##Features

- Upload data as single cell hdf5 format
- Currently support hdf5 transformation from .hic datasets
- View data as interactive plots and matrices
- Search for specific genomic regions
- Containerized using docker compose for fast deployment

##Prerequisites

- Python 3.10.0
- Django 3.2
- h5py 3.6.0
- Numpy 1.23.4

##Installation

#### 1. Clone the repository

```
git clone https://github.com/ChouYunShuo/scHiC_server.git
```

Newest development is on the _Development_ branch.

#### 2. Change into the project directory

```
cd scHiC_server
```

#### 3. Create a virtual environment and activate it or use conda

```
python3 -m venv env
source env/bin/activate
```

####4. Install the required packages

```
pip install -r requirements.txt
```

##Usage
####1.Place your h5, .hic files under the hic_data respository

```
mkdir hic_data
```

#### 2.Create a .env file to store environment variables

```
cd src
vim .env
```

Note that PostgreSQL needs to be installed and configured and set the database related environment variables in .env

_Or you can use the docker-compose to start a PostgreSQL instance with docker_

#### 3. Start the Django development server

```
python manage.py runserver
```

The default port is 8000

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](https://github.com/ChouYunShuo/scHiC_server/LICENSE) file for details.
