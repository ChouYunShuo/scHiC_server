# Django Server for Hosting Single-Cell HiC Data

## Features

- Upload data in single-cell HDF5 format
- Supports HDF5 transformation from .hic datasets
- View data as interactive plots and matrices
- Search for specific genomic regions
- Containerized with Docker Compose for easy deployment

## Prerequisites
- Docker `>=19.03.13` (tested on 19.03.13, 26.1.0)
- Docker-compose `>=v2.29.2` (tested on v2.29.2, v1.26.0)

## Running the Data Server Locally

## Installation

#### 1. Clone the repository

```
git clone https://github.com/ChouYunShuo/scHiC_server
```


#### 2. Change into the project directory

```
cd scHiC_server
```

#### 3. Create a Docker Network

```bash
docker network create cellscope_network

```

### 4. Run the Data server Using Docker Compose
```
docker-compose up --build 
```

## Usage

#### 1. Place Your H5 Files in the hic_data Folder

```bash
cp <your-data-path> ./hic_data
```
 
#### 2. Create the Corresponding Dataset and Session in Django Admin
1. Open your web browser and navigate to the Django Admin interface:
```bash
http://127.0.0.1:8020/admin/
```
(Replace the URL and port with your data server URL and port if they are different.)

2. **Add a Session**:
- In the **Sessions** section, add a new session.
- For the `filepath`, enter the session filename (e.g., `Ramani_et_al.json`).

3. **Add a Dataset**:
- In the **Datasets** section, add a new dataset.
- For the `filepath`, enter the dataset filename (e.g., `Ramani_et_al.h5`).
- For `resolutions`, enter a single resolution (e.g., `500000`).
- Copy the corresponding `session_uuid` and paste it into the `Session_uuid` field.

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](/LICENSE) file for details.
