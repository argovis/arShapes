FROM python:3.11

RUN apt-get update -y ; apt-get install -y nano netcdf-bin
RUN pip install xarray netcdf4 numpy scipy
