FROM python:3.9

RUN apt-get update -y ; apt-get install -y nano netcdf-bin
RUN pip install xarray netcdf4 numpy scipy pymongo shapely geopy argovisHelpers
WORKDIR /app
COPY populate-AR.py populate-AR.py
COPY loaddata.sh loaddata.sh
COPY compute-summaries.py compute-summaries.py
COPY parameters/basinmask_01.nc parameters/basinmask_01.nc
RUN chown -R 1000660000 /app

