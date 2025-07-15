# given one of the upstream files, output a netcdf file for a single timestep with AR mask on [min_lat, max_lat, min_lon, max_lon].
import xarray, numpy

ds = xarray.open_dataset("data/Rutz_ARCatalog_MERRA2_2022.nc")
lat_min, lat_max = 20, 50
lon_min, lon_max = -180, -110

lat_vals = ds['latitude'].values
lon_vals = ds['longitude'].values
lat_inds = numpy.where((lat_vals >= lat_min) & (lat_vals <= lat_max))[0]
lon_inds = numpy.where((lon_vals >= lon_min) & (lon_vals <= lon_max))[0]

for i in range(ds.dims['ntim']):
    year = int(ds['cal_year'][i])
    month = int(ds['cal_mon'][i])
    day = int(ds['cal_day'][i])
    hour = int(ds['cal_hour'][i])
    timestamp_str = f"{year:04d}{month:02d}{day:02d}_{hour:02d}"

    timestep_slice = ds.isel(
        ntim=i,
        nlat=lat_inds,
        nlon=lon_inds
    )

    out_filename = f"AR_timestep_{timestamp_str}.nc"
    timestep_slice.to_netcdf(out_filename)
