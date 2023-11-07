# usage: python populate-AR.py <AR filename, like Rutz_ARCatalog_MERRA2_2000.nc>

import xarray, numpy, scipy, sys, datetime, copy, random, math
from pymongo import MongoClient
from collections import defaultdict
import gridtools
from functools import partial
from geopy import distance

numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(linewidth=200000)

client = MongoClient('mongodb://database/argo')
db = client.argo

xar = xarray.open_dataset(sys.argv[1])

# helper functions

basins = xarray.open_dataset('parameters/basinmask_01.nc')
def find_basin(basins, lon, lat):
    # for a given lon, lat,
    # identify the basin from the lookup table.
    # choose the nearest non-nan grid point.

    gridspacing = 0.5

    basin = basins['BASIN_TAG'].sel(LONGITUDE=lon, LATITUDE=lat, method="nearest").to_dict()['data']
    if math.isnan(basin):
        # nearest point was on land - find the nearest non nan instead.
        lonplus = math.ceil(lon / gridspacing)*gridspacing
        lonminus = math.floor(lon / gridspacing)*gridspacing
        latplus = math.ceil(lat / gridspacing)*gridspacing
        latminus = math.floor(lat / gridspacing)*gridspacing
        grids = [(basins['BASIN_TAG'].sel(LONGITUDE=lonminus, LATITUDE=latminus, method="nearest").to_dict()['data'], distance.distance((lat, lon), (latminus, lonminus)).miles),
                 (basins['BASIN_TAG'].sel(LONGITUDE=lonminus, LATITUDE=latplus, method="nearest").to_dict()['data'], distance.distance((lat, lon), (latplus, lonminus)).miles),
                 (basins['BASIN_TAG'].sel(LONGITUDE=lonplus, LATITUDE=latplus, method="nearest").to_dict()['data'], distance.distance((lat, lon), (latplus, lonplus)).miles),
                 (basins['BASIN_TAG'].sel(LONGITUDE=lonplus, LATITUDE=latminus, method="nearest").to_dict()['data'], distance.distance((lat, lon), (latminus, lonplus)).miles)]

        grids = [x for x in grids if not math.isnan(x[0])]
        if len(grids) == 0:
            # all points on land
            #print('warning: all surrounding basin grid points are NaN')
            basin = -1
        else:
            grids.sort(key=lambda tup: tup[1])
            basin = grids[0][0]
    return int(basin)


def index2coords(longitudes, latitudes, index):
    # index [lat_idx, lon_idx]; return [lon, lat]
    lon = longitudes[index[1]] - 5./16.
    if lon < -180:
        lon += 360.

    if index[0] == 0:
        lat = -90
    elif index[0] == 361:
        lat = 90
    else:
        lat = latitudes[index[0]] - 0.25

    return [lon, lat]

def convert_hour(time):

    hh = int(time)
    mm = (time*60) % 60
    ss = (time*3600) % 60

    return "%d:%02d.%02d" % (hh, mm, ss)

def add_noise(mp):
    # take a multipolygon and add a tiny amount of noise to each vertex
    # to try and dodge https://jira.mongodb.org/browse/SERVER-52928

    MP = copy.deepcopy(mp)
    noise = (random.random()-0.5)/1000000000

    for i, blob in enumerate(mp['coordinates']):
        for j, loop in enumerate(blob):
            for k, vertex in enumerate(loop):
                MP['coordinates'][i][j][k][0] = mp['coordinates'][i][j][k][0] + noise
                MP['coordinates'][i][j][k][1] = mp['coordinates'][i][j][k][1] + noise
            MP['coordinates'][i][j][0] = MP['coordinates'][i][j][-1] # close the loop

    return MP

# generate a metadata document
meta = {
    '_id': 'ar',
    'data_type': 'atmospheric_rivers',
    'data_info': [
        ['ivt'],
        ['units', 'long_name'],
        [['kg/m/s', 'integrated water vapor transport']]
    ],
    'date_updated_argovis': datetime.datetime.now(),
    'source':[
        {
            'doi': 'https://doi.org/10.1175/MWR-D-13-00168.1'
        }
    ]
}
db.extendedMeta.replace_one({'_id': 'ar'}, meta, upsert=True)

# unpack the netcdf
longitudes = xar['longitude'].to_dict()['data']
latitudes = xar['latitude'].to_dict()['data']
cal_mons = xar['cal_mon'].to_dict()['data']
cal_days = xar['cal_day'].to_dict()['data']
cal_years = xar['cal_year'].to_dict()['data']
cal_hours = xar['cal_hour'].to_dict()['data']
ars = xar['ARs']
ivts = xar['IVT']

for timestep in range(len(cal_years)):
#for timestep in [2197]:

    cal_mon = cal_mons[timestep]
    cal_day = cal_days[timestep]
    cal_year = cal_years[timestep]
    cal_hour = cal_hours[timestep]
    hhmmss = convert_hour(cal_hour)
    ar = ars[timestep].to_numpy()
    ivt = ivts[timestep].to_numpy()

    # label the ARs and make periodic on longitude boundary
    labeled_map = gridtools.label_features(ar)
    labels = numpy.unique(labeled_map)
    for label in labels:
    #for label in [1]:
        if label == 0:
            continue
        else:
            cells = numpy.where(labeled_map == label)
            lats = [latitudes[x] for x in cells[0]]
            lons = [longitudes[x] for x in cells[1]]
            basin_set = set()
            for i in range(len(lats)):
                basin_set.add(find_basin(basins, lons[i], lats[i]))
            vapors = [ [ivt[cells[0][i]][cells[1][i]]] for i in range(len(cells[0])) ]
            raster = list(zip(lons, lats, vapors))

            geo, flags = gridtools.generate_geojson(labeled_map, label, partial(index2coords, longitudes, latitudes), reverse_winding=True)
            flags = list(flags)
            flags = ['south_pole' if x=='first_pole' else x for x in flags]
            flags = ['north_pole' if x=='last_pole' else x for x in flags]
            AR = {      
                '_id': f'{cal_year}{cal_mon}{cal_day}{cal_hour}_{label}' , 
                'timestamp': datetime.datetime(year=int(cal_year), month=int(cal_mon), day=int(cal_day), hour=int(cal_hour), minute=int((cal_hour*60) % 60), second=int((cal_hour*3600) % 60) ), 
                'basins': list(basin_set),
                'raster': raster, 
                'flags': flags,
                'geolocation': geo,
                'metadata': ['ar']
            }

        try:
            db.ar.insert_one(AR)
        except Exception as e:
            retries = 0
            while retries < 10:
                try: 
                    # failing due to https://jira.mongodb.org/browse/SERVER-52928
                    print(f'trying to write AR {AR["_id"]} with noise')
                    AR_noise = copy.deepcopy(AR)
                    AR_noise['true_geolocation'] = AR['geolocation']
                    AR_noise['geolocation'] = add_noise(AR_noise['geolocation'])
                    if 'flags' in AR_noise:
                        AR_noise['flags'].append('noise_added')
                    else:
                        AR_noise['flags'] = ['noise_added']
                    db.ar.insert_one(AR_noise)
                    retries = 9999
                except:
                    retries += 1
            if retries < 9999:
                print(f'unable to write AR {AR["_id"]}')


  