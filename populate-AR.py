# usage: python populate-AR.py <AR filename, like Rutz_ARCatalog_MERRA2_2000.nc>

import xarray, numpy, scipy, sys, datetime, copy, random, math
from pymongo import MongoClient
from collections import defaultdict
import gridtools

numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(linewidth=200000)

client = MongoClient('mongodb://database/argo')
db = client.argo

xar = xarray.open_dataset(sys.argv[1])

# helper functions

def index2coords(index, longitudes, latitudes):
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

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[str(item)].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)

def loopsort(elt):
    return elt[1][1] - elt[1][0]

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

# # hack
# def add_noise(mp):
#     # take a multipolygon and add a tiny amount of noise to each vertex
#     # to try and dodge https://jira.mongodb.org/browse/SERVER-52928

#     MP = copy.deepcopy(mp)
#     noise = (random.random()-0.5)/1000000000

#     for i, blob in enumerate(mp['coordinates']):
#         for j, loop in enumerate(blob):
#             for k, vertex in enumerate(loop):
#             #     if MP['coordinates'][i][j][k][0] < 0:
#             #         MP['coordinates'][i][j][k][0] += 360
#             #     MP['coordinates'][i][j][k][0] -= 100
#                 MP['coordinates'][i][j][k][1] += 10
#             MP['coordinates'][i][j][0] = MP['coordinates'][i][j][-1] # close the loop

#     return MP
# # end hack

# # hack
# def add_noise(mp):
#     # take a multipolygon and add a tiny amount of noise to each vertex
#     # to try and dodge https://jira.mongodb.org/browse/SERVER-52928

#     MP = copy.deepcopy(mp)
#     noise = (random.random()-0.5)/1000000000

#     for i, blob in enumerate(mp['coordinates']):
#         for j, loop in enumerate(blob):
#             pass
#             MP['coordinates'][i][j].reverse()
#     return MP
# # end hack


# unpack the netcdf
longitudes = xar['longitude'].to_dict()['data']
latitudes = xar['latitude'].to_dict()['data']
cal_mons = xar['cal_mon'].to_dict()['data']
cal_days = xar['cal_day'].to_dict()['data']
cal_years = xar['cal_year'].to_dict()['data']
cal_hours = xar['cal_hour'].to_dict()['data']
ars = xar['ARs']
ivts = xar['IVT']

#for timestep in range(len(cal_years)):
for timestep in [2197]:

    cal_mon = cal_mons[timestep]
    cal_day = cal_days[timestep]
    cal_year = cal_years[timestep]
    cal_hour = cal_hours[timestep]
    hhmmss = convert_hour(cal_hour)
    ar = ars[timestep].to_numpy()
    ivt = ivts[timestep].to_numpy()
    # if numpy.any(ar[0]) or numpy.any(ar[-1]):
    #     print(timestep, '-----------')
    #     print(ar[0])
    #     print(ar[-1])
    # continue
    #ar[-1] = numpy.zeros(len(ar[-1])) # last latitude bin is strictly the north pole, doesn't make sense to bin in longitude; suppress
    #ar[-10] = numpy.ones(len(ar[-10])) # hack to check for weird behavior
    #ar[-11] = numpy.ones(len(ar[-11]))
    # ar[-2] = numpy.ones(len(ar[-2]))
    # print(timestep, cal_year, cal_mon, cal_day, cal_hour)
    # continue

    # label the ARs and make periodic on longitude boundary
    labeled_map = gridtools.label_features(ar)
    # create non-diagonally connected sublabels for each AR that are globally unique
    sublabel_map = gridtools.label_features(ar, structure=[[0,1,0],[1,1,1],[0,1,0]])
    #print(labeled_map)

    # get distinct AR labels
    labels = numpy.unique(labeled_map)
    sublabels = numpy.unique(sublabel_map)

    # identify blobs
    ARs = []
    #for label in labels:
    for label in [1]:
        flags = set(())
        if label == 0:
            continue
        else:
            cells = numpy.where(labeled_map == label)
            lats = [latitudes[x] for x in cells[0]]
            lons = [longitudes[x] for x in cells[1]]
            vapors = [ [ivt[cells[0][i]][cells[1][i]]] for i in range(len(cells[0])) ]
            raster = list(zip(lons, lats, vapors))
            
            # generate a loop for every sublabel in this AR
            local_map = ars[timestep].to_numpy()
            #local_map[-10] = numpy.ones(len(local_map[-10])) # hack to check for weird behavior
            #local_map[-11] = numpy.ones(len(local_map[-11]))
            # local_map[-2] = numpy.ones(len(local_map[-2]))
            local_map[labeled_map != label] = 0 # gets rid of other ARs in a the AR binary flag map
            local_sublabels_map = numpy.copy(sublabel_map)
            local_sublabels_map[local_map == 0] = 0 # sublabel map only including this AR
            local_sublabels = [x for x in numpy.unique(local_sublabels_map) if x != 0] # these are the rings that belong as top level objects in this AR
            loops = [gridtools.trace_shape(local_sublabels_map, sublabel) for sublabel in local_sublabels]
            for geo in loops:
                if len(geo) == 2:
                    flags.add('annulus')

            #print(loops)
            #print([index2coords(vertex, longitudes, latitudes) for vertex in loops[0][1] ])

            # print(loops) # appears to be exactly correct, tracing the bottom edge in a ring containing the pole

            # flag ARs that touch the poles
            l = [vertex[0] for loop in loops for vertexes in loop for vertex in vertexes]
            if 0 in l:
                flags.add('south_pole')
            if len(latitudes)-1 in l:
                flags.add('north_pole')
            # flag ARs that touch the dateline
            l = [vertex[1] for loop in loops for vertexes in loop for vertex in vertexes]
            if 0 in l or len(longitudes)-1 in l:
                flags.add('dateline')

            ARs.append({"coordinates": loops, "raster": raster, "label": label, "sublabels": local_sublabels, "flags": flags})

    # hole identification
    for i, AR in enumerate(ARs):
        for j, sublabel in enumerate(AR['sublabels']): # note loops and local_sublabels correspond by index
            ## mask off everything that isnt this subregion
            local_sublabels_map = numpy.copy(sublabel_map)
            local_sublabels_map[local_sublabels_map != sublabel] = 0 # suppress all other subregions
            local_sublabels_map[local_sublabels_map == sublabel] = 1 # change to binary flag indicating this lone subregion
            
            ## invert
            b = [[1-y for y in x] for x in local_sublabels_map]
            
            ## identify the exterior of the subregion
            ocean_mask = scipy.ndimage.label(b, structure=[[0,1,0],[1,1,1],[0,1,0]])[0]
            if numpy.array_equal(numpy.unique(ocean_mask), [0,1]):
                # no holes, bail out now
                continue
            for y in range(ocean_mask.shape[0]):
                if ocean_mask[y, 0] > 0 and ocean_mask[y, -1] > 0:
                    ocean_mask[ocean_mask == ocean_mask[y, -1]] = ocean_mask[y, 0]
            values, counts = numpy.unique(ocean_mask, return_counts=True)
            most_common = values[numpy.argmax(counts)]
            ocean_mask[ocean_mask != most_common] = 0
            ocean_mask[ocean_mask == most_common] = 1

            ## label the holes periodically, mask off exterior
            holes = scipy.ndimage.label(b, structure=[[0,1,0],[1,1,1],[0,1,0]])[0] # no diagonal contiguity == don't need to pick apart nested loops
            holes[ocean_mask == 1] = 0
            for y in range(holes.shape[0]):
                if holes[y, 0] > 0 and holes[y, -1] > 0:
                    holes[holes == holes[y, -1]] = holes[y, 0]

            ## 'holes' adjacent to the poles arent holes, they're boundaries
            southholes = numpy.unique(holes[0])
            for s in southholes:
                holes[holes==s] = 0
            northholes = numpy.unique(holes[-1])
            for n in northholes:
                holes[holes==n] = 0

            ## trace boundaries of each hole and add to geojson
            for h in numpy.unique(holes):
                if h == 0:
                    continue
                else:
                    vertexes = gridtools.trace_shape(holes, h)[0]
                    ARs[i]['coordinates'][j].append(vertexes)
                    ARs[i]['flags'].add('holes')

    # straight runs need only the first and last point
    for i, AR in enumerate(ARs):
        for j, loop in enumerate(AR['coordinates']):
            for k, poly in enumerate(loop):
                reduced_poly = [poly[0]]
                for v in range(1,len(poly)-1):
                    if poly[v][0] == reduced_poly[-1][0] and poly[v][0] == poly[v+1][0]:
                        continue
                    elif poly[v][1] == reduced_poly[-1][1] and poly[v][1] == poly[v+1][1]:
                        continue
                    else:
                        reduced_poly.append(poly[v])
                reduced_poly.append(poly[-1])
                # a polygon that loops the entire planet in a straight line will at this point be reduced to just the 'starting' point,
                # need to reinject another point from the line to make valid geojson
                if len(reduced_poly) == 2:
                    reduced_poly.insert(1, poly[math.floor(len(poly)/2)])
                    reduced_poly.insert(2, poly[math.floor(len(poly)/2)+1])
                ARs[i]['coordinates'][j][k] = reduced_poly

    # map indexes back onto real locations
    ARs = [ {   '_id': f'{cal_year}{cal_mon}{cal_day}{cal_hour}_{AR["label"]}' , 
                'timestamp': datetime.datetime(year=int(cal_year), month=int(cal_mon), day=int(cal_day), hour=int(cal_hour), minute=int((cal_hour*60) % 60), second=int((cal_hour*3600) % 60) ), 
                'raster': AR['raster'], 
                'flags': list(AR['flags']),
                'geolocation': {"type": "MultiPolygon", "coordinates": [[[index2coords(index, longitudes, latitudes) for index in poly] for poly in loop] for loop in AR['coordinates']]}
            } for AR in ARs]

    #print(ARs)

    for AR in ARs:
        print(AR)


        # try:
        #     db.arShapes.insert_one(AR)
        # except:
        #     del AR['raster']
        #     del AR['timestamp']
        #     print(AR)


        # try:
        #     db.arShapes.insert_one(AR)
        # except Exception as e:
        #     retries = 0
        #     while retries < 10:
        #         try: 
        #             # failing due to https://jira.mongodb.org/browse/SERVER-52928
        #             print(f'trying to write AR {AR["_id"]} with noise')
        #             AR_noise = copy.deepcopy(AR)
        #             AR_noise['geolocation'] = add_noise(AR_noise['geolocation'])
        #             if 'flags' in AR_noise:
        #                 AR_noise['flags'].append('noise_added')
        #             else:
        #                 AR_noise['flags'] = ['noise_added']
        #             db.arShapes.insert_one(AR_noise)
        #             retries = 9999
        #         except:
        #             retries += 1
        #     if retries < 9999:
        #         print(f'unable to write AR {AR["_id"]}')

    # #xx = {'geolocation': {'type': 'MultiPolygon', 'coordinates': [[[[-170.3125, -45.75], [-177.8125, -45.75], [-170.3125, -45.75], [-170.3125, -46.25], [-167.1875, -46.25], [-167.1875, -46.75], [-164.0625, -46.75], [-164.0625, -47.25], [-162.1875, -47.25], [-162.1875, -47.75], [-160.9375, -47.75], [-160.9375, -48.25], [-160.3125, -48.25], [-160.3125, -48.75], [-158.4375, -48.75], [-158.4375, -49.25], [-157.1875, -49.25], [-157.1875, -49.75], [-155.9375, -49.75], [-155.9375, -50.25], [-154.6875, -50.25], [-154.6875, -50.75], [-153.4375, -50.75], [-153.4375, -51.25], [-152.8125, -51.25], [-152.8125, -51.75], [-150.9375, -51.75], [-150.9375, -52.25], [-149.6875, -52.25], [-149.6875, -52.75], [-148.4375, -52.75], [-148.4375, -53.25], [-147.8125, -53.25], [-147.8125, -53.75], [-146.5625, -53.75], [-146.5625, -54.25], [-145.3125, -54.25], [-145.3125, -54.75], [-144.0625, -54.75], [-144.0625, -55.25], [-142.8125, -55.25], [-142.8125, -55.75], [-140.9375, -55.75], [-140.9375, -56.25], [-139.6875, -56.25], [-139.6875, -56.75], [-137.8125, -56.75], [-137.8125, -57.25], [-136.5625, -57.25], [-136.5625, -57.75], [-134.0625, -57.75], [-134.0625, -58.25], [-132.8125, -58.25], [-132.8125, -58.75], [-131.5625, -58.75], [-131.5625, -59.25], [-130.3125, -59.25], [-130.3125, -59.75], [-128.4375, -59.75], [-128.4375, -60.25], [-127.8125, -60.25], [-127.8125, -60.75], [-126.5625, -60.75], [-126.5625, -61.25], [-125.9375, -61.25], [-125.9375, -61.75], [-125.3125, -61.75], [-125.3125, -62.25], [-124.6875, -62.25], [-124.6875, -62.75], [-124.0625, -62.75], [-124.0625, -63.75], [-123.4375, -63.75], [-123.4375, -65.25], [-122.8125, -65.25], [-122.8125, -67.75], [-122.1875, -67.75], [-122.1875, -68.25], [-119.6875, -68.25], [-119.6875, -67.75], [-117.1875, -67.75], [-117.1875, -66.75], [-116.5625, -66.75], [-116.5625, -64.25], [-117.1875, -64.25], [-117.1875, -62.75], [-117.8125, -62.75], [-117.8125, -62.25], [-118.4375, -62.25], [-118.4375, -61.25], [-119.0625, -61.25], [-119.0625, -60.25], [-119.6875, -60.25], [-119.6875, -58.25], [-120.3125, -58.25], [-120.3125, -57.75], [-119.6875, -57.75], [-119.6875, -55.75], [-120.3125, -55.75], [-120.3125, -54.25], [-120.9375, -54.25], [-120.9375, -53.25], [-121.5625, -53.25], [-121.5625, -52.25], [-122.1875, -52.25], [-122.1875, -51.75], [-122.8125, -51.75], [-122.8125, -50.75], [-123.4375, -50.75], [-123.4375, -50.25], [-124.0625, -50.25], [-124.0625, -49.25], [-124.6875, -49.25], [-124.6875, -48.25], [-125.3125, -48.25], [-125.3125, -47.75], [-125.9375, -47.75], [-125.9375, -47.25], [-127.8125, -47.25], [-127.8125, -46.75], [-129.0625, -46.75], [-129.0625, -46.25], [-130.9375, -46.25], [-130.9375, -45.75], [-134.0625, -45.75], [-134.0625, -45.25], [-134.6875, -45.25], [-134.6875, -45.75], [-135.3125, -45.75], [-135.3125, -45.25], [-137.8125, -45.25], [-137.8125, -45.75], [-147.8125, -45.75], [-147.8125, -46.25], [-150.9375, -46.25], [-150.9375, -46.75], [-152.8125, -46.75], [-152.8125, -46.25], [-155.9375, -46.25], [-155.9375, -45.75], [-158.4375, -45.75], [-158.4375, -45.25], [-175, -45.25], [-175, -45.75], [-170.3125, -45.75]], [[176.5625, -45.25], [176.5625, -45.75], [174.6875, -45.75], [174.6875, -46.25], [172.8125, -46.25], [172.8125, -46.75], [170.9375, -46.75], [170.9375, -47.25], [169.0625, -47.25], [169.0625, -46.75], [167.1875, -46.75], [167.1875, -45.75], [166.5625, -45.75], [166.5625, -45.25], [165.9375, -45.25], [165.9375, -44.75], [165.3125, -44.75], [165.3125, -43.75], [164.0625, -43.75], [164.0625, -43.25], [162.8125, -43.25], [162.8125, -41.75], [162.1875, -41.75], [162.1875, -41.25], [161.5625, -41.25], [161.5625, -40.75], [160.9375, -40.75], [160.9375, -40.25], [160.3125, -40.25], [160.3125, -39.75], [159.6875, -39.75], [159.6875, -38.25], [159.0625, -38.25], [159.0625, -34.25], [157.8125, -34.25], [157.8125, -33.75], [157.1875, -33.75], [157.1875, -32.25], [156.5625, -32.25], [156.5625, -30.75], [157.1875, -30.75], [157.1875, -30.25], [156.5625, -30.25], [156.5625, -29.75], [155.3125, -29.75], [155.3125, -30.25], [154.6875, -30.25], [154.6875, -30.75], [153.4375, -30.75], [153.4375, -30.25], [152.8125, -30.25], [152.8125, -31.25], [151.5625, -31.25], [151.5625, -29.75], [150.9375, -29.75], [150.9375, -29.25], [151.5625, -29.25], [151.5625, -28.75], [150.9375, -28.75], [150.9375, -27.75], [149.6875, -27.75], [149.6875, -27.25], [148.4375, -27.25], [148.4375, -26.75], [147.8125, -26.75], [147.8125, -26.25], [147.1875, -26.25], [147.1875, -25.75], [145.9375, -25.75], [145.9375, -25.25], [145.3125, -25.25], [145.3125, -26.25], [144.6875, -26.25], [144.6875, -26.75], [143.4375, -26.75], [143.4375, -26.25], [142.8125, -26.25], [142.8125, -25.75], [141.5625, -25.75], [141.5625, -25.25], [140.3125, -25.25], [140.3125, -24.75], [137.8125, -24.75], [137.8125, -24.25], [135.3125, -24.25], [135.3125, -24.75], [135.9375, -24.75], [135.9375, -26.25], [138.4375, -26.25], [138.4375, -26.75], [139.6875, -26.75], [139.6875, -27.25], [140.3125, -27.25], [140.3125, -27.75], [142.1875, -27.75], [142.1875, -28.25], [143.4375, -28.25], [143.4375, -28.75], [144.6875, -28.75], [144.6875, -29.25], [145.3125, -29.25], [145.3125, -29.75], [146.5625, -29.75], [146.5625, -30.25], [147.1875, -30.25], [147.1875, -30.75], [147.8125, -30.75], [147.8125, -32.25], [148.4375, -32.25], [148.4375, -34.75], [149.0625, -34.75], [149.0625, -35.25], [149.6875, -35.25], [149.6875, -37.75], [150.9375, -37.75], [150.9375, -39.25], [151.5625, -39.25], [151.5625, -39.75], [152.1875, -39.75], [152.1875, -40.25], [152.8125, -40.25], [152.8125, -40.75], [153.4375, -40.75], [153.4375, -41.75], [154.0625, -41.75], [154.0625, -42.25], [154.6875, -42.25], [154.6875, -43.25], [155.3125, -43.25], [155.3125, -43.75], [154.6875, -43.75], [154.6875, -45.25], [155.3125, -45.25], [155.3125, -46.25], [155.9375, -46.25], [155.9375, -46.75], [156.5625, -46.75], [156.5625, -47.25], [157.1875, -47.25], [157.1875, -49.25], [157.8125, -49.25], [157.8125, -49.75], [158.4375, -49.75], [158.4375, -50.25], [159.0625, -50.25], [159.0625, -50.75], [160.3125, -50.75], [160.3125, -51.25], [162.8125, -51.25], [162.8125, -50.75], [166.5625, -50.75], [166.5625, -50.25], [169.6875, -50.25], [169.6875, -49.75], [170.9375, -49.75], [170.9375, -49.25], [172.8125, -49.25], [172.8125, -48.75], [174.0625, -48.75], [174.0625, -48.25], [174.6875, -48.25], [174.6875, -47.75], [175.9375, -47.75], [175.9375, -47.25], [177.8125, -47.25], [177.8125, -46.75], [179.0625, -46.75], [179.0625, -46.25], [-177.8125, -46.25], [-177.8125, -45.75], [-175, -45.75], [-175, -45.25], [176.5625, -45.25]]]]}}
    # # west half works
    # xx = {'_id': 'west', 'timestamp': None, 'raster': [], 'geolocation': {'type': 'MultiPolygon', 'coordinates': [[[[176.5625, -45.25], [176.5625, -45.75], [174.6875, -45.75], [174.6875, -46.25], [172.8125, -46.25], [172.8125, -46.75], [170.9375, -46.75], [170.9375, -47.25], [169.0625, -47.25], [169.0625, -46.75], [167.1875, -46.75], [167.1875, -45.75], [166.5625, -45.75], [166.5625, -45.25], [165.9375, -45.25], [165.9375, -44.75], [165.3125, -44.75], [165.3125, -43.75], [164.0625, -43.75], [164.0625, -43.25], [162.8125, -43.25], [162.8125, -41.75], [162.1875, -41.75], [162.1875, -41.25], [161.5625, -41.25], [161.5625, -40.75], [160.9375, -40.75], [160.9375, -40.25], [160.3125, -40.25], [160.3125, -39.75], [159.6875, -39.75], [159.6875, -38.25], [159.0625, -38.25], [159.0625, -34.25], [157.8125, -34.25], [157.8125, -33.75], [157.1875, -33.75], [157.1875, -32.25], [156.5625, -32.25], [156.5625, -30.75], [157.1875, -30.75], [157.1875, -30.25], [156.5625, -30.25], [156.5625, -29.75], [155.3125, -29.75], [155.3125, -30.25], [154.6875, -30.25], [154.6875, -30.75], [153.4375, -30.75], [153.4375, -30.25], [152.8125, -30.25], [152.8125, -31.25], [151.5625, -31.25], [151.5625, -29.75], [150.9375, -29.75], [150.9375, -29.25], [151.5625, -29.25], [151.5625, -28.75], [150.9375, -28.75], [150.9375, -27.75], [149.6875, -27.75], [149.6875, -27.25], [148.4375, -27.25], [148.4375, -26.75], [147.8125, -26.75], [147.8125, -26.25], [147.1875, -26.25], [147.1875, -25.75], [145.9375, -25.75], [145.9375, -25.25], [145.3125, -25.25], [145.3125, -26.25], [144.6875, -26.25], [144.6875, -26.75], [143.4375, -26.75], [143.4375, -26.25], [142.8125, -26.25], [142.8125, -25.75], [141.5625, -25.75], [141.5625, -25.25], [140.3125, -25.25], [140.3125, -24.75], [137.8125, -24.75], [137.8125, -24.25], [135.3125, -24.25], [135.3125, -24.75], [135.9375, -24.75], [135.9375, -26.25], [138.4375, -26.25], [138.4375, -26.75], [139.6875, -26.75], [139.6875, -27.25], [140.3125, -27.25], [140.3125, -27.75], [142.1875, -27.75], [142.1875, -28.25], [143.4375, -28.25], [143.4375, -28.75], [144.6875, -28.75], [144.6875, -29.25], [145.3125, -29.25], [145.3125, -29.75], [146.5625, -29.75], [146.5625, -30.25], [147.1875, -30.25], [147.1875, -30.75], [147.8125, -30.75], [147.8125, -32.25], [148.4375, -32.25], [148.4375, -34.75], [149.0625, -34.75], [149.0625, -35.25], [149.6875, -35.25], [149.6875, -37.75], [150.9375, -37.75], [150.9375, -39.25], [151.5625, -39.25], [151.5625, -39.75], [152.1875, -39.75], [152.1875, -40.25], [152.8125, -40.25], [152.8125, -40.75], [153.4375, -40.75], [153.4375, -41.75], [154.0625, -41.75], [154.0625, -42.25], [154.6875, -42.25], [154.6875, -43.25], [155.3125, -43.25], [155.3125, -43.75], [154.6875, -43.75], [154.6875, -45.25], [155.3125, -45.25], [155.3125, -46.25], [155.9375, -46.25], [155.9375, -46.75], [156.5625, -46.75], [156.5625, -47.25], [157.1875, -47.25], [157.1875, -49.25], [157.8125, -49.25], [157.8125, -49.75], [158.4375, -49.75], [158.4375, -50.25], [159.0625, -50.25], [159.0625, -50.75], [160.3125, -50.75], [160.3125, -51.25], [162.8125, -51.25], [162.8125, -50.75], [166.5625, -50.75], [166.5625, -50.25], [169.6875, -50.25], [169.6875, -49.75], [170.9375, -49.75], [170.9375, -49.25], [172.8125, -49.25], [172.8125, -48.75], [174.0625, -48.75], [174.0625, -48.25], [174.6875, -48.25], [174.6875, -47.75], [175.9375, -47.75], [175.9375, -47.25], [177.8125, -47.25], [177.8125, -46.75], [179.0625, -46.75], [179.0625, -46.25], [-177.8125, -46.25], [-177.8125, -45.75], [-175, -45.75], [-175, -45.25], [176.5625, -45.25]]]]}}
    # # east half works
    # xx = {'_id': 'east', 'timestamp': None, 'raster': [], 'geolocation': {'type': 'MultiPolygon', 'coordinates': [[[[-170.3125, -45.75], [-170.3125, -46.25], [-167.1875, -46.25], [-167.1875, -46.75], [-164.0625, -46.75], [-164.0625, -47.25], [-162.1875, -47.25], [-162.1875, -47.75], [-160.9375, -47.75], [-160.9375, -48.25], [-160.3125, -48.25], [-160.3125, -48.75], [-158.4375, -48.75], [-158.4375, -49.25], [-157.1875, -49.25], [-157.1875, -49.75], [-155.9375, -49.75], [-155.9375, -50.25], [-154.6875, -50.25], [-154.6875, -50.75], [-153.4375, -50.75], [-153.4375, -51.25], [-152.8125, -51.25], [-152.8125, -51.75], [-150.9375, -51.75], [-150.9375, -52.25], [-149.6875, -52.25], [-149.6875, -52.75], [-148.4375, -52.75], [-148.4375, -53.25], [-147.8125, -53.25], [-147.8125, -53.75], [-146.5625, -53.75], [-146.5625, -54.25], [-145.3125, -54.25], [-145.3125, -54.75], [-144.0625, -54.75], [-144.0625, -55.25], [-142.8125, -55.25], [-142.8125, -55.75], [-140.9375, -55.75], [-140.9375, -56.25], [-139.6875, -56.25], [-139.6875, -56.75], [-137.8125, -56.75], [-137.8125, -57.25], [-136.5625, -57.25], [-136.5625, -57.75], [-134.0625, -57.75], [-134.0625, -58.25], [-132.8125, -58.25], [-132.8125, -58.75], [-131.5625, -58.75], [-131.5625, -59.25], [-130.3125, -59.25], [-130.3125, -59.75], [-128.4375, -59.75], [-128.4375, -60.25], [-127.8125, -60.25], [-127.8125, -60.75], [-126.5625, -60.75], [-126.5625, -61.25], [-125.9375, -61.25], [-125.9375, -61.75], [-125.3125, -61.75], [-125.3125, -62.25], [-124.6875, -62.25], [-124.6875, -62.75], [-124.0625, -62.75], [-124.0625, -63.75], [-123.4375, -63.75], [-123.4375, -65.25], [-122.8125, -65.25], [-122.8125, -67.75], [-122.1875, -67.75], [-122.1875, -68.25], [-119.6875, -68.25], [-119.6875, -67.75], [-117.1875, -67.75], [-117.1875, -66.75], [-116.5625, -66.75], [-116.5625, -64.25], [-117.1875, -64.25], [-117.1875, -62.75], [-117.8125, -62.75], [-117.8125, -62.25], [-118.4375, -62.25], [-118.4375, -61.25], [-119.0625, -61.25], [-119.0625, -60.25], [-119.6875, -60.25], [-119.6875, -58.25], [-120.3125, -58.25], [-120.3125, -57.75], [-119.6875, -57.75], [-119.6875, -55.75], [-120.3125, -55.75], [-120.3125, -54.25], [-120.9375, -54.25], [-120.9375, -53.25], [-121.5625, -53.25], [-121.5625, -52.25], [-122.1875, -52.25], [-122.1875, -51.75], [-122.8125, -51.75], [-122.8125, -50.75], [-123.4375, -50.75], [-123.4375, -50.25], [-124.0625, -50.25], [-124.0625, -49.25], [-124.6875, -49.25], [-124.6875, -48.25], [-125.3125, -48.25], [-125.3125, -47.75], [-125.9375, -47.75], [-125.9375, -47.25], [-127.8125, -47.25], [-127.8125, -46.75], [-129.0625, -46.75], [-129.0625, -46.25], [-130.9375, -46.25], [-130.9375, -45.75], [-134.0625, -45.75], [-134.0625, -45.25], [-134.6875, -45.25], [-134.6875, -45.75], [-135.3125, -45.75], [-135.3125, -45.25], [-137.8125, -45.25], [-137.8125, -45.75], [-147.8125, -45.75], [-147.8125, -46.25], [-150.9375, -46.25], [-150.9375, -46.75], [-152.8125, -46.75], [-152.8125, -46.25], [-155.9375, -46.25], [-155.9375, -45.75], [-158.4375, -45.75], [-158.4375, -45.25], [-175, -45.25], [-175, -45.75], [-170.3125, -45.75]]]]}}
    # # their powers combined
    # xx = {'_id': 'xxx', 'timestamp': None, 'raster': [], 'geolocation': {'type': 'MultiPolygon', 'coordinates': [[[[-170.3125, -45.75], [-170.3125, -46.25], [-167.1875, -46.25], [-167.1875, -46.75], [-164.0625, -46.75], [-164.0625, -47.25], [-162.1875, -47.25], [-162.1875, -47.75], [-160.9375, -47.75], [-160.9375, -48.25], [-160.3125, -48.25], [-160.3125, -48.75], [-158.4375, -48.75], [-158.4375, -49.25], [-157.1875, -49.25], [-157.1875, -49.75], [-155.9375, -49.75], [-155.9375, -50.25], [-154.6875, -50.25], [-154.6875, -50.75], [-153.4375, -50.75], [-153.4375, -51.25], [-152.8125, -51.25], [-152.8125, -51.75], [-150.9375, -51.75], [-150.9375, -52.25], [-149.6875, -52.25], [-149.6875, -52.75], [-148.4375, -52.75], [-148.4375, -53.25], [-147.8125, -53.25], [-147.8125, -53.75], [-146.5625, -53.75], [-146.5625, -54.25], [-145.3125, -54.25], [-145.3125, -54.75], [-144.0625, -54.75], [-144.0625, -55.25], [-142.8125, -55.25], [-142.8125, -55.75], [-140.9375, -55.75], [-140.9375, -56.25], [-139.6875, -56.25], [-139.6875, -56.75], [-137.8125, -56.75], [-137.8125, -57.25], [-136.5625, -57.25], [-136.5625, -57.75], [-134.0625, -57.75], [-134.0625, -58.25], [-132.8125, -58.25], [-132.8125, -58.75], [-131.5625, -58.75], [-131.5625, -59.25], [-130.3125, -59.25], [-130.3125, -59.75], [-128.4375, -59.75], [-128.4375, -60.25], [-127.8125, -60.25], [-127.8125, -60.75], [-126.5625, -60.75], [-126.5625, -61.25], [-125.9375, -61.25], [-125.9375, -61.75], [-125.3125, -61.75], [-125.3125, -62.25], [-124.6875, -62.25], [-124.6875, -62.75], [-124.0625, -62.75], [-124.0625, -63.75], [-123.4375, -63.75], [-123.4375, -65.25], [-122.8125, -65.25], [-122.8125, -67.75], [-122.1875, -67.75], [-122.1875, -68.25], [-119.6875, -68.25], [-119.6875, -67.75], [-117.1875, -67.75], [-117.1875, -66.75], [-116.5625, -66.75], [-116.5625, -64.25], [-117.1875, -64.25], [-117.1875, -62.75], [-117.8125, -62.75], [-117.8125, -62.25], [-118.4375, -62.25], [-118.4375, -61.25], [-119.0625, -61.25], [-119.0625, -60.25], [-119.6875, -60.25], [-119.6875, -58.25], [-120.3125, -58.25], [-120.3125, -57.75], [-119.6875, -57.75], [-119.6875, -55.75], [-120.3125, -55.75], [-120.3125, -54.25], [-120.9375, -54.25], [-120.9375, -53.25], [-121.5625, -53.25], [-121.5625, -52.25], [-122.1875, -52.25], [-122.1875, -51.75], [-122.8125, -51.75], [-122.8125, -50.75], [-123.4375, -50.75], [-123.4375, -50.25], [-124.0625, -50.25], [-124.0625, -49.25], [-124.6875, -49.25], [-124.6875, -48.25], [-125.3125, -48.25], [-125.3125, -47.75], [-125.9375, -47.75], [-125.9375, -47.25], [-127.8125, -47.25], [-127.8125, -46.75], [-129.0625, -46.75], [-129.0625, -46.25], [-130.9375, -46.25], [-130.9375, -45.75], [-134.0625, -45.75], [-134.0625, -45.25], [-134.6875, -45.25], [-134.6875, -45.75], [-135.3125, -45.75], [-135.3125, -45.25], [-137.8125, -45.25], [-137.8125, -45.75], [-147.8125, -45.75], [-147.8125, -46.25], [-150.9375, -46.25], [-150.9375, -46.75], [-152.8125, -46.75], [-152.8125, -46.25], [-155.9375, -46.25], [-155.9375, -45.75], [-158.4375, -45.75], [-158.4375, -45.25], [-174.99999999, -45.25], [-174.99999999, -45.75], [-170.3125, -45.75]]], [[[176.5625, -45.25], [176.5625, -45.75], [174.6875, -45.75], [174.6875, -46.25], [172.8125, -46.25], [172.8125, -46.75], [170.9375, -46.75], [170.9375, -47.25], [169.0625, -47.25], [169.0625, -46.75], [167.1875, -46.75], [167.1875, -45.75], [166.5625, -45.75], [166.5625, -45.25], [165.9375, -45.25], [165.9375, -44.75], [165.3125, -44.75], [165.3125, -43.75], [164.0625, -43.75], [164.0625, -43.25], [162.8125, -43.25], [162.8125, -41.75], [162.1875, -41.75], [162.1875, -41.25], [161.5625, -41.25], [161.5625, -40.75], [160.9375, -40.75], [160.9375, -40.25], [160.3125, -40.25], [160.3125, -39.75], [159.6875, -39.75], [159.6875, -38.25], [159.0625, -38.25], [159.0625, -34.25], [157.8125, -34.25], [157.8125, -33.75], [157.1875, -33.75], [157.1875, -32.25], [156.5625, -32.25], [156.5625, -30.75], [157.1875, -30.75], [157.1875, -30.25], [156.5625, -30.25], [156.5625, -29.75], [155.3125, -29.75], [155.3125, -30.25], [154.6875, -30.25], [154.6875, -30.75], [153.4375, -30.75], [153.4375, -30.25], [152.8125, -30.25], [152.8125, -31.25], [151.5625, -31.25], [151.5625, -29.75], [150.9375, -29.75], [150.9375, -29.25], [151.5625, -29.25], [151.5625, -28.75], [150.9375, -28.75], [150.9375, -27.75], [149.6875, -27.75], [149.6875, -27.25], [148.4375, -27.25], [148.4375, -26.75], [147.8125, -26.75], [147.8125, -26.25], [147.1875, -26.25], [147.1875, -25.75], [145.9375, -25.75], [145.9375, -25.25], [145.3125, -25.25], [145.3125, -26.25], [144.6875, -26.25], [144.6875, -26.75], [143.4375, -26.75], [143.4375, -26.25], [142.8125, -26.25], [142.8125, -25.75], [141.5625, -25.75], [141.5625, -25.25], [140.3125, -25.25], [140.3125, -24.75], [137.8125, -24.75], [137.8125, -24.25], [135.3125, -24.25], [135.3125, -24.75], [135.9375, -24.75], [135.9375, -26.25], [138.4375, -26.25], [138.4375, -26.75], [139.6875, -26.75], [139.6875, -27.25], [140.3125, -27.25], [140.3125, -27.75], [142.1875, -27.75], [142.1875, -28.25], [143.4375, -28.25], [143.4375, -28.75], [144.6875, -28.75], [144.6875, -29.25], [145.3125, -29.25], [145.3125, -29.75], [146.5625, -29.75], [146.5625, -30.25], [147.1875, -30.25], [147.1875, -30.75], [147.8125, -30.75], [147.8125, -32.25], [148.4375, -32.25], [148.4375, -34.75], [149.0625, -34.75], [149.0625, -35.25], [149.6875, -35.25], [149.6875, -37.75], [150.9375, -37.75], [150.9375, -39.25], [151.5625, -39.25], [151.5625, -39.75], [152.1875, -39.75], [152.1875, -40.25], [152.8125, -40.25], [152.8125, -40.75], [153.4375, -40.75], [153.4375, -41.75], [154.0625, -41.75], [154.0625, -42.25], [154.6875, -42.25], [154.6875, -43.25], [155.3125, -43.25], [155.3125, -43.75], [154.6875, -43.75], [154.6875, -45.25], [155.3125, -45.25], [155.3125, -46.25], [155.9375, -46.25], [155.9375, -46.75], [156.5625, -46.75], [156.5625, -47.25], [157.1875, -47.25], [157.1875, -49.25], [157.8125, -49.25], [157.8125, -49.75], [158.4375, -49.75], [158.4375, -50.25], [159.0625, -50.25], [159.0625, -50.75], [160.3125, -50.75], [160.3125, -51.25], [162.8125, -51.25], [162.8125, -50.75], [166.5625, -50.75], [166.5625, -50.25], [169.6875, -50.25], [169.6875, -49.75], [170.9375, -49.75], [170.9375, -49.25], [172.8125, -49.25], [172.8125, -48.75], [174.0625, -48.75], [174.0625, -48.25], [174.6875, -48.25], [174.6875, -47.75], [175.9375, -47.75], [175.9375, -47.25], [177.8125, -47.25], [177.8125, -46.75], [179.0625, -46.75], [179.0625, -46.25], [-177.8125, -46.25], [-177.8125, -45.75], [-175.00000001, -45.75], [-175.00000001, -45.25], [176.5625, -45.25]]]]}}

    # db.arShapes.insert_one(xx)