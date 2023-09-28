import xarray, numpy, scipy, sys, datetime
from pymongo import MongoClient
from collections import defaultdict

client = MongoClient('mongodb://database/argo')
db = client.argo

xar = xarray.open_dataset('Rutz_ARCatalog_MERRA2_2000.nc')

numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(linewidth=200)

# helper functions
def transform_facing_and_position(currentFacing, change):
    if change == 'P':
        # proceed
        if currentFacing == 'N':
            return 'N', -1, 0
        elif currentFacing == 'E':
            return 'E', 0, 1
        elif currentFacing == 'S':
            return 'S', 1, 0
        elif currentFacing == 'W':
            return 'W', 0, -1
    elif change == 'L':
        # turn left
        if currentFacing == 'N':
            return 'W', 0, -1
        elif currentFacing == 'E':
            return 'N', -1, 0
        elif currentFacing == 'S':
            return 'E', 0, 1
        elif currentFacing == 'W':
            return 'S', 1, 0
    elif change == 'R':
        # turn right
        if currentFacing == 'N':
            return 'E', 0, 1
        elif currentFacing == 'E':
            return 'S', 1, 0
        elif currentFacing == 'S':
            return 'W', 0, -1
        elif currentFacing == 'W':
            return 'N', -1, 0
    else:
        raise Exception(f'no valid change found {currentFacing}, {change}')

def choose_move(label, map, current_iLat, current_iLon, currentFacing):
    # A B C D are top left, top right, bottom left, bottom right cells around current vertex, oriented to true north
    A_iLat = current_iLat - 1
    if A_iLat < 0:
        A = False
    else:
        A_iLon = (current_iLon - 1)%len(map[0])
        A = map[A_iLat][A_iLon] == label

    B_iLat = current_iLat - 1
    if B_iLat < 0:
        B = False
    else:
        B_iLon = current_iLon
        B = map[B_iLat][B_iLon] == label

    C_iLat = current_iLat
    C_iLon = (current_iLon - 1)%len(map[0])
    if C_iLat < len(map):
        C = map[C_iLat][C_iLon] == label
    else:
        C = False

    D_iLat = current_iLat
    D_iLon = current_iLon
    if D_iLat < len(map):
        D = map[D_iLat][D_iLon] == label
    else:
        D = False

    # transform A B C D to match current facing
    if currentFacing == 'N':
        pass
    elif currentFacing == 'E':
        X = A
        A = B
        B = D
        D = C
        C = X
    elif currentFacing == 'S':
        X = A
        A = D
        D = X
        X = B
        B = C
        C = X
    elif currentFacing == 'W':
        X = A
        A = C
        C = D
        D = B
        B = X

    # determine new center vertex and facing based on A B C D
    if C and not A and not B and not D:
        return transform_facing_and_position(currentFacing, 'L')
    elif D and not A and not B and not C:
        return transform_facing_and_position(currentFacing, 'R')
    elif A and C and not B and not D:
        return transform_facing_and_position(currentFacing, 'P')
    elif B and D and not A and not C:
        return transform_facing_and_position(currentFacing, 'P')
    elif A and D and not B and not C:
        return transform_facing_and_position(currentFacing, 'L')
    elif B and C and not A and not D:
        return transform_facing_and_position(currentFacing, 'R')
    elif A and B and C and not D:
        return transform_facing_and_position(currentFacing, 'R')
    elif A and B and D and not C:
        return transform_facing_and_position(currentFacing, 'L')
    else:
        raise Exception(f'unconsidered option {A} {B} {C} {D}')

def index2coords(index, longitudes, latitudes):
    # index [lat_idx, lon_idx]; return [lon, lat]
    return [longitudes[index[1]], latitudes[index[0]]]

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[str(item)].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)

def loopsort(elt):
    return elt[1][1] - elt[1][0]

def label_features(image):
    # given a 2D grid image[latitude][longitude] labeling features with 1 and voids with 0,
    # label distinct isolated features with a periodic boundary on the inner index 
    labeled_map = scipy.ndimage.label(image, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
    for y in range(labeled_map.shape[0]):
        if labeled_map[y, 0] > 0 and labeled_map[y, -1] > 0:
            labeled_map[labeled_map == labeled_map[y, -1]] = labeled_map[y, 0]

    return labeled_map

def trace_shape(labeled_map, label):
    # trace the shape labeled with label in labeled_map

    nlon = len(labeled_map[0])

    # find a northern edge, and take its two top vertexes as the first two boundary vertexes, in order
    cells = numpy.where(labeled_map == label)
    print(cells)
    vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
    facing = 'E'
    nrun = 0
    while not numpy.array_equal(vertexes[0], vertexes[-1]):
        # determine which pattern we're in as a function of present vertex and direction
        # make the appropriate move to generate nextvertex, and append it to vertexes
        oldfacing = facing
        facing, delta_iLat, delta_iLon = choose_move(label, labeled_map, vertexes[-1][0], vertexes[-1][1], facing)
        # straight runs only need the first and last point
        if facing == oldfacing:
            nrun += 1
        elif nrun > 2:
            del vertexes[-1*(nrun-2):-1]
            nrun = 0
        vertexes.append([vertexes[-1][0]+delta_iLat, (vertexes[-1][1]+delta_iLon)%nlon])

    return vertexes

def convert_hour(time):

    hh = int(time)
    mm = (time*60) % 60
    ss = (time*3600) % 60

    return "%d:%02d.%02d" % (hh, mm, ss)

# unpack the netcdf
longitudes = xar['longitude'].to_dict()['data']
latitudes = xar['latitude'].to_dict()['data']
print(len(latitudes))
cal_mons = xar['cal_mon'].to_dict()['data']
cal_days = xar['cal_day'].to_dict()['data']
cal_years = xar['cal_year'].to_dict()['data']
cal_hours = xar['cal_hour'].to_dict()['data']
ars = xar['ARs']
ivts = xar['IVT']

#for timestep in range(len(cal_years)):
for timestep in [1618]:
    cal_mon = cal_mons[timestep]
    cal_day = cal_days[timestep]
    cal_year = cal_years[timestep]
    cal_hour = cal_hours[timestep]
    hhmmss = convert_hour(cal_hour)
    ar = ars[timestep].to_numpy()
    ivt = ivts[timestep].to_numpy()

    # label the clusters and make periodic on longitude boundary
    labeled_map = label_features(ar)

    # get distinct labels
    labels = numpy.unique(labeled_map)

    # identify blobs
    ARs = []
    # for label in labels:
    for label in [21]:
        flags = set(())
        #if label == 0:
        if label == 0:
            continue
        else:
            cells = numpy.where(labeled_map == label)
            lats = [latitudes[x] for x in cells[0]]
            lons = [longitudes[x] for x in cells[1]]
            vapors = [ [ivt[cells[0][i]][cells[1][i]]] for i in range(len(cells[0])) ]
            raster = list(zip(lons, lats, vapors))
            
            print(min(cells[0]), max(cells[0]))
            print(min(cells[1]), max(cells[1]))
            #print(labeled_map.shape)
            #if timestep == 0 and label == 11:
            if timestep == 1618 and label == 21:
                print(labeled_map[267:,450:514])

            # flag ARs that touch the poles
            vertexes = trace_shape(labeled_map, label)
            l = [x[0] for x in vertexes]
            if 0 in l:
                flags.add('south_pole')
            if len(latitudes)-1 in l:
                flags.add('north_pole')
            # flag ARs that touch the dateline
            l = [x[1] for x in vertexes]
            if 0 in l or len(longitudes)-1 in l:
                flags.add('dateline')

            # ARs with multiple regions joined only by a single vertex must be written as a list of polygons to be indexed in mongo
            dupes = list(list_duplicates(vertexes))
            dupes.sort(key=loopsort)
            shapes = []
            while len(dupes) > 1:
                shapes.append([vertexes[dupes[0][1][0] : dupes[0][1][1]+1]])
                del vertexes[dupes[0][1][0] : dupes[0][1][1]+1]
                vertexes.insert(dupes[0][1][0], [int(x) for x in dupes[0][0].strip('[]').split(', ')])
                dupes = list(list_duplicates(vertexes))
                dupes.sort(key=loopsort)
            shapes.append([vertexes])

            ARs.append({"coordinates": shapes, "raster": raster, "label": label, "flags": flags})

    # invert the map and do it over again, looking for holes
    b = [[1-y for y in x] for x in ar]

    # label the holes, drop the most common since that's just the open ocean, and make periodic on longitude boundary
    holes = scipy.ndimage.label(b, structure=[[0,1,0],[1,1,1],[0,1,0]])[0] # no diagonal contiguity == don't need to pick apart nested loops
    values, counts = numpy.unique(holes, return_counts=True)
    most_common = values[numpy.argmax(counts)]
    holes[holes == most_common] = 0
    for y in range(holes.shape[0]):
        if holes[y, 0] > 0 and holes[y, -1] > 0:
            holes[holes == holes[y, -1]] = holes[y, 0]

    # get distinct labels
    labels = numpy.unique(holes)
    h = []
    # identify blobs
    for label in labels:
        if label == 0:
            continue
        else:
            vertexes = trace_shape(holes, label)

            # identify which AR this hole belongs to
            hole_min_lat = min([val[0] for val in vertexes])
            hole_max_lat = max([val[0] for val in vertexes])
            hole_min_lon = min([val[1] for val in vertexes])
            hole_max_lon = max([val[1] for val in vertexes])
            for i, AR in enumerate(ARs):
                for j, loop in enumerate(AR['coordinates']):
                    AR_min_lat = min([val[0] for val in loop[0]])
                    AR_max_lat = max([val[0] for val in loop[0]])
                    AR_min_lon = min([val[1] for val in loop[0]])
                    AR_max_lon = max([val[1] for val in loop[0]])
                    if AR_min_lat < hole_min_lat and AR_max_lat > hole_max_lat and AR_min_lon < hole_min_lon and AR_max_lon > hole_max_lon:
                        ARs[i]['coordinates'][j].append(vertexes)
                        ARs[i]['flags'].add('holes')

    # map indexes back onto real locations
    ARs = [ {   '_id': f'{cal_year}{cal_mon}{cal_day}{cal_hour}_{AR["label"]}' , 
                'timestamp': datetime.datetime(year=int(cal_year), month=int(cal_mon), day=int(cal_day), hour=int(cal_hour), minute=int((cal_hour*60) % 60), second=int((cal_hour*3600) % 60) ), 
                'raster': AR['raster'], 
                'flags': list(AR['flags']),
                'geolocation': {"type": "MultiPolygon", "coordinates": [[[index2coords(index, longitudes, latitudes) for index in poly] for poly in loop] for loop in AR['coordinates']]}
            } for AR in ARs]

    db.blobs.insert_many(ARs)
