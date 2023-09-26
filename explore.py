import xarray, numpy, scipy, sys
from pymongo import MongoClient
from collections import defaultdict

client = MongoClient('mongodb://database/argo')
db = client.argo

xar = xarray.open_dataset('Rutz_ARCatalog_MERRA2_2000.nc')

ntim = 0
nlon = 0
nlat = 0

#print(xar['cal_year'][ntim].to_dict()['data'], xar['cal_mon'][ntim].to_dict()['data'], xar['cal_day'][ntim].to_dict()['data'], xar['cal_hour'][ntim].to_dict()['data'])
#print(xar['longitude'].to_dict()['data'], xar['latitude'].to_dict()['data'])
#print(xar['ARs'][ntim][nlat][nlon].to_dict()['data'])
#print(xar['IVT'][ntim][nlat][nlon].to_dict()['data'])

print(numpy.unique(xar['ARs'][ntim].to_dict()['data']))

#print(len(xar['longitude'].to_dict()['data']))
#print(len(xar['latitude'].to_dict()['data']))

ars = xar['ARs'][ntim].to_numpy()
ars = scipy.ndimage.label(ars, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
#print(ars)
#print(ars.shape)
#print(ars[360][575])

label_image = scipy.ndimage.label(ars, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
for y in range(label_image.shape[0]):
    if label_image[y, 0] > 0 and label_image[y, -1] > 0:
        label_image[label_image == label_image[y, -1]] = label_image[y, 0]
#print(label_image)
#print(numpy.unique(label_image))
#print((label_image == 1).nonzero())
#print('xxxx')
numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%.3g" % x))
print(label_image[30:70, 10:82])

#-------------------------------------

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

# toy [lat][lon] grid for a single timestamp direct from the netcdf
#a = [[1, 1, 1, 1, 0, 0, 0, 0],[0, 0, 1, 1, 0, 0, 0, 0],[0, 0, 1, 1, 0, 0, 0, 0],[0,0,0,1,1,0,1,0],[0,0,0,0,0,1,1,0],[0,0,1,1,0,0,0,0],[1,0,1,1,0,0,0,1],[1,1,1,1,0,0,1,1]]
#a = [[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,1,1,1,0,0,0],[0,0,1,0,1,0,0,0],[0,0,1,0,1,0,0,0],[0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1],[0,0,0,0,0,0,1,1]]
a = label_image

# some constants
nlon = len(a[0])
longitudes = xar['longitude'].to_dict()['data']
latitudes = xar['latitude'].to_dict()['data']

# label the clusters and make periodic on longitude boundary
map = scipy.ndimage.label(a, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
for y in range(map.shape[0]):
    if map[y, 0] > 0 and map[y, -1] > 0:
        map[map == map[y, -1]] = map[y, 0]
#print(map)
# get distinct labels
labels = numpy.unique(map)
# identify blobs
ARs = []
for label in labels:
    if label == 0:
        continue
    else:
        # find a northern edge, and take its two top vertexes as the first two boundary vertexes, in order
        cells = numpy.where(map == label)
        vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
        facing = 'E'
        while not numpy.array_equal(vertexes[0], vertexes[-1]):
            # determine which pattern we're in as a function of present vertex and direction
            # make the appropriate move to generate nextvertex, and append it to vertexes
            facing, delta_iLat, delta_iLon = choose_move(label, map, vertexes[-1][0], vertexes[-1][1], facing)
            vertexes.append([vertexes[-1][0]+delta_iLat, (vertexes[-1][1]+delta_iLon)%nlon])
        #print(vertexes)
        # ARs with multiple regions joined only by a single vertex must be written as a list of polygons to be indexed in mongo
        dupes = list(list_duplicates(vertexes))
        dupes.sort(key=loopsort)
        shapes = []
        while len(dupes) > 1:
            shapes.append([vertexes[dupes[0][1][0] : dupes[0][1][1]+1]])
            del vertexes[dupes[0][1][0] : dupes[0][1][1]+1]
            dupes = list(list_duplicates(vertexes))
            dupes.sort(key=loopsort)
        shapes.append([vertexes])

        ARs.append({"type": "MultiPolygon", "coordinates": shapes})

print(ARs)


# # invert the map and do it over again, looking for holes
# b = [[1-y for y in x] for x in a]

# # label the holes, drop the most common since that's just the open ocean, and make periodic on longitude boundary
# holes = scipy.ndimage.label(b, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
# values, counts = numpy.unique(holes, return_counts=True)
# most_common = values[numpy.argmax(counts)]
# holes[holes == most_common] = 0
# for y in range(holes.shape[0]):
#     if holes[y, 0] > 0 and holes[y, -1] > 0:
#         holes[holes == holes[y, -1]] = holes[y, 0]
# #print(holes)
# # get distinct labels
# labels = numpy.unique(holes)
# h = []
# # identify blobs
# for label in labels:
#     if label == 0:
#         continue
#     else:
#         # find a northern edge, and take its two top vertexes as the first two boundary vertexes, in order
#         cells = numpy.where(holes == label)
#         vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
#         facing = 'E'
#         while not numpy.array_equal(vertexes[0], vertexes[-1]):
#             # determine which pattern we're in as a function of present vertex and direction
#             # make the appropriate move to generate nextvertex, and append it to vertexes
#             facing, delta_iLat, delta_iLon = choose_move(label, holes, vertexes[-1][0], vertexes[-1][1], facing)
#             vertexes.append([vertexes[-1][0]+delta_iLat, (vertexes[-1][1]+delta_iLon)%nlon])
#         # identify which AR this hole belongs to
#         hole_min_lat = min([val[0] for val in vertexes])
#         hole_max_lat = max([val[0] for val in vertexes])
#         hole_min_lon = min([val[1] for val in vertexes])
#         hole_max_lon = max([val[1] for val in vertexes])
#         for i, AR in enumerate(ARs):
#             AR_min_lat = min([val[0] for val in AR['coordinates'][0][0][0]])
#             AR_max_lat = max([val[0] for val in AR['coordinates'][0][0][0]])
#             AR_min_lon = min([val[1] for val in AR['coordinates'][0][0][0]])
#             AR_max_lon = max([val[1] for val in AR['coordinates'][0][0][0]])
#             if AR_min_lat < hole_min_lat and AR_max_lat > hole_max_lat and AR_min_lon < hole_min_lon and AR_max_lon > hole_max_lon:
#                 ARs[i]['coordinates'][0][0].append()

# map indexes back onto real locations
# ARs = [ {"type": "MultiPolygon", "coordinates": [[[index2coords(index, longitudes, latitudes) for index in poly] for poly in AR["coordinates"][0]]]} for AR in ARs]

#for AR in ARs:
#    db.blobs.insert_one({'geolocation': AR})
