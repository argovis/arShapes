import xarray, numpy, scipy, sys, datetime
from pymongo import MongoClient
from collections import defaultdict

client = MongoClient('mongodb://database/argo')
db = client.argo

xar = xarray.open_dataset('Rutz_ARCatalog_MERRA2_2000.nc')

numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(linewidth=200000)

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
    # A B C D are top left, top right, bottom left, bottom right cells around current vertex, oriented upwards (ie to smaller first index) in the matrix
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

    # transform A B C D to match current facing (N to smaller first index, W to smaller second index)
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
    return [longitudes[index[1]], min(-90.0 + 0.5*index[0], 90) ] # half degree bins -90 to 90, associate value with left edge; values in the 90 bin are by definition right on the left edge of the bin, and so are included with the same border as the previous bin

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[str(item)].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)

def loopsort(elt):
    return elt[1][1] - elt[1][0]

def label_features(image, structure=[[1,1,1],[1,1,1],[1,1,1]]):
    # given a 2D grid image[latitude][longitude] labeling features with 1 and voids with 0,
    # label distinct isolated features with a periodic boundary on the inner index 
    labeled_map = scipy.ndimage.label(image, structure=structure)[0]
    for y in range(labeled_map.shape[0]):
        if labeled_map[y, 0] > 0 and labeled_map[y, -1] > 0:
            labeled_map[labeled_map == labeled_map[y, -1]] = labeled_map[y, 0]

    return labeled_map

def trace_shape(labeled_map, label):
    # trace the shape labeled with label in labeled_map

    nlon = len(labeled_map[0])

    # find a northern edge, and take its two top vertexes as the first two boundary vertexes, in order
    cells = numpy.where(labeled_map == label)
    vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
    facing = 'E'
    nrun = 0
    while not numpy.array_equal(vertexes[0], vertexes[-1]):
        # determine which pattern we're in as a function of present vertex and direction
        # make the appropriate move to generate nextvertex, and append it to vertexes
        oldfacing = facing
        facing, delta_iLat, delta_iLon = choose_move(label, labeled_map, vertexes[-1][0], vertexes[-1][1], facing)
        # # straight runs only need the first and last point
        # if facing == oldfacing:
        #     nrun += 1
        # elif nrun > 2:
        #     del vertexes[-1*(nrun-2):-1]
        #     nrun = 0
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
    #print(labeled_map)

    # get distinct labels
    labels = numpy.unique(labeled_map)

    # identify blobs
    ARs = []
    for label in labels:
    #for label in [12]:
        flags = set(())
        if label == 0:
            continue
        else:
            cells = numpy.where(labeled_map == label)
            lats = [latitudes[x] for x in cells[0]]
            lons = [longitudes[x] for x in cells[1]]
            vapors = [ [ivt[cells[0][i]][cells[1][i]]] for i in range(len(cells[0])) ]
            raster = list(zip(lons, lats, vapors))
            
            # re-label just this region, this time not connecting diagonally so diagonally connected regions make separate loops
            local_map = ars[timestep].to_numpy()
            local_map[labeled_map != label] = 0
            local_label_map = label_features(local_map, [[0,1,0],[1,1,1],[0,1,0]])
            local_labels = numpy.unique(local_label_map)
            loops = [[trace_shape(local_label_map, local_label)] for local_label in local_labels if local_label != 0]

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

            # # ARs with multiple regions joined only by a single vertex must be written as a list of polygons to be indexed in mongo
            # dupes = list(list_duplicates(vertexes))
            # dupes.sort(key=loopsort)
            # shapes = []
            # while len(dupes) > 1:
            #     shapes.append([vertexes[dupes[0][1][0] : dupes[0][1][1]+1]])
            #     del vertexes[dupes[0][1][0] : dupes[0][1][1]+1]
            #     vertexes.insert(dupes[0][1][0], [int(x) for x in dupes[0][0].strip('[]').split(', ')])
            #     dupes = list(list_duplicates(vertexes))
            #     dupes.sort(key=loopsort)
            # shapes.append([vertexes])

            ARs.append({"coordinates": loops, "raster": raster, "label": label, "flags": flags})

    # invert the map and do it over again, looking for holes
    b = [[1-y for y in x] for x in ar]
    ## start by identifying the open ocean, including non-AR points diagonally connected to it, and mask them all out
    ocean_mask = scipy.ndimage.label(b, structure=[[1,1,1],[1,1,1],[1,1,1]])[0] 
    values, counts = numpy.unique(ocean_mask, return_counts=True)
    most_common = values[numpy.argmax(counts)]
    ocean_mask[ocean_mask != most_common] = 0
    ocean_mask[ocean_mask == most_common] = 1

    # label the holes, mask off the open ocean, and make periodic on longitude boundary
    holes = scipy.ndimage.label(b, structure=[[0,1,0],[1,1,1],[0,1,0]])[0] # no diagonal contiguity == don't need to pick apart nested loops
    holes[ocean_mask == 1] = 0
    for y in range(holes.shape[0]):
        if holes[y, 0] > 0 and holes[y, -1] > 0:
            holes[holes == holes[y, -1]] = holes[y, 0]

    # print('-------------------')
    # print(holes)

    # get distinct labels
    labels = numpy.unique(holes)
    h = []
    # identify blobs
    for label in labels:
        if label == 0:
            continue
        else:
            vertexes = trace_shape(holes, label)
            placedhole = False
            #print(label)
            # identify which AR this hole belongs to
            for i, AR in enumerate(ARs):
                for j, loop in enumerate(AR['coordinates']):
                    # identify which loop in this AR the hole might belong to
                    ## is there a loop bound directly north of the first hole vertex?
                    northbound = [x for x in loop[0] if x[0] > vertexes[0][0] and x[1] == vertexes[0][1]]
                    #print(AR['label'], northbound)
                    if len(northbound) == 0:
                        continue
                    ## directly south?
                    southbound = [x for x in loop[0] if x[0] < vertexes[0][0] and x[1] == vertexes[0][1]]
                    #print(AR['label'], southbound)
                    if len(southbound) == 0:
                        continue
                    ARs[i]['coordinates'][j].append(vertexes)
                    ARs[i]['flags'].add('holes')
                    placedhole = True
            if not placedhole:
                print(f'warning: didnt place hole in timestamp {timestep}, hole label {label}, scanning on vertex {vertexes[0]}')
                print([index2coords(index, longitudes, latitudes) for index in vertexes])

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
                ARs[i]['coordinates'][j][k] = reduced_poly


    # map indexes back onto real locations
    ARs = [ {   '_id': f'{cal_year}{cal_mon}{cal_day}{cal_hour}_{AR["label"]}' , 
                'timestamp': datetime.datetime(year=int(cal_year), month=int(cal_mon), day=int(cal_day), hour=int(cal_hour), minute=int((cal_hour*60) % 60), second=int((cal_hour*3600) % 60) ), 
                'raster': AR['raster'], 
                'flags': list(AR['flags']),
                'geolocation': {"type": "MultiPolygon", "coordinates": [[[index2coords(index, longitudes, latitudes) for index in poly] for poly in loop] for loop in AR['coordinates']]}
            } for AR in ARs]

    db.blobs.insert_many(ARs)
