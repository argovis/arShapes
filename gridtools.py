import scipy, numpy, copy

def label_features(feature_map, structure=[[1,1,1],[1,1,1],[1,1,1]], connected_poles=True, periodic_dateline=True):
    # given a 2D numpy array feature_map[latitude][longitude] labeling features with 1 and voids with 0,
    # label distinct isolated features. 
    # periodic_dateline = True makes a periodic boundary on the inner index
    # connected_poles = True makes all features touching a pole connected
    labeled_map = scipy.ndimage.label(feature_map, structure=structure)[0]
    # periodic boundary
    if periodic_dateline:
        for y in range(labeled_map.shape[0]):
            if labeled_map[y, 0] > 0 and labeled_map[y, -1] > 0:
                labeled_map[labeled_map == labeled_map[y, -1]] = labeled_map[y, 0]
    # connected poles
    if connected_poles:
        spole_features = [x for x in numpy.unique(labeled_map[0]) if x > 0]
        if len(spole_features) > 1:
            for feature in spole_features[1:]:
                labeled_map[labeled_map == feature] = spole_features[0]
        npole_features = [x for x in numpy.unique(labeled_map[-1]) if x > 0]
        if len(npole_features) > 1:
            for feature in npole_features[1:]:
                labeled_map[labeled_map == feature] = npole_features[0]

    return labeled_map

def trace_shape(labeled_map, label, nlatsteps=361):
    # trace the shape labeled with <label> in the map returned by label_features

    nlon = len(labeled_map[0])

    # find a top edge, and take its two top vertexes as the first two boundary vertexes, in order
    cells = numpy.where(labeled_map == label)
    vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
    facing = 'R'
    while not numpy.array_equal(vertexes[0], vertexes[-1]):
        # determine which pattern we're in as a function of present vertex and direction
        # make the appropriate move to generate nextvertex, and append it to vertexes
        oldfacing = facing
        facing, delta_iLat, delta_iLon = choose_move(label, labeled_map, vertexes[-1][0], vertexes[-1][1], facing)
        vertexes.append([vertexes[-1][0]+delta_iLat, (vertexes[-1][1]+delta_iLon)%nlon])

    # if AR wraps around the whole planet in an annulus, the above will only be the top edge; find the bottom edge.
    if [cells[0][-1]+1,cells[1][-1]] not in vertexes:
        n_vertexes = [[cells[0][-1]+1,cells[1][-1]], [cells[0][-1]+1,(cells[1][-1]+1) % nlon]]
        facing = 'R'
        while not numpy.array_equal(n_vertexes[0], n_vertexes[-1]):
            # determine which pattern we're in as a function of present vertex and direction
            # make the appropriate move to generate nextvertex, and append it to vertexes
            oldfacing = facing
            facing, delta_iLat, delta_iLon = choose_move(label, labeled_map, n_vertexes[-1][0], n_vertexes[-1][1], facing)
            n_vertexes.append([n_vertexes[-1][0]+delta_iLat, (n_vertexes[-1][1]+delta_iLon)%nlon])
        # decide who is the blob and who is the hole
        n_vertexes_north = len([x for x in n_vertexes if x[0] > nlatsteps/2])
        if n_vertexes_north / len(n_vertexes):
            # most points are in the northern hemisphere, n_vertexes is the hole
            return [vertexes, n_vertexes]
        else:
            return [n_vertexes, vertexes]
    else:
        return [vertexes]

    return [vertexes]

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

    # transform A B C D to match current facing (U(p) to smaller first index, L(eft) to smaller second index)
    if currentFacing == 'U':
        pass
    elif currentFacing == 'R':
        X = A
        A = B
        B = D
        D = C
        C = X
    elif currentFacing == 'D':
        X = A
        A = D
        D = X
        X = B
        B = C
        C = X
    elif currentFacing == 'L':
        X = A
        A = C
        C = D
        D = B
        B = X

    # determine new center vertex and facing based on A B C D
    if C and not A and not B and not D:
        return transform_facing_and_position(currentFacing, 'turnleft')
    elif D and not A and not B and not C:
        return transform_facing_and_position(currentFacing, 'turnright')
    elif A and C and not B and not D:
        return transform_facing_and_position(currentFacing, 'proceed')
    elif B and D and not A and not C:
        return transform_facing_and_position(currentFacing, 'proceed')
    elif A and D and not B and not C:
        return transform_facing_and_position(currentFacing, 'turnleft')
    elif B and C and not A and not D:
        return transform_facing_and_position(currentFacing, 'turnright')
    elif A and B and C and not D:
        return transform_facing_and_position(currentFacing, 'turnright')
    elif A and B and D and not C:
        return transform_facing_and_position(currentFacing, 'turnleft')
    else:
        raise Exception(f'unconsidered option {A} {B} {C} {D}')

def transform_facing_and_position(currentFacing, change):
    if change == 'proceed':
        # proceed
        if currentFacing == 'U':
            return 'U', -1, 0
        elif currentFacing == 'R':
            return 'R', 0, 1
        elif currentFacing == 'D':
            return 'D', 1, 0
        elif currentFacing == 'L':
            return 'L', 0, -1
    elif change == 'turnleft':
        # turn left
        if currentFacing == 'U':
            return 'L', 0, -1
        elif currentFacing == 'R':
            return 'U', -1, 0
        elif currentFacing == 'D':
            return 'R', 0, 1
        elif currentFacing == 'L':
            return 'D', 1, 0
    elif change == 'turnright':
        # turn right
        if currentFacing == 'U':
            return 'R', 0, 1
        elif currentFacing == 'R':
            return 'D', 1, 0
        elif currentFacing == 'D':
            return 'L', 0, -1
        elif currentFacing == 'L':
            return 'U', -1, 0
    else:
        raise Exception(f'no valid change found {currentFacing}, {change}')

def generate_geojson(labeled_map, label, index2coords, connected_poles=True, periodic_dateline=True):
    # given a map <labeled_map> returned by label_features and the <label> of interest,
    # and a function index2coords that takes [lat_idx, lon_idx] and returns [lon, lat]
    # return a valid geojson MultiPolygon representing the labeled feature.
    # <connected_poles> and <periodic_dateline> should be the same as used when creating <labeled_map>

    flags = set(())
    local_map = copy.deepcopy(labeled_map)
    local_map[labeled_map != label] = 0 # gets rid of other ARs in a the AR binary flag map
    local_sublabels_map = label_features(local_map, structure=[[0,1,0],[1,1,1],[0,1,0]], connected_poles=True, periodic_dateline=True)
    local_sublabels = [x for x in numpy.unique(local_sublabels_map) if x != 0] # these are the rings that belong as top level objects in this blob

    # get the outer loops
    loops = [trace_shape(local_sublabels_map, sublabel) for sublabel in local_sublabels]
    for geo in loops:
        if len(geo) == 2:
            flags.add('annulus')

    # flag blobs that touch the poles
    l = [vertex[0] for loop in loops for vertexes in loop for vertex in vertexes]
    if 0 in l:
        flags.add('first_pole')
    if len(labeled_map)-1 in l:
        flags.add('last_pole')
    # flag ARs that touch the dateline
    l = [vertex[1] for loop in loops for vertexes in loop for vertex in vertexes]
    if 0 in l or len(labeled_map[0])-1 in l:
        flags.add('dateline')

    # put holes in the outer loops
    for j, sublabel in enumerate(local_sublabels): # note loops and local_sublabels correspond by index
        ## mask off everything that isnt this subregion
        subregion_map = numpy.copy(local_sublabels_map)
        subregion_map[local_sublabels_map != sublabel] = 0 # suppress all other subregions
        subregion_map[local_sublabels_map == sublabel] = 1 # change to binary flag indicating this lone subregion
        
        ## invert
        b = [[1-y for y in x] for x in subregion_map]
        
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
                vertexes = trace_shape(holes, h)[0]
                loops[j].append(vertexes)
                flags.add('holes')

    # straight runs need only the first and last point
    for j, loop in enumerate(loops):
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
            loops[j][k] = reduced_poly

    # map indexes back onto real locations and return the geojson
    return {"type": "MultiPolygon", "coordinates": [[[index2coords(index) for index in poly] for poly in loop] for loop in loops]}, flags






