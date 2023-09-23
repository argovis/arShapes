import xarray, numpy, scipy

xar = xarray.open_dataset('Rutz_ARCatalog_MERRA2_2000.nc')

ntim = 0
nlon = 0
nlat = 0

print(xar['cal_year'][ntim].to_dict()['data'], xar['cal_mon'][ntim].to_dict()['data'], xar['cal_day'][ntim].to_dict()['data'], xar['cal_hour'][ntim].to_dict()['data'])
print(xar['longitude'][nlon].to_dict()['data'], xar['latitude'][nlat].to_dict()['data'])
print(xar['ARs'][ntim][nlat][nlon].to_dict()['data'])
print(xar['IVT'][ntim][nlat][nlon].to_dict()['data'])

#print(numpy.unique(xar['ARs'][ntim].to_dict()['data']))

#print(len(xar['longitude'].to_dict()['data']))
#print(len(xar['latitude'].to_dict()['data']))

ars = xar['ARs'][0].to_numpy()
ars = scipy.ndimage.label(ars, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
#print(ars)
#print(ars.shape)
#print(ars[360][575])

label_image = scipy.ndimage.label(ars, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
for y in range(label_image.shape[0]):
    if label_image[y, 0] > 0 and label_image[y, -1] > 0:
        label_image[label_image == label_image[y, -1]] = label_image[y, 0]
print(label_image)
print(numpy.unique(label_image))

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
    A_iLat = max(current_iLat - 1, 0)
    A_iLon = (current_iLon - 1)%len(map[0])
    A = map[A_iLat][A_iLon] == label
    print(f'A, {A_iLat}, {A_iLon}')

    B_iLat = max(current_iLat - 1, 0)
    B_iLon = current_iLon
    B = map[B_iLat][B_iLon] == label
    print(f'B, {B_iLat}, {B_iLon}')

    C_iLat = current_iLat
    C_iLon = (current_iLon - 1)%len(map[0])
    C = map[C_iLat][C_iLon] == label
    print(f'C, {C_iLat}, {C_iLon}')

    D_iLat = current_iLat
    D_iLon = current_iLon
    D = map[D_iLat][D_iLon] == label
    print(f'D, {D_iLat}, {D_iLon}')

    # transform A B C D to match current facing
    print(A, B, C, D)
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
    print(A, B, C, D)

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


# some constants
nlon = 8 # number of grid steps in longitude

# toy [lat][lon] grid for a single timestamp direct from the netcdf
a = [[0, 0, 0, 0, 0, 0, 0, 0],[0, 0, 1, 1, 0, 0, 0, 0],[0, 0, 1, 1, 0, 0, 0, 0],[0,0,0,1,1,0,1,0],[0,0,0,0,0,1,1,0],[0,0,1,1,0,0,0,0],[0,0,1,1,0,0,0,0],[0,0,0,0,0,0,0,0]]

# label the clusters and make periodic on longitude boundary
map = scipy.ndimage.label(a, structure=[[1,1,1],[1,1,1],[1,1,1]])[0]
for y in range(map.shape[0]):
    if map[y, 0] > 0 and map[y, -1] > 0:
        map[map == map[y, -1]] = map[y, 0]
print(map)
# get distinct labels
labels = numpy.unique(map)
print(labels)
for label in labels:
    if label == 0:
        continue
    else:
        # find a northern edge, and take its two top vertexes as the first two boundary vertexes, in order
        cells = numpy.where(map == label)
        vertexes = [[cells[0][0],cells[1][0]], [cells[0][0],(cells[1][0]+1) % nlon]]
        facing = 'E'
        print(f'starting {vertexes}, facing {facing}')
        while not numpy.array_equal(vertexes[0], vertexes[-1]):
            # determine which pattern we're in as a function of present vertex and direction
            # make the appropriate move to generate nextvertex, and append it to vertexes
            print('---------------')
            print(facing, vertexes[-1])
            facing, delta_iLat, delta_iLon = choose_move(label, map, vertexes[-1][0], vertexes[-1][1], facing)
            vertexes.append([vertexes[-1][0]+delta_iLat, vertexes[-1][1]+delta_iLon])
        print(label, vertexes)


print(-1%10)
