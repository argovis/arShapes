from pymongo import MongoClient
import time

client = MongoClient('mongodb://database/argo')
db = client.argo

while True:
	# get a random AR
	original = list(db.ar.aggregate([{"$sample": {"size": 1}}]))[0]
	munged = list(db.ar_new.find({"_id":original['_id']}))[0]
	print('Checking AR id ' + str(original['_id']))

	# compare raster to data
	longitudes = [x[0] for x in original['raster']]
	latitudes = [x[1] for x in original['raster']]
	ivts = [x[2][0] for x in original['raster']]

	if [longitudes, latitudes, ivts] != munged['data']:
		print('mismatch on data versus raster')

	if original['timestamp'] != munged['timestamp']:
		print('mismatch on timestamp')

	if original['geolocation'] != munged['geolocation']:
		print('mismatch on geolocation')

	if original['basins'] != munged['basins']:
		print('mismatch on basins')

	if original['flags'] != munged['flags']:
		print('mismatch on flags')
	
	if original['metadata'] != munged['metadata']:
		print('mismatch on metadata')

	time.sleep(1)

