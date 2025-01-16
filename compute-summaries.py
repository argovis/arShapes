from pymongo import MongoClient
from bson.son import SON
import datetime, json, copy, re
    
client = MongoClient('mongodb://database/argo')
db = client.argo

def get_timestamp_range(db, collection_name):
    collection = db[collection_name]
             
    # Find the earliest timestamp
    filter = {}
    earliest_doc = collection.find_one(filter, sort=[("timestamp", 1)])
    if earliest_doc and "timestamp" in earliest_doc:
        earliest_timestamp = earliest_doc["timestamp"]
    else:
        return None, None  # Return None if no timestamps are found
 
    # Find the latest timestamp or current time, whichever is earlier
    filter = {}
    latest_doc = collection.find_one(filter, sort=[("timestamp", -1)])
    current_time = datetime.datetime.utcnow()
    
    if latest_doc and "timestamp" in latest_doc:
        latest_timestamp = min(latest_doc["timestamp"], current_time)
    else:
        latest_timestamp = current_time  # If no documents, default to current time
    
    # Convert timestamps to ISO 8601 format
    try:
        earliest_iso = earliest_timestamp.isoformat() + "Z"
        latest_iso = latest_timestamp.isoformat() + "Z"
        return earliest_iso, latest_iso
    except:
        return None, None

startDate, endDate = get_timestamp_range(db, 'ar')
entry = {"metagroups": ["id"], "startDate": startDate, "endDate": endDate}

rldoc = db.summaries.find_one({"_id": 'ratelimiter'})
if rldoc:
    rldoc['metadata']['ar'] = entry
else:
    rldoc = {"_id": "ratelimiter", "metadata": {"ar": entry}}

try:        
    db.summaries.replace_one({"_id": 'ratelimiter'}, rldoc, upsert=True)
except BaseException as err:
    print('error: db write failure')
    print(err)