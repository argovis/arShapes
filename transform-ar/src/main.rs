use mongodb::{bson::doc, Client, options::{FindOptions, ClientOptions, ResolverConfig}};
use std::env;
use tokio_stream::StreamExt;
use serde::{Deserialize, Serialize};
use chrono::Utc;
use mongodb::bson::{Bson, DateTime};
use chrono::prelude::*;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {

    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct GeoJSONPolygon {
        #[serde(rename = "type")]
        polygon_type: String,
        #[serde(rename = "coordinates")]
        polygon_coordinates: Vec<Vec<Vec<Vec<f64>>>>,
    }

    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct DataSchema {
        _id: String,
        timestamp: DateTime,
        basins: Vec<f64>,
        raster: Option<Vec<(f64, f64, Vec<f64>)>>,
        flags: Vec<String>,
        geolocation: GeoJSONPolygon,
        metadata: Vec<String>,
    }

    // Connect to MongoDB
    let client_uri = env::var("MONGODB_URI").expect("You must set the MONGODB_URI environment var!"); 
    // A Client is needed to connect to MongoDB:
    // An extra line of code to work around a DNS issue on Windows:
    let options =
       ClientOptions::parse_with_resolver_config(&client_uri, ResolverConfig::cloudflare())
          .await?;
    let client = Client::with_options(options)?; 

    let db = client.database("argo");
    let collection = db.collection::<DataSchema>("ar");

    // Find all documents in the collection
    let find_options = FindOptions::builder().build();
    let mut cursor = collection.find(None, find_options).await?;

    // Iterate over the documents
    while let Some(result) = cursor.next().await {
        match result {
            Ok(document) => {
                // Get the 'raster' array from the document
                if let Some(raster) = document.raster {
                    // Extract longitudes, latitudes, and ivt from the 'raster' array

                    let longitudes: Vec<f64> = raster
                        .iter()
                        .map(|(longitude, _, _)| *longitude)
                        .collect();
                    
                    let latitudes: Vec<f64> = raster
                        .iter()
                        .map(|(_, latitude, _)| *latitude)
                        .collect();
                    
                    let ivt_first_elements: Vec<f64> = raster
                        .iter()
                        .filter_map(|(_, _, ivt)| ivt.first().cloned())
                        .collect();                    

                    let mut new_document = doc! {
                        "data": [longitudes, latitudes, ivt_first_elements],
                        "geolocation": {
                            "type": &document.geolocation.polygon_type,
                            "coordinates": &document.geolocation.polygon_coordinates
                        },
                        "_id": document._id.clone(),
                        "basins": document.basins.clone(),
                        "flags": document.flags.clone(),
                        "metadata": document.metadata.clone(),
                        "timestamp": document.timestamp.clone()
                    };

                    // Insert the new document into a new collection
                    let new_collection = db.collection("ar_new");
                    new_collection.insert_one(new_document, None).await?;
                }
            }
            Err(e) => {
                eprintln!("Error: {}", e);
            }
        }
    }

    Ok(())
    }