import pandas as pd
import xarray as xr 
import cdsapi 

dataset = "reanalysis-era5-single-levels-timeseries"
request = {
    "variable": ["mean_sea_level_pressure"],#significant_height_of_combined_wind_waves_and_swell
    "location": {"longitude": 166.7820, "latitude": -19.3096},
    "date": ["2023-09-21/2024-12-16"],
    "data_format": "csv"
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()
