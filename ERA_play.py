import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
ds=pd.read_csv("/Users/seansantellanes/Documents/water_level_codes/reanalysis-era5-single-levels-timeseries-wavfc_wwfvs.csv")
plt.plot(ds["swh"])
plt.savefig("swh.png")
print(ds)