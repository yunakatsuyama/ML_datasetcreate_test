import xarray as xr
from ea_grid import EAGridRegridder
from cams_grid import CAMSExtractor

# -------------------------
# Read datasets here
# -------------------------
ds_eagrid = xr.open_dataset("Emission_grid_all.nc")

path2cams = "CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_nox_v6.2_yearly_2010.nc"
camsds = xr.open_dataset(path2cams)
cams_extract = CAMSExtractor(camsds) 

target_lon, target_lat = 135.5009 ,34.6913
cams_lon_array, cams_lat_array, cams_emission = cams_extract.extract_window(lon_center=target_lon, lat_center=target_lat, pix = 5)  # WHAT WILL BE HERE ?????)


# -------------------------
# Run regridding
# -------------------------
regridder = EAGridRegridder(
    ds_eagrid,
    cams_lon_array,
    cams_lat_array,
    cams_emission,
    k=10
)

regridder.build_fine_grid()
regridder.extract_eagrid()
regridder.regrid_all()

print("CAMS conserv sum:", regridder.cams_conserv.sum())
print("EAGrid conserv sum:", regridder.eagrid_conserv.sum())