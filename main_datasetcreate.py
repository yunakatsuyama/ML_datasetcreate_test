import xarray as xr
from ea_grid import EAGridRegridder
from cams_grid import CAMSExtractor
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cmcrameri import cm
import matplotlib.colors as mcolors

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

print("CAMS before regrid :", regridder.cams_emission.sum())
print("CAMS conserv sum:", regridder.cams_conserv.sum())
print("EAGrid before regrid: ", regridder.emission_grid.sum())
print("EAGrid conserv sum:", regridder.eagrid_conserv.sum())

# ------------------------
# Plot
# ------------------------

mpl.rcParams.update({
    "font.size": 16,          # base size
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 20,
})
log_ticks = np.arange(2.5, 6.6, 1.0)  # 2.5, 3.5, ..., 6.5
tick_labels = [rf"$10^{{{t:.1f}}}$" for t in log_ticks]


fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # 横並び2列、各7x7相当

# ------------------- CAMS before -------------------
ax = axes[0,0]
cmap = cm.roma_r
norm = mcolors.Normalize(vmin=2.5, vmax=6.5)

sc = ax.pcolormesh(cams_lon_array, cams_lat_array, np.log10(cams_emission),
                   cmap=cmap, norm=norm,
                   edgecolors='lightgray', linewidth=0.2)
ax.set_aspect('equal')

ax.set(xlim=[cams_lon_array[0], cams_lon_array[-1]], ylim=[cams_lat_array[0], cams_lat_array[-1]])
cb = plt.colorbar(sc, ax=ax, shrink=0.7, aspect=20, fraction=0.046)
#cb.ax.tick_params(labelsize=7)
cb.set_ticks(log_ticks)
cb.set_ticklabels(tick_labels)
cb.set_label(r'Annual NO$_x$ Emissions (kg km$^{-2}$)')
ax.set_title('(a) Original CAMS')

# ------------------- CAMS after -------------------
ax = axes[0,1]
sc = ax.pcolormesh(regridder.eagrid_lon, regridder.eagrid_lat, np.log10(regridder.cams_conserv),
                   cmap=cmap, norm=norm,
                   edgecolors='lightgray', linewidth=0.2)
ax.set_aspect('equal')
ax.set(xlim=[regridder.eagrid_lon[0], regridder.eagrid_lon[-1]], ylim=[regridder.eagrid_lat[0], regridder.eagrid_lat[-1]])
cb = plt.colorbar(sc, ax=ax, shrink=0.7, aspect=20, fraction=0.046)
#cb.ax.tick_params(labelsize=7)
cb.set_ticks(log_ticks)
cb.set_ticklabels(tick_labels)
cb.set_label(r'Annual NO$_x$ Emissions (kg km$^{-2}$)')
ax.set_title('(b) Regridded CAMS')
# plt.title('CAMS')
# plt.tight_layout()
# plt.show()



# == EAGrids before and after regridding =======
#fig, axes = plt.subplots(1, 2, figsize=(14, 7))  # 横並び2列、各7x7相当

# -------------------EAGRid before -------------------
ax = axes[1,0]
cmap = cm.roma_r
norm = mcolors.Normalize(vmin=2.5, vmax=6.5)

sc = ax.pcolormesh(regridder.lon_e_edges, regridder.lat_e_edges, np.log10(regridder.emission_grid),
                   cmap=cmap, norm=norm,
                   edgecolors='lightgray', linewidth=0.2)
ax.set(xlim=[regridder.lon_e_edges[0], regridder.lon_e_edges[-1]], ylim=[regridder.lat_e_edges[0], regridder.lat_e_edges[-1]])
ax.set_aspect('equal')

cb = plt.colorbar(sc, ax=ax, shrink=0.7, aspect=20, fraction=0.046)
#cb.ax.tick_params(labelsize=7)
cb.set_ticks(log_ticks)
cb.set_ticklabels(tick_labels)
cb.set_label(r'Annual NO$_x$ Emissions (kg km$^{-2}$)')
ax.set_title('(c) Original EAGrid')

# -------------------EAgrid after -------------------
ax = axes[1,1]
sc = ax.pcolormesh(regridder.eagrid_lon, regridder.eagrid_lat, np.log10(regridder.eagrid_conserv),
                   cmap=cmap, norm=norm,
                   edgecolors='lightgray', linewidth=0.2)
ax.set(xlim=[regridder.eagrid_lon[0], regridder.eagrid_lon[-1]], ylim=[regridder.eagrid_lat[0], regridder.eagrid_lat[-1]])
ax.set_aspect('equal')

cb = plt.colorbar(sc, ax=ax, shrink=0.7, aspect=20, fraction=0.046)
cb.set_ticks(log_ticks)
cb.set_ticklabels(tick_labels)
#cb.ax.tick_params(labelsize=7)
cb.set_label(r'Annual NO$_x$ Emissions (kg km$^{-2}$)')
ax.set_title('(d) Regrid EAGrid')

# # == ==== fill =======
# ax = axes[2]
# sc = ax.pcolormesh(eagrid_lon_array, eagrid_lat_array, value_grid,
#                    cmap=cmap, 
#                    edgecolors='lightgray', linewidth=0.2)
# ax.set(xlim=[eagrid_lon_array[0], eagrid_lon_array[-1]], ylim=[eagrid_lat_array[0], eagrid_lat_array[-1]])
# cb = plt.colorbar(sc, ax=ax, shrink=0.7, aspect=20, fraction=0.046)
# cb.ax.tick_params(labelsize=7)
# cb.set_label('after', fontsize=7)
# ax.set_title('After', fontsize=10)

# plt.title('After fill')
plt.tight_layout()
plt.savefig('plots/regridded_tokyo.png')
plt.show()
# print(f'before fill {np.isnan(value_grid).sum()}')
#print(f'after fill {np.isnan(value_grid_fill).sum()}')