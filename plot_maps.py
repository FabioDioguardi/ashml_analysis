import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import scipy.ndimage
import os
import shutil

def plot_fall3d_result(test_folder, test, time_step_in):
	map = Basemap(projection='merc',llcrnrlon=-30.,llcrnrlat=40.,urcrnrlon=10.,urcrnrlat=70.,resolution='i') 
	map.drawcoastlines()
	fall3d_output = os.path.join(test_folder, 'eyja2010.res.nc')
	nc_fall3d = Dataset(fall3d_output, 'r')
	lat = nc_fall3d.variables['lat'][:]
	lon = nc_fall3d.variables['lon'][:]
	ash_loading = nc_fall3d.variables['tephra_col_mass'][time_step_in - 1]
	ash_loading_masked = np.ma.masked_where(ash_loading <= 0, ash_loading)
	lons, lats = np.meshgrid(lon, lat)
	x,y = map(lons,lats)
	clevs = np.arange(0,20,2)
	cs_filled_fall3d = map.contourf(x, y, ash_loading_masked, clevs, cmap=mpl.colormaps['Greys'], extend='max')
	cs_fall3d = map.contour(x, y, ash_loading, clevs, colors='black', linewidths=0.5)
	plt.clabel(cs_fall3d, fontsize=6, inline=True, inline_spacing=0) # contour labels
	plt.title('Ash mass loading')
	clb = plt.colorbar(cs_filled_fall3d)
	clb.ax.tick_params(labelsize=8) 
	clb.set_label('ML (g m⁻²)', rotation=90)
	plt.savefig(os.path.join(cwd, 'output_maps', 'fall3d_' + test + '.png'), dpi=600)
	plt.savefig(os.path.join(cwd, 'output_maps', 'fall3d_' + test + '.svg'))
	plt.close('all')


def plot_sat(time_step_in):
	from datetime import datetime, timedelta
	start_time_datetime = datetime.strptime(start_time, '%Y%m%d-%H:%M')
	actual_time = start_time_datetime + timedelta(hours=time_step_in)
	sat_outputs_folder = os.path.join(cwd, 'sat_retrievals_maps')
	sat_outputs = os.listdir(sat_outputs_folder)
	min_time_difference = 1000000
	for output in sat_outputs:
		if 'Aqua' in output and 'LUT' in output:
			time_validity = output.split('EYJAFJALLAJOKULL-')[1]
			time_validity = time_validity.split('-fv1.nc')[0]
			time_validity_datetime = datetime.strptime(time_validity, '%Y%m%d-%H%M%S')
			time_difference = time_validity_datetime - actual_time
			time_difference_in_min = abs(divmod(time_difference.total_seconds(), 60)[0])
			if time_difference_in_min <= min_time_difference:
				min_time_difference = time_difference_in_min
				sat_output_actual_time = os.path.join(sat_outputs_folder, output)
	map = Basemap(projection='merc',llcrnrlon=-30.,llcrnrlat=40.,urcrnrlon=10.,urcrnrlat=70.,resolution='i') 
	map.drawcoastlines()
	nc_sat = Dataset(sat_output_actual_time,'r')
	lat = nc_sat.variables['lat'][:]
	lon = nc_sat.variables['lon'][:]
	ash_loading = nc_sat.variables['ash_loading'][:]
	ash_loading_masked = np.ma.masked_where(ash_loading <= 0, ash_loading)
	lons, lats = np.meshgrid(lon, lat)
	x,y = map(lons,lats)
	clevs = np.arange(0,20,2)
	cs_filled_sat = map.contourf(x, y, ash_loading_masked, clevs, cmap=mpl.colormaps['Greys'], extend='max')
	cs_sat = map.contour(x, y, ash_loading, clevs, colors='black', linewidths=0.5)
	plt.clabel(cs_sat, fontsize=6, inline=True, inline_spacing=0) # contour labels
	plt.title('Ash mass loading')
	clb = plt.colorbar(cs_filled_sat)
	clb.ax.tick_params(labelsize=8) 
	clb.set_label('ML (g m⁻²)', rotation=90)
	plt.savefig(os.path.join(cwd, 'output_maps','sat.png'), dpi=600)
	plt.savefig(os.path.join(cwd, 'output_maps','sat.svg'))
	plt.close('all')


time_step = 638
start_time = '20100414-00:00'
cwd = os.getcwd()
test_folders = ['test1', 'test2']
test_cases = []
test_cases_folders = []
for folder in test_folders:
	for subfolder in os.listdir(os.path.join(cwd, folder)):
		if 'wood_rcorr_0.3' not in subfolder:
			test_cases_folders.append(os.path.join(cwd, folder, subfolder))
			test_cases.append(folder + '_' + subfolder)
try:
	os.mkdir(os.path.join(cwd, 'output_maps'))
except FileExistsError:
	shutil.rmtree(os.path.join(cwd, 'output_maps'))
	os.mkdir(os.path.join(cwd, 'output_maps'))
for i in range(len(test_cases_folders)):
	plot_fall3d_result(test_cases_folders[i], test_cases[i], time_step)
plot_sat(time_step)