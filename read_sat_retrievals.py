import numpy as np
import netCDF4
from datetime import datetime, timedelta
import os
from multiprocessing import Pool
import shutil


def analyse_data(run):
   def get_sat_time_validity(sat_file):
      sat_retrieval_data = netCDF4.Dataset(sat_file)
      start_time = sat_retrieval_data.getncattr('time_coverage_start')
      datetime_start_time = datetime.strptime(start_time, '%Y%m%dT%H%M%SZ')
      end_time = sat_retrieval_data.getncattr('time_coverage_end')
      datetime_end_time = datetime.strptime(end_time, '%Y%m%dT%H%M%SZ')
      time_coverage_length = (datetime_end_time - datetime_start_time).seconds
      sat_retrieval_time = datetime_start_time + timedelta(seconds=time_coverage_length/2)
      return sat_retrieval_time
      

   def calc_rmse(sat_file):
      sat_retrieval_data = netCDF4.Dataset(sat_file)
      decoded_time_steps = []
      datetime_time_steps = []
      for time_step in time_steps_fall3d:
         decoded_time_step = ''
         i = 0
         for character in time_step:
            try:
               if i == 2:
                  decoded_time_step += character.decode('utf-8').upper()
               else:
                  decoded_time_step += character.decode('utf-8')
            except AttributeError:
               continue
            i += 1
         decoded_time_steps.append(decoded_time_step.strip(' '))
         datetime_time_steps.append(datetime.strptime(decoded_time_step.strip(' '), '%d%b%Y_%H:%M'))
      # Get FALL3D output time step that is needed
      for i in range(0, len(datetime_time_steps) - 1):
         if datetime_time_steps[i] <= sat_retrieval_time_validity <= datetime_time_steps[i + 1]:
            if datetime_time_steps[i + 1] - sat_retrieval_time_validity < sat_retrieval_time_validity - datetime_time_steps[i]:
               fall3d_output_time_validity = datetime_time_steps[i + 1]
               i_time = i + 1
            else:
               fall3d_output_time_validity = datetime_time_steps[i]
               i_time = i
            break
      ash_mass_loading_fall3d = fall3d_output.variables['tephra_col_mass'][i_time]
      # Sat outputs
      min_sat_loading = 0.5
      lats = sat_retrieval_data.variables['latitude'][:]
      lons = sat_retrieval_data.variables['longitude'][:]
      n_cells_lats = lats.shape[0]
      n_cells_lons = lons.shape[1]
      ash_mass_loading_sat = sat_retrieval_data.variables['ash_mass'][:]
      try:
         ash_mass_loading_unc_sat = sat_retrieval_data.variables['ash_mass_uncertainty'][:] #
         ash_mass_loading_unc_sat_fall3d_domain = np.zeros((ash_mass_loading_fall3d.shape[0], ash_mass_loading_fall3d.shape[1])) #
         uncertainty_available = True
      except KeyError:
         uncertainty_available = False
      lat_min_sat = sat_retrieval_data.getncattr('geospatial_lat_min')
      lat_max_sat = sat_retrieval_data.getncattr('geospatial_lat_max')
      lon_min_sat = sat_retrieval_data.getncattr('geospatial_lon_min')
      lon_max_sat = sat_retrieval_data.getncattr('geospatial_lon_max')
      lat_resolution_sat = (lat_max_sat - lat_min_sat) / n_cells_lats
      lon_resolution_sat = (lon_max_sat - lon_min_sat) / n_cells_lons
      cell_area_sat = lat_resolution_sat * lon_resolution_sat
      conversion_factor = cell_area_sat / cell_area_fall3d
      ash_mass_loading_sat_above_0_ijs = (ash_mass_loading_sat > min_sat_loading).nonzero()
      # Aggiungere estrapolazione degli output di FALL3D nei punti dove il mass loading del satellite è > 0. Calcolare la media su questi punti per FALL3D e sat (a questa aggiungere/sottrarre l'incertezza) e salvarle in un file di time series
      ash_mass_loading_sat_fall3d_domain = np.zeros((ash_mass_loading_fall3d.shape[0], ash_mass_loading_fall3d.shape[1]))
      store_avg_loadings = [] 
      store_unc_avg_loadings = [] #
      for j_lon in range(len(lons_fall3d)):
         for j_lat in range(len(lats_fall3d)):
            store_avg_loadings.append([j_lat, j_lon, 0, 0])
            if uncertainty_available:
               store_unc_avg_loadings.append([j_lat, j_lon, 0, 0]) #
      for i in range(0, len(ash_mass_loading_sat_above_0_ijs[0])):
         lat_point = lats[ash_mass_loading_sat_above_0_ijs[0][i]][ash_mass_loading_sat_above_0_ijs[1][i]]
         lon_point = lons[ash_mass_loading_sat_above_0_ijs[0][i]][ash_mass_loading_sat_above_0_ijs[1][i]]
         j_lat = -1000
         j_lon = -1000
         for j_lat_temp in range(0, len(lats_fall3d) - 1):
            if lats_fall3d[j_lat_temp] <= lat_point <= lats_fall3d[j_lat_temp + 1]:
               j_lat = j_lat_temp
               break
         for j_lon_temp in range(0, len(lons_fall3d) - 1):
            if lons_fall3d[j_lon_temp] <= lon_point <= lons_fall3d[j_lon_temp + 1]:
               j_lon = j_lon_temp
               break   
         if j_lat >= 0 and j_lon >= 0:
            for k in range(0, len(store_avg_loadings)):
               if store_avg_loadings[k][0] == j_lat and store_avg_loadings[k][1] == j_lon:
                  store_avg_loadings[k][2] += ash_mass_loading_sat[ash_mass_loading_sat_above_0_ijs[0][i]][ash_mass_loading_sat_above_0_ijs[1][i]]
                  store_avg_loadings[k][3] += 1
                  if uncertainty_available:
                     store_unc_avg_loadings[k][2] += ash_mass_loading_unc_sat[ash_mass_loading_sat_above_0_ijs[0][i]][ash_mass_loading_sat_above_0_ijs[1][i]]
                     store_unc_avg_loadings[k][3] += 1
      for i in range(0, len(store_avg_loadings)):
         if store_avg_loadings[i][2] > 0:
            ash_mass_loading_sat_fall3d_domain[store_avg_loadings[i][0]][store_avg_loadings[i][1]] = store_avg_loadings[i][2] / store_avg_loadings[i][3]
            if uncertainty_available:
               ash_mass_loading_unc_sat_fall3d_domain[store_unc_avg_loadings[i][0]][store_unc_avg_loadings[i][1]] = store_unc_avg_loadings[i][2] / store_unc_avg_loadings[i][3]
      deviations = (ash_mass_loading_sat_fall3d_domain - ash_mass_loading_fall3d) ** 2 
      deviations = np.where(deviations < 1.e-10, 0, deviations)
      rmse_file = (np.sum(np.where(deviations < 1.e-10, 0, deviations)) / np.where(deviations > 0)[0].shape[0]) ** 0.5   
      ash_mass_loading_sat_fall3d_domain_above_0_ijs = (ash_mass_loading_sat_fall3d_domain > min_sat_loading).nonzero() 
      avg_ash_fall3d_output = np.ma.average(np.ma.masked_where(ash_mass_loading_sat_fall3d_domain < min_sat_loading, ash_mass_loading_fall3d))
      avg_ash_sat_output = np.ma.average(np.ma.masked_where(ash_mass_loading_sat_fall3d_domain < min_sat_loading, ash_mass_loading_sat_fall3d_domain))
      if uncertainty_available:
         avg_ash_sat_output_plus_uncertainty = np.ma.average(np.ma.masked_where(ash_mass_loading_sat_fall3d_domain < min_sat_loading, np.add(ash_mass_loading_sat_fall3d_domain, ash_mass_loading_unc_sat_fall3d_domain)))
         avg_ash_sat_output_minus_uncertainty = np.ma.average(np.ma.masked_where(ash_mass_loading_sat_fall3d_domain < min_sat_loading, np.subtract(ash_mass_loading_sat_fall3d_domain, ash_mass_loading_unc_sat_fall3d_domain)))
      else:
         avg_ash_sat_output_plus_uncertainty = 0.0
         avg_ash_sat_output_minus_uncertainty = 0.0 
      sat_retrieval_to_plot = sat_file.split(os.sep)[-1]
      if 'Aqua' in sat_file.split(os.sep):
         sat_retrieval_to_plot = 'Aqua-' + sat_retrieval_to_plot
      else:
         sat_retrieval_to_plot = 'Terra-' + sat_retrieval_to_plot
      sat_retrieval_to_plot = os.path.join(sat_retrievals_maps_folder, sat_retrieval_to_plot)
      with netCDF4.Dataset(sat_retrieval_to_plot, 'w', 'NETCDF4') as ncfile:
         ncfile.createDimension('lat', ash_mass_loading_sat_fall3d_domain.shape[0])
         lat = ncfile.createVariable('lat', np.float32, ('lat'))
         lat[:] = lats_fall3d
         ncfile.createDimension('lon', ash_mass_loading_sat_fall3d_domain.shape[1])
         lon = ncfile.createVariable('lon', np.float32, ('lon'))
         lon[:] = lons_fall3d
         data_var = ncfile.createVariable('ash_loading', np.float64, ('lat', 'lon'))
         data_var[:] = ash_mass_loading_sat_fall3d_domain
         devs = ncfile.createVariable('deviations', np.float64, ('lat', 'lon'))
         devs[:] = deviations
      test_output = netCDF4.Dataset(sat_retrieval_to_plot)
      return rmse_file, avg_ash_fall3d_output, avg_ash_sat_output, avg_ash_sat_output_minus_uncertainty, avg_ash_sat_output_plus_uncertainty

   os.chdir(run)
   for item in os.listdir(run):
      if item.endswith('.csv'):
         os.remove(os.path.join(run, item))
   output_aqua_lut = open(os.path.join(run, 'comparison_aqua_lut.csv'), 'a')
   output_aqua_vpr = open(os.path.join(run, 'comparison_aqua_vpr.csv'), 'a')
   output_terra_lut = open(os.path.join(run, 'comparison_terra_lut.csv'), 'a')
   output_terra_vpr = open(os.path.join(run, 'comparison_terra_vpr.csv'), 'a')
   fall3d_output = netCDF4.Dataset(os.path.join(run, 'eyja2010.res.nc'))
   lats_fall3d = fall3d_output.variables['lat'][:]
   lons_fall3d = fall3d_output.variables['lon'][:]
   time_steps_fall3d = fall3d_output.variables['date'][:]
   lat_resolution_fall3d = 0.25
   lon_resolution_fall3d = 0.25
   cell_area_fall3d = lat_resolution_fall3d * lon_resolution_fall3d
   for sat_retrieval_file in aqua_lut_files:
      sat_retrieval_time_validity = get_sat_time_validity(sat_retrieval_file)
      sat_retrieval_time_validity_str = datetime.strftime(sat_retrieval_time_validity, '%Y%m%d-%H:%M:%S')
      rmse, avg_fall3d_output, avg_sat_output, avg_sat_output_minus_uncertainty, avg_sat_output_plus_uncertainty = calc_rmse(sat_retrieval_file)
      output_aqua_lut.write(sat_retrieval_time_validity_str + ',' + '{:.3f}'.format(rmse)+ ',' + '{:.3f}'.format(avg_fall3d_output) + ',' + '{:.3f}'.format(avg_sat_output) + ',' + '{:.3f}'.format(avg_sat_output_minus_uncertainty) + ',' + '{:.3f}'.format(avg_sat_output_plus_uncertainty) + '\n')
   for sat_retrieval_file in aqua_vpr_files:
      sat_retrieval_time_validity = get_sat_time_validity(sat_retrieval_file)
      sat_retrieval_time_validity_str = datetime.strftime(sat_retrieval_time_validity, '%Y%m%d-%H:%M:%S')
      rmse, avg_fall3d_output, avg_sat_output, avg_sat_output_minus_uncertainty, avg_sat_output_plus_uncertainty = calc_rmse(sat_retrieval_file)
      output_aqua_vpr.write(sat_retrieval_time_validity_str + ',' + '{:.3f}'.format(rmse)+ ',' + '{:.3f}'.format(avg_fall3d_output) + ',' + '{:.3f}'.format(avg_sat_output) + ',' + '{:.3f}'.format(avg_sat_output_minus_uncertainty) + ',' + '{:.3f}'.format(avg_sat_output_plus_uncertainty) + '\n')
   for sat_retrieval_file in terra_lut_files:
      sat_retrieval_time_validity = get_sat_time_validity(sat_retrieval_file)
      sat_retrieval_time_validity_str = datetime.strftime(sat_retrieval_time_validity, '%Y%m%d-%H:%M:%S')
      rmse, avg_fall3d_output, avg_sat_output, avg_sat_output_minus_uncertainty, avg_sat_output_plus_uncertainty = calc_rmse(sat_retrieval_file)
      output_terra_lut.write(sat_retrieval_time_validity_str + ',' + '{:.3f}'.format(rmse)+ ',' + '{:.3f}'.format(avg_fall3d_output) + ',' + '{:.3f}'.format(avg_sat_output) + ',' + '{:.3f}'.format(avg_sat_output_minus_uncertainty) + ',' + '{:.3f}'.format(avg_sat_output_plus_uncertainty) + '\n')
   for sat_retrieval_file in terra_vpr_files:
      sat_retrieval_time_validity = get_sat_time_validity(sat_retrieval_file)
      sat_retrieval_time_validity_str = datetime.strftime(sat_retrieval_time_validity, '%Y%m%d-%H:%M:%S')
      rmse, avg_fall3d_output, avg_sat_output, avg_sat_output_minus_uncertainty, avg_sat_output_plus_uncertainty = calc_rmse(sat_retrieval_file)
      output_terra_vpr.write(sat_retrieval_time_validity_str + ',' + '{:.3f}'.format(rmse)+ ',' + '{:.3f}'.format(avg_fall3d_output) + ',' + '{:.3f}'.format(avg_sat_output) + ',' + '{:.3f}'.format(avg_sat_output_minus_uncertainty) + ',' + '{:.3f}'.format(avg_sat_output_plus_uncertainty) + '\n')
   output_aqua_lut.close()
   output_aqua_vpr.close()
   output_terra_lut.close()
   output_terra_vpr.close()

aqua_lut_files = []
aqua_vpr_files = []
terra_lut_files = []
terra_vpr_files = []
fall3d_run_folders = []
root = os.getcwd()
sat_retrievals_maps_folder = os.path.join(root, 'sat_retrievals_maps')
try:
   os.mkdir(sat_retrievals_maps_folder)
except FileExistsError:
   shutil.rmtree(sat_retrievals_maps_folder)
   os.mkdir(sat_retrievals_maps_folder)
for path, subdirs, files in os.walk(root):
   for name in files:
      if '.nc' in name:
         if 'Aqua' in path:
            if 'lut' in path:
               aqua_lut_files.append(os.path.join(path, name))
            else:
               aqua_vpr_files.append(os.path.join(path, name))
         elif 'Terra' in path:
            if 'lut' in path:
               terra_lut_files.append(os.path.join(path, name))
            else:
               terra_vpr_files.append(os.path.join(path, name))
         elif 'eyja2010' in path:
            if 'res.nc' in name:
               fall3d_run_folders.append(path)

for run_folder in fall3d_run_folders:
    analyse_data(run_folder)