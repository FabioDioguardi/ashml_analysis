import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import shutil

def plot_results(test_number, sat_type, algorithm_in):
	test_folder = 'test' + str(test_number)
	reference_simulation_folder = 'reference'
	if test_number == 1:
		test_simulation_folder = 'db_rcorr_0.3'
	else:
		test_simulation_folder = 'wood_rcorr_0.5'
	output_file = 'comparison_' + sat_type + '_' + algorithm_in + '.csv'
	# ref_solution
	output_file_ref = os.path.join(os.getcwd(), test_folder, reference_simulation_folder, 'with_v2', output_file)
	output_ref_temp = np.genfromtxt(output_file_ref, delimiter=',', dtype=[('S17'),('f8'),('f8'),('f8'),('f8'),('f8')])
	output_ref_temp = np.sort(output_ref_temp,0)
	rows_with_nans = []
	for i in range(len(output_ref_temp)):
		if any(each!=each for each in output_ref_temp[i]):
			rows_with_nans.append(i)
	output_ref = np.delete(output_ref_temp, rows_with_nans, axis=0)
	time_stamps_ref = [output_ref_step[0] for output_ref_step in output_ref]
	fall3d_ml_ref = [output_ref_step[2] for output_ref_step in output_ref]
	sat_ml = [output_ref_step[3] for output_ref_step in output_ref]
	sat_ml_minus_unc = [output_ref_step[4] for output_ref_step in output_ref]
	sat_ml_plus_unc = [output_ref_step[5] for output_ref_step in output_ref]
	# test_solution
	output_file_test = os.path.join(os.getcwd(), test_folder, test_simulation_folder, 'with_v2', output_file)
	output_test_temp = np.genfromtxt(output_file_test, delimiter=',', dtype=[('S17'),('f8'),('f8'),('f8'),('f8'),('f8')])
	output_test_temp = np.sort(output_test_temp,0)
	rows_with_nans = []
	for i in range(len(output_test_temp)):
		if any(each!=each for each in output_test_temp[i]):
			rows_with_nans.append(i)
	output_test = np.delete(output_test_temp, rows_with_nans, axis=0)
	time_stamps_test = [output_test_step[0] for output_test_step in output_test]
	fall3d_ml_test = [output_test_step[2] for output_test_step in output_test]
	fig, ax = plt.subplots()
	ax.plot(fall3d_ml_ref, color='black', label='Reference')
	ax.plot(fall3d_ml_test, color='black', linestyle='dashed', label='Test')
	ax.plot(sat_ml, color='grey', label='Sat average')
	if algorithm_in == 'lut':
		ax.plot(sat_ml_plus_unc, color='grey', linestyle='dotted', label='Sat average - uncertainty')
		ax.plot(sat_ml_minus_unc, color='grey', linestyle='dashed', label='Sat average + uncertainty')
	ax.set(xlabel='time stamp', ylabel='ML (g m⁻²)',
		   title='Test' + str(test_number) + ' ' + sat_type.title() + ' ' + algorithm_in.upper())
	ax.set_xlim([0, 82])
	ax.set_ylim([0, 120])
	ax.legend(loc='best')
	fig.savefig(os.path.join(cwd, 'output_time_series', 'test' + str(test_number) + '_' + sat_type.title() + '_' + algorithm_in.upper() + '.png'))
	fig.savefig(os.path.join(cwd, 'output_time_series', 'test' + str(test_number) + '_' + sat_type.title() + '_' + algorithm_in.upper() + '.svg'))
	with open(os.path.join(cwd, 'output_time_series', 'test' + str(test_number) + '_' + sat_type + '_' + algorithm_in + '_time_stamps_test.csv'), 'w') as time_stamps_file:
		for time_stamp in time_stamps_test:
			time_stamps_file.write(time_stamp.decode() + '\n')
	with open(os.path.join(cwd, 'output_time_series', 'test' + str(test_number) + '_' + sat_type + '_' + algorithm_in + '_time_stamps_ref.csv'), 'w') as time_stamps_file:
		for time_stamp in time_stamps_ref:
			time_stamps_file.write(time_stamp.decode() + '\n')


def plot_rmse(sat_type, algorithm_in):
	output_file = 'comparison_' + sat_type + '_' + algorithm_in + '.csv'
	rmse_test = []
	rmse_ref = []
	for test_number in range(1,3):
		test_folder = 'test' + str(test_number)
		reference_simulation_folder = 'reference'
		if test_number == 1:
			test_simulation_folder = 'db_rcorr_0.3'
		else:
			test_simulation_folder = 'wood_rcorr_0.5'
		# ref_solution
		output_file_ref = os.path.join(os.getcwd(), test_folder, reference_simulation_folder, 'with_v2', output_file)
		output_ref_temp = np.genfromtxt(output_file_ref, delimiter=',', dtype=[('S17'),('f8'),('f8'),('f8'),('f8'),('f8')])
		output_ref_temp = np.sort(output_ref_temp,0)
		rows_with_nans = []
		for i in range(len(output_ref_temp)):
			if any(each!=each for each in output_ref_temp[i]):
				rows_with_nans.append(i)
		output_ref = np.delete(output_ref_temp, rows_with_nans, axis=0)
		time_stamps_ref = [output_ref_step[0] for output_ref_step in output_ref]
		rmse_ref.append([output_ref_step[1] for output_ref_step in output_ref])
		# test_solution
		output_file_test = os.path.join(os.getcwd(), test_folder, test_simulation_folder, 'with_v2', output_file)
		output_test_temp = np.genfromtxt(output_file_test, delimiter=',', dtype=[('S17'),('f8'),('f8'),('f8'),('f8'),('f8')])
		output_test_temp = np.sort(output_test_temp,0)
		rows_with_nans = []
		for i in range(len(output_test_temp)):
			if any(each!=each for each in output_test_temp[i]):
				rows_with_nans.append(i)
		output_test = np.delete(output_test_temp, rows_with_nans, axis=0)
		time_stamps_test = [output_test_step[0] for output_test_step in output_test]
		rmse_test.append([output_test_step[1] for output_test_step in output_test])
	fig, ax = plt.subplots()
	ax.plot(rmse_ref[0], color='grey', label='Test 1 - Reference')
	ax.plot(rmse_test[0], color='grey', linestyle='dashed', label='Test 1 - Test')
	ax.plot(rmse_ref[1], color='black', label='Test 2 - Reference')
	ax.plot(rmse_test[1], color='black', linestyle='dashed', label='Test 2 - Test')
	ax.set(xlabel='time stamp', ylabel='ML (g m⁻²)',
		   title='RMSE' + ' ' + sat_type.title())
	ax.set_xlim([0, 80])
	ax.set_ylim([0, 80])
	ax.legend(loc='best')
	fig.savefig(os.path.join(cwd, 'output_time_series', 'rmse' + '_' + sat_type.title() + '.png'))
	fig.savefig(os.path.join(cwd, 'output_time_series', 'rmse' + '_' + sat_type.title() + '.svg'))
	return np.mean(rmse_ref[0]), np.mean(rmse_test[0]), np.mean(rmse_ref[1]), np.mean(rmse_test[1])

cwd = os.getcwd()
try:
	os.mkdir(os.path.join(cwd, 'output_time_series'))
except FileExistsError:
	shutil.rmtree(os.path.join(cwd, 'output_time_series'))
	os.mkdir(os.path.join(cwd, 'output_time_series'))
sats = ['aqua', 'terra']
algorithms = ['lut', 'vpr']
for i in range(1,3):
	for sat in sats:
		for algorithm in algorithms:
			plot_results(i, sat, algorithm)
avg_rmse_file = open(os.path.join(cwd, 'output_time_series', 'avg_rmse.csv'), 'a')
avg_rmse_file.write('Sat, Algorithm, Test1_ref, Test1_test, Test2_ref, Test2_test\n')
for sat in sats:
	for algorithm in algorithms:	
		avg_rmse_test1_reference, avg_rmse_test1_test, avg_rmse_test2_reference, avg_rmse_test2_test = plot_rmse(sat, algorithm)
		avg_rmse_file.write('{},{},{:.2f},{:.2f},{:.2f},{:.2f}\n'.format(sat, algorithm, avg_rmse_test1_reference, avg_rmse_test1_test, avg_rmse_test2_reference, avg_rmse_test2_test))
		

