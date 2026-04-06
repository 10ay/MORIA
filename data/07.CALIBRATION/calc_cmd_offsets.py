import pandas as pd
import numpy as np
import pdb

def get_average_mag_offsets(n_calib_stars=5):

	data = []
	with open("VI_HST_ogle_Cal_matches4.dat") as f:
	    for line in f:
	        # Skip empty or comment lines
	        if line.strip() == "" or line.startswith("#"):
	            continue
	        parts = line.strip().split()  # splits by any whitespace
	        data.append(parts)

	# Check the max number of columns
	max_len = max(len(row) for row in data)
	#print(f"Max columns in any row: {max_len}")

	# Pad short rows (if needed)
	for row in data:
	    if len(row) < max_len:
	        row.extend([None] * (max_len - len(row)))

	initial_mags = pd.DataFrame(data)
	initial_mags.columns = ['id', 'x',        'y',     'nmat', 'nbmat', 'I_ogle',  'V_ogle',  'I_hst1',  'I_hst',   'I_hfs',   'Io-Ihfs', 'lg_c2Imx', 'I_hstB',  'V_hst1',  'V_hst',   'V_hfs',   'Vo-Vhfs', 'lg_c2Vmx', 'V_hstB',  'Io-Ih1',  'Io-Ih',   'Vo-Vh1',  'Vo-Vh']


	#print(initial_mags)
	#print("------------")
	#print("HST V and I offsets:")
	#print("------------")


	# Step 1: Read the last 10 non-empty, non-comment lines
	with open("fit_HST_IV_ogle_col.log") as f:
	    lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

	last_10_lines = lines[-n_calib_stars:]

	# Step 2: Split each line on whitespace
	data = [line.split() for line in last_10_lines]

	# Step 3: Pad rows if they have inconsistent lengths (optional but recommended)
	max_len = max(len(row) for row in data)
	for row in data:
	    if len(row) < max_len:
	        row.extend([None] * (max_len - len(row)))

	# Step 4: Create a DataFrame
	calib_mags = pd.DataFrame(data)
	pdb.set_trace()

	calib_mags.columns = ['id',      'V_ogle',   'V_oglem',  'V_og_err', 'chi^2',   'I_hst',   'V_hst',   'V-I']

	#print(calib_mags)


	#Next step is to match between the two arrays, then calculate the magnitude offset for each calibration star:

	# Make sure 'id' is the same type in both DataFrames
	initial_mags['id'] = initial_mags['id'].astype(str)
	calib_mags['id'] = calib_mags['id'].astype(str)

	# Build lookup for I_hst1 from initial_mags
	ihst1_lookup = initial_mags.set_index('id')['I_hst1'].to_dict()

	# Lists to store values
	i_hst_diff = []
	i_hst1_vals = []

	# Loop through calib_mags
	for _, row in calib_mags.iterrows():
	    row_id = row['id']
	    if row['V_hst']!=None:
	    	i_hst_calib = float(row['I_hst'])

	    if row_id in ihst1_lookup:
	        i_hst1 = float(ihst1_lookup[row_id])
	        diff = i_hst_calib - i_hst1
	    else:
	        i_hst1 = None
	        diff = None

	    i_hst1_vals.append(i_hst1)
	    i_hst_diff.append(diff)

	# Add to calib_mags DataFrame
	calib_mags['I_hst1'] = i_hst1_vals
	calib_mags['I_hst_diff'] = i_hst_diff

	# Print desired columns
	#print(calib_mags[['id', 'I_hst1', 'I_hst', 'I_hst_diff']])
	avg_offset_I = np.mean(calib_mags['I_hst_diff'])
	print("avg offset I = ", avg_offset_I)
	# Compute RMS difference
	I_rms_diff = np.sqrt(np.mean(calib_mags['I_hst_diff']))
	print("RMS of I offset = ", I_rms_diff)
	#print("------------")

	# Build lookup for V_hst1 from initial_mags
	vhst1_lookup = initial_mags.set_index('id')['V_hst1'].to_dict()

	# Lists to store values
	v_hst_diff = []
	v_hst1_vals = []

	# Loop through calib_mags
	for _, row in calib_mags.iterrows():
	    row_id = row['id']
	    #pdb.set_trace()
	    if row['V_hst']!=None:
	    	v_hst_calib = float(row['V_hst'])

	    if row_id in vhst1_lookup:
	        v_hst1 = float(vhst1_lookup[row_id])
	        diff = v_hst_calib - v_hst1
	    else:
	        v_hst1 = None
	        diff = None
	    pdb.set_trace()
	    v_hst1_vals.append(v_hst1)
	    v_hst_diff.append(diff)

	# Add to calib_mags DataFrame
	calib_mags['V_hst1'] = v_hst1_vals
	calib_mags['V_hst_diff'] = v_hst_diff

	# Print desired columns
	#print(calib_mags[['id', 'V_hst1', 'V_hst', 'V_hst_diff']])
	avg_offset_V = np.mean(calib_mags['V_hst_diff'])
	print("avg offset V = ", avg_offset_V)
	# Compute RMS difference
	V_rms_diff = np.sqrt(np.mean(calib_mags['V_hst_diff']))
	print("RMS of V offset = ", V_rms_diff)

	return avg_offset_I,avg_offset_V

