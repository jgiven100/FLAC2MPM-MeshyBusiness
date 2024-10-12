import pandas as pd

# --- INPUTS --- #
in_dir = '~/FLAC2MPM/InputGeneration/LSFD/Cal3_St2p0/'
input_csv = 'LSF_Cal3_8per_zone.csv'
# in_dir = '~/FLAC2MPM/InputGeneration/USFD/Cal3_St2p0/'
# input_csv = 'USF_Cal3_8per_zone.csv'

output_material = 'material.txt'
output_stress = 'stress.txt'
col_names = ['ZoneNum', 'cent_x', 'cent_y', 'esyy', 'esxx', 'sxy', 'pp', 'density', 'drnphi', 'drnc', 'cursu', 'ressu', 'str2rem', 'shear', 'bulk', 'descrip']

# --- MAIN --- #
# Load the CSV file with expected column names
input_path = in_dir + input_csv
data = pd.read_csv(input_path, names=col_names, skiprows=1)

# Strip any leading/trailing spaces from the column names
data.columns = data.columns.str.strip()

# Extract the required columns for material.txt
# zid  shear_modulus  bulk_modulus  density  drained_phi drained_cohesion current_su residual_su remaining_pdstrain
material_columns = ['ZoneNum', 'shear', 'bulk', 'density', 'drnphi', 'drnc', 'cursu', 'ressu', 'str2rem']
material_data = data[material_columns]

# Extract the required columns for stress.txt and reverse the sign of the 'pp' column
# zid  sigma'_xx  sigma'_yy  sigma_xy  -u
stress_columns = ['ZoneNum', 'esxx', 'esyy', 'sxy', 'pp']
stress_data = data[stress_columns].copy()
stress_data['pp'] = -stress_data['pp']

# Save to text files without headers and space delimiter
material_path = in_dir + output_material
stress_path = in_dir + output_stress
material_data.to_csv(material_path, sep=' ', index=False, header=False)
stress_data.to_csv(stress_path, sep=' ', index=False, header=False)
