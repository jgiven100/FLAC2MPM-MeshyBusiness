import pandas as pd

# --- INPUTS --- #
in_dir = '~/FLAC2MPM/CoreInputFiles/'
input_csv = 'Example_Dam_FR_SS_Ru_CoreSu_0p05_zone.csv'
output_material = '05_material.txt'
output_stress = '05_stress.txt'

# --- MAIN --- #
# Load the CSV file
input_path = in_dir + input_csv
data = pd.read_csv(input_path)

# Strip any leading/trailing spaces from the column names
data.columns = data.columns.str.strip()

# Extract the required columns for material.txt
# zid  shear_modulus  bulk_modulus  density  c_liq  c_dry
material_columns = ['ZoneNum', 'shear', 'bulk', 'density', 'liqstr', 'drnstr']
material_data = data[material_columns]

# Extract the required columns for stress.txt and reverse the sign of the 'pp' column
#zid  sigma'_xx  sigma'_yy  sigma_xy  -u
stress_columns = ['ZoneNum', 'esxx', 'esyy', 'sxy', 'pp']
stress_data = data[stress_columns].copy()
stress_data['pp'] = -stress_data['pp']

# Save to text files without headers and space delimiter
material_path = in_dir + output_material
stress_path = in_dir + output_stress
material_data.to_csv(material_path, sep=' ', index=False, header=False)
stress_data.to_csv(stress_path, sep=' ', index=False, header=False)
