import re
import numpy as np


####################################################  README  ###########################################################################
""" 
Purpose: Counts water atoms in MD frames to determine ligand hydration shell

Usage: Run <python3 LigandHydration_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.		
																																		
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb"

	- This script expects MD frames to only contain ligand and spherical zone of water molecules exported from e.g. Maestro (Schrodinger LCC)

Example: Example with 5 MD frames is provided



Written in geany by André Fischer (2018/2019)	
"""

###############################################  USER DEFINED VARIABLES  ################################################################

nfiles = 5								# Number of files/frames for analysis	
water_type = 'T3P'						# Residue Type of water molecules (e.g. T3P, HOH, H2O, ...)
included_atoms = 'OH'					# included water atoms: write 'OH' for all atoms, write 'O' for only oxygen
water_atom_type = 'HETATM'				# water atom type in PDB file (usually does not have to be changed)
jobname = 'hydration_example1'			# jobname

################################################  EXECUTION  #############################################################################

# Opening Output files
valuefile_name = jobname + '.val'
valuefile = open(valuefile_name,'w')

resfile_name = jobname + '.res'
resfile = open(resfile_name, 'w')

# Array with water ResIDs
water_res_ID = []

# Frameloop
for frame in range(1,(nfiles+1)):
	
	print('--processing file ' + str(frame) + ' --')
	# Opening frame
	curfile_name = str(frame) + '.pdb'
	curfile_open = open(curfile_name, 'r')
	curfile = curfile_open.read()
	#HETATM   58  H1  T3P  AW86      -0.135  11.256  21.378  1.00  0.00           H  
	#
	waterpattern = r'' + re.escape(water_atom_type) + r'(\s+)?\S+\s+[' + re.escape(included_atoms) + r']\d?\s+' + re.escape(water_type) + r'\s+(\S+)'
	waterpattern_compile = re.compile(waterpattern)
	water_matches = waterpattern_compile.finditer(curfile)
	
	for matchcounter1 in water_matches:
		water_res_ID.append(matchcounter1.group(2))
	
	water_res_ID_clean = list(set(water_res_ID))	# Cleans list from duplicates
	print('ligand meets the following water residue IDs: ')
	print(water_res_ID_clean)
	waters_in_frame = len(water_res_ID_clean)		# waters in frame = number of unique ResIDs after cleaning duplicates from array = length of the ResID array
	print('waters in frame: ' + str(waters_in_frame))
	valuefile.write(str(waters_in_frame))			# Appending current waters to file
	valuefile.write('\n')
	
	# Cleaning array after file
	water_res_ID_clean = []
	water_res_ID = []
	
valuefile.close()

## Average calculation from value file

# Opening valuefile
valuefile2_open = open(valuefile_name, 'r')
valuefile2 = valuefile2_open.read().split('\n')

# Array modification
value_array = valuefile2[:(len(valuefile2))-1]	# Cleaning array from last, empty entry
value_array = [ int(val) for val in value_array ]

# Average+SD calculation
sum_values = sum(value_array)
n_values = len(value_array)
average_wat = sum_values / n_values
stdev_wat = np.std(value_array, ddof=1)


#### OUTPUT

# CLI
print('___CALCULATION FINISHED____')
print('In average, the ligand was followed by ' + str(average_wat) + ' water molecules')
print('The standard deviation was: ' + str(stdev_wat))

# Outfile
resfile.write('Calculation for: ' + jobname)
resfile.write('In average, the ligand was followed by ' + str(average_wat) + ' water molecules')
resfile.write('The standard deviation was: ' + str(stdev_wat))


