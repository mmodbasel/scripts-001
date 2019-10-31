import re
import numpy as np

####################################################  README  ###########################################################################
""" 
Purpose: This script is designed to the burying depth of CYP2D6 in the membrane. Calculates the distance between 		

Usage: Run <python3 BuryingDepth_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.
																																		
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb" that are aligned (protein in center, membrane below)	

Example: Example with 5 MD frames is provided



Written in gedit by André Fischer (2018/2019)										
																													
"""

################################################  EXAMPLE PDB INPUT  ####################################################################
""" 
ATOM    403  CA  MET A   1      -9.447 -19.308  52.259  1.00  0.00           C 
ATOM   6409  CA  ARG A 380     -18.375 -12.590   8.149  1.00 54.99           C  													
								......																									
HETATM 8841   C1 POPC    5      23.981  23.955  48.926  1.00  0.00           C    
HETATM11655   C1 POPC   26      19.894  39.747  54.458  1.00  0.00           C  
"""                                    				
#########################################################################################################################################


###############################################  USER DEFINED VARIABLES  ################################################################

nfiles = 5						# Number of files to process (frames)
memtype = "POPC"					# Membrane residue type (in PDB File)
jobname = "bd_example1"					# Your jobname
chain = "A"						# Chain of protein (e.g. "A")

################################################  EXECUTION  #############################################################################

print('Performing burying depth calculation for ' + jobname)
print('---')

# Opening output files
outputfile_basename = jobname + '.out'
outputfile=open(outputfile_basename, "w")

values_basename = jobname + '.val'
valuefile = open(values_basename, "w")

# Creating array with values for average calculation
bd_array = []

# File loop
for i in range(1,(nfiles+1)):
		
	# Reading current PDB file
	filename = str(i) + '.pdb'
	curfile_open = open(filename, 'r')
	curfile = curfile_open.read()
	print('Reading file ' + str(i) + '.pdb')
	
	#  Creating arrays to fill with coordinates and results
	CA_z_array = []
	C1_z_array = []
	
	# Adding POPC C1 atoms to array and calculate Mass center
	C1_pattern = r'HETATM\s?\d+\s+' + re.escape('C1') + r'\s+' + re.escape(memtype) + r'\s+\d+\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'
	C1_compile = re.compile(C1_pattern)
	C1_matches = C1_compile.finditer(curfile)
	
	for C1_count in C1_matches:
		C1_z_array.append(float(C1_count.group(1)))			# Appending z-coords of protein CA atoms into array
		
	sum_C1_array = sum(C1_z_array)
	len_C1_array = len(C1_z_array)
	average_C1_array = sum_C1_array / len_C1_array			# Mass center of Membrane C1 Atoms
	print('The mass center of the POPC C1 atoms is ' + str(average_C1_array))
	
	# Adding Protein CA atoms to array and calculate Mass center
	CA_pattern = r'ATOM\s+\d+\s+CA\s+[A-Z-]+\s+' + re.escape(chain) + r'\s+\d+\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'
	CA_compile = re.compile(CA_pattern)
	CA_matches = CA_compile.finditer(curfile)
	
	for CA_count in CA_matches:
		CA_z_array.append(float(CA_count.group(1)))			# Appending z-coords of protein CA atoms into array
	
	sum_CA_array = sum(CA_z_array)	
	len_C2_array = len(CA_z_array)
	average_CA_array = sum_CA_array / len_C1_array			# Mass center of Protein CA atoms
	print('The mass center of the protein CA atoms is ' + str(average_CA_array))
	
	# Measuring distance between mass centers
	burying_depth = average_C1_array - average_CA_array
	valuefile.write(str(burying_depth))
	valuefile.write('\n')
	bd_array.append(float(burying_depth))					# Appending results into array
	print('The burying depth in the file' + filename + ' is ' + str(burying_depth) + ' Angstroms')
	print('----')
	
# Average
sum_bdarray = sum(bd_array)
len_bdarray = len(bd_array)
average_bd = sum_bdarray / len_bdarray						# Average over all frames

# standard deviation
bd_stdev = np.std(bd_array,  ddof=1)


# Check
if len(bd_array) == nfiles:
	pass
else:
	print("Something went wrong in your calculation (maybe not all files processed)")

################################################  OUTPUT PROCESSING ##########################################################################

# Output Terminal
print(2*'\n')
print('-----RESULTS-----')
print('The calculation included ' + str(len_C1_array) + ' POPC residues and ' + str(len_C2_array) + ' protein residues')
print('The average burying depth in ' + str(nfiles) + ' files was ' + str(average_bd) + ' Angstroms')
print('\n')
print('The following files have been written: ')
print('   - ' + values_basename + ' (values of burying depth for all frames)')
print('   - ' + outputfile_basename + ' (logfile with average and general information)')

# Output file
outputfile.write('-----RESULTS-----')
outputfile.write(3*'\n')
outputfile.write('The calculation included ' + str(len_C1_array) + ' POPC residues and ' + str(len_C2_array) + ' protein residues')
outputfile.write('\n')
outputfile.write('The average burying depth in ' + str(nfiles) + ' files was ' + str(average_bd))
outputfile.write('\n')
outputfile.write('The standard deviation was ' + str(bd_stdev))

