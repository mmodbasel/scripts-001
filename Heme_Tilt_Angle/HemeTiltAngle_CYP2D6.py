import re
import numpy as np

####################################################  README  ###########################################################################
""" 
Purpose: This script is designed to determine the heme tilt angle in a CYP Protein in MD Frames	
																																		
Usage: Run <python3 HemeTiltAngle_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.
						
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb" that are aligned with protein in center (membrane perpendicular to z-axis)

Example: Example with 5 MD frames is provided



Written in gedit by André Fischer (2018/2019)				
																																		
"""

################################################  EXAMPLE INPUT  ########################################################################
""" 
HETATM57644  NA  HEM A 800      -2.646  -0.822  -2.803  1.00 20.73           N1-
HETATM57645  NB  HEM A 800      -0.692   0.068  -4.890  1.00 20.27           N  
HETATM57646  NC  HEM A 800      -2.575   2.371  -5.588  1.00 21.77           N1-
HETATM57647  ND  HEM A 800      -4.636   1.268  -3.716  1.00 21.26           N  
HETATM57648 FE   HEM A 800      -2.626   0.664  -4.347  1.00 18.76          Fe2+
"""                                    				
#########################################################################################################################################


###############################################  USER DEFINED VARIABLES  ################################################################

nfiles = 5					# Number of files to process (frames)
jobname = "tilt_example1"			# Jobname
chain = "A"					# Chain of protein (e.g. "A")
hemnum = 800					# Residue number of HEM

################################################  EXECUTION  #############################################################################

# Terminal Intro
print('Performing heme tilt angle calculation for ' + jobname)
print('---')

# Creating output files
valuefile_basename = jobname + '.val'
valuefile = open(valuefile_basename, "w")

logfile_basename = jobname + '.log'
logfile = open(logfile_basename, "w")

# Creating array for result values (average calculation)
tilt_array = []

# Loop for all files
for i in range(1,(nfiles+1)):
	# Opening file
	filename = str(i) + '.pdb'
	currentfile_open = open(filename,"r")
	curfile = currentfile_open.read()
	
	print('processing file ' + str(i))

	# Saving HEM NA coordinates
	NA_pattern = r'HETATM\s?\d+\s+' + re.escape("NA") +  r'\s+HEM\s+' + re.escape(chain) + r'\s+' + re.escape(str(hemnum)) + r'\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
	NA_pattern_compile = re.compile(NA_pattern)
	NA_matches = NA_pattern_compile.finditer(curfile)
	for NA_counter in NA_matches:
		NA_x = float(NA_counter.group(1))
		NA_y = float(NA_counter.group(2))
		NA_z = float(NA_counter.group(3))
		print('NA Coords: ')
		print(NA_x, NA_y, NA_z)
	
	# Saving HEM NB coordinates
	NB_pattern = r'HETATM\s?\d+\s+' + re.escape("NB") +  r'\s+HEM\s+' + re.escape(chain) + r'\s+' + re.escape(str(hemnum)) + r'\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
	NB_pattern_compile = re.compile(NB_pattern)
	NB_matches = NB_pattern_compile.finditer(curfile)
	for NB_counter in NB_matches:
		NB_x = float(NB_counter.group(1))
		NB_y = float(NB_counter.group(2))
		NB_z = float(NB_counter.group(3))	
		print('NB Coords: ')
		print(NB_x, NB_y, NB_z)
	
	# Saving HEM NC coordinates
	NC_pattern = r'HETATM\s?\d+\s+' + re.escape("NC") +  r'\s+HEM\s+' + re.escape(chain) + r'\s+' + re.escape(str(hemnum)) + r'\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
	NC_pattern_compile = re.compile(NC_pattern)
	NC_matches = NC_pattern_compile.finditer(curfile)
	for NC_counter in NC_matches:
		NC_x = float(NC_counter.group(1))
		NC_y = float(NC_counter.group(2))
		NC_z = float(NC_counter.group(3))
		print('NC Coords: ')
		print(NC_x, NC_y, NC_z)
		
	# Initial vectors
	vAB_x = NB_x - NA_x
	vAB_y = NB_y - NA_y
	vAB_z = NB_z - NA_z
	
	print('first vector: ')
	print(vAB_x, vAB_y, vAB_z)
	
	vAC_x = NC_x - NA_x
	vAC_y = NC_y - NA_y
	vAC_z = NC_z - NA_z
	
	print('second vector: ')
	print(vAC_x, vAC_y, vAC_z)
	
	# Crossproduct vector
	cp_x = (vAB_y * vAC_z) - (vAB_z * vAC_y)
	cp_y = (vAB_z * vAC_x) - (vAB_x * vAC_z)
	cp_z = (vAB_x * vAC_y) - (vAB_y * vAC_x)
	
	print('crossP vector: ')
	print(cp_x,cp_y,cp_z)
	
	# Squaring Crossproduct vector
	cp_x_sq = cp_x * cp_x
	cp_y_sq = cp_y * cp_y
	cp_z_sq = cp_z * cp_z
	print('squared crossP vector: ')
	print(cp_x_sq, cp_y_sq, cp_z_sq)
	
	# Betrag of Crossproduct (SQRT of added squares)
	cp_sum = cp_x_sq + cp_y_sq + cp_z_sq
	cp_length = np.sqrt(cp_sum)
	
	print('cp length: ' + str(cp_length))	
	
	# Q5 calculation
	Q5 = cp_z / cp_length
	
	# anti-tilt calculation
	anti_tilt_rad = np.arccos(Q5)				# calculate arccos of Q5
	anti_tilt = np.rad2deg(anti_tilt_rad)			# convert radian value to degrees
	
	print('anti tilt in rad is: ' + str(anti_tilt_rad))
	print('anti tilt is: ' + str(anti_tilt))
	
	# tilt calculation
	tilt = 90 - anti_tilt
	tilt_array.append(float(tilt))						# writing result to array for average calculation
	valuefile.write(str(tilt))
	valuefile.write('\n')
	print('The heme tilt angle in file ' + str(i) + ' (' + filename + ') was ' + str(tilt))
	
# average tilt angle
sum_tiltarray = sum(tilt_array)
len_tiltarray = len(tilt_array)
tilt_average  = sum_tiltarray / len_tiltarray

# standard deviation
tilt_stdev = np.std(tilt_array,  ddof=1)

# security check
if len_tiltarray == nfiles:
	pass
else:
	print('Error: total values do not match total files')
	
################################################  OUTPUT  #############################################################################

# Terminal
print('-----')
print('--RESULT--')
print('This run called ' + jobname + ' processed ' + str(nfiles) + ' files')
print('The average heme tilt angle in these frames was ' + str(tilt_average) + ' degrees')
print('The following files have been written: ')
print('   - ' + valuefile_basename + ' (values of tilt angle for all frames)')
print('   - ' + logfile_basename + ' (logfile with average and general information)')

# Logfile
logfile.write('--RESULT--')
logfile.write(2*'\n')
logfile.write('This run called ' + jobname + ' processed ' + str(nfiles) + ' files')
logfile.write('\n')
logfile.write('The average heme tilt angle in these frames was ' + str(tilt_average) + ' degrees')
logfile.write('\n')
logfile.write('The standard deviation was: ' + str(tilt_stdev) + ' degrees')
logfile.write('\n')
	
	
	
	
	
	
