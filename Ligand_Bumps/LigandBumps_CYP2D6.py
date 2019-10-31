import re
import numpy as np

####################################################  README  ###########################################################################
""" 
Purpose: This script is designed to determine the number of contacts between a ligand of interest (LOF) to other ligands in MD frames

Usage: Run <python3 LigandBumps_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.		
																																		
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb"

	- This script expects your ligands to be sequentially numbered	and to be of the same type
	- It is enough if the frames used for the calculation only contain the ligand coordinates
	- Only heavy atoms are included for the calculation

Example: Example with 5 MD frames is provided



Written in gedit by André Fischer (2018/2019)	
"""

################################################  EXAMPLE INPUT  ########################################################################
""" 
MODEL        1
ATOM      1  O1  ACE   901      41.112   8.305  -7.215  1.00  0.00           O  
ATOM      2  O2  ACE   901      35.471   6.880 -10.040  1.00  0.00           O  
ATOM      3  N1  ACE   901      37.501   6.893 -11.116  1.00  0.00           N  
ATOM      4  C1  ACE   901      38.375   7.126  -9.982  1.00  0.00           C  
"""                                    				
#########################################################################################################################################


###############################################  USER DEFINED VARIABLES  ################################################################

nfiles = 5					# number of individual files / MD frames to analyze
jobname = 'ligbump_example1'			# your jobname
ligand_type = 'ACE'				# ligand name in PDB files
ligands_of_interest = [918]			# PDB numbers of ligands of interest
ligand_atomtype = 'ATOM'			# PDB atomtype for ligands to be compared (e.g. ATOM, HETATM)
startlig = 901					# first ligand number in PDB File
nlig = 20					# total number of ligands to be compared
heavy_atom_types = 'CNO' 			# heavy atom types of ligand of interest as in PDB (e.g. N for N1), alphabetic order required


################################################  EXECUTION  #############################################################################

# Opening output files
logfile_basename = jobname + '.log'
logfile = open(logfile_basename, "w")

# Getting endlig nr
endlig = startlig + nlig
	
# File loop (Loop for MD Frames)
for i in range(1,(nfiles+1)):
	print('Processing file: ' + str(i))
	# Opening file
	curfile_basename = str(i) + '.pdb'
	curfile_open = open(curfile_basename, "r")
	curfile = curfile_open.read()
		
### Adding Ligand of interest heavy atom coordinates into array##
	
	# Loop over ligands of interest
	
	for lig_of_int in ligands_of_interest:
				
		# Arrays to fill for current ligand
		lig_x_array = []
		lig_y_array = []
		lig_z_array = []
		alllig_x_array = []
		alllig_y_array = []
		alllig_z_array = []
		dif_x_array = []
		dif_y_array = []
		dif_z_array = []
		vector_length_array = []
	
		# finditer for LOF
		lig_interest_pattern = r'' + re.escape(ligand_atomtype) + r'\s+\d+\s+[' + re.escape(heavy_atom_types) + r']\d\s+' + re.escape(ligand_type) + r'\s+' + re.escape(str(lig_of_int)) + r'\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
		lig_interest_pattern_compile = re.compile(lig_interest_pattern)
		lig_interest_matches = lig_interest_pattern_compile.finditer(curfile)
		
		# adding LOF coords to arrays
		for ligcounter in lig_interest_matches:

			lig_xn = ligcounter.group(1)
			lig_x_array.append(float(lig_xn))	
			
			lig_yn = ligcounter.group(2)
			lig_y_array.append(float(lig_yn))
			
			lig_zn = ligcounter.group(3)
			lig_z_array.append(float(lig_zn))
				
		
		# Creating current all_ligand array without LOF e.g. for ACE 907 LOF  ->  [901, 902, ..., 906, 908, ..., lignum]
		all_lig_nums = []
		for u2 in range(startlig,endlig):
			if int(u2) != int(lig_of_int):
				all_lig_nums.append(u2)
			else:
				pass
		
		# adding all, but LOF into arrays
		for all_lig_c in all_lig_nums:
			
			# finditer for all ligands
			all_but_lig_pattern = r'' + re.escape(ligand_atomtype) + r'\s+\d+\s+[' + re.escape(heavy_atom_types) + r']\d\s+' + re.escape(ligand_type) + r'\s+' + re.escape(str(all_lig_c)) + r'\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
			all_but_lig_pattern_compile = re.compile(all_but_lig_pattern)
			all_but_lig_matches = all_but_lig_pattern_compile.finditer(curfile)
			
			# adding all ligand coordinates (exctept for current LOF)
			for allligcounter in all_but_lig_matches:
				
				alllig_xn = allligcounter.group(1)
				alllig_x_array.append(float(alllig_xn))
				
				alllig_yn = allligcounter.group(2)
				alllig_y_array.append(float(alllig_yn))
				
				alllig_zn = allligcounter.group(3)
				alllig_z_array.append(float(alllig_zn))
		

		# x-coord differences LOF vs. all ligands 
		for cur_x_LOF in lig_x_array:
			for cur_x_all in alllig_x_array:
				dif_x = cur_x_all - cur_x_LOF
				dif_x_array.append(dif_x)
		# y-coord differences LOF vs. all ligands 
		for cur_y_LOF in lig_y_array:
			for cur_y_all in alllig_y_array:
				dif_y = cur_y_all - cur_y_LOF
				dif_y_array.append(dif_y)

		# y-coord differences LOF vs. all ligands 
		for cur_z_LOF in lig_z_array:
			for cur_z_all in alllig_z_array:
				dif_z = cur_z_all - cur_z_LOF
				dif_z_array.append(dif_z)
		# Calculating vector lengths and appending to array
		for difcounter in range(0,len(dif_x_array)):
			vector_x_sqaured = dif_x_array[difcounter] * dif_x_array[difcounter]
			vector_y_squared = dif_y_array[difcounter] * dif_y_array[difcounter]
			vector_z_squared = dif_z_array[difcounter] * dif_z_array[difcounter]
			vector_length = np.sqrt((vector_x_sqaured + vector_y_squared + vector_z_squared))
			vector_length_array.append(vector_length)


		#### Looping through length array and checking distances  (here you also need to decode the position and refer to the atoms)
		
			# for decoding: get an array in the form of ['O1', 'O2', 'N1', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8'] at the example of ACE
		heavy_atom_codes = []
		heavyatompattern = r'' + re.escape(ligand_atomtype) + r'\s+\d+\s+([' + re.escape(heavy_atom_types) + r']\d)\s+' + re.escape(ligand_type) + r'\s+' + re.escape(str(startlig)) + r'\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}'
		heavyatompattern_compile = re.compile(heavyatompattern)
		heavyatom_matches = heavyatompattern_compile.finditer(curfile)
		
		for heavyatomcounter in heavyatom_matches:	
			heavy_atom_codes.append(heavyatomcounter.group(1))
				
		for v_length_counter in range(0,len(vector_length_array)):
			# All hits below 5A will be appended to logfile
			if vector_length_array[v_length_counter] <= 5:			 # hit if vector shorter or = 5A
				
				position = v_length_counter
				LOF_HA_index = int((position / ((nlig-1)*(len(heavy_atom_codes))))) + 1  #changed from 209 to   (nlig-1) * n_ha =  (nlig-1) * len(heavy_atoms_codes) = 19 * 11
				hit_LOF_HA = heavy_atom_codes[(LOF_HA_index-1)]  # Outputs e.g. 'C3'
				hit_current_lig = ligand_type + str(lig_of_int) # Outputs e.g. 'ACE907'
							
				helper1 = (position / ((nlig-1)*(len(heavy_atom_codes)))) - (int((position / ((nlig-1)*(len(heavy_atom_codes))))))
				helper2 = helper1 * ((nlig-1)*(len(heavy_atom_codes)))							#Example1: 166 missing positions for this LOF HA (for hit at position 1002 in big array)
				helper3 = helper2 / (len(heavy_atom_codes))										#Example1: 166/11 = 15.09
				acenr_hit = int(helper3) + 1 + 900  											#Example1: outputs e.g. 916			
							
				helper4 = helper3 - int(helper3) 												#Example1: 0.09			
				helper5 = int(helper4*(len(heavy_atom_codes)))	
				aceatom_hit = heavy_atom_codes[(helper5)]		# Outputs e.g. O1  # helper5-1 changed 15:22
				
							# Creating a log: Frame 1 | C1 ACE907 | O1 ACE902 | 2.252 Angstroms
				logstring = 'Frame ' + str(i) + ' | ' + hit_LOF_HA + ' ' + hit_current_lig + ' | '  + aceatom_hit + ' ' + ligand_type + str(acenr_hit) + ' | ' + str(vector_length_array[v_length_counter]) + ' Angstroms'
				logfile.write(logstring)
				logfile.write('\n')									
			else:
				pass
# Status: File with all heavy atom hits was written

# Closing Logfile
logfile.close()

print('\n')
print('----------------------------------------------------')
print('----------------------------------------------------')
print('Processing Hits')
print('----------------------------------------------------')	
print('----------------------------------------------------')	
print('\n')

# Opening input logfile
inputfile_name = jobname + '.log'
inputfile_open = open(inputfile_name)
inputfile_content = inputfile_open.read()

cleanme2_array = []
cleanme3_array = []

# Opening a clean outout file
outputfile_name = jobname + '.clog'
outputfile = open(outputfile_name,'w')

# Opening a results file
resultfile_basename = jobname + '.res'
resultfile = open(resultfile_basename,'w')

for frame in range(1,(nfiles+1)):
	for otherlig in range(startlig, (startlig+19)):
		
		mainpattern = r'Frame\s+' + re.escape(str(frame)) + r'\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'' + re.escape(str(otherlig)) + r'\s+\|\s+(\d\.\d+)'
		mainpattern_compile = re.compile(mainpattern)
		main_matches = mainpattern_compile.finditer(inputfile_content)
		testarray = []  # Array to test if some matches are there
		
		for maincounter in main_matches:
			testarray.append(maincounter.group(0))
			
		if len(testarray) == 0:   # if no match -> pass (there is no data matching the criteria)
			pass
		
		else:					  # if match -> append array with all lines for 1 hit to testarray  [ [...],[....],[...] ]
				cleanme2_array.append(testarray)

			
newstrings_all=''
for i_1 in range(0,len(cleanme2_array)):
	
	current_array = cleanme2_array[i_1]   # current array is i_1'th index of the big conactenated array
	
	if len(current_array) == 1:
		cleanme3_array.append(current_array[0]) # appending array with only 1 entry to clean arrray (level 3)
	
	else:
		current_string=''
		for i_2 in range(0,len(current_array)):
			 current_string = current_string + '\n' + current_array[i_2]
		second_pattern = r'Frame\s+\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+(\d\.\d+)'
		second_pattern_compile = re.compile(second_pattern)
		second_matches = second_pattern_compile.finditer(current_string)
		dist_array = []
		for i_3 in second_matches:
			dist_array.append(float(i_3.group(1)))	# Appending distances of a single array entry to a separate array
		
		# Searching the line that shows minimal distance in the distances_array by using the particular distance in the pattern
		third_pattern = r'Frame\s+\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+' + re.escape(str(min(dist_array)))
		third_pattern_compile = re.compile(third_pattern)
		third_pattern_matches = third_pattern_compile.finditer(current_string)			# matches  1 correct line
		for i_4 in third_pattern_matches:
			localhit = 	i_4.group(0)
		cleanme3_array.append(localhit) # Writing array with clean lines to clean output file, in which the hits are finally counted


for i_5 in range(0,(len(cleanme3_array))):
	outputfile.write(cleanme3_array[i_5])		
	outputfile.write('\n')
	
# Clean logfile has now been generated and will be closed for further usage
outputfile.close()

####### Processing Clean logfile

# Getting length of clean logfile (for sanity check below)
clean_log_file = open(outputfile_name,'r')
clean_log_file_array = clean_log_file.read().split('\n')	# Length is  measured by array length below

# Reading logfile for usage
clean_log_file.seek(0,0) # Getting reading-pointer back to start, # can maybe be deleted ( if stuff is deleted above)
clean_log_file_content = clean_log_file.read()

final_distance_array = []		# Array with distances of "pure" hits

fourth_pattern=r'Frame\s+\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+(\d\.\d+)'
fourth_pattern_compile = re.compile(fourth_pattern)
fourth_pattern_matches = fourth_pattern_compile.finditer(clean_log_file_content)
for i_6 in fourth_pattern_matches:
	final_distance_array.append(float(i_6.group(1)))

###### Status: Now Array with distances of only hits is ready to be counted through

# Global variables to increment
hits_between_0A2A = 0
hits_between_2A3A = 0
hits_between_3A4A = 0
hits_between_4A5A = 0


# Check for ranges e.g. between 2 and 3 A

for i_7 in range(0,(len(final_distance_array))):
	
	if final_distance_array[i_7] <= 2:	# below 2A radius check
		hits_between_0A2A = hits_between_0A2A + 1
		
	elif final_distance_array[i_7] > 2 and final_distance_array[i_7] <= 3: 	# between 2 and 3A radius check
		hits_between_2A3A = hits_between_2A3A + 1
		
	elif final_distance_array[i_7] > 3 and final_distance_array[i_7] <= 4:   	# between 3 and 4A radius check
		hits_between_3A4A = hits_between_3A4A + 1
		
	elif final_distance_array[i_7] > 4 and final_distance_array[i_7] <= 5:   	# between 4 and 5A radius check
		hits_between_4A5A = hits_between_4A5A  + 1
		
	else:
		print('ERROR: No condition fullfilled at distance comparison')
			
# Calculate cumulative results based on range results
total_below_5A = hits_between_0A2A + hits_between_2A3A + hits_between_3A4A + hits_between_4A5A
total_below_4A = hits_between_0A2A + hits_between_2A3A + hits_between_3A4A
total_below_3A = hits_between_0A2A + hits_between_2A3A
total_below_2A = hits_between_0A2A


# Sanity check: checking if all hits were processed and no false hits were in array
if (len(clean_log_file_array)-1) != total_below_5A:
	print('ALERT: Are there hits above 5A? Are there different lines in the file?')
else:
	pass


# Closing clean logfile to reopen
clean_log_file.close()


# Calculating percentage of frames with bumps below 5A 
myfile_open = open(outputfile_name,'r')
myfile_open.seek(0,0)
myfile = myfile_open.read()

frame_single_array = []
for frame_2 in range(0,(nfiles+1)):
	percentage_pattern = r'Frame\s+' + re.escape(str(frame_2)) + r'\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\d\.\d+'
	percentage_pattern_compile = re.compile(percentage_pattern)
	percentage_matches = percentage_pattern_compile.finditer(myfile)
	
	dummy_array = []
	for percentage_counter in percentage_matches:
		dummy_array.append((percentage_counter.group(0)))
	
	if len(dummy_array) == 0:
		pass
	else:
		frame_single_array.append(int(frame_2))
		
total_frame_hits = len(frame_single_array)
percentage_5A = (total_frame_hits / nfiles) * 100

myfile_open.close()

### Creating a time-evolved datasheet

# Reading Clean log file (again)
timeres_file_open = open(outputfile_name,'r')
timeres_file_open.seek(0,0)
timeres_file = timeres_file_open.read()

# Opening timeresolution output file
timeres_outfile_name = jobname + '_time-evolved.res'
timeres_outfile = open(timeres_outfile_name, 'w')


data_frame_array = []
for frame_3 in range(0,(nfiles+1)):
	timeres_pattern = r'Frame\s+' + re.escape(str(frame_3)) + r'\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\w\d\s+' + re.escape(ligand_type) + r'\d+\s+\|\s+\d\.\d+'
	timeres_pattern_compile = re.compile(timeres_pattern)
	timeres_matches = timeres_pattern_compile.finditer(timeres_file)
	for timeres_counter in timeres_matches:
		data_frame_array.append(int(frame_3))

for i_10 in range(0,(nfiles+1)):
	matchcount = data_frame_array.count(i_10)
	
	if matchcount !=0:
		timeres_outfile.write('Frame  ' + str(i_10) + '  ' + str(matchcount) + '   Bumps')
		timeres_outfile.write('\n')
	elif matchcount == 0:
		timeres_outfile.write('Frame  ' + str(i_10) + '  0   Bumps')
		timeres_outfile.write('\n')
	else:
		print('something went wrong counting your matches')
		
		

################################################  OUTPUT GENERATION  ######################################################################


# OUTPUT: output Console
print('----------------------------------------------------')
print('----------------------------------------------------')
print('Calculation finished')
print('----------------------------------------------------')	
print('----------------------------------------------------')
print('\n')	
print('---- Run info -----')	
print('Number of files processed: ' + str(nfiles))
print('Type of ligand in the files: ' + str(ligand_type))
print('Total number of ligands in files: ' + str(nlig))
print('Bumps checked for ligand numbers: ')
print(ligands_of_interest)
print('\n')
print('-----Results (ranges)-----')
print('Hits between 0A and 2A: ' + str(hits_between_0A2A))
print('Hits between 2A and 3A: ' + str(hits_between_2A3A))
print('Hits between 3A and 4A: ' + str(hits_between_3A4A))
print('Hits between 4A and 5A: ' + str(hits_between_4A5A))
print('\n')
print('-----Reults (cumulative)-----')
print('Hits below 2A: ' + str(total_below_2A))
print('Hits below 3A: ' + str(total_below_3A))
print('Hits below 4A: ' + str(total_below_4A))
print('Hits below 5A: ' + str(total_below_5A))
print('In ' + str(percentage_5A) + '% of the frames, a bump below 5A was registered')
print('\n')
print('-----File Info-----')
print('The following files were generated: ')
print('      - ' + resultfile_basename + ' | contains result summary')
print('      - ' + inputfile_name + ' | contains all hits for all atoms below 5A')
print('      - ' + outputfile_name + ' | contains filtered hits (only 1 per molecule that had bumps)')
print('      - ' + timeres_outfile_name + ' | cotains time-evolved data on bumps (frames = time)')


# OUTPUT: Result output File (.res)
resultfile.write('-------RESULTS FOR JOB: ' + jobname + ' --------')
resultfile.write(2*'\n')
resultfile.write('---- Run info -----')
resultfile.write('\n')
resultfile.write('Number of files processed: ' + str(nfiles))
resultfile.write('\n')
resultfile.write('Type of ligand in the files: ' + str(ligand_type))
resultfile.write('\n')
resultfile.write('Total number of ligands in files: ' + str(nlig))
resultfile.write('\n')
resultfile.write('Bumps checked for ' + str(len(ligands_of_interest)) + ' ligand(s)')
resultfile.write('\n')
resultfile.write('Total bumps below 5A: ' + str(total_below_5A))
resultfile.write('\n')
resultfile.write('In ' + str(percentage_5A) + '% of the frames (total: ' + str(total_frame_hits) + '), a bump below 5A (with one or multiple ligands) was registered')
resultfile.write(2*'\n')
resultfile.write('-----Results (ranges)-----')
resultfile.write('\n')
resultfile.write('Hits between 0A and 2A: ' + str(hits_between_0A2A))
resultfile.write('\n')
resultfile.write('Hits between 2A and 3A: ' + str(hits_between_2A3A))
resultfile.write('\n')
resultfile.write('Hits between 3A and 4A: ' + str(hits_between_3A4A))
resultfile.write('\n')
resultfile.write('Hits between 4A and 5A: ' + str(hits_between_4A5A))
resultfile.write(2*'\n')
resultfile.write('-----Results (cumulative)-----')
resultfile.write('\n')
resultfile.write('Hits below 2A: ' + str(total_below_2A))
resultfile.write('\n')
resultfile.write('Hits below 3A: ' + str(total_below_3A))
resultfile.write('\n')
resultfile.write('Hits below 4A: ' + str(total_below_4A))
resultfile.write('\n')
resultfile.write('Hits below 5A: ' + str(total_below_5A))
resultfile.write(3*'\n')
resultfile.write('The following files were generated: ')
resultfile.write('\n')
resultfile.write('      - ' + resultfile_basename + ' | current file')
resultfile.write('\n')
resultfile.write('      - ' + inputfile_name + ' | contains all hits for all atoms below 5A')
resultfile.write('\n')
resultfile.write('      - ' + outputfile_name + ' | contains filtered hits (only 1 per molecule with bumps)')
resultfile.write('\n')
resultfile.write('      - ' + timeres_outfile_name + ' | cotains time-evolved data on bumps (frames = time)')

print('-----------------------------------------------------------------------------')
			

			
