import re


####################################################  README  ###########################################################################
""" 
Purpose: This script is designed to determine the position of ligand molecules relative to the membrane in an CYP2D6 system			
																																		
Usage: Run <python3 LigandPositionMembrane_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.		
																																		
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb"	
												
 		 - Frames need to be aligned, so that protein is centered and membrane is below it
		 - This script expects ligands to be sequentially numbered in the PDB File (e.g. 901, 902, 903, ...)						
																													
Output: Percentages of the ligand in the compartements: Membrane, Headgroups, Solvent/Protein

Example: Example with 5 MD frames is provided



Written in gedit by André Fischer (2018/2019)	
"""

################################################  EXAMPLE INPUT  ########################################################################
""" 
ATOM      1  O1  ACE   901     -15.429 -16.024  -4.709  1.00  0.00           O  
ATOM      2  O2  ACE   901      -9.316 -14.596  -6.952  1.00  0.00           O  
ATOM      3  N1  ACE   901     -10.008 -14.694  -4.794  1.00  0.00           N  
ATOM      4  C1  ACE   901     -11.397 -15.083  -4.914  1.00  0.00           C  
ATOM      5  C2  ACE   901     -11.901 -15.633  -3.677  1.00  0.00           C  
ATOM      6  C3  ACE   901     -12.216 -14.820  -6.019  1.00  0.00           C  														
HETATM 9996 H16R POPC   13       8.000   8.051  44.784  1.00  0.00           H  
HETATM 9997 H16S POPC   13       8.984   8.616  43.366  1.00  0.00           H  
HETATM 9998   HX POPC   13      -4.672  13.847  49.031  1.00  0.00           H  
HETATM 9999   HY POPC   13      -3.624  13.002  50.254  1.00  0.00           H  
HETATM10000  H2X POPC   13      -5.548   9.426  47.843  1.00  0.00           H  
HETATM10001  H2Y POPC   13      -3.821   9.698  47.446  1.00  0.00           H  
HETATM10002  H3X POPC   13      -4.888  11.868  46.203  1.00  0.00           H  
HETATM10003  H3Y POPC   13      -6.404  11.080  46.357  1.00  0.00           H    
"""                                    				
#########################################################################################################################################


###############################################  USER DEFINED VARIABLES  ################################################################

maxframes = 5							# Number of aligned PDB files to analyze (e.g. 500)
resname = "ACE"    						# Residue of your ligand to be monitored (e.g. ACE)
totalres = 20 							# Total numbers of residues of interest	 (e.g. 20)	
lowresrebound = 901 						# Starting number of your residues (e.g. 901)
memtype = "POPC"						# Type of membrane molecules used (e.g. POPC)
output_basename = 'ligpos_example1' 				# Your outputfile basename	


################################################  EXECUTION  #############################################################################

# Opening Output file
outputfilename = output_basename + '.out'
outputfile = open(outputfilename,"w")

# Preparing global output counter variable
HG_count=0
M_count=0
S_count=0

# Looping over the files 1.pdb, 2.pdb, ..., n.pdb
for i in range(1, (maxframes+1)):
	print("processing file " + str(i))		
		
	# Reading the file
	filestring = str(i) + '.pdb'
	curfile_open = open(filestring, "r")
	curfile = curfile_open.read()
	
	# N/C2-z-coordiate array to fill
	C2_z_POPC_coords=[]
	N_z_POPC_coords=[]
	
	# Determine total membrane residues
	allpattern = r'HETATM\s?\d+\s+' + re.escape("C2") + r'\s+' + re.escape(memtype) + r'\s+\d+\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'
	whole_allpatern = re.compile(allpattern)
	matches_all = whole_allpatern.findall(curfile) 							# Coords of all C2-atoms in array
	total_POPC = len(matches_all)									# Length of this array = number of membrane residues
	
	# Adding C2-z-coordinates to array
	for current_resnum in range(1, (total_POPC+1)):
		residuepatternC2 = r'HETATM\s?\d+\s+' + re.escape("C2") + r'\s+' + re.escape(memtype) + r'\s+' + re.escape(str(current_resnum)) + r'\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'	
		wholepatternC2 = re.compile(residuepatternC2)		
		matchesC2_2 = wholepatternC2.finditer(curfile)
		
		for lcount_C2 in matchesC2_2:
			z_coord = lcount_C2.group(1)
			C2_z_POPC_coords.append(float(z_coord))						# appends z-coordinates of POPC-N-Atoms into Array		
							
	# Sorting z-coordinates of POPC-C2-Atoms according to size
	C2_z_POPC_coords.sort()												
	
	# Determine upper POPC from C2 Array; this value will also be used for N-Atoms later
	allpopc=total_POPC
	arraycounter=0
	for arraycounter in range(0,(allpopc-1)):
		lowerbound = C2_z_POPC_coords[arraycounter]
		upperbound=C2_z_POPC_coords[arraycounter+1]
		C2dif = upperbound - lowerbound
		if C2dif > 12: 										# checking for at least 12A difference between upper and lower membrane leaflet points
			upper_POPC = arraycounter + 2
		else: 
			pass
			
	# Slicing array with C2-POPC z-coordinates into two arrays with upper and lower leaflet	
	C2_top_array = C2_z_POPC_coords[0:(upper_POPC-1)]
	C2_low_array = C2_z_POPC_coords[(upper_POPC-1):len(C2_z_POPC_coords)]
	
	# Calculating Average value of top / lower C2 coords
	sum_C2_top = sum(C2_top_array)
	len_C2_top = len(C2_top_array)
	sum_C2_lower = sum(C2_low_array)
	len_C2_lower = len(C2_low_array)
	average_C2_top = sum_C2_top / len_C2_top 							# Mass center of C2-POPC atoms of upper leaflet
	average_C2_lower = sum_C2_lower / len_C2_lower							# Mass center of C2-POPC atoms of lower leaflet
	
	# Adding N-z-coordinates to array
	for current_resnum in range(1, (total_POPC+1)):
		residuepatternN = r'HETATM\s?\d+\s+' + re.escape("N") + r'\s+' + re.escape(memtype) + r'\s+' + re.escape(str(current_resnum)) + r'\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'	
		wholepatternN = re.compile(residuepatternN)		
		matchesN = wholepatternN.finditer(curfile)
		
		for lcount_N in matchesN:
			z_coord = lcount_N.group(1)
			N_z_POPC_coords.append(float(z_coord))						# appends z-coordinates of POPC-N-Atoms into Array
	
	# Sorting z-coordinates of POPC-C2-Atoms according to size
	N_z_POPC_coords.sort()
	
	# Slicing array with N-POPC z-coordinates into two arrays with upper and lower leaflet														
	N_top_array = N_z_POPC_coords[0:(upper_POPC-1)]						
	N_low_array = N_z_POPC_coords[(upper_POPC-1):len(N_z_POPC_coords)]		
	
	# Calculating Average value of top / lower N coords
	sum_N_top = sum(N_top_array)
	len_N_top = len(N_top_array)
	sum_N_lower = sum(N_low_array)
	len_N_lower = len(N_low_array)
	average_N_top = sum_N_top / len_N_top									# Mass center of N-POPC atoms of upper leaflet
	average_N_lower = sum_N_lower / len_N_lower								# Mass center of N-POPC atoms of lower leaflet
	
	
	# Correction factors for POPC C2/N Mass centers
	average_C2_top = average_C2_top + 1.5
	average_C2_lower = average_C2_lower - 1.5
	average_N_top = average_N_top - 1.5
	average_N_lower = average_N_lower + 1.5
	
	
	# Looping over all the Ligands
	lignr_upper_end = int(lowresrebound + totalres)
	for lignr in range(lowresrebound,lignr_upper_end):
		ACE_mc_array=[]
		residuepatternACE = r'ATOM\s+\d+\s+[A-Z0-9]+\s+' + re.escape(resname) + r'\s+' + re.escape(str(lignr)) + r'\s+[0-9-]+\.\d{3}\s+[0-9-]+\.\d{3}\s+([0-9-]+\.\d{3})'
		wholepatternACE = re.compile(residuepatternACE)		
		matchesACE = wholepatternACE.finditer(curfile)
		for lcount_ACE in matchesACE:
			z_coord = lcount_ACE.group(1)
			ACE_mc_array.append(float(z_coord))							# appends z-coordinates of ligand_n atoms into Array
		
		# Calculating mass center of ligand_n
		sum_ACE=sum(ACE_mc_array)
		len_ACE=len(ACE_mc_array)
		ACE_mc = sum_ACE / len_ACE
	
		
		# Checking Ligand position Relative to the membrane
		if ACE_mc < average_C2_top and ACE_mc > average_N_top:
			HG_count=HG_count + 1									# incrementing M-counter if hit
			
		elif ACE_mc > average_C2_top and ACE_mc < average_C2_lower:
			M_count = M_count + 1									# incrementing HG-counter if hit
			
		elif ACE_mc > average_C2_lower and ACE_mc < average_N_lower:		
			HG_count=HG_count + 1									# incrementing HG-counter if hit
			
		else:
			S_count=S_count + 1									# incrementing S-counter if hit
			
			
################################################  OUTPUT PROCESSING ##########################################################################

# Calculation Percentages
total_count = HG_count + S_count + M_count
percent_HG=((HG_count / total_count)*100)
percent_S=((S_count / total_count)*100)
percent_M=((M_count / total_count)*100)

# Result Output for Terminal
print('::::RESULT::::')
print('HG: ' + str(HG_count))
print('S: ' + str(S_count))
print('M: ' + str(M_count))
print('Total: ' + str((HG_count + S_count + M_count)))

		
# Result output File
outputfile.write('RESULT for ' + output_basename + 'with Ligand ' + resname + ' and numbers from ' + str(lowresrebound) + ' to ' + str(lignr_upper_end-1))
outputfile.write(3*'\n')
outputfile.write('HG: ' + str(HG_count))
outputfile.write('\n')
outputfile.write('HG Equals ' + str(percent_HG) + '% ')
outputfile.write(2*'\n')
outputfile.write('S: ' + str(S_count))
outputfile.write('\n')
outputfile.write('S Equals ' + str(percent_S) + '% ')
outputfile.write(2*'\n')
outputfile.write('M: ' + str(M_count))
outputfile.write('\n')
outputfile.write('M Equals ' + str(percent_M) + '% ')
outputfile.write(2*'\n')
outputfile.write('Total: ' + str((HG_count + S_count + M_count)))

# Check expected length of output
explength = maxframes * totalres
if total_count != explength:
	print("Something might have went wrong in your calculation")
	outputfile.write("Something might have went wrong in your calculation")
	outputfile.write(2*'\n')

else:
	print("Your calculation went as panned")
	outputfile.write(2*'\n')	
