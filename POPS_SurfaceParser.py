# coding=utf-8
import os
import sys
import re
import pickle
import numpy as np

#filter parametres
threshold_percent_sasa= 0.25  #amino-acids having a lower fraction of their nominal SASA exposed to the solvent, with respect to this threshold, will be discarded
file_intro_surface="Results for each class.\nOne row for each aminoacid. Legend:\n\tR : residue code.\n\t# : occurrence count on the surface of proteins belonging to this class.\n\tM : average surface occurrences per protein (average over the class).\n\t% : percentage of surface residues in this class correspodning to this aminoacid."
file_intro_protein="Results for each class.\nOne row for each aminoacid. Legend:\n\tR : residue code.\n\t# : occurrence count in the structure of proteins belonging to this class.\n\tM : average occurrences per protein (average over the class).\n\t% : percentage of residues in this class correspodning to this aminoacid."

#path parametres
path_input = "Data/pdb_list.txt"
path_enzymes = "Data/lists/all_enzymes.txt"
path_downloaded="Data/downloaded/"
path_state="Data/counting_progress.pkl"
path_results = "Results/"
url_base = "https://files.rcsb.org/download/"
classes = ["enzymes","non_enzymes"]

#class that keeps track of the counting progress
class CountingProgress:
	
	#constructor	 
	def __init__(self):
		self.last_protein = None
		self.class_index = 0
		self.analyzed=0
		#maps that keep track of the raw number of amino-acids of each type
		self.surface_map = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}
		self.protein_map = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}

	#analyzes the given protein and update the state variables
	def analyzeProtein(self, pdb_code):
		#download the given protein's structure file
		os.system("wget "+url_base+pdb_code.upper()+".pdb -O Data/"+pdb_code+".pdb")
		#call POPS (Parameter OPtimized Surface of proteins and nucleic acids) on the given protein
		try:
			os.system("pops --pdb Data/"+pdb_code+".pdb --popsOut Data/"+pdb_code+"_area.txt --residueOut")
		except:
			os.remove("Data/"+pdb_code+".pdb")
			return False			
		#scan the POPS output file
		try:
			pops_file = open("Data/"+pdb_code+"_area.txt","r")
			pops_text = pops_file.read()
			pops_file.close()
			pops_lines = pops_text.splitlines(pops_text.count("\n"))
			del pops_text
		except:
			os.remove("Data/"+pdb_code+".pdb")
			return False
		#scan the output table
		for pl in pops_lines:
			cells = re.split("[\s]+",pl)
			#skip molecule lines and empty lines
			if len(cells) < 9:
				continue
			#skip headers
			if cells[0] == "ResidNe":
				continue
			#analyze current line
			res_name = cells[1].lower()
			res_qsasa = float(cells[8])
			#skip heteroatoms
			if res_name == "het":
				continue
			#update surface stats if this residue is located on the surface
			if res_qsasa >= threshold_percent_sasa:
				if res_name in self.surface_map.keys():
					self.surface_map[res_name] += 1
				else:
					self.surface_map[res_name] = 1
			#always update protein occurrence stats, in any case
			if res_name in self.protein_map.keys():
				self.protein_map[res_name] += 1
			else:
				self.protein_map[res_name] = 1
		#update general stats
		self.analyzed+=1
		#remove the temporary files (structure and POPS output)
		os.remove("Data/"+pdb_code+".pdb")
		os.remove("Data/"+pdb_code+"_area.txt")
		os.remove("popsb.out")
		os.remove("sigma.out")
		return True

	#saves the data relative to the surface occurrence in the current class to the output file
	def finalizeSurfaceMap(self, out_path, class_name):
		output_file = open(out_path,'a')
		#calculating additional stats
		keylist = sorted(self.surface_map.keys())
		total_residues = 0
		for k in range(len(keylist)):
			total_residues += self.surface_map[keylist[k]]
		output_file.write("\n\n\n*************************************************************************************\n\n\n")
		output_file.write("\n"+class_name+"\n\n")
		output_file.write("Proteins analyzed: "+str(self.analyzed)+"\n")
		output_file.write("Total surface residues: "+str(total_residues)+"\n\n")
		output_file.write("Map of surface residues: \n\n")
		output_file.write("\tR\t#\tM\t%\n\n")
		for k in range(len(keylist)):
			c = self.surface_map[keylist[k]]
			M = float(c)/total_residues
			local_percent = 100*(float(c)/total_residues)
			output_file.write("\t"+keylist[k]+"\t"+str(c)+"\t"+formatNumber(M)+"\t"+formatNumber(local_percent)+"\t\n\n")
		output_file.close()

	#saves the data relative to the protein occurrence in the current class to the output file
	def finalizeProteinMap(self, out_path, class_name):
		output_file = open(out_path,'a')
		#calculating additional stats
		keylist = sorted(self.protein_map.keys())
		total_residues = 0
		for k in range(len(keylist)):
			total_residues += self.protein_map[keylist[k]]
		output_file.write("\n\n\n*************************************************************************************\n\n\n")
		output_file.write("\n"+class_name+"\n\n")
		output_file.write("Proteins analyzed: "+str(self.analyzed)+"\n")
		output_file.write("Total surface residues: "+str(total_residues)+"\n\n")
		output_file.write("Map of surface residues: \n\n")
		output_file.write("\tR\t#\tM\t%\n\n")
		for k in range(len(keylist)):
			c = self.protein_map[keylist[k]]
			M = float(c)/total_residues
			local_percent = 100*(float(c)/total_residues)
			output_file.write("\t"+keylist[k]+"\t"+str(c)+"\t"+formatNumber(M)+"\t"+formatNumber(local_percent)+"\t\n\n")
		output_file.close()
		
	#goes to the next class (use "finalizeClass" before this method, in order to save your data relative to the previous class)
	def nextClass(self):
		self.last_protein = None
		self.class_index += 1
		self.analyzed = 0
		self.surface_map = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}
		self.protein_map = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}

def formatNumber(p):
	return ('%.2f' % p)

#script starts here
print("Initializing script")
if not os.path.exists(os.path.dirname(path_results)):
	os.makedirs(os.path.dirname(path_results))
path_output_surface = path_results+"POPS_surface_occurrence.txt"
path_output_protein = path_results+"POPS_protein_occurrence.txt"
print("Acquiring codes from non redundant list")
non_redundant_codes = np.loadtxt(path_input, delimiter = ",", dtype = np.unicode_)
for i in range(len(non_redundant_codes)):
	non_redundant_codes[i] = non_redundant_codes[i].lower().strip()
#delete empty lines or truncated codes
del_indices = list()
for i in range(len(non_redundant_codes)):
	if len(non_redundant_codes[i]) < 4:
		del_indices.append(i)
del_indices.sort(reverse=True)
for di in del_indices:
	print("Deleting line : |"+non_redundant_codes[di]+"|")
	del non_redundant_codes[di]
#check if a counting is in progress
CP = CountingProgress()
if os.path.exists(path_state):
	#resume counting progress
	state_file = open(path_state, "rb")
	CP = pickle.load(state_file)
	state_file.close()
else:
	#clear output files before starting a new count
	outfile = open(path_output_surface,'w')
	outfile.write(file_intro_surface)
	outfile.close()
	outfile = open(path_output_protein,'w')
	outfile.write(file_intro_protein)
	outfile.close()
#count surface residues
print("Initializing counting procedure for surface and total residues")
codes_list = non_redundant_codes
percent_complete=1
start_j = 0
#if the class counting was interrupted, the start index corresponds to the index following that of the last analyzed protein
if CP.last_protein is not None:
	for j in range(len(codes_list)):
		if codes_list[j]==CP.last_protein:
			start_j = j+1
print("Parsing protein structures... "+str(0)+"%")
for j in range(start_j, len(codes_list)):
	while((float(j)/len(codes_list))>(float(percent_complete)/100)):
		print("Parsing protein structures... "+str(percent_complete)+"%")
		percent_complete+=1
	#try and analyze the structure, skip the structures POPS does not work on
	if not CP.analyzeProtein(codes_list[j]):
		continue
	CP.last_protein = codes_list[j]
	#save counting progress
	state_file = open(path_state, "wb")
	pickle.dump(CP, state_file)
	state_file.close()
print("Parsing protein structures... "+str(100)+"%")
print("Printing results")
CP.finalizeSurfaceMap(path_output_surface, "All Proteins")
CP.finalizeProteinMap(path_output_protein, "All Proteins")

#delete the counting progress file
os.remove(path_state)
#script ends here
print("Script terminated succesfully!")

