# coding=utf-8
import os
import sys
import signal
from lxml import etree
import numpy as np

#filter parametres
min_surface=150.00; #Å^2 (Square Ångström)
sym_op="x,y,z";
file_intro="Results for each class.\nOnly protein-ligand interfaces fulfilling all the following criteria are eligible: area greater than "+str(min_surface)+"Å^2,\nsymmetry operator \"x, y, z\" for both molecules.\n\nEach table row corresponds to an aminoacid. Legend:\n\tR : residue code.\n\t# : number of occurrences of this residue in the interfaces of proteins in this class.\n\tM : average number of residue  occurrences per protein of this class.\n\t% : percentage of residues in interfaces of this class corresponding to this aminoacid."

#tag parametres
tag_root="pisa_interfaces"
tag_subroot="pdb_entry"
tag_pdb="pdb_code"
tag_status="status"
tag_interface="interface"
tag_area="int_area"
tag_molecule="molecule"
tag_symmetry="symop"
tag_class="class"
tag_chain="residues"
tag_energy="solv_en"
tag_code="name"
tag_residue="residue"
tag_serial="ser_no"
tag_sequential = "seq_num"
class_protein="Protein"
class_ligand="Ligand"

#path parametres
path_input = "Data/pdb_list.txt"
path_enzymes = "Data/lists/all_enzymes.txt"
path_downloaded = "Data/downloaded/"
path_results = "Results/"
classes = ["enzymes","non_enzymes"]

#returns the first child with tag "tag" of node "node"
def getChild(node, tag):
	for i in range(len(node)):
		if(node[i].tag == tag):
			return node[i]
	return None

#returns a list which contains all the children with tag "tag" of node "node"
def getChildren(node, tag):
	res=[]
	for i in range(len(node)):
		if(node[i].tag == tag):
			res.append(node[i])
	return res

#gets a "Molecule" node's class
def getMoleculeClass(molecule):
	class_node = getChild(molecule,tag_class)
	if(class_node == None):
		#tag not found, xml acquisition error
		print("Tag class not found for molecule.")
		return "" 
	else:
		return str(class_node.text)

#checks the symmetry operator applied to this interface (to avoid interfaces which are a side result of the crystallization process)
def checkSymOp(molecule):
	s = str(getChild(molecule,tag_symmetry).text).lower()
	return (s == sym_op)

#returns true if an interface met the parametres and was analyzed, false otherwise
def analyzeEnzyme(enzyme, activesmap, allostermap, stats_dict):
	interfaces = getChildren(enzyme,tag_interface)
	valid_interfaces_list = list() 
	valid_interfaces = 0
	largest_area = 0
	active_site_index = -1
	#iterate over all the interfaces
	for k in range(len(interfaces)):
		area=0
		molecules=[]
		current = interfaces[k]
		#search for the interface's surface value and for the two molecules involved
		area=float(getChild(current,tag_area).text)
		molecules=getChildren(current,tag_molecule)
		if(len(molecules)!=2):
			print("Error: number of interface molecules not equal to 2.")
			continue
		#protein-ligand interface filter
		if( not (((getMoleculeClass(molecules[0])==class_ligand) and (getMoleculeClass(molecules[1])==class_protein)) or ((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_ligand)))):
			stats_dict["protein-protein"] += 1
			continue
		#symmetry operator filter (both molecules must have symmetry operator equal to "sym_op")
		if( not ((checkSymOp(molecules[0])) and (checkSymOp(molecules[1])))):
			stats_dict["symmetry"] += 1
			continue
		#minimum surface filter
		if(area<min_surface):
			stats_dict["too_small"] += 1
			continue
		#all conditions satisfied
		valid_interfaces += 1
		valid_interfaces_list.append(current)
		#check if this interface has a larger surface area than the current best
		if area>largest_area:
			largest_area = area
			active_site_index = valid_interfaces - 1 
	#if no valid interfaces were found, return False
	if valid_interfaces == 0:
		return False
	#analyze the valid interfaces
	for k in range(len(valid_interfaces_list)):
		current = valid_interfaces_list[k]
		#find the molecule entry which corresponds to the protein
		protein = None
		molecules = getChildren(current, tag_molecule)	
		if getMoleculeClass(molecules[0])==class_protein:
			protein = molecules[0]
		else:
			protein = molecules[1]
		#extract the protein chain
		chain = None
		pcn = list(protein)
		for l in range(len(pcn)):
			if(pcn[l].tag==tag_chain):
				chain=pcn[l]
				break
		#if this is the active site, update the active site map 
		if k == active_site_index:
			#update stats dict
			stats_dict["active_sites"] += 1
			#iterate over the chain sequence and count the residues (grouped by name)
			for r in range(len(chain)):
				#skip text nodes between residues
				if(chain[r].tag!=tag_residue): continue
				#analyze residue
				residue=chain[r]
				if(float(getChild(residue,tag_energy).text)!=0):
					activesmap[getChild(residue,tag_code).text.lower()] = activesmap.get(getChild(residue,tag_code).text.lower(), 0) + 1
		#otherwise update the allosteric site map
		else:
			#update stats dict
			stats_dict["allosteric_sites"] += 1
			#iterate over the chain sequence and count the residues (grouped by name)
			for r in range(len(chain)):
				#skip text nodes between residues
				if(chain[r].tag!=tag_residue): continue
				#analyze residue
				residue=chain[r]
				if(float(getChild(residue,tag_energy).text)!=0):
					allostermap[getChild(residue,tag_code).text.lower()] = allostermap.get(getChild(residue,tag_code).text.lower(), 0) + 1	
	#valid interfaces were found and analyzed, then return True	
	return True

#returns true if an interface met the parametres and was analyzed, false otherwise
def analyzeNonEnzyme(non_enzyme, classmap, stats_dict):
	interfaces = getChildren(non_enzyme,tag_interface)
	valid_interfaces = False
	#iterate over all the interfaces
	for k in range(len(interfaces)):
		area=0
		molecules=[]
		current = interfaces[k]
		#search for the interface's surface value and for the two molecules involved
		area=float(getChild(current,tag_area).text)
		molecules=getChildren(current,tag_molecule)
		if(len(molecules)!=2):
			print("Error: number of interface molecules not equal to 2.")
			continue
		#protein-ligand interface filter
		if( not (((getMoleculeClass(molecules[0])==class_ligand) and (getMoleculeClass(molecules[1])==class_protein)) or ((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_ligand)))):
			stats_dict["protein-protein"] += 1
			continue
		#symmetry operator filter (both molecules must have symmetry operator equal to "sym_op")
		if( not ((checkSymOp(molecules[0])) and (checkSymOp(molecules[1])))):
			stats_dict["symmetry"] += 1
			continue
		#minimum surface filter
		if(area<min_surface):
			stats_dict["too_small"] += 1
			continue
		#all conditions satisfied
		valid_interfaces=True
		stats_dict["valid"] += 1
		protein = None
		if getMoleculeClass(molecules[0])==class_protein:
			protein = molecules[0]
		else:
			protein = molecules[1]
		#carry on with the amino-acid occurrence count
		chain = None
		pcn = list(protein)
		for l in range(len(pcn)):
			if(pcn[l].tag==tag_chain):
				chain=pcn[l]
				break
		#iterate over the chain sequence and counts the residues (grouped by name)
		for r in range(len(chain)):
			#skip text nodes between residues
			if(chain[r].tag!=tag_residue): continue
			#analyze residue
			residue=chain[r]
			if(float(getChild(residue,tag_energy).text)!=0):
				classmap[getChild(residue,tag_code).text.lower()] = classmap.get(getChild(residue,tag_code).text.lower(), 0) + 1
	return valid_interfaces

#function to standardize string format over all the output files
def formatNumber(p):
	return ('%.2f' % p)

#script starts here
print("Initializing script")
if not os.path.exists(os.path.dirname(path_results)):
	os.makedirs(os.path.dirname(path_results))
#write output file header
path_output=path_results+"all_"+str(min_surface)+"A.txt"
outfile = open(path_output,'w')
outfile.write(file_intro)
outfile.close()
#acquire global list of codes
print("Acquiring codes from non redundant list")
non_redundant_codes = np.loadtxt(path_input, delimiter = ",", dtype = np.unicode_)
#trim white spaces
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
#build lists of enzymes and non_enzymes
codes_dict = {"enzymes":[] , "non_enzymes":[]}
print("Reading codes for enzymes")
list_of_enzymes = np.loadtxt(path_enzymes, delimiter = ",", dtype = np.unicode_)
for i in range(len(list_of_enzymes)):
	list_of_enzymes[i] = list_of_enzymes[i].lower().strip()
print("Creating non-redundant parsing lists")
parse_dict = {"enzymes":[] , "non_enzymes":[]}
for i in range(len(non_redundant_codes)):
	if non_redundant_codes[i] in list_of_enzymes:
		parse_dict["enzymes"].append(non_redundant_codes[i])
	else:
		parse_dict["non_enzymes"].append(non_redundant_codes[i])

#parse enzymes
print("Initializing count for class \"Enzymes\"")
codes_list = parse_dict["enzymes"]
analyzed=0
not_found=0
no_interfaces=0
valid_proteins=0
stats_dictionary = {"active_sites":0, "allosteric_sites":0, "symmetry":0, "protein-protein":0, "too_small":0}
activesmap = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}
allostermap = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}
percent_complete=1
print("Parsing enzymes... "+str(0)+"%")
#for each pdb code in class
for j in range(len(codes_list)):
	if((float(j)/len(codes_list))>(float(percent_complete)/100)):
		print("Parsing enzymes... "+str(percent_complete)+"%")
		percent_complete+=1
	#reading and parsing interface data from the interfaces file
	try:
		enzymeDOM = etree.parse(path_downloaded+codes_list[j]+".xml")
	except KeyboardInterrupt:
		print("Execution terminated by keyboard")
		sys.exit()
	except:
		print("File "+codes_list[j]+".xml was acquired incorrectly")
		os.remove(path_downloaded+codes_list[j]+".xml")
		print("Removed file: "+codes_list[j]+".xml")
		print("Please run the downloader again before proceeding")
		sys.exit()
	root = enzymeDOM.getroot()
	pdb_entry = getChild(root, tag_subroot)
	#checking PDB entry status
	if getChild(pdb_entry,tag_status).text == "Ok":
		res=analyzeEnzyme(pdb_entry, activesmap, allostermap, stats_dictionary)
		if not res: no_interfaces+=1
		else: valid_proteins+=1
	else:
		not_found+=1
	analyzed+=1
print("Parsing enzymes... "+str(100)+"%")

#process map of active sites
print("Printing results for enzymatic active sites")
output_file = open(path_output,'a')
#calculate additional stats
percentage = 100*(float(valid_proteins)/analyzed)
keylist = sorted(activesmap.keys())
total_residues = 0
for k in range(len(keylist)):
	total_residues += activesmap[keylist[k]]
#write to output file
output_file.write("\n\n\n*************************************************************************************\n\n\n")
output_file.write("\nEnzimi\n\n")
output_file.write("Proteine analizzate: "+str(analyzed)+"\n")
output_file.write("Non trovate in PDBePISA: "+str(not_found)+"\n")
output_file.write("Senza interfacce eleggibili: "+str(no_interfaces)+"\n")
output_file.write("Con interfacce valide: "+str(valid_proteins)+" ( "+formatNumber(percentage)+"% )\n\n")
output_file.write("Interfacce proteina-proteina: "+str(stats_dictionary["protein-protein"])+"\n")
output_file.write("Interfacce non x,y,z: "+str(stats_dictionary["symmetry"])+"\n")
output_file.write("Interfacce troppo piccole: "+str(stats_dictionary["too_small"])+"\n\n")
output_file.write("Residui totali nelle interfacce analizzate: "+str(total_residues)+"\n\n")
output_file.write("Siti Attivi: "+str(stats_dictionary["active_sites"])+"\n")
output_file.write("Mappa dei residui: \n\n")
output_file.write("\tR\t#\tM\t%\n\n\n")
for k in range(len(keylist)):
	c = activesmap[keylist[k]];
	average = (float(c))/float(stats_dictionary["active_sites"])
	local_percent = 100*(float(c)/total_residues)
	output_file.write("\t"+keylist[k]+"\t"+str(c)+"\t"+formatNumber(average)+"\t"+formatNumber(local_percent)+"\t\n\n")
output_file.close()

#process map of allosteric sites
print("Printing results for enzymatic allosteric sites")
output_file = open(path_output,'a')
#calculate additional stats
percentage = 100*(float(valid_proteins)/analyzed)
keylist = sorted(allostermap.keys())
total_residues = 0
for k in range(len(keylist)):
	total_residues += allostermap[keylist[k]]
#write to output file
output_file.write("Siti Allosterici: "+str(stats_dictionary["allosteric_sites"])+"\n")
output_file.write("Mappa dei residui: \n\n")
output_file.write("\tR\t#\tM\t%\n\n\n")
for k in range(len(keylist)):
	c = allostermap[keylist[k]];
	average = (float(c))/float(stats_dictionary["allosteric_sites"])
	local_percent = 100*(float(c)/total_residues)
	output_file.write("\t"+keylist[k]+"\t"+str(c)+"\t"+formatNumber(average)+"\t"+formatNumber(local_percent)+"\t\n\n")
output_file.close()

#parse non-enzymes
print("Initializing count for non-enzymes")
codes_list = parse_dict["non_enzymes"]
analyzed=0
not_found=0
no_interfaces=0
valid_proteins=0
stats_dictionary = {"valid":0, "symmetry":0, "protein-protein":0, "too_small":0}
classmap = {"ala":0, "arg":0, "asn":0, "asp":0, "cys":0, "gln":0, "glu":0, "gly":0, "his":0, "ile":0, "leu":0, "lys":0, "met":0, "phe":0, "pro":0, "ser":0, "thr":0, "trp":0, "tyr":0, "val":0}
percent_complete=1
print("Parsing non-enzymes... "+str(0)+"%")
#for each pdb code in class
for j in range(len(codes_list)):
	if((float(j)/len(codes_list))>(float(percent_complete)/100)):
		print("Parsing non_enzymes... "+str(percent_complete)+"%")
		percent_complete+=1
	#reading and parsing interface data from the interfaces file
	try:
		enzymeDOM = etree.parse(path_downloaded+codes_list[j]+".xml")
	except KeyboardInterrupt:
		print("Execution terminated by keyboard")
		sys.exit()
	except:
		print("File "+codes_list[j]+".xml was acquired incorrectly")
		os.remove(path_downloaded+codes_list[j]+".xml")
		print("Removed file: "+codes_list[j]+".xml")
		print("Please run the downloader again before proceeding")
		sys.exit()
	root = enzymeDOM.getroot()
	pdb_entry = getChild(root, tag_subroot)
	#checking PDB entry status
	if getChild(pdb_entry,tag_status).text == "Ok":
		res=analyzeNonEnzyme(pdb_entry, classmap, stats_dictionary)
		if not res: no_interfaces+=1
		else: valid_proteins+=1
	else:
		not_found+=1
	analyzed+=1
print("Parsing non-enzymes... "+str(100)+"%")

#process map of non-enzymes
print("Printing results for class \"Non-enzymes\"")
output_file = open(path_output,'a')
#calculating additional stats
percentage = 100*(float(valid_proteins)/analyzed)
keylist = sorted(classmap.keys())
total_residues = 0
for k in range(len(keylist)):
	total_residues += classmap[keylist[k]]
#write to output file	
output_file.write("\n\n\n*************************************************************************************\n\n\n")
output_file.write("\nNon-Enzimi\n\n")
output_file.write("Proteine analizzate: "+str(analyzed)+"\n")
output_file.write("Non trovate in PDBePISA: "+str(not_found)+"\n")
output_file.write("Senza interfacce eleggibili: "+str(no_interfaces)+"\n")
output_file.write("Con interfacce valide: "+str(valid_proteins)+" ( "+formatNumber(percentage)+"% )\n\n")
output_file.write("Interfacce valide analizzate: "+str(stats_dictionary["valid"])+"\n")
output_file.write("Interfacce proteina-proteina: "+str(stats_dictionary["protein-protein"])+"\n")
output_file.write("Interfacce non x,y,z: "+str(stats_dictionary["symmetry"])+"\n")
output_file.write("Interfacce troppo piccole: "+str(stats_dictionary["too_small"])+"\n\n")
output_file.write("Residui totali nelle interfacce analizzate: "+str(total_residues)+"\n\n")
output_file.write("Mappa dei residui: \n\n")
output_file.write("\tR\t#\tM\t%\n\n")
for k in range(len(keylist)):
	c = classmap[keylist[k]];
	average = (float(c))/float(stats_dictionary["valid"])
	local_percent = 100*(float(c)/total_residues)
	output_file.write("\t"+keylist[k]+"\t"+str(c)+"\t"+formatNumber(average)+"\t"+formatNumber(local_percent)+"\t\n\n")
output_file.close()

