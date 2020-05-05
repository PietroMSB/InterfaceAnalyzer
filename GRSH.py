# coding=utf-8
import os
import sys
import signal
from lxml import etree
import numpy as np

#filter parametres
min_surface=150.00 #Å^2 (Square Ångström)
min_glycine_ratio = 0.3
sym_op="x,y,z"
file_intro="Results for each class.\nOnly protein-ligand interfaces with a surface greater than "+str(min_surface)+"Å^2,\n in which both molecules have symmetry operator \"x, y, z\" are considered.\n\nEach row in the table describes a Glycine Rich Site. Legend:\n\tPDB : PDB code of the protein.\n\tG : number of glycines in the site.\n\tT : total number of residues in the site.\n\t% : percentage of glycine residues in the site."

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

#returns true if an interface met the parametres and the protein was analyzed, false otherwise
def analyzeProtein(protein_object, classmap, stats_dict, pdb_code):
	interfaces = getChildren(protein_object,tag_interface)
	valid_interfaces = False
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
		#count glycines and other aminoacids in the active site
		count_glycines = 0
		count_total = 0
		chain = None
		pcn = list(protein)
		for l in range(len(pcn)):
			if(pcn[l].tag==tag_chain):
				chain=pcn[l]
				break
		#iterate over the chain sequence and count the residues
		for r in range(len(chain)):
			#skip text nodes between residues
			if(chain[r].tag!=tag_residue): continue
			#analyze residue
			residue=chain[r]
			#check if the residue belongs to the interface
			if(float(getChild(residue,tag_energy).text)!=0):
				if getChild(residue,tag_code).text.lower() == 'gly':
					count_glycines += 1
				count_total += 1
		#normalize the local classmap, dividing each count by the total number of residues in the map
		glycine_ratio = float(count_glycines)/float(count_total)
		#filter out sites which do not reach the minimum glycine ratio
		if glycine_ratio >= min_glycine_ratio:
			#append the entry corresponding to this interface to the classmap list
			classmap.append([pdb_code, count_glycines, count_total, glycine_ratio])
	#return True if at least one valid interface was found, False otherwise
	return valid_interfaces

#formats floating point values as strings with a given number of digits
def formatNumber(p):
	return ('%.6f' % p)

#script starts here
print("Initializing script")
if not os.path.exists(os.path.dirname(path_results)):
	os.makedirs(os.path.dirname(path_results))
#write ouutput file header
path_output=path_results+"glycine_rich_sites.txt"
outfile = open(path_output,'w')
outfile.write(file_intro)
outfile.close()
#load list of PDB codes
print("Acquiring codes from non redundant list")
non_redundant_codes = np.loadtxt(path_input, delimiter = ",", dtype = np.unicode_)
#trim white spaces
for i in range(len(non_redundant_codes)):
	non_redundant_codes[i] = non_redundant_codes[i].lower()
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
input_file = open(path_enzymes,'r')
input_list = input_file.read()
list_of_enzymes = input_list.split("\n",input_list.count("\n"))
input_file.close()
print("Creating non-redundant parsing lists")
parse_dict = {"enzymes":[] , "non_enzymes":[]}
for i in range(len(non_redundant_codes)):
	if non_redundant_codes[i] in list_of_enzymes:
		parse_dict["enzymes"].append(non_redundant_codes[i])
	else:
		parse_dict["non_enzymes"].append(non_redundant_codes[i])

#parse enzymes
print("Initializing search for class \"Enzymes\"")
codes_list = parse_dict["enzymes"]
analyzed=0
not_found=0
no_interfaces=0
valid_proteins=0
stats_dictionary = {"valid":0, "symmetry":0, "protein-protein":0, "too_small":0}
classmap = list()
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
		res=analyzeProtein(pdb_entry, classmap, stats_dictionary, codes_list[j])
		if not res: no_interfaces+=1
		else: valid_proteins+=1
	else:
		not_found+=1
	analyzed+=1
print("Parsing enzymes... "+str(100)+"%")

#process countmap of enzymes
print("Printing results for class \"Enzymes\"")
#percentage of valid proteins over the total number of proteins in the class
percentage = 100*(float(valid_proteins)/analyzed)
#write to output file
output_file = open(path_output,'a')
output_file.write("\n\n\n*************************************************************************************\n\n\n")
output_file.write("\nEnzimi\n\n")
output_file.write("Proteine analizzate: "+str(analyzed)+"\n")
output_file.write("Non trovate in PDBePISA: "+str(not_found)+"\n")
output_file.write("Senza interfacce eleggibili: "+str(no_interfaces)+"\n")
output_file.write("Con interfacce valide: "+str(valid_proteins)+" ( "+formatNumber(percentage)+"% )\n\n")
output_file.write("Interfacce proteina-proteina: "+str(stats_dictionary["protein-protein"])+"\n")
output_file.write("Interfacce non x,y,z: "+str(stats_dictionary["symmetry"])+"\n")
output_file.write("Interfacce troppo piccole: "+str(stats_dictionary["too_small"])+"\n\n")
output_file.write("Tabella dei siti ricchi di glicina: \n\n")
output_file.write("\tPDB\tG\tT\t%\n\n")
for k in range(len(classmap)):
	glycine_percent = 100*classmap[k][3]
	output_file.write("\t"+classmap[k][0]+"\t"+str(classmap[k][1])+"\t"+str(classmap[k][2])+"\t"+formatNumber(glycine_percent)+"\t\n\n")
output_file.close()

#parse non-enzymes
print("Initializing count for non-enzymes")
codes_list = parse_dict["non_enzymes"]
analyzed=0
not_found=0
no_interfaces=0
valid_proteins=0
stats_dictionary = {"valid":0, "symmetry":0, "protein-protein":0, "too_small":0}
classmap = list()
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
		res=analyzeProtein(pdb_entry, classmap, stats_dictionary, codes_list[j])
		if not res: no_interfaces+=1
		else: valid_proteins+=1
	else:
		not_found+=1
	analyzed+=1
print("Parsing non-enzymes... "+str(100)+"%")

#process countmap of non-enzymes
print("Printing results for class \"Non-enzymes\"")
#percentage of valid proteins over the total number of proteins in the class
percentage = 100*(float(valid_proteins)/analyzed)
#write to output file
output_file = open(path_output,'a')	
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
output_file.write("Tabella dei siti ricchi di glicina: \n\n")
output_file.write("\tPDB\tG\tT\t%\n\n")
for k in range(len(classmap)):
	glycine_percent = 100*classmap[k][3]
	output_file.write("\t"+classmap[k][0]+"\t"+str(classmap[k][1])+"\t"+str(classmap[k][2])+"\t"+formatNumber(glycine_percent)+"\t\n\n")
output_file.close()

