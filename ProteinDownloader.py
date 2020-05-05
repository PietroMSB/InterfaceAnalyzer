# coding=utf-8
import os
import sys
import requests
import numpy as np
from lxml import etree
import numpy as np

#path parametres
path_input = "Data/pdb_list.txt"
path_lists = "Data/lists/"
path_output = "Data/downloaded/"
classes = ["enzymes", "non_enzymes"]
base_url="http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?"

#start execution
print("Initializing download script")
print("Acquiring codes from non redundant list")
#load global list of PDB codes
non_redundant_codes = np.loadtxt(path_input, delimiter = ",", dtype = np.unicode_)
#trim white spaces
for i in range(len(non_redundant_codes)):
	non_redundant_codes[i] = non_redundant_codes[i].strip()
#start downloading
print("Starting download procedure")
for i in range(len(non_redundant_codes)):
		#build path string
		filepath = path_output+non_redundant_codes[i].lower()+".xml"
		#check if the file already exists and
		if os.path.exists(filepath):
			print("Checking file: "+non_redundant_codes[i]+".xml ( "+str(i+1)+" of "+str(len(non_redundant_codes))+" )", end='\r')
			#check file integrity and, in a positive case, continue to the next PDB code
			try:
				dom_check = etree.parse(filepath)
				continue
			#catch keyboard interrupts
			except KeyboardInterrupt:
				print("Execution terminated by keyboard")
				sys.exit()
			#catch integrity check failures and, in case, just download the file again
			except:
				print("Integrity check failed")
		#download the file (if it was missing or failed the integrity check)
		print("Retrieving interface file for "+non_redundant_codes[i]+" ( "+str(i+1)+" of "+str(len(non_redundant_codes))+" )", end='\r')
		response=requests.get(base_url+non_redundant_codes[i].lower())
		if(response.status_code==200):
			outfile = open(filepath,'wb')
			outfile.write(response.content)
			outfile.close()
print("Download script terminated succesfully!")


