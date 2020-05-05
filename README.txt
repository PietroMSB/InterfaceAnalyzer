ProteinDownloader.py downloads the xml interface files needed for the site analysis.

ProteinParser.py calculates the occurrence of each amino-acid in protein-ligand interfaces, keeping separate stats for: active sites of enzymes, allosteric sites of enzymes, allosteric sites of non-enzymes.

POPS_SurfaceParser.py calculates the occurrence of each amino-acid on the surface of the proteins in the dataset, and their occurrence in the whole structure of the proteins in the dataset. Instead of making use of the xml interface files stored in "Data/downloaded", for each protein, it downloads a PDB structure file on the go. The strictire file is analyzed with POPS, in order to determine which aminoacids are located on the surface.

GRSH.py (Glycine-Rich Site Hunter) explores the sites of the proteins in the dataset, in search of glycine rich sites.

Running ProteinDownloader.py is needed before running ProteinParser.py and GRSH.py, as the latter work on the interface files downloaded by the former script. POPS_SurfaceParser.py operates independently, as it works on PDB structure files. To minimize disk space occupation, PDB structure files are downloaded on the go by POPS_SurfaceParser.py. If you need to speed-up the execution and disk space is not your concern, it is relatively easy to change the code, and split POPS_SurfaceParser.py in two modules, so as one will just download the files once, and the other will be able to analyze them multiple times without having to download them again.
