     ##############################
     ### --- TFIIH complex --- ####
     ##############################


These scripts demonstrate the use of IMP and PMI in the modeling of the TFIIH complex using chemical cross-links and an electron microscopy (EM) density maps.

For the publication, the following versions of IMP and PMI were used:
- IMP: git develop branch, git hash [c93138f](https://github.com/salilab/imp/tree/c93138f863070db03646bfb110b61e2547332b40)
- PMI: git develop branch, git hash [c522e0a](https://github.com/salilab/pmi/tree/c522e0abc690afc49d33d34dded70d84ad4ec3a7)


### --- List of files and directories:

 - imp-21Oct13: contains a version of IMP that works with the provided scripts

 - yeast:
   - inputs: contains comparative models, structures, cross-linking data, and EM data
   - outputs:
     - contains an example of output files (*.rmf and *.dat) from sampling
     - models:
       - subunits: contains a chimera session file and solution density maps for all subunits
       - domain_decomposition: contains a chimera session file and solution density maps for all domains
   - scripts: 
     - sample_multires.py: sampling script
     - analyze.py: script for analysis of the solutions

 - human:
   - inputs: contains comparative models, structures, cross-linking data, and EM data
   - outputs:
     - contains an example of output files (*.rmf and *.dat) from sampling
     - models:
       - subunits: contains a chimera session file and solution density maps for all subunits
       - domain_decomposition: contains a chimera session file and solution density maps for all domains
   - scripts: 
     - sample_multires.py: sampling script
     - analyze.py: script for analysis of the solutions


### --- Running the IMP/PMI scripts for the TFIIH

 -yeast:
 	- cd yeast/scripts
	- $PATH_TO_IMP/setup_environment.sh python sample_multires.py 
	- $PATH_TO_IMP/setup_environment.sh python analyze.py
 
 -human:
	- cd human/scriptss
	- $PATH_TO_IMP/setup_environment.sh python sample_multires.py 
	- $PATH_TO_IMP/setup_environment.sh python analyze.py
        - $PATH_TO_IMP/setup_environment.sh python clustering.py

Note that in both yeast and human cases, the script analyze.py has many different uses: calculating the top 10% of scores, getting density maps, calculating the percentage
of violated constraints, getting PDB files for the models, clustering the top scoring models, getting contact maps etc.

To use each of these functions, one can just comment out the remaining part of the code and retain only the neccessary function each time. 

Note that clustering is a 2 step process for the human modeling: the cluster subroutine in analyze.py outputs the distance matrix and the clustering.py script run on that
distance matrix produces the actual clusters.


### --- Information

Authors: Peter Cimermancic, Shruthi Viswanath, Riccardo Pellarin, Charles Greenberg
       

