These scripts demonstrate the use of [IMP](https://integrativemodeling.org)
and [PMI](https://github.com/salilab/pmi) in the modeling of the TFIIH complex
using chemical cross-links and electron microscopy (EM) density maps.

For the publication, the following versions of IMP and PMI were used:
- IMP: git develop branch, git hash [c93138f](https://github.com/salilab/imp/tree/c93138f863070db03646bfb110b61e2547332b40)
- PMI: git develop branch, git hash [c522e0a](https://github.com/salilab/pmi/tree/c522e0abc690afc49d33d34dded70d84ad4ec3a7)


## List of files and directories:

 - `yeast`:
   - `inputs`: contains comparative models, structures, cross-linking data, and EM data
   - `outputs`:
     - contains an example of output files (*.rmf and *.dat) from sampling
     - `models`:
       - `subunits`: contains a chimera session file and solution density maps for all subunits
       - `domain_decomposition`: contains a chimera session file and solution density maps for all domains
   - `scripts`: 
     - `sample_multires.py`: sampling script
     - `analyze.py`: script for analysis of the solutions

 - `human`:
   - `inputs`: contains comparative models, structures, cross-linking data, and EM data
   - `outputs`:
     - contains an example of output files (*.rmf and *.dat) from sampling
     - `models`:
       - `subunits`: contains a chimera session file and solution density maps for all subunits
       - `domain_decomposition`: contains a chimera session file and solution density maps for all domains
   - `scripts`: 
     - `sample_multires.py`: sampling script
     - `analyze.py`: script for analysis of the solutions


## Running the IMP/PMI scripts for the TFIIH

 - `yeast`:
   - `cd yeast/scripts`
   - `python sample_multires.py`
   - `python analyze.py`
 
 - `human`:
   - `cd human/scripts`
   - `python sample_multires.py`
   - `python analyze.py`
   - `python clustering.py`

Note that in both yeast and human cases, the script `analyze.py` has many different uses: calculating the top 10% of scores, getting density maps, calculating the percentage
of violated constraints, getting PDB files for the models, clustering the top scoring models, getting contact maps etc.

To use each of these functions, one can just comment out the remaining part of the code and retain only the neccessary function each time. 

Note that clustering is a 2 step process for the human modeling: the cluster subroutine in `analyze.py` outputs the distance matrix and the `clustering.py` script run on that
distance matrix produces the actual clusters.

## Information

_Author(s)_: Peter Cimermančič, Shruthi Viswanath, Riccardo Pellarin, Charles Greenberg

_Date_: October 6th, 2014

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
