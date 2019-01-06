# ProPU

Inka Leroy  
January 2019
## Presentation

ProPU is a software that analyses a protein structure and suggests different subunits within the protein. Those subunits, named protein units (PU), are independent from the rest of the protein and compact.

## Perequisites :  

* [Python3](https://www.python.org/downloads/) with packages :
  * sys
  * re
  * math
  * numpy : `conda install numpy` or `sudo apt-get install python3-numpy`
  * matplotlib : 
  * statistics
  * scipy
  * progress
  
To install python3 packages it is recommended to install pip with the command :  
  `sudo apt install python3-pip`  
and then to use pip to install packages with the command :  
  `pip install package`

Once ProPU is downloaded, follow this command to make the program executable :  
  `chmod +x bin/dssp-2.0.4-linux-amd64`  
  `chmod +x ProPU`  
  
Now you can run ProPU :  
  `./ProPU -i [PDB_FILES_PATH]`

Options are available :

-h, --help                   Displays help
-i, --input                  Directory where pdb files are or path to a pdb file
--min                        Minimum size for a PU (10 by default)
--max                        Maximum size for a PU (40 by default)
--delta                      Parameter of the logistic probability function (1.5 by default)
--dist                       Distance cut-off for interactions (8.0 by default)

An example is provided in the directory example/  
Run this example with the command line : `./ProPU -i example/`

## Description of the program

This program analyses a pdb file given in entry and extract alpha carbon atoms on the precised chain. It creates a contacts matrix based on distances between atoms with, for each pair of atoms, a probability of contact based on a ligistic probability equation. Then, it cuts the protein into protein units (PUs) and calculates three criteria:
  * partition index (PI)
  * separation criterion (sigma)
  * compactness criterion (k)
  
Based on those values, it defines best PUs along the protein. 
The best PU has :
 * a PI value near 1
 * a sigma value near 0
 * the highest k value as possible 

PUs are cut according to minimum and maximum sizes defined in options. Then, based on a normalization (z-scores) of criteria values and p-value calculation, criteria are defined as significant or not with a threshold of 0.05. It is whether z-scores are negative of positive that defined if a PU can be considered as good. PUs with significantly too high sigma or too low PI and k are left.
PI values prevail on the selection of significant PUs as it assesses spliting quality quantifying the PUs independence based on contacts probabilities. Sigma and k values provide complementary information but a PU without a significant PI will not be kept, even if sigma or k were significant. 
To choose best PUs among significant ones, ProPU searches for PUs with all three significant criteria first, then with two, and then with only PI values. Best PUs do not overlap on each other. However, ProPU suggests other PUs that could be interesting to study as it is hard to well defined protein units. 

At the end, two .txt files are created in the directory resultPU/ :
* chain_name2.txt : which contains all significant PUs
* chain_name.txt : which contains best significant PUs that do not overlap on each other
  
Also, as many .png files are created as there are best PUs found by ProPU. Those files show boundaries of PUs on the contacts matrix of the protein. 

This program uses formulas and derived formulas from: 
Gelly, J. C., C. Etchebest, S. Hazout, and A. G. de Brevern. 2006. “Protein Peeling 2: A Web Server to Convert Protein Structures into Series of Protein Units.” Nucleic Acids Research 34 (WEB. SERV. ISS.).  
and  
Najmanovich, R.J., Torrance, J.W., and Thornton, J.M. (2005). Prediction of protein function from structure: Insights from methods for the detection of local structural similarities. BioTechniques 38, 847–851.


## Thanks for using ProPU !
