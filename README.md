# pKa analysis scripts

Author: Eric Lang 
<eric.jm.lang@gmail.com>
University of Canterbury

Licensed under the GNU General Public License

If you use these scripts, please cite:

  >Calculated pKa Variations Expose Dynamic Allosteric Communication Networks
  >Eric J. M. Lang, Logan C. Heyes, Geoffrey B. Jameson, and Emily J. Parker
  >Journal of the American Chemical Society 2016 138 (6), 2036-2045
  >DOI: 10.1021/jacs.5b13134 


Table of contents:

[1. Introduction](#1. Introduction)
[2. Installation](#2. Installation)
[3. How to perform the analysis](#3. How to perform the analysis)
[4. Examples of custom modification of the scripts](#4. Examples of custom modification of the scripts)

-----------------------------------------------------------------------------
1. Introduction
-----------------------------------------------------------------------------

This is a collection of scripts used in the paper _"Calculated pKa Variations 
Expose Dynamic Allosteric Communication Networks"_ published in JACS (see 
citation above).

Most of these scripts are custom made for the analysis of the data presented 
in the paper. However, with minimum modifications, they can be applied to 
other systems to identify allosteric communication pathways.

In the following I provide a description of how to typically run the 
analysis.In the header of each scripts, the part that need to or can be 
modified according to your needs are indicated.

If you have any problems or comments, you can to contact me.

-----------------------------------------------------------------------------
2. Installation
-----------------------------------------------------------------------------

This is a collection of scripts so no installation of scripts are required,
they just need to be in your working directory or in your path.
However, these scripts rely on a number of programs or libraries that need
to be installed.

- PROPKA 3.1 (https://github.com/jensengroup/propka-3.1)

- Python 2: the following python modules needs to be installed (e.g. using `pip`):
 * glob
 * numpy 

- Batcher. This software enables to parallelise multiple PROPKA runs.  
  Batcher is part of the ADLB library and can be downloaded from the 
  following site:
  https://www.cs.mtsu.edu/~rbutler/adlb/ and follow the instructions of the 
  authors to install it. Batcher is advantageously installed on a cluster so 
  it can be used  with a large number of cores, meaning that a large number of
  PROPKA runs can be performed at the same time.

- VMD (http://www.ks.uiuc.edu/Research/vmd/)


-----------------------------------------------------------------------------
3. How to perform the analysis
-----------------------------------------------------------------------------

You need to start in your working directory with your trajectory file(s) and 
create a directory

	mkdir pdbfiles

First convert the trajectory into multiple PDB files that do not include 
hydrogen atoms (one PDB file for each frame). e.g. in vmd starting with 
a dcd and psf file, modify `DCD2PDBs.tcl` to add the file name and number 
of frames and run:

	vmd -dispdev text -e DCD2PDBs.tcl 

If the PDB are not "standard" PDB files, convert them into standard ones. 
e.g. for CHARMM or NAMD generated PDB use `convert_PDB_files.py` included.

	python convert_PDB_files.py > convert_PDB_files.log

Then prepare the input file for Batcher using the 
`create_input_command_batcher.sh` batch script:

	bash create_input_command_batcher.sh

This will generate the file `batcher_input_command` which can then be passed on 
to Batcher, e.g.:

	mpiexec -n X batcher_input_command

where X is an integer corresponding to the number of CPUs.
If you are processing a large number of PDB files, Batcher will need to be 
installed on a local cluster and run via your usual job scheduler.I 
processed for example 40000 PDB files, from one trajectory, so PROPKA needs 
to run 40000 times which was made possible by using
batcher on 128 CPUs.

This will generate a number of .pka files (PROPKA output), one for each PDB 
file. These pka files
contain a lot of information that needs to be parsed. For this you need to 
use the parse_propka_output.py script which is run with:

	python parse_propka_output.py

this script need a file structure.pdb which is for example the first frame 
of the trajectory and will process all the .pka files to generate a number 
of files. In detail it does the following:

- Read structure.pdb and retrieve the identifier of each residues (residue 
  name, number and chain), and various pieces of information (number of 
  residues, chains, etc.) to generate a dictionary that will associate a 
  unique number to each unique residue and generate a list of residue: 
  'residues.dat'

- Extract the total pKa value and the total contribution of the desolvation 
  effect (in unit of pKa) and save them in two files: `pka_all.dat` and 
  `desolvation_all.dat` in which each row corresponds to a different 
  ionizable residue and each column corresponds to a trajectory frame.

- For one pka file, create 3 matrices of dimension number of residues x 
  number of residues, one for the coulombic interactions, 
  one for the side chain H-bonds, one for the backbone H-bonds.

 * Extract the pKa change induced by a residue via side chain H-
   bonding as well as the identifier of this interacting residue. 
   this value is saved in the matrix associated with side chain H-
   bonds and for which the row number corresponds to the unique value 
   of the ionisable residue generated with the dictionary and the 
   column number corresponds to the unique value of the interacting 
   residue generated with the dictionary.

 * Do the same for side-chain H-bonds and coulombic interactions.

 * The generated matrices are not symmetrical but take the form:
	0 X 0
	-X 0 Y
	Z 0 0
   In order to be able to average the value by chain during the 
   analysis, a number of operations need to be done 
   on the matrices in order to have them in the form:
	0 X -Z
	0 0 Y
	0 0 0
   This if for one pKa file. This is repeated for all pKa files and each new 
   set of three matrices is added to the first set
   of matrices generated and the resulting set of matrices are saved 
   in `sidechain_Hbond_mat.dat` `backbone_Hbond_mat.dat`, `coulombic_mat.dat`.



Once the parsing is done, the resulting files can be analysed with other 
scripts. For the total pKa value and the desolvation effect, the values are 
averaged per chain and the standard deviation calculated with the 
`analyse_pka_global.py` script:

	python analyse_pka_global.py

which will generate one file for the averaged pKa value and one for the std 
dev.

The same is done for the desolvation running:

	python analyse_desolvation.py

The files `sidechain_Hbond_mat.dat`,`backbone_Hbond_mat.dat` and 
`coulombic_mat.dat`can then be used to identify the H-bond
and coulombic interactions networks. This is done with 
`analyse_pka_interactions.py` script:

	python analyse_pka_interactions.py

The interest being to compare a system in the absence or presence of an 
allosteric effector, the three set of files is required for both the pKa 
analysis in the presence and abscence of effector. it also requires the 
`residues.dat` file generated during the parsing of the multiple .pka files, 
the `AVERAGE_pKa.dat` files calculated with `analyse_pka_global.py` both in the 
presence and abscence of effector,and an Identifier.dat file which is 
extracted from one .pka file (the summary section where the residues are 
listed i.e. the `SUMMARY OF THIS PREDICTION` section - an example is given 
'Example_Identifier.dat'). Finally a number of parameters need to be added: 
 pH, number of frames in the initial trajectories, etc.

This script mainly calculates the difference in average interaction between 
residues in the presence and abscence of effector and then generate 
automatically the tcl scripts that will draw the interaction maps in vmd.

You can then run the tcl scripts in VMD.


-----------------------------------------------------------------------------
4. Examples of custom modifications of the scripts
-----------------------------------------------------------------------------

As mentioned previously the scripts need to be modified according to your
needs. This requires some work, but in the Example directory you will find
how the scripts have been modified for a particular structure
(`Example_structure.pdb`)

An example of a single .pka file is also provided `Example.pka`

Modification of most of the scripts should not be too much a source of 
complications, except for `parse_propka_output.py` and 
`analyse_pka_interactions.py`


 4.1. Modifying `parse_propka_output.py`
-----------------------------------------------------------------------------

Please refer to `Example_parse_propka_output.py` to see the modifications of the 
original script

First you need to modify the `SearchHETATM` (line  of the original script)
string to identify `HETATM` entries in your PDB file. 
This is simply used to identify each residue. You should only 
need to modify the part `[CAMN]` in:

	SearchHETATM='^[HETATM\s\d]{11}\s+[CAMN]{2}\s+\w{2}[\w\s]\s(\w)\s+\d+'

Here the HETATM entries correspond to the Mn<sup>2+</sup> ions and the PHE ligand (which 
as a CA) so we look for CA or MN in the regular expression string. 
Here I used CA but it could have been any other unique PHE atoms such as CB, 
in which case it would have been `[CBMN]` in the `SearchHETATM` string. If we had 
Zn<sup>2+</sup> instead of Mn<sup>2+</sup> it could have been `[CBZN]`.

The second thing to modify is the function 'dictionnaries'. This function is 
used to asign a unique ID to each unique residue that can interact with a 
ionizable residue. You only need to worry about the section on the 
ligands (line 161 to 166 of the original script). 

Here one of our ligand is PHE, so it has both an amine and an acid carboxylic
group. Therefore its pKa will be calculated by PROPKA. if you look at 
`Example.pka` you can see how the amine and carboxylate functions are defined 
in PROPKA. Where you usually have the residue number, you find the function: N 
or C (for amine and carboxylate resp.). This is what you need to identify 
them, and you can modify the original script to have:

        if chain in Chains_het: 
            if resid == 'N': 
                Id_2_Numb[ID]= (len(Chains) - len(Chains_het))*Resid_count + (Chain_2_Numb[chain] - (len(Chains) - len(Chains_het)))
            if resid == 'C':
                Id_2_Numb[ID]= len(Identifier) + (Chain_2_Numb[chain] - (len(Chains) - len(Chains_het)))

if you cannot easily find a way to calculate the number associated with each 
ligand ID, using for example:

	len(Chains) - len(Chains_het))*Resid_count + (Chain_2_Numb[chain] - (len(Chains) - len(Chains_het)) 

you can always find out the number manually. 

Because we have two different charged function on the ligand it is easier to 
work as if there were two ligands each carrying only one charge. This is shown 
at the end of `Example_residues.dat` where PH1 corresponds to the amine 
function of PHE and PH2 to the carboxylate function of PHE. This part can be 
reused "as is" for any amino acid ligands

Sometimes your ligands might not be ionizable but interact with the ionizable
residues nontheless and thus need to be added to the dictionnary.
for example on _line 123_ of `Example.pka` there is the following:

	ASP 267 A   2.72*  100 %    2.90  611   0.54    0    0.00 XXX   0 X   -0.77 SER 269 A   -1.10 MN   MN I

So the Manganese ion forms coulombic interactions with residue `ASP 267`, 
it needs therefore to be added to the dictionnary `Id_2_Numb`. 
There are 4 Mn<sup>2+</sup> in the structure file, to add them to the dictionnary 

	Id_2_Numb["MN  MN I"] = Id_2_Numb["MN    1 I"]
	Id_2_Numb["MN  MN J"] = Id_2_Numb["MN    1 J"]
	Id_2_Numb["MN  MN K"] = Id_2_Numb["MN    1 K"]
	Id_2_Numb["MN  MN L"] = Id_2_Numb["MN    1 L"]

where `MN    1 I` is the id extracted from the PDB by the script

Finally you can modify the filenames in _line 178 to 186_.


 4.2. Modifying the `analyse_pka_interactions.py`
-----------------------------------------------------------------------------

Please refer to `Example_analyse_pka_interactions.py` to see the modifications 
of the original script

First we need to specify the charged ligand in the residuepKa function
Here we have PH1 (amine of PHE) and MN which are positively charged and PH2 is 
negativelly charged so we add them to the lists:

	positive_charge = ['PH1   1', 'MN   1' ] 
	negative_charge = ['PH2   1'] 


This example script includes a function, transMat, that makes possible to 
average the interactions per chain, which requires a lot of operations that 
are specific to the organisation of protein and its ligands in space. some 
explanations are given in the script. But this part can be avoided is you are
not interrested in averaging per chain or if you have a monomeric protein.
If this funtion is used, it needs to be called in the matPrep function. Note 
the modifications in this function in the `Example_analyse_pka_interactions.py`
scripts to obtain the averaged interaction maps for both the dimer and the 
tetramer.

Finally the name of the files can be changed

