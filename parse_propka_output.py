from __future__ import print_function
import re, glob, os, gc, sys
import numpy as np

#############################################################################
#
# USAGE 
#        python parse_propka_output.py
#              
# DESCRIPTION
#        This python script parses the propka output files and generate a number of
#		 files storing all the data of interest and which can then be further processed  
#
# USER PARAMETERS
#
#		 line 47: modify [CAMN] depending on your ligand residues.
# 				  for example if ZN2+ is present and an organic ligand that does 
#				  not include CA but other atom types, including for example C1, then you can 
#				  simply replace [CAMN] with [C1ZN]
#
#		 line 161 to 166: Modify according to your ligands 
#
#
#    	 line 178: replace 'structure.pdb' with the name of your structure file 
#    			  (typically xray structure or first frame)
#
# 		 lines 181 to 186: adjust the name of the output files according to your needs
#
# AUTHOR
#        Eric Lang <eric.jm.lang@gmail.com>
#
# Version 1.6 - 27/05/2015
#
#############################################################################



def listResid(pdb_file):
    """A function to identify all the residues in a PDB file
    and return a list containing residue name, residue number and chain"""
    Identifier = []
    Chains =[]
    Chains_het = []
    Resid_count = 0
	# Here the SearchHETATM searches for CA and MN: needs to be modified if needed
    SearchPDB='^[ATOM\s\d]{11}\s+[CA]{2}\s+(\w{2}[\w\s])\s(\w)\s+(\d+)'
    SearchHETATM='^[HETATM\s\d]{11}\s+[CAMN]{2}\s+\w{2}[\w\s]\s(\w)\s+\d+'
    for Line in pdb_file:
        Line = Line.strip('\n')
        if re.search(SearchPDB, Line):
            Result = re.search(SearchPDB, Line)
            Resname = str(Result.group(1))
            Chain = str(Result.group(2))
            Resid = int(Result.group(3))
            Id2 = "%s %3s %s" % (Resname, Resid, Chain) # %3s is needed for the formating
            Identifier.append(Id2)
            if Chain not in Chains:
                Chains.append(Chain)
            if Resid > Resid_count:
                Resid_count = Resid
        if re.search(SearchHETATM, Line):
            Result_het = re.search(SearchHETATM, Line)
            Chain_het = str(Result_het.group(1))
            if Chain_het not in Chains_het :
                Chains_het .append(Chain_het)
    pdb_file.close()
    return(Identifier, Chains,  Resid_count, Chains_het)


def parsePropka1(propka_output):
    """Parse the output from propka and store the results of interest in lists"""
    result_pka_file = open(propka_output, "r")
    list_results = []
    for l in result_pka_file:
        if not l.strip():
            continue
        else:
            if len(l.strip().split()) == 22:
                list_results.append([l.strip().split()[0], l.strip().split()[1], l.strip().split()[2], l.strip().split()[3], l.strip().split()[6], l.strip().split()[8]])
    result_pka_file.close()
    return(list_results)


def parsePropka2(propka_output):
    """Parse the output from propka and store the results of interest in lists"""
    result_pka_file = open(propka_output, "r")
    list_interactions = []
    for l in result_pka_file:
        if not l.strip():
            continue
        else:
            if len(l.strip().split()) == 22:
                list_interactions.append([l.strip().split()[0], l.strip().split()[1], l.strip().split()[2], l.strip().split()[10], l.strip().split()[11], l.strip().split()[12], l.strip().split()[13], l.strip().split()[14], l.strip().split()[15], l.strip().split()[16], l.strip().split()[17], l.strip().split()[18], l.strip().split()[19], l.strip().split()[20], l.strip().split()[21]])
            if len(l.strip().split()) == 15 and l.strip().split()[0] != 'Could':
                list_interactions.append([l.strip().split()[0], l.strip().split()[1], l.strip().split()[2], l.strip().split()[3], l.strip().split()[4], l.strip().split()[5], l.strip().split()[6], l.strip().split()[7], l.strip().split()[8], l.strip().split()[9], l.strip().split()[10], l.strip().split()[11], l.strip().split()[12], l.strip().split()[13], l.strip().split()[14]])
    result_pka_file.close()
    return(list_interactions)



def extractValuepKa(pkalist, desolvlist, FileList,OutFile1, OutFile2):
    for InfileName in FileList:
        list_results = parsePropka1(InfileName)
        RecordNum = 0
        for l in list_results:
            resname = l[0]
            resid = l[1]
            chain = l[2]
            Id = "%s %3s %s" % (resname, resid, chain)
            pka = l[3]
            desolv = float(l[4]) + float(l[5])
            if "*" in pka:
                pka = pka.replace("*", "")
            #pka = float(pka)
            #pkavalue[RecordNum,FileNum] = float(pka)
            pkalist[RecordNum] = str(pkalist[RecordNum]) 
            pkalist[RecordNum] += '\t' + pka
            desolvlist[RecordNum] = str(desolvlist[RecordNum]) 
            desolvlist[RecordNum] += '\t' + str(desolv)
            RecordNum += 1 
    for Item in pkalist:
        OutFile1.write(Item + '\n')
    for Item in desolvlist:
        OutFile2.write(Item + '\n')
    OutFile1.close()
    OutFile2.close()
    del list_results
    del desolvlist
    del pkalist
    return

def dictionnaries(Identifier, Chains, Resid_count, Chains_het, List, OutFile3):
    pkalist = []
    desolvlist = []
    # create a dictionnary to relate residue identifier to a number
    Chain_2_Numb={}
    for i in range(0, len(Chains)):
        Chain_2_Numb[Chains[i]] = i
    Id_2_Numb = {}
    for j in range(0, len(Identifier)):
        Id_2_Numb[Identifier[j]] = j
    # create a reverse dictionnary
    Numb_2_Id={}
    for item in Id_2_Numb:
        Numb = Id_2_Numb[item]
        Numb_2_Id[Numb]=item
    list_results = parsePropka1(List)
    for l in list_results:
        resname = l[0]
        resid = str(l[1])
        chain = str(l[2])
        ID = "%s %3s %s" % (resname, resid, chain)
        pkalist.append(ID)
        desolvlist.append(ID)
        # deal with N terminal
        if resname == 'N+':
            Id_2_Numb[ID]= Chain_2_Numb[chain]*Resid_count 
        # deal with C terminal
        if resname == 'C-':
            Id_2_Numb[ID]= Chain_2_Numb[chain]*Resid_count + (Resid_count-1)
        if chain in Chains_het: # these chains are for the ligand can be removed if ligand does not have a pKa calculated by PROPKA
            if resid == 'N': # In PROPKA for a ligand, the region normally corresponding to the residue number is replaced with the ionisable atom of the ligand. Here an amine group for example (N)
                Id_2_Numb[ID]= (len(Chains) - len(Chains_het))*Resid_count + (Chain_2_Numb[chain] - (len(Chains) - len(Chains_het)))

    # If the ligand does not have a pKa but interact with ionizable residues nonetheless, it needs to be added to the list:
    Id_2_Numb["MN  MN I"] = Id_2_Numb["MN    1 I"] # Here the manganese ion in chain I which appear like this in .pka file: "MN  MN I", is given an ID similar to the others "MN    1 I", note the number of spaces: 3 as usual + 1 beacuse the residue name only include three characteres instead of three 

    keys = list(Numb_2_Id.keys())
    keys.sort()
    for i in range(0,len(keys)):
        OutFile3.write("%s %5s\n" % (keys[i], Numb_2_Id[keys[i]]))
    OutFile3.close()
    return(pkalist, desolvlist, Id_2_Numb, Numb_2_Id)


# Open input files
FileList = glob.glob('*.pka')
pdb_file = open('structure.pdb', 'r')

# Open Output files
OutMat1 = open('sidechain_Hbond_mat.dat', 'w')
OutMat2 = open('backbone_Hbond_mat.dat', 'w')
OutMat3 = open('coulombic_mat.dat', 'w')
OutFile1 = open('pka_all.dat', 'w')
OutFile2 = open('desolvation_all.dat', 'w')
OutFile3 = open('residues.dat', 'w')

# Read PDB file
(Identifier, Chains, Resid_count, Chains_het)=listResid(pdb_file)

#calculate the number of file to process
number_pkafiles = 0
for file in os.listdir('.'):    
    if file.endswith(".pka"):
        number_pkafiles +=1
print("Number of .pka files to process: %d" % (number_pkafiles))


(pkalist, desolvlist, Id_2_Numb, Numb_2_Id)=dictionnaries(Identifier, Chains, Resid_count, Chains_het, FileList[0], OutFile3)

Mat = np.zeros((3,len(Numb_2_Id),len(Numb_2_Id)))
for InfileName in FileList:
    list_interactions = parsePropka2(InfileName)
    Matrix = np.zeros((3,len(Numb_2_Id),len(Numb_2_Id)))
    for l in list_interactions:
        resname_int = l[0]
        resid_int = l[1]
        chain_int = l[2]
        Id_int = "%s %3s %s" % (resname_int, resid_int, chain_int)
        SH_pka = float(l[3])
        SH_resname = l[4]
        SH_resid = l[5]
        SH_chain = l[6]
        Id_SH = "%s %3s %s" % (SH_resname, SH_resid, SH_chain)
        BH_pka = float(l[7])
        BH_resname = l[8]
        BH_resid = l[9]
        BH_chain = l[10]
        Id_BH = "%s %3s %s" % (BH_resname, BH_resid, BH_chain)
        CI_pka = float(l[11])
        CI_resname = l[12]
        CI_resid = l[13]
        CI_chain = l[14]
        Id_CI = "%s %3s %s" % (CI_resname, CI_resid, CI_chain)
        ID1 = Id_2_Numb[Id_int]
        if Id_SH != 'XXX   0 X':
            ID2 = Id_2_Numb[Id_SH]
            Matrix[0,ID1,ID2] = SH_pka
        if Id_BH != 'XXX   0 X':
            ID3 = Id_2_Numb[Id_BH]
            Matrix[1,ID1,ID3] = BH_pka
        if Id_CI != 'XXX   0 X':
            ID4 = Id_2_Numb[Id_CI]
            Matrix[2,ID1,ID4] = CI_pka
    for i in range(0,3):
        triu = np.triu(Matrix[i])
        tril = np.tril(Matrix[i])
        tril = tril.T
        mask = np.isclose(triu,-tril) & ~np.isclose(triu,0)
        Matrix[i] = np.abs(Matrix[i])
        Matrix[i,mask]=0
        Matrix[i] = np.tril(Matrix[i],-1).transpose() + np.triu(Matrix[i],0) # check if diagonal can be non 0
        Mat[i] +=Matrix[i]
np.savetxt(OutMat1, Mat[0], fmt='%.2f', delimiter=' ', newline='\n')
np.savetxt(OutMat2, Mat[1], fmt='%.2f', delimiter=' ', newline='\n')
np.savetxt(OutMat3, Mat[2], fmt='%.2f', delimiter=' ', newline='\n')
OutMat1.close()
OutMat2.close()
OutMat3.close()


extractValuepKa(pkalist, desolvlist, FileList,OutFile1, OutFile2)



