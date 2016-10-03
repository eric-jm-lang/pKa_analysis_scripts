#! /usr/bin/python

from __future__ import print_function
import numpy as np

#############################################################################
#
# USAGE 
#        python analyse_pka_interactions.py
#              
# DESCRIPTION
#        This python script generate the interaction maps
#
# USER PARAMETERS
#
#		 line 65 & 66: Modify with your ligand residue names and residue number
#
#		 line 218: change structure.pdb filename if needed
#
#		 lines 278-301: change filenames if needed
#
# 		 lines 304-308: change with your own values 
#
# AUTHOR
#        Eric Lang <eric.jm.lang@gmail.com>
#
# Version 1.6 - 27/05/2015
#
#############################################################################

def TriangOp(matrix):
    """A function to perform matrix operations"""
    matrix = np.tril(matrix,-1).transpose() + np.triu(matrix,0) + np.triu(matrix,1).transpose() + np.tril(matrix,0)
    return(matrix)


def dictionaries(residue_name):
    """A function that creates a dictionary of residues"""
    residue_file = open(residue_name, 'r')
    list_residues = []
    for l in residue_file:
        list_residues.append([l.strip().split()[0], l.strip().split()[1], l.strip().split()[2], l.strip().split()[3]])
    residue_file.close()
    Number=[]
    Identifier=[]
    for l in list_residues:
        Numb = int(l[0])
        resname = l[1]
        resid = l[2]
        chain = l[3]
        Id = "%s %3s %s" % (resname, resid, chain)
        Identifier.append(Id)
        Number.append(Numb)
    Numb_2_Id={}
    for i in range(0, len(list_residues)):
        Numb_2_Id[Number[i]]=Identifier[i]
    return(Numb_2_Id)


def residuepKa(File1, File2):
    """Identify which residues are positively charged, which ones are negatively charged and which ones are neutral based on the pH"""
    list_residu = []
    list_pka = []
    Identif=[]
    Name = []
    positive_charge = ['LIG   1', 'MN   1' ] # charged ligand -  needs to be adapted to user's needs: residue name, three spaces and residue number if nothing to add keep empty list []
    negative_charge = [] # charged ligands - needs to be adapted to user's needs: residue name, three spaces and residue number. if nothing to add keep empty list []
    neutral = [] # neutral ligand 
    pos_res = ['ARG', 'LYS', 'HIS']
    neg_res = ['GLU', 'ASP', 'CYS', 'TYR']

    lines1 = File1.readlines()
    for l in lines1:
        list_residu.append([l.strip().split()[0], l.strip().split()[1]])
    #File1.close()
    for l in list_residu:
        resname = l[0]
        resid = l[1]
        Id = "%s %3s" % (resname, resid)
        Identif.append(Id)
        Name.append(resname)
    #print(Identif)

    lines2 = File2.readlines()
    for l in lines2:
        list_pka.append(float(l.strip().split()[0]))
    #File2.close()
    #print(list_pka)
    for i in range(0, len(Identif)):
        if Name[i] in neg_res:  
            if list_pka[i] > ph:
                neutral.append(Identif[i])
            if list_pka[i] < ph:
                negative_charge.append(Identif[i])
        if Name[i] in pos_res:
            if list_pka[i] > ph:
                positive_charge.append(Identif[i])
            if list_pka[i] < ph:
                neutral.append(Identif[i])
    print("Neutral\n")
    print(neutral)
    print("\n")
    print("Negative\n")
    print(negative_charge)
    print("\n")
    print("Positive\n")
    print(positive_charge)
    print("\n")
    return(neutral, negative_charge, positive_charge)

def writeTCL(Matrix,type, color1, color2, residue_name, cutoff, neutral, negative_charge, positive_charge):
    """Write the TCL script for displaying interaction in VMD"""
    Numb_2_Id=dictionaries(residue_name)
    str  = "%s\n" % ( getTCLHeader() ) # write TCL header
    for i in range(len(Matrix)):
        for j in range(len(Matrix[i])):
            if Matrix[i][j]>= cutoff:
                weight = abs(Matrix[i][j])
                id1 = Numb_2_Id[i]
                id2 = Numb_2_Id[j]
                id1 = id1.split(' ')
                id1 = [x.rstrip('\n') for x in id1 if x != '']
                name1 = id1[0]
                chain1 = id1[2]
                resid1 = int(id1[1])
                id2 = id2.split(' ')
                id2 = [x.rstrip('\n') for x in id2 if x != '']
                name2 = id2[0]
                chain2 = id2[2]
                resid2 = int(id2[1])
                Id1 = "%s %3s" % (name1, resid1)
                Id2 = "%s %3s" % (name2, resid2)
                #print(Id1, Id2)
                if type == 'SH':
                    if Id1 in negative_charge and Id2 in negative_charge: # no possible side chain hbond between two deprotonated acids
                        pass
                    else:
                        str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color1)
                if type == 'BH':
                    str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color1)
                if type == 'CI':
                    if Id1 in positive_charge and Id2 in positive_charge:
                        str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color2)
                    if Id1 in negative_charge and Id2 in negative_charge:
                        str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color2)
                    if Id1 in positive_charge and Id2 in negative_charge:
                        str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color1)
                    if Id2 in positive_charge and Id1 in negative_charge:
                        str += getNodesEdges(chain1, resid1, 'CA', chain2, resid2, 'CA', weight, color1)
                    
                    else:
                        pass
    return(str)


def getTCLHeader():
    """Returns the header of the VMD visualisation script"""
    str  = "# This is an automatically generated tcl script for vmd.\n"
    str += "#\n"
    str += "# The purpose of this script is to visualise the H-bonds\n"
    str += "# and coulombic interactions responsible for the pKa variations\n"
    str += "# of each ionisable residues over the course of a MD trajectory\n"
    str += "proc vmdrestoremymaterials {} {\n"
    str += "set mlist { Opaque Transparent }\n"
    str += "set mymlist [material list]\n"
    str += "foreach mat $mlist {\n"
    str += "if { [lsearch $mymlist $mat] == -1 } {\n" 
    str += "material add $mat\n"
    str += "  }\n"
    str += "}\n"
    str += "material change ambient Opaque 0.050000\n"
    str += "material change diffuse Opaque 0.850000\n"
    str += "material change specular Opaque 0.000000\n"
    str += "material change shininess Opaque 0.000000\n"
    str += "material change mirror Opaque 0.000000\n"
    str += "material change opacity Opaque 1.000000\n"
    str += "material change outline Opaque 0.000000\n"
    str += "material change outlinewidth Opaque 0.000000\n"
    str += "material change transmode Opaque 0.000000\n"
    str += "material change ambient Transparent 0.370000\n"
    str += "material change diffuse Transparent 0.330000\n"
    str += "material change specular Transparent 0.350000\n"
    str += "material change shininess Transparent 0.200000\n"
    str += "material change mirror Transparent 0.000000\n"
    str += "material change opacity Transparent 0.250000\n"
    str += "material change outline Transparent 2.850000\n"
    str += "material change outlinewidth Transparent 0.720000\n"
    str += "material change transmode Transparent 0.000000\n"
    str += "}\n"
    str += "vmdrestoremymaterials\n"
    str += "display eyesep       0.065000\n"
    str += "display focallength  2.000000\n"
    str += "display height       6.000000\n"
    str += "display distance     -2.000000\n"
    str += "display projection   Orthographic\n"
    str += "display nearclip set 0.010000\n"
    str += "display farclip  set 10.000000\n"
    str += "display depthcue   on\n"
    str += "display cuestart   0.500000\n"
    str += "display cueend     10.000000\n"
    str += "display cuedensity 0.200000\n"
    str += "display cuemode    Exp2\n"
    str += "proc vmdrestoremycolors {} {\n"
    str += "  set colorcmds {\n"
    str += "    {color Display {Background} white}\n"
    str += "  }\n"
    str += "  foreach colcmd $colorcmds {\n"
    str += "    set val [catch {eval $colcmd}]\n"
    str += "  }\n"
    str += "  color change rgb 19 0.67 0.86 0.16\n" # My green: 19
    str += "  color change rgb 23 0.29 0.51 0.71\n" # My blue: 23
    str += "  color change rgb 26 0.67 0.60 0.75\n" # My purple: 26
    str += "  color change rgb 29 0.91 0.20 0.19\n" # My red: 29
    str += "  color change rgb 31 0.96 0.62 0.34\n" # My orange: 31
    str += "}\n"
    str += "vmdrestoremycolors\n"
    str += "mol new structure.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"  # change structure.pdb filename if needed
    str += "mol delrep 0 top\n"
    str += "mol representation NewCartoon 0.300000 30.000000 4.100000 0\n"
    str += "mol color colorID 8\n"
    str += "mol selection {chain A B}\n"
    str += "mol material Transparent\n"
    str += "mol addrep top\n"
    return(str)


def getNodesEdges(chain1, resid1, atname1, chain2, resid2, atname2, weight, colorid):
    """Returns the section of the visualisation script dedicated to plotting nodes and edges"""
    str  = "mol representation VDW 0.800000 30.000000\n"
    str += "mol color ColorID %d\n" % (colorid)
    str += "mol selection {chain %s and resid %d and name %s}\n" % (chain1, resid1, atname1)
    str += "mol material Opaque\n"
    str += "mol addrep top\n"
    str += "mol representation VDW 0.800000 30.000000\n"
    str += "mol color ColorID %d\n" % (colorid)
    str += "mol selection {chain %s and resid %d and name %s}\n" % (chain2, resid2, atname2)
    str += "mol material Opaque\n"
    str += "mol addrep top\n"
    str += "set sel1 [atomselect top \"chain %s and resid %d and name %s\"]\n" % (chain1, resid1, atname1)
    str += "set sel2 [atomselect top \"chain %s and resid %d and name %s\"]\n" % (chain2, resid2, atname2)
    str += "set coord1 [lindex [$sel1 get {x y z}] 0]\n"
    str += "set coord2 [lindex [$sel2 get {x y z}] 0]\n"
    str += "draw color %d\n" % (colorid)
    str += "graphics top cylinder $coord1 $coord2 radius %f resolution 30\n" % (weight)
    str += "$sel1 delete\n"
    str += "$sel2 delete\n"
    return(str)


def matPrep(InMat_holo, n_holo, InMat_apo, n_apo, type, color1, color2, color3, color4, OutVizA, OutVizB, x, residue_name1, cutoff):

    Matrix_holo = np.array([[float(i)/n_holo for i in line.split()] for line in InMat_holo])
    Matrix_apo = np.array([[float(i)/n_apo for i in line.split()] for line in InMat_apo])

    size_diff =  len(Matrix_apo) - len(Matrix_holo)

    if size_diff < 0:
        Matrix_apo=np.hstack((Matrix_apo,np.zeros((Matrix_apo.shape[0], abs(size_diff)))))
        Matrix_apo=np.vstack((Matrix_apo,np.zeros((abs(size_diff), Matrix_apo.shape[1]))))

    if size_diff > 0:
        Matrix_holo=np.hstack((Matrix_holo,np.zeros((Matrix_holo.shape[0], abs(size_diff)))))
        Matrix_holo=np.vstack((Matrix_holo,np.zeros((abs(size_diff), Matrix_holo.shape[1]))))

    Matrix_apo_holo = Matrix_apo - Matrix_holo
    Matrix_holo_apo = Matrix_holo - Matrix_apo

    Viz1 = writeTCL(Matrix_apo_holo, type, color1, color2, residue_name1, cutoff, neutral_apo, negative_charge_apo, positive_charge_apo)
    OutVizA.write(Viz1)
    Viz2 = writeTCL(Matrix_holo_apo, type, color3, color4, residue_name1, cutoff, neutral_holo, negative_charge_holo, positive_charge_holo)
    OutVizB.write(Viz2)

    return


# Open input files
InMat1 = open('sidechain_Hbond_mat_holo.dat', 'r')
InMat2 = open('backbone_Hbond_mat_holo.dat', 'r')
InMat3 = open('coulombic_mat_holo.dat', 'r')
InMat4 = open('sidechain_Hbond_mat_apo.dat', 'r')
InMat5 = open('backbone_Hbond_mat_apo.dat', 'r')
InMat6 = open('coulombic_mat_apo.dat', 'r')
InFile1=open('Identifier_holo.dat', 'r')
InFile2=open('AVERAGE_pka_holo_.dat', 'r')
InFile3=open('Identifier_apo.dat', 'r')
InFile4=open('AVERAGE_pka_apo.dat', 'r')


# Open Output files
OutViz1 = open('sidechain_Hbond_apo.tcl', 'w')
OutViz2 = open('sidechain_Hbond_holo.tcl', 'w')

OutViz5 = open('backbone_Hbond_apo.tcl', 'w')
OutViz6 = open('backbone_Hbond_holo.tcl', 'w')

OutViz9 = open('coulombic_apo.tcl', 'w')
OutViz10 = open('coulombic_holo.tcl', 'w')

residue_name1 = 'residues.dat'


n_holo = 23999 # number of pdb files  processed
n_apo = 23999 # number of pdb files processed 
x = 351 # number of protein residues in one chain
ph = 7.30 # working pH
cutoff = 0.1 # cutoff for displaying interactions (interaction maps not displayed if the difference in interactions between apo and holo is lower than the cutoff

(neutral_holo, negative_charge_holo, positive_charge_holo) = residuepKa(InFile1, InFile2)
(neutral_apo, negative_charge_apo, positive_charge_apo) = residuepKa(InFile3, InFile4)

# the number 19, 29, 23, 29 correspond to the colorID in VMD and can be modified

matPrep(InMat1, n_holo, InMat4, n_apo, 'SH', 19, 29, 23, 29, OutViz1, OutViz2, x, residue_name1,  cutoff)
matPrep(InMat2, n_holo, InMat5, n_apo, 'BH', 19, 29, 23, 29, OutViz5, OutViz6, x, residue_name1, cutoff)
matPrep(InMat3, n_holo, InMat6, n_apo, 'CI', 19, 29, 23, 29, OutViz9, OutViz10, x, residue_name1, cutoff)


