#! /usr/bin/python

from __future__ import print_function
import numpy as np

#############################################################################
#
# USAGE 
#        python analyse_pka_interactions.py
#              
# DESCRIPTION
#        This python script parses the propka output files and generate a number of
#		 files storing all the data of interest and which can then be further processed  
#
# CUSTOM MODIFICATIONS
#
#		 line 45: 
#
#
#		 line 184: 
#
# 		 lines 
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
    list_residu = []
    list_pka = []
    Identif=[]
    Name = []
    positive_charge = ['PH1   1', 'MN   1' ] # needs to be adapted to user's needs
    negative_charge = ['PH2   1'] # needs to be adapted to user's needs
    neutral = []
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
                    if Id2 == "PH2   1" or Id1 == "PH2   1" or Id2 == "PH1   1" or Id1 == "PH1   1" : # avoid attraction reppultion on the same residue
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
    str += "mol new structure.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
    str += "mol delrep 0 top\n"
    str += "mol representation NewCartoon 0.300000 30.000000 4.100000 0\n"
    str += "mol color colorID 8\n"
    str += "mol selection {chain A B}\n"
    str += "mol material Transparent\n"
    str += "mol addrep top\n"
    return(str)


def getNodesEdges(chain1, resid1, atname1, chain2, resid2, atname2, weight, colorid):
    """Returns the section of the visualisation script dedicated to plotting nodes and edges"""
    # to do: depending on residue select the atoms wich will form the interactions (ex N or O for backbone, etc.)
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


def transMat(Matrix, x):
    """This function perform the averaging per chain"""
    # matrix has the form:
    # #
    # AA AB AC AD 0   AF1 AG1 0   AI 0  0  0  0    AF2  AG2  0
    # 0  BB BC BD 0   BF1 BG1 0   0  0  0  BL 0    BF2  BG2  0
    # 0  0  CC CD CE1 0   0   CH1 0  0  CK 0  CE2  0    0    CH2
    # 0  0  0  DD DE1 0   0   DH1 0  DJ 0  0  DE2  0    0    DH2
    # 0  0  0  0  0   0   0   0   0  0  0  0  E1E2 0    0    0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    F1F2 0    0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    G1G2 0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    0    H1H2
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    0    0
    # .....
    # 0  0  0  0  0   0   0   0   0  0  0  0   0   0    0    0
    # #
    # define each submatrices. Very crude at the momeent.
    # might simplified with a dictionary for example.
    # Here AA means interaction of chain A with Chain A, AB means interaction of chain A with Chain B, etc.
    # A, B, C, D are the 4 chains of the tetrameric protein
    # E1, F1, G1, H1, correspond to the amine funtion of allosteric PHE
    # I, J, K, L, correspond to the MN ions
    # E2, F2, G2, H2, correspond to the carboxylate funtion of allosteric PHE
    # Chain I (Mn), only interacts with chain A so AI exists but the other (BI, CI, DI) don't
    AA = Matrix[0:x,0:x]
    BB = Matrix[x:2*x,x:2*x]
    CC = Matrix[2*x:3*x,2*x:3*x]
    DD = Matrix[3*x:4*x,3*x:4*x]
    # for the following matrices, we need to average by chain
    # so e.g. interactions betw. residue 1 chain A and residue 20 chain B = interactions betw. residue 20 chain A and residue 1 chain B
    AB = TriangOp(Matrix[0:x,x:2*x])
    AC = TriangOp(Matrix[0:x,2*x:3*x])
    AD = TriangOp(Matrix[0:x,3*x:4*x])
    BC = TriangOp(Matrix[x:2*x,2*x:3*x])
    BD = TriangOp(Matrix[x:2*x,3*x:4*x])
    CD = TriangOp(Matrix[2*x:3*x,3*x:4*x])
    # the following are single values
    CE1 = Matrix[2*x:3*x,4*x:4*x+1]
    DE1 = Matrix[3*x:4*x,4*x:4*x+1]
    AF1 = Matrix[0:x,4*x+1:4*x+2]
    BF1 = Matrix[x:2*x,4*x+1:4*x+2]
    AG1 = Matrix[0:x,4*x+2:4*x+3]
    BG1 = Matrix[x:2*x,4*x+2:4*x+3]
    CH1 = Matrix[2*x:3*x,4*x+3:4*x+4]
    DH1 = Matrix[3*x:4*x,4*x+3:4*x+4]
    AI = Matrix[0:x,4*x+4:4*x+5]
    BL = Matrix[x:2*x,4*x+7:4*x+8]
    CK = Matrix[2*x:3*x,4*x+6:4*x+7]
    DJ = Matrix[3*x:4*x,4*x+5:4*x+6]
    CE2 = Matrix[2*x:3*x,4*x+8:4*x+9]
    DE2 = Matrix[3*x:4*x,4*x+8:4*x+9]
    AF2 = Matrix[0:x,4*x+9:4*x+10]
    BF2 = Matrix[x:2*x,4*x+9:4*x+10]
    AG2 = Matrix[0:x,4*x+10:4*x+11]
    BG2 = Matrix[x:2*x,4*x+10:4*x+11]
    CH2 = Matrix[2*x:3*x,4*x+11:4*x+12]
    DH2 = Matrix[3*x:4*x,4*x+11:4*x+12]
    E1E2 = Matrix[4*x:4*x+1,4*x+8:4*x+9]
    F1F2 = Matrix[4*x+1:4*x+2,4*x+9:4*x+10]
    G1G2 = Matrix[4*x+2:4*x+3,4*x+10:4*x+11]
    H1H2 = Matrix[4*x+3:4*x+4,4*x+11:4*x+12]

    XX = AA + BB + CC + DD
    XY = AB + CD
    XZ = AC + BD
    YZ = AD + BC
    XL1 = AG1 + BF1 + CH1 + DE1
    YL1 = AF1 + BG1 + CE1 + DH1
    XL2 = AG2 + BF2 + CH2 + DE2
    YL2 = AF2 + BG2 + CE2 + DH2
    XM = AI + BL + CK + DJ
    L1L2 = E1E2 + F1F2 + G1G2 + H1H2

    o = np.zeros((1,1))
    O = np.zeros((x,1))
    OO = np.zeros((x,x))

    M1=np.concatenate((XX, XY, XZ, YZ, O, YL1, XL1, O, XM, O, O, O, O, YL2, XL2, O), axis=1)
    M2=np.concatenate((OO, XX, YZ, XZ, O, XL1, YL1, O, O, O, O, XM, O, XL2, YL2, O), axis=1)
    M3=np.concatenate((OO, OO, XX, XY, YL1, O, O, XL1, O, O, XM, O, YL2, O, O, XL2), axis=1)
    M4=np.concatenate((OO, OO, OO, XX, XL1, O, O, YL1, O, XM, O, O, XL2, O, O, YL2), axis=1)
    M5=np.concatenate((np.zeros((1,4*x+8)), L1L2, o, o, o), axis=1)
    M6=np.concatenate((np.zeros((1,4*x+8)), o, L1L2, o, o), axis=1)
    M7=np.concatenate((np.zeros((1,4*x+8)), o, o, L1L2, o), axis=1)
    M8=np.concatenate((np.zeros((1,4*x+8)), o, o, o, L1L2), axis=1)

    Tetramer = np.concatenate((M1, M2, M3, M4, M5, M6, M7, M8, np.zeros((8,4*x+12))), axis=0)/4

    # matrix Tetramer has  the form:
    # #
    # XX XY XZ YZ 0   YL1 XL1 0   XM 0  0  0  0    YL2  XL2  0
    # 0  XX YZ XZ 0   XL1 YL1 0   0  0  0  XM 0    XL2  YL2  0
    # 0  0  XX XY YL1 0   0   XL1 0  0  XM 0  YL2  0    0    XL2
    # 0  0  0  XX XL1 0   0   YL1 0  XM 0  0  XL2  0    0    YL2
    # 0  0  0  0  0   0   0   0   0  0  0  0  L1L2 0    0    0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    L1L2 0    0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    L1L2 0
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    0    L1L2
    # 0  0  0  0  0   0   0   0   0  0  0  0  0    0    0    0
    # .....
    # 0  0  0  0  0   0   0   0   0  0  0  0   0   0    0    0
    # #

    N1=np.concatenate((XX, XY, YL1, XL1, XM, O, YL2, XL2), axis=1)
    N2=np.concatenate((OO, XX, XL1, YL1, O, XM, XL2, YL2), axis=1)
    N3=np.concatenate((np.zeros((1,2*x+4)), L1L2, o), axis=1)
    N4=np.concatenate((np.zeros((1,2*x+4)), o, L1L2), axis=1)


    Dimer = np.concatenate((N1, N2, N3, N4, np.zeros((4,2*x+6))), axis=0)/4

    # matrix Dimer has  the form:
    # #
    # XX XY YL1 XL1 XM 0  YL2  L2
    # 0  XX XL1 YL1 0  XM XL2  YL2
    # 0  0  0   0   0  0  L1L2 0
    # 0  0  0   0   0  0  0    L1L2
    # 0  0  0   0   0  0  0    0
    # .....
    # 0  0  0   0   0  0  0    0
    # #

    return(Tetramer, Dimer)

def matPrep(InMat_phe, n_phe, InMat_apo, n_apo, type, color1, color2, color3, color4, OutVizA, OutVizB, OutVizC, OutVizD, x, residue_name1, residue_name2,  cutoff):

    Matrix_phe = np.array([[float(i)/n_phe for i in line.split()] for line in InMat_phe])
    Matrix_apo = np.array([[float(i)/n_apo for i in line.split()] for line in InMat_apo])

    size_diff =  len(Matrix_apo) - len(Matrix_phe)

    if size_diff < 0:
        Matrix_apo=np.hstack((Matrix_apo,np.zeros((Matrix_apo.shape[0], abs(size_diff)))))
        Matrix_apo=np.vstack((Matrix_apo,np.zeros((abs(size_diff), Matrix_apo.shape[1]))))

    if size_diff > 0:
        Matrix_phe=np.hstack((Matrix_phe,np.zeros((Matrix_phe.shape[0], abs(size_diff)))))
        Matrix_phe=np.vstack((Matrix_phe,np.zeros((abs(size_diff), Matrix_phe.shape[1]))))

    # The following two lines can be commented so no averaging per chain is conducted
    (Tetramer_phe, Dimer_phe) = transMat(Matrix_phe, x)
    (Tetramer_apo, Dimer_apo) = transMat(Matrix_apo, x)

    Tetramer_apo_phe = Tetramer_apo - Tetramer_phe
    Tetramer_phe_apo = Tetramer_phe - Tetramer_apo

    Dimer_apo_phe = Dimer_apo - Dimer_phe
    Dimer_phe_apo = Dimer_phe - Dimer_apo

    Viz1 = writeTCL(Tetramer_apo_phe, type, color1, color2, residue_name1, cutoff, neutral_apo, negative_charge_apo, positive_charge_apo)
    OutVizA.write(Viz1)
    Viz2 = writeTCL(Tetramer_phe_apo, type, color3, color4, residue_name1, cutoff, neutral_phe, negative_charge_phe, positive_charge_phe)
    OutVizB.write(Viz2)
    Viz3 = writeTCL(Dimer_apo_phe, type, color1, color2, residue_name2,cutoff, neutral_apo, negative_charge_apo, positive_charge_apo)
    OutVizC.write(Viz3)
    Viz4 = writeTCL(Dimer_phe_apo, type, color3, color4, residue_name2, cutoff, neutral_phe, negative_charge_phe, positive_charge_phe)
    OutVizD.write(Viz4)

    #np.savetxt('apo-phe_tetra_%s.dat' %(type), Tetramer_apo_phe,  fmt='%.6f', delimiter=' ', newline='\n')
    #np.savetxt('apo-phe_dimer_%s.dat' %(type), Dimer_apo_phe,  fmt='%.6f', delimiter=' ', newline='\n')

    return



# Open input files
InMat1 = open('sidechain_Hbond_mat_phes.dat', 'r')
InMat2 = open('backbone_Hbond_mat_phe.dat', 'r')
InMat3 = open('coulombic_mat_phe.dat', 'r')
InMat4 = open('sidechain_Hbond_mat_apo.dat', 'r')
InMat5 = open('backbone_Hbond_mat_apo.dat', 'r')
InMat6 = open('coulombic_mat_apo.dat', 'r')
InFile1=open('Identifier_phe.dat', 'r')
InFile2=open('AVERAGE_pka_phe_.dat', 'r')
InFile3=open('Identifier_apo.dat', 'r')
InFile4=open('AVERAGE_pka_apo.dat', 'r')


# Open Output files
OutViz1 = open('sidechain_Hbond_tetramer_apo.tcl', 'w')
OutViz2 = open('sidechain_Hbond_tetramer_phe.tcl', 'w')
OutViz3 = open('sidechain_Hbond_dimer_apo.tcl', 'w')
OutViz4 = open('sidechain_Hbond_dimer_phe.tcl', 'w')

OutViz5 = open('backbone_Hbond_tetramer_apo.tcl', 'w')
OutViz6 = open('backbone_Hbond_tetramer_phe.tcl', 'w')
OutViz7 = open('backbone_Hbond_dimer_apo.tcl', 'w')
OutViz8 = open('backbone_Hbond_dimer_phe.tcl', 'w')

OutViz9 = open('coulombic_tetramer_apo.tcl', 'w')
OutViz10 = open('coulombic_tetramer_phe.tcl', 'w')
OutViz11 = open('coulombic_dimer_apo.tcl', 'w')
OutViz12 = open('coulombic_dimer_phe.tcl', 'w')

residue_name1 = 'residues.dat'
residue_name2 = 'residues_dimer.dat'



n_phe = 23999
n_apo = 23999
x = 351
ph = 7.30
cutoff = 0.1

(neutral_phe, negative_charge_phe, positive_charge_phe) = residuepKa(InFile1, InFile2)
(neutral_apo, negative_charge_apo, positive_charge_apo) = residuepKa(InFile3, InFile4)

matPrep(InMat1, n_phe, InMat4, n_apo, 'SH', 19, 29, 23, 29, OutViz1, OutViz2, OutViz3, OutViz4, x, residue_name1, residue_name2, cutoff)
matPrep(InMat2, n_phe, InMat5, n_apo, 'BH', 19, 29, 23, 29, OutViz5, OutViz6, OutViz7, OutViz8, x, residue_name1, residue_name2, cutoff)
matPrep(InMat3, n_phe, InMat6, n_apo, 'CI', 19, 29, 23, 29, OutViz9, OutViz10, OutViz11, OutViz12, x, residue_name1, residue_name2, cutoff)


