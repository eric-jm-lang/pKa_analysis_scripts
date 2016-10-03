#! /usr/bin/python

import re, glob

#############################################################################
#
# USAGE 
#        python convert_PDB_files.py
#              
# DESCRIPTION
#        This python script convert a NAMD/CHARMM type of PDB file into a 
#		 more standard PDB format that can be understood by PROPKA
#		 The script process all the PDB files in a folder and generate 
#		 the modified PDB files with New_ in front of the file name
#
# AUTHOR
#        Eric Lang <eric.jm.lang@gmail.com>
#
# Version 1.2 - 26/11/2014
#
#############################################################################


# read all the PDB files in the current directory
FileList = glob.glob('*.pdb')

# Differences CHARMM vs PDB
SearchStr_Ile = ' CD  ILE'
SearchStr_HSD = ' HSD '
SearchStr_HSE = ' HSE '
SearchStr_HSP = ' HSP '
SearchStr_OT1 = ' OT1 '
SearchStr_OT2 = ' OT2 '
# Replacement strings
New_ILE = ' CD1 ILE'
New_HIS = ' HIS '
New_OT1 = ' O   '
New_OT2 = ' OXT '

for PDB_File  in FileList:
    InFile = open(PDB_File , 'r')
    OutFile = open('New_' + PDB_File, 'w')
    for Line in InFile:
        Line = Line.strip('\n')
        if re.search(SearchStr_Ile, Line):
            Line = re.sub(SearchStr_Ile, New_ILE, Line)
        if re.search(SearchStr_HSD, Line):
            Line = re.sub(SearchStr_HSD, New_HIS, Line)
        if re.search(SearchStr_HSE, Line):
            Line = re.sub(SearchStr_HSE, New_HIS, Line)
        if re.search(SearchStr_HSP, Line):
            Line = re.sub(SearchStr_HSP, New_HIS, Line)
        if re.search(SearchStr_OT1, Line):
            Line = re.sub(SearchStr_OT1, New_OT1, Line)
        if re.search(SearchStr_OT2, Line):
            Line = re.sub(SearchStr_OT2, New_OT2, Line)

        else:
            Line = Line
        OutFile.write(Line + '\n')
    InFile.close()
    OutFile.close()

