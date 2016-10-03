#! /usr/bin/python

from __future__ import print_function
import re, math, os
import numpy as np

#############################################################################
#
# USAGE 
#        python analyse_desolvation.py
#              
# DESCRIPTION
#        This python script calculates the average total pKa value per chain and
#	     the standard deviation
#
# AUTHOR
#        Eric Lang <eric.jm.lang@gmail.com>
#
# Version 1.0 - 27/05/2015
#
#############################################################################

filename = 'desolvation_all.dat'
InFile1 = open(filename, 'r')
#InFile2 = open('Identifier.dat', 'r')
OutFile1 = open('AVERAGE_' + filename, 'w')
OutFile2 = open('RMSD_' + filename, 'w')

source = np.genfromtxt(InFile1, delimiter='\t')

x=(len(source))/4
A = source[0:x,1:len(source[0])]
B = source[x:2*x,1:len(source[0])]
C = source[2*x:3*x,1:len(source[0])]
D = source[3*x:4*x,1:len(source[0])]

AB=np.append(A,B, axis=1)
CD=np.append(C,D, axis=1)
Results=np.append(AB,CD, axis=1)



#A = source[0:x,1:len(source[0])]
#B = source[x:2*x,1:len(source[0])]
#C = source[2*x:3*x,1:len(source[0])]
#D = source[3*x:4*x,1:len(source[0])]
#E1 = source[4*x:4*x+1,1:len(source[0])]
#E2 = source[4*x+1:4*x+2,1:len(source[0])]
#F1 = source[4*x+2:4*x+3,1:len(source[0])]
#F2 = source[4*x+3:4*x+4,1:len(source[0])]
#G1 = source[4*x+4:4*x+5,1:len(source[0])]
#G2 = source[4*x+5:4*x+6,1:len(source[0])]
#H1 = source[4*x+6:4*x+7,1:len(source[0])]
#H2 = source[4*x+7:4*x+8,1:len(source[0])]#

#AB=np.append(A,B, axis=1)
#CD=np.append(C,D, axis=1)
#ABCD=np.append(AB,CD, axis=1)#

#EF1=np.append(E1,F1, axis=1)
#EF2=np.append(E2,F2, axis=1)
#GH1=np.append(G1,H1, axis=1)
#GH2=np.append(G2,H2, axis=1)
#EFGH1=np.append(EF1,GH1, axis=1)
#EFGH2=np.append(EF2,GH2, axis=1)#

#Results = np.append(ABCD, EFGH1, axis=0)
#Results = np.append(Results, EFGH2, axis=0)

Averagepka = np.mean(Results, axis=1)
RMSDpka = np.std(Results, axis=1)


np.savetxt(OutFile1, Averagepka, fmt='%.2f', delimiter=' ', newline='\n')
np.savetxt(OutFile2, RMSDpka, fmt='%.2f', delimiter=' ', newline='\n')
InFile1.close()
OutFile1.close()
OutFile2.close()


