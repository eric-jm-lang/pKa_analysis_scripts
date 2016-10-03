#! /usr/bin/python

from __future__ import print_function
import re, math, os
import numpy as np

#############################################################################
#
# USAGE 
#        python analyse_pka_global.py
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


pkafile = 'pka_all.dat'
InFile1 = open(pkafile, 'r')
#InFile2 = open('Identifier.dat', 'r')
OutFile1 = open('AVERAGE_' + pkafile, 'w')
OutFile2 = open('RMSD_' + pkafile, 'w')

source = np.genfromtxt(InFile1, delimiter='\t')

x=(len(source))/4
A = source[0:x,1:len(source[0])]
B = source[x:2*x,1:len(source[0])]
C = source[2*x:3*x,1:len(source[0])]
D = source[3*x:4*x,1:len(source[0])]

AB=np.append(A,B, axis=1)
CD=np.append(C,D, axis=1)
Results=np.append(AB,CD, axis=1)

Averagepka = np.mean(Results, axis=1)
RMSDpka = np.std(Results, axis=1)


np.savetxt(OutFile1, Averagepka, fmt='%.2f', delimiter=' ', newline='\n')
np.savetxt(OutFile2, RMSDpka, fmt='%.2f', delimiter=' ', newline='\n')
InFile1.close()
OutFile1.close()
OutFile2.close()


