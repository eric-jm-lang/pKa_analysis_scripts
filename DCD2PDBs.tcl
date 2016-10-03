#############################################################################
#
# USAGE 
#       vmd -dispdev text -e DCD2PDBs.tcl 
#              
# DESCRIPTION
#        This VMD script convert a DCD trajectory into multiple PDB files,
#		 one for each frame which will be saved in a dedicated directory pdbfiles
#		 The name of the dcd and psf files (here protein.psf and protein.dcd 
#        need to be changed by the user as well as the total number of frames (here 2000)
#
#############################################################################

mol new protein.psf

  for {set i 0} {$i < 2000 } {incr i} {
     animate read dcd protein.dcd beg $i end $i
     set sel [atomselect top noh]
     $sel writepdb ./pdbfiles/$i.pdb

     $sel delete
     animate delete all
  }

