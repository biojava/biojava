#!/bin/bash

# example:

# show help:
#  bash runCECP.sh -h 

# show alignment GUI (and download the PDB files automatically if they don't exist)
#  bash runCECP.sh -pdb1 3cna.A -pdb2 2pel.A -pdbFilePath /tmp/ -autoFetch -show3d 

# print output as XML
#  bash runCECP.sh -pdb1 3cna.A -pdb2 2pel.A -pdbFilePath /tmp/ -printXML 

# print output in CE style
#  bash runCECP.sh -pdb1 3cna.A -pdb2 2pel.A -pdbFilePath /tmp/ -printCE

# print output in FatCat style
#  bash runCECP.sh -pdb1 3cna.A -pdb2 2pel.A -pdbFilePath /tmp/ -printFatCat

# load files from local file system. Note: aligned the whole file, i.e. all chains. If you want to break this up into regions, you need to manipulate the files first manually.
# bash runCECP.sh -file1 /tmp/cn/pdb3cna.ent.gz -file2  file:///tmp/pe/pdb2pel.ent.gz -show3d

# send the arguments to the java app
# allows to specify a different config file
args="$*"

java -Xmx500M -cp "$PWD/jars/*" org.biojava.bio.structure.align.ce.CeCPMain $args
