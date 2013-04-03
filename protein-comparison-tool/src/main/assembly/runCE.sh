#!/bin/bash

# example:

# show help:
#  bash runCE.sh -h 

# show alignent GUI (and download the PDB files automatically if they don't exist)
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -autoFetch -show3d 

# print output as XML
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printXML 

# print output in CE style
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printCE

# print output in FatCat style
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printFatCat

# load files from a URL. Note: alignes the whole file, i.e. all chains. 
# bash runCE.sh -file1 ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/cd/pdb1cdg.ent.gz -file2 ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ti/pdb1tim.ent.gz -printCE

# load files from local file system. Note: aligned the whole file, i.e. all chains. If you want to break this up into regions, you need to manipulate the files first manually.
# bash runCE.sh -file1 /tmp/cd/pdb1cdg.ent.gz -file2  file:///tmp/ti/pdb1tim.ent.gz -printCE

# print more information about the aligned fragment pairs (before optimization)
# bash runCE.sh -file1 /tmp/cd/pdb1cdg.ent.gz -file2  file:///tmp/ti/pdb1tim.ent.gz -printCE -showAFPRanges

# run a DB search
# bash runCE.sh -pdbFilePath /tmp/ -alignPairs ./example.lst -outFile db.out

# view DB search results
# bash runCE.sh -pdbFilePath /tmp/ -showDBresult db.out


# send the arguments to the java app
# allows to specify a different config file
args="$*"

java -Xmx500M -cp "$PWD/jars/*" org.biojava.bio.structure.align.ce.CeMain $args

# To run CE-CP (detection of circular permutations:
#java -Xmx500M -cp "$PWD/jars/*" org.biojava.bio.structure.align.ce.CeCPMain $args
# test with something like this:
# bash runCE.sh -pdb1 d1qdma1 -pdb2 d1nkla_ -pdbFilePath /tmp/ -autoFetch -show3d