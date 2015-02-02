#!/bin/bash

# example:

# show help:
#  bash runCE.sh -h 

# show alignment GUI (and download the PDB files automatically if they don't exist)
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -autoFetch -show3d 

# print output as XML
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printXML 

# print output in CE style
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printCE

# print output in FatCat style
#  bash runCE.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printFatCat

# load files from a URL. Note: aligns the whole file, i.e. all chains. 
# bash runCE.sh -file1 ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/cd/pdb1cdg.ent.gz -file2 ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ti/pdb1tim.ent.gz -printCE

# load files from local file system. Note: aligned the whole file, i.e. all chains. If you want to break this up into regions, you need to manipulate the files first manually.
# bash runCE.sh -file1 /tmp/cd/pdb1cdg.ent.gz -file2  file:///tmp/ti/pdb1tim.ent.gz -printCE

# print more information about the aligned fragment pairs (before optimization)
# bash runCE.sh -file1 /tmp/cd/pdb1cdg.ent.gz -file2  file:///tmp/ti/pdb1tim.ent.gz -printCE -showAFPRanges

# run a DB search
# bash runCE.sh -pdbFilePath /tmp/ -alignPairs ./example.lst -outFile db.out

# view DB search results
# bash runCE.sh -pdbFilePath /tmp/ -showDBresult db.out


### Execute jar ###

# Get the base directory of the argument.
# Can resolve single symlinks if readlink is installed
function scriptdir {
    cd "$(dirname "$1")"
    cd "$(dirname "$(readlink "$1" 2>/dev/null || basename "$1" )")"
    pwd
}
DIR="$(scriptdir "$0" )"
# send the arguments to the java app
java -Xmx500M -cp "$DIR/${project.build.finalName}.jar" org.biojava.nbio.structure.align.ce.CeMain "$@"
