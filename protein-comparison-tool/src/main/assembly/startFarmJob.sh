#!/bin/bash

# example:
# startFarmJob.sh -pdbFilePath /Users/ap3/WORK/PDB -nrAlignments 10

# send the arguments to the java app
# allows to specify a different config file
args="$*"

if [ -f $OSG_APP/engage/jdk1.6.0_16/bin/java ]; then
    $OSG_APP/engage/jdk1.6.0_16/bin/java -Xmx1G -cp "$PWD/jars/*" org.biojava.bio.structure.align.FarmJob $args
else
    if [ -f /osg/osg-app/engage/jdk1.6.0_03/bin/java ]; then
        /osg/osg-app/engage/jdk1.6.0_03/bin/java  -Xmx1G -cp "$PWD/jars/*" org.biojava.bio.structure.align.FarmJob $args
     else
        which java
        java -version
        java  -Xmx1G -cp "$PWD/jars/*" org.biojava.bio.structure.align.FarmJob $args
    fi
fi

exit $?