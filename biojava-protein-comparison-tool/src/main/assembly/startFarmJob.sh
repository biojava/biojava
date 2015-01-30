#!/bin/bash

# example:
# startFarmJob.sh -pdbFilePath /Users/ap3/WORK/PDB -nrAlignments 10

# send the arguments to the java app
# allows to specify a different config file

if [ -f $OSG_APP/engage/jdk1.6.0_16/bin/java ]; then
    $OSG_APP/engage/jdk1.6.0_16/bin/java -Xmx1G -cp "$PWD/${project.build.finalName}.jar" org.biojava.nbio.structure.align.FarmJob "$@"
else
    if [ -f /osg/osg-app/engage/jdk1.6.0_03/bin/java ]; then
        /osg/osg-app/engage/jdk1.6.0_03/bin/java  -Xmx1G -cp "$PWD/${project.build.finalName}.jar" org.biojava.nbio.structure.align.FarmJob "$@"
    else
        which java
        java -version
        java  -Xmx1G -cp "$PWD/${project.build.finalName}.jar" org.biojava.nbio.structure.align.FarmJob "$@"
    fi
fi

exit $?

