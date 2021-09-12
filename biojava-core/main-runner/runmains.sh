#!/bin/bash
CLASSPATH_FILE=./main-runner/cp.txt
#mvn clean compile -f ../pom.xml
#mvn dependency:build-classpath -Dmdep.outputFile=$CLASSPATH_FILE -f ../pom.xml

CLASSPATH="../target/classes:"$(<cp.txt)
echo $CLASSPATH
FAILED_COUNT=0
while  read -r line; do
   java -cp $CLASSPATH  $line
   if [[ $? != 0 ]]; then
	echo "$line failed"
	FAILED_COUNT=$(($FAILED_COUNT + 1 ))
   fi
done < invocations.txt
echo "Total failures = $FAILED_COUNT"
