#!/bin/bash
#### tests out main methods 
#### run this script from within main-runner directory
CLASSPATH_FILE=./main-runner/cp.txt
mvn clean compile -f ../pom.xml

## generate classpath from maven dependencies
mvn dependency:build-classpath -Dmdep.outputFile=$CLASSPATH_FILE -f ../pom.xml

## add compiled classes to classpath
CLASSPATH="../target/classes:"$(<cp.txt)

FAILED_COUNT=0
## read and execute java command lines, check non-zero exit code
while  read -r line; do
   java -cp $CLASSPATH  $line
   if [[ $? != 0 ]]; then
	echo "$line failed"
	FAILED_COUNT=$(($FAILED_COUNT + 1 ))
   fi
done < invocations.txt
echo "Total failures = $FAILED_COUNT"
