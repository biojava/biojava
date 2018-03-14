#!/bin/bash
# Adds the BioJava LGPL license statement to the top of every java file

BASEDIR=$(dirname "$0")

find . -iname '*.java' -exec grep -L 'http://www.gnu.org/copyleft/lesser.html' '{}' '+'|
xargs grep -Li 'copyright' |
while read file; do
    echo "$file"
    cat $BASEDIR/../HEADER.txt > tmp.java
    echo >> tmp.java
    cat "$file" >> tmp.java
    mv tmp.java "$file"
done
