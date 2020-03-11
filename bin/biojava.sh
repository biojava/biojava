#!/usr/bin/env bash

# this runs any main function in the source tree.
# to change the class to be run
# EXECMAIN=org.biojava.not.demo.main bin/biojava.sh
mkdir -p target/lib/                                          && {
ln -s 2>/dev/null    $PWD/*/target{,/lib}/*jar  target/lib
cp 2>/dev/null  -as  $PWD/*/src/main/resources/* target/
}

set -fx
JDIR=$PWD/$(dirname $0)/../
pushd target
exec java  -classpath "$JDIR/target/*:$JDIR/target/lib/*"  ${EXECMAIN:=demo.DemoSixFrameTranslation} "$@"
popd