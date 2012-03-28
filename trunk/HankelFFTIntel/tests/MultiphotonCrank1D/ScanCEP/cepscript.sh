#!/bin/bash
export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
n=22
name=CEPMultiPhoton${n}
mkdir $name
cp multi $name
cd $name
time -p ./multi $n &