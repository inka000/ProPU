#!/bin/bash

#####Repositories creation#####
rm -Rf Query DSSP suivi
mkdir Query DSSP suivi

#####PDB queries will be relocalized in Query repository#####

for query in "$1"* ; do
    cp $query ./Query
done

######Secondary structure assignation with DSSP#####
for f in Query/* ; do
	pdb=`echo $f | cut -d/ -f2`
	name=`echo ${f##*/} | cut -d. -f1`
	./dssp-2.0.4-linux-amd64 -i ./$f -o DSSP/$name.out &>> suivi/tmp.log
	#python3 ./parse_dssp.py DSSP/$name.out
done

