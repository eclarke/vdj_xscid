#!/usr/bin/env bash

##
# Usage: blastn_parallel.sh <fasta file>
# Returns: <original file name>.blast (hopefully in same dir as original file)

# Splits a large fasta file into parts and runs blast on them. Specify the 
# number of parts in the $CHUNKS variable below.
#
# Faster than just feeding the whole file directly to blast;
# also preserves some options without having to type them. 
#
# By Erik Clarke (ecl@mail.med.upenn.edu)
##

set -o errexit

CHUNKS=16
DB="$HOME/vdj_xscid/unfiltered.all.db"
BLAST="blastn -db ${DB} -outfmt 7 -query "

fname=`basename $1`
dname=`dirname $1`
starting_dir=`pwd`
echo "Working on file ${dname}/${fname} in ${starting_dir}."
echo "Copying file to /tmp..."
cp $1 /tmp/
cd /tmp/
echo "Splitting file into ${CHUNKS} pieces..."
split -n l/${CHUNKS} /tmp/${fname} ${fname}.part.
for f in ${fname}.part.*; do
	echo "Blasting ${f}..."
	eval $blast $f -out ${f}.out &
	pidlist="$pidlist $!"
done
echo "Working... this may take a while..."
for job in $pidlist; do
	wait $job
done

echo "Done blasting, concatenating and moving back to ${dname}"
cat ${fname}.part.*.out > ${fname}.blast
cd $starting_dir
mv /tmp/${fname}.blast ${dname}
echo "File outputted to ${dname}/${fname}.blast\n"
