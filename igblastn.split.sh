#!/usr/bin/env bash

##
# Splits a large fasta file into parts and runs igblastn on them
# simultaneously, with no concern for the number of cores on your machine
# or how choked up everything will get spawning 10-15 heavy jobs :)
#
# Marginally faster than just feeding the whole file directly to igblastn;
# also preserves some options without having to type them. 
# Very hackish.
#
# Usage: igblastn.split.sh <fasta file>
# Returns: <original file name>.blast (hopefully in same dir as original file)
#
# Erik Clarke (ecl@mail.med.upenn.edu)
##

set -o errexit

igdir="/usr/local/ncbi/igblast"
igdb="${igdir}/database"
igblast="igblastn -germline_db_V ${igdb}/human_gl_V -germline_db_D ${igdb}/human_gl_D -germline_db_J ${igdb}/human_gl_J -organism human -domain_system kabat -auxiliary_data ${igdir}/data/optional_file/human_gl.aux -outfmt 7 -query "

fname=`basename $1`
dname=`dirname $1`
starting_dir=`pwd`
echo "Working on file ${dname}/${fname} in ${starting_dir}."
echo "Copying file to /tmp..."
cp $1 /tmp/
cd /tmp/
echo "Splitting file into 50,000-line pieces..."
split -l 50000 -a 2 /tmp/${fname} ${fname}.part.
for f in ${fname}.part.*; do
	echo "Igblasting ${f}..."
	eval $igblast $f -out ${f}.out &
	pidlist="$pidlist $!"
done
for job in $pidlist; do
	echo $job
	wait $job
done

echo "Done igblasting, concatenating and moving back to ${dname}"
cat ${fname}.part.*.out > ${fname}.blast
cd $starting_dir
mv /tmp/${fname}.blast ${dname}
