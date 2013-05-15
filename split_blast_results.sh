#!/usr/bin/env bash
nlines=84153931
chunk=$(( $nlines/16 ))
line=1
for i in $(seq 1 $chunk $nlines); do
    oldline=$(( $line+1 ))
    accn=`sed "${i}q;d" unfiltered.all.blast.tdf | cut -f1`
    line=`grep -n "^${accn}" unfiltered.all.blast.tdf | head -1 | cut -d: -f1`
    echo "${accn}    ${oldline}:${line}"
    sed -n "${oldline},${line}w chunk.${line}.out" unfiltered.all.blast.tdf
done