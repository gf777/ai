#!/bin/bash

cd QV

ln -s ../${2}_${1}/${2}_${1}.contigs.fasta ${2}_${1}.contigs.fasta

../../spectra-cn.sh raw_data.meryl 18 ${2}_${1}.contigs.fasta ${2}_${1}

until [ -f "${2}_${1}.qv" ] ; do sleep 60; done

gsize=$(perl ../../countFasta.pl ../${2}_${1}/${2}_${1}.contigs.fasta | sed -n 's/^.*Total length of sequence:\t//p' | sed 's/ bp//g')
ctgn=$(perl ../../countFasta.pl ../${2}_${1}/${2}_${1}.contigs.fasta | sed -n 's/^.*Total number of sequences:\t//p' | sed 's/ bp//g')

stats50=($(perl ../../countFasta.pl ../${2}_${1}/${2}_${1}.contigs.fasta | sed -n 's/^N50 stats.*in the //p' | grep -o "[0-9]*"))
L50=${stats50[0]}
N50=${stats50[1]}


printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${1%-*}" "${1#*-}" $(awk '{print $4}' ${2}_${1}.qv) ${gsize} ${N50} ${ctgn} ${L50} >> ../report.log
