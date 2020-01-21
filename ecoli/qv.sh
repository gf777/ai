#!/bin/bash

cd QV

ln -s ../${2}_${1}/${2}_${1}.contigs.fasta ${2}_${1}.contigs.fasta

../../spectra-cn.sh raw_data.meryl 18 ${2}_${1}.contigs.fasta ${2}_${1}
