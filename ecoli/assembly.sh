#!/bin/bash

/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools/canu-1.9/Linux-amd64/bin/canu \
	-p ${4}_${1} -d ${4}_${1} \
	genomeSize=${2} \
	corErrorRate=${3} \
	-pacbio-raw ${5} \
	useGrid=${6} \
	corMhapMemory=1 obtOvlMemory=1 utgOvlMemory=1 oeaMemory=1 corPartitions=200 corPartitionMin=1 cnsPartitions=200 cnsPartitionMin=1 batThreads=8 \
	gridOptions=" --partition=vgl,vgl_bigmem,hpc --time=00:00:00 --nice=0"
