#!/bin/bash

asm=${1}
data=${2}
par1=${3}
par2=$(printf "%.15f" ${4})
pas=${5}

mkdir logs QV

wait_file() {
  local file="$1"; shift

  until [ -f $file ] ; do sleep 60; done
  
}

if ! [[ -e ../raw_data.meryl ]]; then

	echo ${data} > input.fofn

	$tools/meryl/scripts/_submit_meryl2_build.sh 18 input.fofn raw_data

fi

ln -s ../../raw_data.meryl QV/raw_data.meryl

rm -f report_ai.log

asm_qv() {

	
	genomeSize=${1}
	corErrorRate=${2}

	if ! [[ -e "${asm}_${genomeSize}_${corErrorRate}" ]]; then

		printf "sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=1 -o logs/asm_${asm}_${genomeSize}_${corErrorRate}.log assembly.sh ${genomeSize}_${corErrorRate} ${genomeSize} ${corErrorRate} ${asm} ${data}\n"
		sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=1 -o logs/asm_${asm}_${genomeSize}_${corErrorRate}.log ../assembly.sh ${genomeSize}_${corErrorRate} ${genomeSize} ${corErrorRate} ${asm} ${data} true | awk '{print $4}' > job.jid
	
	fi
	
	wait_file "${asm}_${genomeSize}_${corErrorRate}/${asm}_${genomeSize}_${corErrorRate}.contigs.fasta"

	if	! [[ -e "QV/${asm}_${genomeSize}_${corErrorRate}.qv" ]]; then

		until [ -d "QV/raw_data.meryl" ] ; do sleep 60; done
	
		printf "\nsbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=8 -o logs/qv_${asm}_${genomeSize}_${corErrorRate}.log --dependency=$WAIT qv.sh ${genomeSize}_${corErrorRate} ${asm} ${data}\n"
		sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=4 -o logs/qv_${asm}_${genomeSize}_${corErrorRate}.log ../qv.sh ${genomeSize}_${corErrorRate} ${asm} ${data}
		
	fi
	
	wait_file "QV/${asm}_${genomeSize}_${corErrorRate}.qv"

	gsize=$(perl ../countFasta.pl ${asm}_${genomeSize}_${corErrorRate}/${asm}_${genomeSize}_${corErrorRate}.contigs.fasta | sed -n 's/^.*Total length of sequence:\t//p' | sed 's/ bp//g')
	ctgn=$(perl ../countFasta.pl ${asm}_${genomeSize}_${corErrorRate}/${asm}_${genomeSize}_${corErrorRate}.contigs.fasta | sed -n 's/^.*Total number of sequences:\t//p' | sed 's/ bp//g')
	stats50=($(perl ../countFasta.pl ${asm}_${genomeSize}_${corErrorRate}/${asm}_${genomeSize}_${corErrorRate}.contigs.fasta | sed -n 's/^N50 stats.*in the //p' | grep -o "[0-9]*"))
	L50=${stats50[0]}
	N50=${stats50[1]}
	
}

asm_qv ${par1} ${par2}

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" ${genomeSize} ${corErrorRate} $(awk '{print $4}' QV/${asm}_${genomeSize}_${corErrorRate}.qv) ${gsize} ${N50} ${ctgn} ${L50} >> report_ai.log

for ((idx=0; idx<500; ++idx));
do

	iQV=45
	igsize=4641652
	ictgn=1

	gsize=$(perl ../countFasta.pl ${asm}_${par1}_${par2}/${asm}_${par1}_${par2}.contigs.fasta | sed -n 's/^.*Total length of sequence:\t//p' | sed 's/ bp//g')
	ctgn=$(perl ../countFasta.pl ${asm}_${par1}_${par2}/${asm}_${par1}_${par2}.contigs.fasta | sed -n 's/^.*Total number of sequences:\t//p' | sed 's/ bp//g')
	QV=$(awk '{print $4}' QV/${asm}_${par1}_${par2}.qv)
	
	printf "\n\nStarting position on x axis:%s\n" ${par1}
	printf "Starting position on y axis:%s\n\n" ${par2}
	printf "Starting QV:%s\n\n" ${QV}
	printf "Starting genome size:%s\n\n" ${gsize}
	printf "Starting contig N:%s\n\n" ${ctgn}
	
	dgsize=$(printf "%.0f" $(echo "scale=0; ( ${gsize} - ${igsize} )^2" | bc))
	dctgn=$(printf "%.0f" $(echo "scale=0; ( ${ctgn} - ${ictgn} )^2" | bc))
	dQV=$(printf "%.15f" $(echo "scale=15; ( ${QV} - ${iQV} )^2" | bc))
	


	new_par1=$(printf "%.0f" $(echo "scale=0; ${par1} + ${par1} * 0.01" | bc))
	new_par2=$(printf "%.15f" $(echo "scale=15; ${par2} + ${par2} * 0.01" | bc))
	
	printf "New position on x axis:%s\n" ${new_par1}
	printf "New position on y axis:%s\n\n" ${new_par2}

	asm_qv ${new_par1} ${par2} &
	asm_qv ${par1} ${new_par2} &
	
	wait
	
	gsize1=$(perl ../countFasta.pl ${asm}_${new_par1}_${par2}/${asm}_${new_par1}_${par2}.contigs.fasta | sed -n 's/^.*Total length of sequence:\t//p' | sed 's/ bp//g')
	gsize2=$(perl ../countFasta.pl ${asm}_${par1}_${new_par2}/${asm}_${par1}_${new_par2}.contigs.fasta | sed -n 's/^.*Total length of sequence:\t//p' | sed 's/ bp//g')

	printf "Genome size in x axis is:%s\n" ${gsize1}
	printf "Genome size in y axis is:%s\n\n" ${gsize2}

	dgsize1=$(printf "%.0f" $(echo "scale=0; ( ${gsize1} - ${igsize} )^2" | bc))
	dgsize2=$(printf "%.0f" $(echo "scale=0; ( ${gsize2} - ${igsize} )^2" | bc))
	
	printf "Cost in genome size on x axis is:%s\n" ${dgsize1}
	printf "Cost in genome size on y axis is:%s\n\n" ${dgsize2}



	ctgn1=$(perl ../countFasta.pl ${asm}_${new_par1}_${par2}/${asm}_${new_par1}_${par2}.contigs.fasta | sed -n 's/^.*Total number of sequences:\t//p' | sed 's/ bp//g')
	ctgn2=$(perl ../countFasta.pl ${asm}_${par1}_${new_par2}/${asm}_${par1}_${new_par2}.contigs.fasta | sed -n 's/^.*Total number of sequences:\t//p' | sed 's/ bp//g')

	printf "Contig N in x axis is:%s\n" ${ctgn1}
	printf "Contig N in y axis is:%s\n\n" ${ctgn2}

	dctgn1=$(printf "%.0f" $(echo "scale=0; ( ${ctgn1} - ${ictgn} )^2" | bc))
	dctgn2=$(printf "%.0f" $(echo "scale=0; ( ${ctgn2} - ${ictgn} )^2" | bc))
	
	printf "Cost in contig N on x axis is:%s\n" ${dctgn1}
	printf "Cost in contig N on y axis is:%s\n\n" ${dctgn2}



	QV1=$(awk '{print $4}' QV/${asm}_${new_par1}_${par2}.qv)
	QV2=$(awk '{print $4}' QV/${asm}_${par1}_${new_par2}.qv)

	printf "QV in x axis is:%s\n" ${QV1}
	printf "QV in y axis is:%s\n\n" ${QV2}

	dQV1=$(printf "%.15f" $(echo "scale=15; ( ${QV1} - ${iQV} )^2" | bc))
	dQV2=$(printf "%.15f" $(echo "scale=15; ( ${QV2} - ${iQV} )^2" | bc))
	
	printf "Cost in QV on x axis is:%s\n" ${dQV1}
	printf "Cost in QV on y axis is:%s\n\n" ${dQV2}



	echo "scale=15; ( (${dQV1} / ${iQV}) + (${dctgn1} / ${ictgn}) + (${dgsize1} / ${igsize}) ) - ( (${dQV} / ${iQV}) + (${dctgn} / ${ictgn}) + (${dgsize} / ${igsize}) ) / (${par1} * 0.01 / ${par1})"
	echo "scale=15; ( (${dQV2} / ${iQV}) + (${dctgn2} / ${ictgn}) + (${dgsize2} / ${igsize}) ) - ( (${dQV} / ${iQV}) + (${dctgn} / ${ictgn}) + (${dgsize} / ${igsize}) ) / (${par2} * 0.01 / ${par2})"

	g1=$(printf "%.15f" $(echo "scale=15; ( (${dQV1} / ${iQV}) + (${dctgn1} / ${ictgn}) + (${dgsize1} / ${igsize}) ) - ( (${dQV} / ${iQV}) + (${dctgn} / ${ictgn}) + (${dgsize} / ${igsize}) ) / (${par1} * 0.01 / ${par1})" | bc))
	g2=$(printf "%.15f" $(echo "scale=15; ( (${dQV2} / ${iQV}) + (${dctgn2} / ${ictgn}) + (${dgsize2} / ${igsize}) ) - ( (${dQV} / ${iQV}) + (${dctgn} / ${ictgn}) + (${dgsize} / ${igsize}) ) / (${par2} * 0.01 / ${par2})" | bc))

	printf "Slope on the x axis is:%s\n" ${g1}
	printf "Slope on the y axis is:%s\n\n" ${g2}

	par1=$(printf "%.0f" $(echo "scale=0; ${par1} - 200 * ${g1}" | bc))
	par2=$(printf "%.15f" $(echo "scale=15; ${par2} - 0.00001 * ${g2}" | bc))
	
	printf "New x is:%s\n" ${par1}
	printf "New y is:%s\n\n" ${par2}

	asm_qv ${par1} ${par2}
	
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" ${genomeSize} ${corErrorRate} $(awk '{print $4}' QV/${asm}_${genomeSize}_${corErrorRate}.qv) ${gsize} ${N50} ${ctgn} ${L50} >> report_ai.log	
	
done

