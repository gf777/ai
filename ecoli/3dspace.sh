#!/bin/bash

genomeSizes=3405000
corErrorRates=0.258

mkdir logs QV

if ! [[ -e raw_data.meryl ]]; then

echo ${2} > input.fofn

	$tools/meryl/scripts/_submit_meryl2_build.sh 18 input.fofn raw_data

fi

ln -s ../raw_data.meryl QV/raw_data.meryl

if  ! [[ -e report.log ]]; then

	printf "genomeSize\terr_rate\tQV\tgsize\tN50\tL50\n" > report.log

fi

until [ -d "raw_data.meryl" ] ; do sleep 60; done

for genomeSize in "${genomeSizes[@]}" 
do 
    for corErrorRate in "${corErrorRates[@]}"
    do
        combinations+=("${genomeSize}-${corErrorRate}") 
    done 
done

for ((idx=0; idx<${#combinations[@]}; ++idx));
do

	IFS=- read genomeSize corErrorRate <<< "${combinations[${idx}]}"

	if ! [[ -e "${1}_${combinations[${idx}]}" ]]; then

		printf "sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=16 -o logs/asm_${1}_${combinations[${idx}]}.log assembly.sh ${combinations[${idx}]} ${genomeSize} ${corErrorRate} ${1} ${2}\n"
		sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=16 -o logs/asm_${1}_${combinations[${idx}]}.log ../assembly.sh ${combinations[${idx}]} ${genomeSize} ${corErrorRate} ${1} ${2} false | awk '{print $4}' > job.jid
	
	fi

	if ! [[ -e "QV/${1}_${combinations[${idx}]}.qv" ]]; then
	
		WAIT="afterok:"
		jid=`cat job.jid`
		WAIT=$WAIT$jid
	
		printf "\nsbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=16 -o logs/qv_${1}_${combinations[${idx}]}.log --dependency=$WAIT qv.sh ${combinations[${idx}]} ${1} ${2}\n"
		sbatch --partition=vgl,vgl_bigmem,hpc --cpus-per-task=16 -o logs/qv_${1}_${combinations[${idx}]}.log --dependency=$WAIT ../qv_3dspace.sh ${combinations[${idx}]} ${1} ${2}
	
	fi
	
	echo ""

done

