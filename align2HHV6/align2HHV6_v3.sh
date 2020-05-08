#!/bin/bash

# Changing Log
# Sep 1, 2019: counting variant number and ignore annotation lines '^#' 
# Sep 4, 2019: add RG infomarion for BWA MEM
# Sep 9, 2019: add patch for ICGC bam files 
# Oct 8, 2019: add align2hhv6_unmap_fq

SAMTOOLS=/my/miniconda2/bin/samtools
BWA=/my/miniconda2/bin/bwa
BCFTOOLS=/my/miniconda2/bin/bcftools
OUTDIR=HHV6

if [ ! -d ${OUTDIR} ];
	then
	mkdir ${OUTDIR}
fi 

# SubtypeName Length DR_sta DR_end U_sta U_end REF_Fasta_path 
readonly HHV6_subtypes=(       
                            'U1102|151288|1|8089|8090|151288|HHV6REF/HHV6A_U1102_wo_DRR.fasta'
							'Z29|152930|1|8793|8911|152930|HHV6REF/HHV6B_z29_wo_DRR.fasta'						
)

###############################################
# Input: Bam file and desired sample name
# Process:  extract unmapped reads as fastq files 
#           and do alignment with ${BWA} mem for 6A and 6B 
# Usage example : align2hhv6_map NA18999.bam NA18999
###############################################
align2hhv6_bam () {
	bam=$1
	sample=$2
	if [ ! -f ${bam}.bai ]; 
	then
	${SAMTOOLS} index ${bam}
	fi

	${SAMTOOLS} fastq -f 12 -1 ${sample}.unmap.R1.fq -2 ${sample}.unmap.R2.fq ${bam}
	# ${SAMTOOLS} fastq -1 ${sample}.unmap.R1.fq -2 ${sample}.unmap.R2.fq ${bam}
	local subtype length dr_sta dr_end U_sta U_end ref_fasta n_var
	for HHV6 in ${HHV6_subtypes[@]}
	do
            IFS=$'|' read -r subtype length dr_sta dr_end u_sta u_end ref_fasta <<< "$HHV6"
	        ${BWA} mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}_library_1" ${ref_fasta} ${sample}.unmap.R1.fq ${sample}.unmap.R2.fq | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam

			# echo "${BWA} mem ${ref_fasta} ${sample}.unmap.R1.fq ${sample}.unmap.R2.fq | ${SAMTOOLS} view -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam"
            
            local cov depth cov_dr dp_dr cov_u dp_u
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'"
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1  -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"
            # echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"

        	IFS=' ' read -r cov depth <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'`
            IFS=' ' read -r cov_dr dp_dr <<<  `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1 -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			IFS=' ' read -r cov_u dp_u <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			
			n_var=`${BCFTOOLS} mpileup -Ou -f ${ref_fasta} ${OUTDIR}/${sample}.${subtype}.bam | ${BCFTOOLS} call -mv -Ov --ploidy 1 | grep -v '^#' | wc -l`
            # echo ${subtype} ${cov} ${depth} 
			echo "subtype=${subtype} cov=${cov} dp=${depth} cov_dr=${cov_dr} dp_dr=${dp_dr} cov_u=${cov_u} dp_u=${dp_u} n_var=${n_var} file=${OUTDIR}/${sample}.${subtype}.bam"
	done
}

###############################################
# Input: unmap.bam file and desired sample name
# Process:  extract unmapped reads to aq fastq files 
#           and do alignment with BWA mem for 6A and 6B 
# Usage example 
###############################################
align2hhv6_unmap_bam () {
	unmapBAM=$1
	sample=$2	
	${SAMTOOLS} view ${unmapBAM} | awk '{printf "@%s\n%s\n+\n%s\n",$1,$10,$11}' > ${OUTDIR}/${sample}.unmap.fq

	local subtype length dr_sta dr_end U_sta U_end ref_fasta n_var
	for HHV6 in ${HHV6_subtypes[@]}
	do
	        
            IFS=$'|' read -r subtype length dr_sta dr_end u_sta u_end ref_fasta <<< "$HHV6"
	        ${BWA} mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}_library_1" ${ref_fasta} -M ${OUTDIR}/${sample}.unmap.fq | ${SAMTOOLS} view -hu -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam

			# echo "${BWA} mem -M ${ref_fasta} ${OUTDIR}/${sample}.unmap.fq | ${SAMTOOLS} view -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam"
            
            local cov depth cov_dr dp_dr cov_u dp_u
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'"
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1  -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"
            # echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"


        	IFS=' ' read -r cov depth <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'`
            IFS=' ' read -r cov_dr dp_dr <<<  `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1 -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			IFS=' ' read -r cov_u dp_u <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`

			n_var=`${BCFTOOLS} mpileup -Ou -f ${ref_fasta} ${OUTDIR}/${sample}.${subtype}.bam | ${BCFTOOLS} call -mv -Ov --ploidy 1 | grep -v '^#' | wc -l`

            # echo ${subtype} ${cov} ${depth} 
			echo "subtype=${subtype} cov=${cov} dp=${depth} cov_dr=${cov_dr} dp_dr=${dp_dr} cov_u=${cov_u} dp_u=${dp_u} n_var=${n_var} file=${OUTDIR}/${sample}.${subtype}.bam"
	done	

}

###############################################
# Input: Pair-end fastq files (R1 and R2) and sample name
# Process: alignment with BWA mem for 6A and 6B 
# Usage example:
# align2hhv6_unmap_fq_pe S1_R1.fq S1_R2.fq S1
###############################################
align2hhv6_unmap_fq_pe () {
	unmap_fq_r1=$1
	unmap_fq_r2=$2
	sample=$3
	local subtype length dr_sta dr_end U_sta U_end ref_fasta n_var
	for HHV6 in ${HHV6_subtypes[@]}
	do
	        
            IFS=$'|' read -r subtype length dr_sta dr_end u_sta u_end ref_fasta <<< "$HHV6"
	        ${BWA} mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}_library_1" ${ref_fasta} ${unmap_fq_r1} ${unmap_fq_r2} | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam

			# echo "${BWA} mem -M ${ref_fasta} ${sample}.unmap.fq | ${SAMTOOLS} view -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam"
            
            local cov depth cov_dr dp_dr cov_u dp_u
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'"
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1  -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"
            # echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"


        	IFS=' ' read -r cov depth <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'`
            IFS=' ' read -r cov_dr dp_dr <<<  `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1 -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			IFS=' ' read -r cov_u dp_u <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`

            n_var=`${BCFTOOLS} mpileup -Ou -f ${ref_fasta} ${OUTDIR}/${sample}.${subtype}.bam | ${BCFTOOLS} call -mv -Ov --ploidy 1 | grep -v '^#' | wc -l`
            # echo ${subtype} ${cov} ${depth} 
			echo "subtype=${subtype} cov=${cov} dp=${depth} cov_dr=${cov_dr} dp_dr=${dp_dr} cov_u=${cov_u} dp_u=${dp_u} n_var=${n_var} file=${OUTDIR}/${sample}.${subtype}.bam"
	done	
}

align2hhv6_unmap_fq () {
	unmap_fq=$1
	sample=$2
	local subtype length dr_sta dr_end U_sta U_end ref_fasta n_var
	for HHV6 in ${HHV6_subtypes[@]}
	do
	        
            IFS=$'|' read -r subtype length dr_sta dr_end u_sta u_end ref_fasta <<< "$HHV6"
	        ${BWA} mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}_library_1" ${ref_fasta} -M ${unmap_fq} | ${SAMTOOLS} view -hu -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam

			# echo "${BWA} mem -M ${ref_fasta} ${OUTDIR}/${sample}.unmap.fq | ${SAMTOOLS} view -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam"
            
            local cov depth cov_dr dp_dr cov_u dp_u
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'"
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1  -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"
            # echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"


        	IFS=' ' read -r cov depth <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'`
            IFS=' ' read -r cov_dr dp_dr <<<  `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1 -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			IFS=' ' read -r cov_u dp_u <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`

			n_var=`${BCFTOOLS} mpileup -Ou -f ${ref_fasta} ${OUTDIR}/${sample}.${subtype}.bam | ${BCFTOOLS} call -mv -Ov --ploidy 1 | grep -v '^#' | wc -l`

            # echo ${subtype} ${cov} ${depth} 
			echo "subtype=${subtype} cov=${cov} dp=${depth} cov_dr=${cov_dr} dp_dr=${dp_dr} cov_u=${cov_u} dp_u=${dp_u} n_var=${n_var} file=${OUTDIR}/${sample}.${subtype}.bam"
	done	

}

###############################################
# Input: Unmap bam from ICGC
# Process: alignment with BWA mem for 6A and 6B 
# Usage example:
# align2hhv6_unmap_bam_icgc S1.unmap.bam S1
###############################################

align2hhv6_unmap_bam_icgc () {
	unmapBAM=$1
	sample=$2	
	#${SAMTOOLS} view ${unmapBAM} | awk '{printf "@%s\n%s\n+\n%s\n",$1,$10,$11}' > ${OUTDIR}/${sample}.unmap.fq
	${SAMTOOLS} view ${unmapBAM} | awk '{split($14, a, ":"); if( and($2,0x4) || and($2,0x8) || (a[1] =="AS" && a[3] < 30) ) printf ">%s\n%s\n", NR,$10 }' > ${OUTDIR}/${sample}.unmap.fq
	local subtype length dr_sta dr_end U_sta U_end ref_fasta n_var
	for HHV6 in ${HHV6_subtypes[@]}
	do
	        
            IFS=$'|' read -r subtype length dr_sta dr_end u_sta u_end ref_fasta <<< "$HHV6"
	        ${BWA} mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}_library_1" ${ref_fasta} -M ${OUTDIR}/${sample}.unmap.fq | ${SAMTOOLS} view -hu -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam

			# echo "${BWA} mem -M ${ref_fasta} ${OUTDIR}/${sample}.unmap.fq | ${SAMTOOLS} view -F 4 | ${SAMTOOLS} sort -o ${OUTDIR}/${sample}.${subtype}.bam; ${SAMTOOLS} index ${OUTDIR}/${sample}.${subtype}.bam"
            
            local cov depth cov_dr dp_dr cov_u dp_u
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'"
			# echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1  -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"
            # echo "${SAMTOOLS} depth ${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}} END{print c/LEN, dp/LEN}'"


        	IFS=' ' read -r cov depth <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${length} '{c++; dp+=\$3}END{print c/LEN, dp/LEN}'`
            IFS=' ' read -r cov_dr dp_dr <<<  `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${dr_end}-${dr_sta}+1 -v S=${dr_sta} -v E=${dr_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`
			IFS=' ' read -r cov_u dp_u <<< `${SAMTOOLS} depth ${OUTDIR}/${sample}.${subtype}.bam | awk -v LEN=${u_end}-${u_sta}+1 -v S=${u_sta} -v E=${u_end} '{if(\$2>=S && \$2<=E){c++; dp+=\$3}}; END{print c/LEN, dp/LEN}'`

			n_var=`${BCFTOOLS} mpileup -Ou -f ${ref_fasta} ${OUTDIR}/${sample}.${subtype}.bam | ${BCFTOOLS} call -mv -Ov --ploidy 1 | grep -v '^#' | wc -l`

            # echo ${subtype} ${cov} ${depth} 
			echo "subtype=${subtype} cov=${cov} dp=${depth} cov_dr=${cov_dr} dp_dr=${dp_dr} cov_u=${cov_u} dp_u=${dp_u} n_var=${n_var} file=${OUTDIR}/${sample}.${subtype}.bam"
	done	

}





