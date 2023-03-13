#!/bin/bash
threads=20
mem=100
wd=~/storage1/Ly49/vivier_2022/
ref=~/genomes/mm10/cellranger/mm10_gencode_klraps/
#######
# custom reference is offered at: 
# https://wangftp.wustl.edu/~cfan/Ly49pub/large_files/mm10_gencode_klraps
# it's not on github due to size
#######
feature_ref=~/storage1/Ly49/vivier_2022/feature_ref.csv
mkdir -p ${wd}/count/

# samples=(Gut_rep1 Gut_rep2)
samples=(`ls $wd/fastq_ori/*fastq.gz | grep -v 'Gut_rep[12]' | sed 's;.*fastq_ori/;;g' | sed 's/_ADT.*$//g' | sed 's/_scRNA.*//g' | uniq`)
echo ${samples[@]}
for sample in ${samples[@]}
do
	csv=${wd}/fastq/${sample}/${sample}.csv
	ls $csv
	LSF_DOCKER_PRESERVE_ENVIRONMENT=false LSF_DOCKER_VOLUMES='/storage1/fs1/hprc:/storage1/fs1/hprc' bsub \
	-q general -n $threads -G compute-hprc \
	-R "span[hosts=1] select[mem>${mem}G] rusage[mem=${mem}G]" -M ${mem}G \
	-a 'docker(cumulusprod/cellranger:6.0.0)' \
	-oo ${wd}/logs/log.${sample}.count.txt \
	/bin/bash -c "cd ${wd}/count && cellranger count --id ${sample} --transcriptome $ref --libraries $csv --feature-ref ${feature_ref} --localcores $((threads-1)) --localmem $((mem-10))"
done

