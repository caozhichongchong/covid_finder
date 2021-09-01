# covid_finder
COVID variants identification in sequence data
# process
* SNPs calling\
`python SNPcalling.py -i /scratch/users/fuqingwu/covid19/fastq_210521/adaptcutted/ -fq _1.fastq.ca -ref /scratch/users/fuqingwu/covid19/fastq_210521/SARS2_MN908947.3.fasta -s /scratch/users/anniz44/scripts/covid/ -o /scratch/users/anniz44/genomes/covid/ -t 40 `
* SNP filtering
`python SNPfilter.py -i /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -vcf .raw.vcf -fq _1.fastq.ca -s /scratch/users/anniz44/scripts/covid/`
