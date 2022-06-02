# covid_finder
COVID variants identification in sequence data
# process clinical variants
* SNPs calling\
`python clinical_variant.py`\
`sh clinical_align.sh`\
`python clinical_sum.py`
# process metagenomes
* SNPs calling\
`python SNPcalling.py -i /scratch/users/fuqingwu/covid19/fastq_210521/adaptcutted/ -fq _1.fastq.ca -ref /scratch/users/anniz44/scripts/covid/trial/reference.fasta -s /scratch/users/anniz44/scripts/covid/ -o /scratch/users/anniz44/genomes/covid/ -t 40 `\
* SNP filtering\
`python SNPfilter.py -i /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -s /scratch/users/anniz44/scripts/covid/`
* Coverage evaluation of each sample\
`python primer_cov.py -o /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -primer /scratch/users/anniz44/genomes/covid/SNPcalling/primer_info.txt`
* Sequencing error evaluation of each sample\
`python sequencing_error.py -o /scratch/users/anniz44/genomes/covid/SNPcalling/merge/`
* SNP sum\
'jobmit /scratch/users/anniz44/scripts/covid/state_sum.sh'
`#python SNPsum.py -clinical /scratch/users/anniz44/genomes/covid/clinical/allearlymut.txt -i /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -ref /scratch/users/anniz44/scripts/covid/trial/reference.fasta -cov /scratch/users/anniz44/genomes/covid/SNPcalling/merge/all.primer.coverage.txt`
`#python SNPsum.py -state TX -clinical /scratch/users/anniz44/genomes/covid/clinical/allearlymutTX.txt -i /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -ref /scratch/users/anniz44/scripts/covid/trial/reference.fasta -cov /scratch/users/anniz44/genomes/covid/SNPcalling/merge/all.primer.coverage.txt`
`#python SNPsum.py -state MA -clinical /scratch/users/anniz44/genomes/covid/clinical/allearlymutMA.txt -i /scratch/users/anniz44/genomes/covid/SNPcalling/merge/ -ref /scratch/users/anniz44/scripts/covid/trial/reference.fasta -cov /scratch/users/anniz44/genomes/covid/SNPcalling/merge/all.primer.coverage.txt`
`python SNPsumallstate.py`
`python dnds.py`
* download /scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/all* and /scratch/users/anniz44/genomes/covid/SNPcalling/all*
* notebook Covid_quality Covid-state