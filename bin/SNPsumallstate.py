# summary state SNPs into one file
import glob
import os
import math
alldepthfiles = glob.glob('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/*.freqsum.txt')
allsumfile = []
allsumfile.append('State\tColor\tClinical\tClinical_time\tSampling_time\tGoodSNP\tPOS\tREF\tALT\tCount\tCount_ALTfreq_cutoff\tPrevalence\tPrevalence_ALTfreq_cutoff\tClinical_prevalence\tAvg_ALT_freq\tGene\tGene_POS\tN_S\tAAchange\n')
allvariant = dict()
for depthfile in alldepthfiles:
    state = os.path.split(depthfile)[-1].split('.')[0]
    for lines in open(depthfile,'r'):
        if not lines.startswith('Color'):
            lines_set = lines.split('\t')[:18]
            allsumfile.append('%s\t%s\n'%(state,'\t'.join(lines_set)))
            variant = '%s\t%s'%(lines_set[5],lines_set[7])
            allvariant[variant] = [allvariant.get(variant,[0,0])[0] + int(lines_set[8]),
                                   allvariant.get(variant,[0,0])[1] + int(lines_set[9])]

allvariantsum = []
allvariantsum.append('POS\tALT\tPrevalence\tPrevalence_ALTfreq_cutoff\n')
for variant in allvariant:
    allvariantsum.append('%s\t%s\t%s\n'%(variant,allvariant[variant][0],allvariant[variant][1]))

f1 = open('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/all.mutsum.txt','w')
f1.write(''.join(allsumfile))
f1.close()

f1 = open('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/all.mutsumpre.txt','w')
f1.write(''.join(allvariantsum))
f1.close()
