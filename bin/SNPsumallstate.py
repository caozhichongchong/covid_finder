# summary state SNPs into one file
import glob
import os
import math
alldepthfiles = glob.glob('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/*.freqsum.txt')
allsumfile = []
allsumfile.append('State\tColor\tClinical\tClinical_time\tSampling_time\tGoodSNP\tPOS\tREF\tALT\tCount\tPrevalence\tClinical_prevalence\tAvg_ALT_freq\tGene\tGene_POS\tN_S\tAAchange\n')
allvariant = dict()

try:
    os.mkdir('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/')
except IOError:
    pass

for depthfile in alldepthfiles:
    state = os.path.split(depthfile)[-1].split('.')[0]
    for lines in open(depthfile,'r'):
        if not lines.startswith('Color'):
            lines_set = lines.split('\t')[:16]
            allsumfile.append('%s\t%s\n'%(state,'\t'.join(lines_set)))
            variant = '%s\t%s'%(lines_set[5],lines_set[7])
            allvariant[variant] = allvariant.get(variant,0) + int(lines_set[8])

allvariantsum = []
allvariantsum.append('POS\tALT\tPrevalence\n')
for variant in allvariant:
    allvariantsum.append('%s\t%s\n'%(variant,allvariant[variant]))

f1 = open('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/all.mutsum.txt','w')
f1.write(''.join(allsumfile))
f1.close()

f1 = open('/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/all.mutsumpre.txt','w')
f1.write(''.join(allvariantsum))
f1.close()

os.system('mv /scratch/users/anniz44/genomes/covid/SNPcalling/merge/all* /scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/')