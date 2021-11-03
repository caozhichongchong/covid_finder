################################################### END ########################################################
################################################### SET PATH ########################################################
# After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-o",
                      help="path to all freq files",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/merge/freqfiles/',
                      metavar='merge/')
required.add_argument("-primer",
                      help="ref file dir",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/merge/primer_info.txt',
                      metavar='primer_info.txt')

################################################## Definition ########################################################
args = parser.parse_args()
allfreqfiles = glob.glob('%s/*.allfreq.txt'%(args.o))
primerfile = args.primer
total_primer = 97
depth_cutoff = 50
################################################### Function ########################################################
def load_primer(primerfile):
    allprimer = dict()
    previous_end = 31
    for lines in open(primerfile,'r'):
        if not lines.startswith('start'):
            start,end,primerID,direction = lines.split('\n')[0].split('\t')
            start = int(start)
            end = int(end)
            primerID = int(primerID)
            if direction == 'Frd':
                if primerID > 1:
                    primerindex = primerID-2
                    for i in range(newstart,start):
                        allprimer.setdefault(i,primerindex)
                newstart = previous_end+1
            else:
                previous_end = end
    return allprimer

def process_primer_cov(freqfile,allprimer):
    primer_cov = [0]*total_primer
    primer_depth = [0] * total_primer
    newoutput = []
    newoutput.append('primer\tavg_depth\n')
    primer_cov_number = 0
    samplename = os.path.split(freqfile)[-1].split('.mapper1.vcf.final.allfreq.txt')[0]
    for lines in open(freqfile, 'r'):
        if not lines.startswith('CHR'):
            CHR,POS,REF,ALT,Depth = lines.split('\n')[0].split('\t')[:5]
            POS = int(POS)
            Depth = float(Depth)
            if POS in allprimer:
                primer_depth[allprimer[POS]]+=Depth
                primer_cov[allprimer[POS]] += 1
    for i in range(0,total_primer):
        if primer_cov[i] == 0:
            avg_depth = 0
        else:
            avg_depth = primer_depth[i]/primer_cov[i]
        newoutput.append('%s\t%.1f\n'%(i + 1,avg_depth))
        if avg_depth >= depth_cutoff:
            primer_cov_number += 1
    alloutput.append('%s\t%s\n'%(samplename,primer_cov_number))
    f1 = open(freqfile + '.primersum.txt','w')
    f1.write(''.join(newoutput))
    f1.close()


################################################### Main ########################################################
# load primer unique region (only covered by one primer)
allprimer = load_primer(primerfile)
# process freq files -> number primers with decent depth (depth_cutoff)
alloutput = []
alloutput.append('sample\tprimer_cover_number\n')
for freqfile in allfreqfiles:
    process_primer_cov(freqfile,allprimer)

f1 = open(args.o + '/../all.primer.coverage.txt','w')
f1.write(''.join(alloutput))
f1.close()