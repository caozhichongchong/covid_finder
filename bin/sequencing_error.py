import glob
import os
import subprocess
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
################################################## Definition ########################################################
args = parser.parse_args()
allsnpfiles = glob.glob('%s/*.snpfreq.txt'%(args.o))

def snp_all(freqfile):
    samplename = os.path.split(freqfile)[-1].split('.mapper1.vcf.final.snpfreq.txt')[0]
    allfreqfile = freqfile.replace('snpfreq','allfreq')
    allcov = int(subprocess.check_output("wc -l %s"%(allfreqfile), shell=True).decode("utf-8").split(' ')[0])-1
    allsnp = int(subprocess.check_output("wc -l %s"%(freqfile), shell=True).decode("utf-8").split(' ')[0])-1
    alloutput.append('%s\t%s\t%s\t%s\n' % (samplename, allsnp/allcov,allsnp,allcov))


################################################### Main ########################################################
alloutput = []
alloutput.append('sample\tsequencing_error\tNo.SNPs\tCoverage\n')
for freqfile in allsnpfiles:
    snp_all(freqfile)

f1 = open(args.o + '/../all.sequencing.error.txt','w')
f1.write(''.join(alloutput))
f1.close()