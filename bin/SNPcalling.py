# start
# round 1-3 run genome assembly and map genomes to a reference genome
import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-fq",
                      help="file extension of WGS fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
required.add_argument("-ref",
                      help="reference genome",
                      type=str, default='reference.fa',
                      metavar='reference.fa')

# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='WGS/',
                      metavar='WGS/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-bw',
                          help="Optional: complete path to bowtie2 if not in PATH",
                          metavar="/usr/local/bin/bowtie2",
                          action='store', default='bowtie2', type=str)
optional.add_argument('-bcf',
                      help="Optional: complete path to bcftools if not in PATH",
                      metavar="/usr/local/bin/bcftools",
                      action='store', default='bcftools', type=str)
optional.add_argument('-sam',
                      help="Optional: complete path to bwa if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default='samtools', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
input_script_vcf = input_script + '/SNPcalling'
fastq_name = args.fq
sample_fastqall = glob.glob('%s/*%s'%(args.i,fastq_name))
output_dir = args.o + '/SNPcalling'

try:
    os.mkdir(args.o)
except IOError:
    pass

try:
    os.mkdir(output_dir)
except IOError:
    pass

reference_genome = output_dir+'/reference.fa'

os.system('cp %s %s'%(args.ref,reference_genome))

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass


################################################## Function ########################################################

def run_vcf(files,files2,database,tempbamoutput):
    # generate bam files
    cmds = ''
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds = args.bw + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
            tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    # merge samples bam files to vcf files
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            args.bcf, min(40, args.t), database,
            ' '.join(allsam), args.bcf, min(40, args.t), vcfoutput)
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += '%s view -H -v snps,indels %s.raw.vcf > %s.flt.snp.vcf \n' % (
            args.bcf, vcfoutput, vcfoutput)
    return cmds

def mapping_fastq(sample_fastqall):
    # mapping all fastq to reference genome
    for sample_fastq in sample_fastqall:
        original_folder, fastq_file_name = os.path.split(sample_fastq)
        sample = fastq_file_name.split(fastq_name)[0]
        filesize = 0
        try:
            filesize = int(os.path.getsize(os.path.join(output_dir + '/merge/', sample + '.all.flt.snp.vcf')))
        except FileNotFoundError:
            pass
        if filesize == 0:
            print('generate mapping code for %s' % (sample))
            cmds = ''
            sample_fastq2 = os.path.join(original_folder,sample + fastq_name.replace('1', '2'))
            outputbwa = os.path.join(output_dir + '/bwa',
                                               sample)
            results = run_vcf(sample_fastq, sample_fastq2,
                                  reference_genome,
                                  outputbwa)
            outputvcf = os.path.join(output_dir + '/merge', sample)
            cmds += results[0]
            cmds += merge_sample(reference_genome, outputvcf, [results[1]])
            f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (sample)), 'w')
            # for blast
            f1.write(
                '#!/bin/bash\nsource ~/.bashrc\n%s' % (
                    ''.join(cmds)))
            f1.close()


################################################## Main ########################################################
# generate code
mapping_fastq(sample_fastqall)

# sum all codes
f1 = open(os.path.join(input_script, 'allSNPcalling.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
# build bowtie2
f1.write(os.path.join(os.path.split('args.bw')[0], 'bowtie2-build') + ' %s %s\n' % (
        reference_genome, reference_genome))
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allSNPcalling.sh'%(input_script))
################################################### END ########################################################
