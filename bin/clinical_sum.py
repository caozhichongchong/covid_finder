################################################### END ########################################################
################################################### SET PATH ########################################################
import glob
import os
import argparse
import datetime
from copy import deepcopy
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/covid/clinical/',
                      metavar='input/')

################################################## Definition ########################################################
args = parser.parse_args()
allvcf = glob.glob('%s/*.fasta.filtered.align.vcf'%(args.i))
allvcf.sort()
print(allvcf)
DNA_ambiguous = {
    'R':['A','G'],
    'Y':['T','C'],
    'K':['G','T'],
    'M':['A','C'],
    'S':['G','C'],
    'W':['A','T'],
    'B':['C','G','T'],
    'D':['A','G','T'],
    'H':['A','C','T'],
    'V':['A','C','G']
}
################################################### Function ########################################################
# set up functions
def checksampletime(sample,sampleID,samplingtime,sample_sampleID,time_sample_dict,time_sample_set):
    if sample.startswith('EPI_ISL_1008203:Achenkirch'):
        sample = 'EPI_ISL_1008203:Achenkirch_2021-01-21'
    time_thissample = time_sample(sample, samplingtime)
    sample_sampleID.setdefault(sample, sampleID)
    time_sample_dict.setdefault(time_thissample, [])
    time_sample_dict[time_thissample].append(sample)
    time_sample_set.add(time_thissample)
    return [sample_sampleID,time_sample_dict,time_sample_set]

def timeorder(time_sample_set,time_sample_dict,sampleID_order,sample_sampleID):
    # order by time
    time_sample_set = list(time_sample_set)
    time_sample_set.sort()
    for time_thissample in time_sample_set:
        for sample in time_sample_dict[time_thissample]:
            sampleID_order.append(sample_sampleID[sample])
    return sampleID_order

def load_vcf(vcffile,variant_sum_set,firstoutput = True):
    samplingtime = os.path.split(vcffile)[-1].split('_')[0].split('-')
    sample_sampleIDTX = dict()
    time_sample_dictTX = dict()
    time_sample_setTX = set()
    sampleID_orderTX = []
    newoutputTX = []
    sample_sampleIDMA = dict()
    time_sample_dictMA = dict()
    time_sample_setMA = set()
    sampleID_orderMA = []
    sample_sampleID = dict()
    time_sample_dict = dict()
    time_sample_set = set()
    sampleID_order = []
    newoutputMA = []
    alloutput = []
    if firstoutput:
        newoutputMA.append('POS\tREF\tALT\tyear\tmonth\tday\tsamples\n')
        newoutputTX.append('POS\tREF\tALT\tyear\tmonth\tday\tsamples\n')
        alloutput.append('POS\tREF\tALT\tdate\tstate\tsamples\n')
        fall = open('%s/allmut.txt' % (args.i), 'w')
        if False:
            fTX = open('%s/allearlymutTX.txt' % (args.i), 'w')
            fMA = open('%s/allearlymutMA.txt' % (args.i), 'w')
    else:
        fall = open('%s/allmut.txt' % (args.i), 'a')
        if False:
            fTX = open('%s/allearlymutTX.txt' % (args.i), 'a')
            fMA = open('%s/allearlymutMA.txt' % (args.i), 'a')
    for lines in open(vcffile,'r'):
        if lines.startswith('#CHROM'):
            allsamplename_line = lines.split('\n')[0].split('\t')
            total_sample = len(allsamplename_line)
            for sampleID in range(10, total_sample):
                sample = allsamplename_line[sampleID]
                if '_USA_' in sample:
                    sample_sampleID, time_sample_dict, time_sample_set = checksampletime(sample, sampleID,
                                                                                               samplingtime,
                                                                                               sample_sampleID,
                                                                                               time_sample_dict,
                                                                                               time_sample_set)
                if False:
                    if '_TX_' in sample:
                        sample_sampleIDTX, time_sample_dictTX, time_sample_setTX = checksampletime(sample, sampleID,samplingtime,
                                                                                             sample_sampleIDTX,
                                                                                             time_sample_dictTX,
                                                                                             time_sample_setTX)
                    if 'MA' in sample:
                        sample_sampleIDMA, time_sample_dictMA, time_sample_setMA = checksampletime(sample, sampleID, samplingtime,
                                                                                             sample_sampleIDMA,
                                                                                             time_sample_dictMA,
                                                                                             time_sample_setMA)
            # order by time
            sampleID_order = timeorder(time_sample_set, time_sample_dict, sampleID_order, sample_sampleID)
            if False:
                sampleID_orderTX = timeorder(time_sample_setTX,time_sample_dictTX,sampleID_orderTX,sample_sampleIDTX)
                sampleID_orderMA = timeorder(time_sample_setMA, time_sample_dictMA,sampleID_orderMA,sample_sampleIDMA)
            if len(sampleID_order) == 0:
                print('no USA samples found in %s'%(vcffile))
                fall.write(''.join(alloutput))
                fall.close()
                return
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS,useless,REF,ALT = lines_set[:5]
            allALT = ALT.split(',')
            allALTset = dict()
            allALTset.setdefault('0',REF)
            if lines_set[9]!= '0':
                print('wrong lines %s %s'%(REF,lines_set[9]), lines)
            for i in range(0,len(allALT)):
                ALT = allALT[i]
                if ALT not in ['*','N']:
                    if ALT in DNA_ambiguous:
                        ALTsubset = deepcopy(DNA_ambiguous[ALT])
                        if REF in ALTsubset:
                            ALTsubset.remove(REF)
                        allALTset.setdefault(str(i+1), ALTsubset)
                    else:
                        allALTset.setdefault(str(i + 1), [ALT])
            if len(allALTset) > 1:
                # with ALT
                # all mut
                for sampleID in sampleID_order:
                    sampleALT = lines_set[sampleID]
                    if sampleALT in allALTset and sampleALT!= '0':
                        ALTsubset = allALTset[sampleALT]
                        for ALT in ALTsubset:
                            mut = '%s\t%s\t%s'%(POS,REF,ALT)
                            sample = allsamplename_line[sampleID]
                            alloutput.append('%s\t%s\t%s\t%s\n' % (
                                mut, sample.split('_')[-1], sample.split(':')[1].split('_')[1],sample))

                if False:
                    # TX
                    for sampleID in sampleID_orderTX:
                        sampleALT = lines_set[sampleID]
                        if sampleALT in allALTset and sampleALT!= '0':
                            ALTsubset = allALTset[sampleALT]
                            for ALT in ALTsubset:
                                mut = '%s\t%s\t%s'%(POS,REF,ALT)
                                sample = allsamplename_line[sampleID]
                                if mut not in variant_sum_set:
                                    variant_sum_set.add(mut)
                                    newoutputTX.append('%s\t%s\t%s\n' % (
                                    mut, '\t'.join(sample.split('_')[-1].split('-')), sample))
                    # MA
                    for sampleID in sampleID_orderMA:
                        sampleALT = lines_set[sampleID]
                        if sampleALT in allALTset and sampleALT != '0':
                            ALTsubset = allALTset[sampleALT]
                            for ALT in ALTsubset:
                                mut = '%s\t%s\t%s' % (POS, REF, ALT)
                                sample = allsamplename_line[sampleID]
                                if mut not in variant_sum_set:
                                    variant_sum_set.add(mut)
                                    newoutputMA.append('%s\t%s\t%s\n' % (
                                        mut, '\t'.join(sample.split('_')[-1].split('-')), sample))
    if False:
        fTX.write(''.join(newoutputTX))
        fTX.close()
        fMA.write(''.join(newoutputMA))
        fMA.close()
    fall.write(''.join(alloutput))
    fall.close()
    return variant_sum_set

def time_sample(sample,samplingtime):
    try:
        settime = [int(i) for i in sample.split('_')[-1].split('-')]
    except ValueError:
        print(sample)
        settime = [int(samplingtime[0]),int(samplingtime[1]),1]
    return datetime.date(settime[0], settime[1], settime[2])

def update_early_time(oldsample,newsample,sample_time): # not used
    date1 = sample_time[oldsample]
    date2 = sample_time[newsample]
    if date1 <= date2:
        return oldsample
    else:
        return newsample

################################################### Main ########################################################
# process all vcf
variant_sum_set = set()
for vcffile in allvcf:
    print('process %s'%(vcffile))
    # summarize variants
    variant_sum_set = load_vcf(vcffile,variant_sum_set,vcffile==allvcf[0])

