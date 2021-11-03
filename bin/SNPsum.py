################################################### END ########################################################
################################################### SET PATH ########################################################
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import math
from scipy.stats import binom
import statistics
# set up path
import argparse
import datetime
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/merge/freqfiles/',
                      metavar='input/')
required.add_argument("-ref",
                      help="ref file dir",
                      type=str, default='/scratch/users/anniz44/scripts/covid/trial/reference.fasta',
                      metavar='reference.fasta')
required.add_argument("-cov",
                      help="primer coverage file",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/merge/all.primer.coverage.txt',
                      metavar='all.primer.coverage.txt (output of primer_cov.py)')
required.add_argument("-clinical",
                      help="summary of clinical samples",
                      type=str, default='/scratch/users/anniz44/genomes/covid/clinical/allearlymut.txt',
                      metavar='allearlymut.txt  (output of clinical_sum.py)')
optional.add_argument("-state",
                      help="samples of which state",
                      type=str, default='None',
                      metavar='for example, MA or TX')

################################################## Definition ########################################################
args = parser.parse_args()

################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
separate_donor_genome = []
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels['-']=4
Allels['+A']=5
Allels['+T']=6
Allels['+G']=7
Allels['+C']=8
Allels_order = ['A','T','G','C','-','+A','+T','+G','+C']
normal_alleles = ['A','T','G','C']
# set up filtering parameter
SNP_prevalence_cutoff = 5
SNP_prevalence_cutoff_strict_ALT_freq = 3
max_sequencing_error_rate = 0.008
################################################### Function ########################################################
# set up functions
def alt_depth(Allels_count,ALT):
    return sum(Allels_count[Allels_order.index(ALT)*4:Allels_order.index(ALT)*4+4])

def load_snp(snpfile,allSNP,total_number_samples,this_sample_number):
    for lines in open(snpfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            POS = lines_set[1]
            ALTset = lines_set[3].split(',')
            REF = lines_set[2]
            ALTset = [ALT for ALT in ALTset if ALT != REF]
            Allels_count = lines_set[-1].split(';')
            Allels_count = [float(i) for i in Allels_count]
            alldepth = sum(Allels_count)
            for i in range(0,len(ALTset)):
                ALT = ALTset[i]
                POS_REF_ALT = '%s\t%s\t%s'%(POS,REF,ALT)
                if POS_REF_ALT not in allSNP:
                    # set up SNP
                    allSNP.setdefault(POS_REF_ALT,[[REF],[0]*total_number_samples,[0]*total_number_samples])# infor, sample
                # store ALT depth and ALT freq of this sample
                ALTdepth = alt_depth(Allels_count,ALT)
                allSNP[POS_REF_ALT][1][this_sample_number] += ALTdepth
                allSNP[POS_REF_ALT][2][this_sample_number] += float(ALTdepth)/float(alldepth)
    return allSNP

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def find_SNP_geneinfor(allSNP):
    CHR = 'MN908947.3'
    for POS_REF_ALT in allSNP:
        POS, REF, ALT = POS_REF_ALT.split('\t')
        gene_info = contig_to_gene(CHR, float(POS))
        REF = allSNP[POS_REF_ALT][0][0]
        temp_snp_line_NS = ['None', 'None', 'None', '']
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '','']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_NS[-1] = Ref_seq_aa + str(int(codon_start/3) + 1)
                        if ALT in normal_alleles:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_NS[-1] += SNP_seq_aa
                            temp_snp_line_NS[-2] = dnORds(Ref_seq_aa, SNP_seq_aa)
                        else:
                            temp_snp_line_NS[-1] = 'Indel'
                            temp_snp_line_NS[-2] = 'N'
        allSNP[POS_REF_ALT][0] = '\t'.join(temp_snp_line_NS)
    return allSNP

def binom_prob(depthsum,freqsum,allfreq_order):
    for i in range(0, total_number_samples):
        if depthsum[i] > 0:
            x = depthsum[i]
            n = int(depthsum[i]/freqsum[i])
            prob_error = binom.pmf(x, n, max_sequencing_error_rate)
            prob_mut = binom.pmf(x, n, allfreq_order[i])
            if prob_mut <= prob_error:
                freqsum[i]=0
                depthsum[i]=0
    return[depthsum,freqsum]

def outputSNP(allSNP,allsamples,allfreq,clinical_variant):
    allfreq_order = []
    pass_num = 0
    sample_time = dict()
    for i in range(0,total_number_samples):
        sample = allsamples[i]
        thissample_time = sample.split('_')[-1]
        if thissample_time == '1.1':
            thissample_time = '11.1'
        thissample_time = thissample_time.split('.')
        sample_time.setdefault(i,datetime.date(2020,int(thissample_time[0]),int(thissample_time[1])))
        if sample in allfreq:
            allfreq_order.append(allfreq[sample])
        else:
            allfreq_order.append(1)
    alloutput = []
    alloutput.append('Color\tClinical\tClinical_time\tSampling_time\tGoodSNP\tPOS\tREF\tALT\tPrevalence\tPrevalence_ALTfreq_cutoff\tGene\tGene_POS\tN_S\tAAchange\t%s\n'%('\t'.join(allsamples)))
    alloutputfreq = []
    alloutputfreq.append('Color\tClinical\tClinical_time\tSampling_time\tGoodSNP\tPOS\tREF\tALT\tPrevalence\tPrevalence_ALTfreq_cutoff\tAvg_ALT_freq\tGene\tGene_POS\tN_S\tAAchange\t%s\n' % ('\t'.join(allsamples)))
    for POS_REF_ALT in allSNP:
        geneinfor, depthsum, freqsum = allSNP[POS_REF_ALT]
        prevalence = len([i for i in depthsum if i > 0])
        depthsum, freqsum = binom_prob(depthsum,freqsum,allfreq_order)
        prevalence_strict = len([i for i in depthsum if i > 0])
        INclinical = POS_REF_ALT in clinical_variant
        clinical_time = clinical_variant.get(POS_REF_ALT,'')
        if prevalence_strict > 0:
            if INclinical:
                clinical_time2 = clinical_time.split('-')
                clinical_time2 = datetime.date(int(clinical_time2[0]), int(clinical_time2[1]), int(clinical_time2[2]))
                allsampletime_withmut = [sample_time[i] for i in range(0, total_number_samples) if depthsum[i] > 0]
                allsampletime_withmut.sort()
                if allsampletime_withmut[0] < clinical_time2:
                    Color = 'red'
                else:
                    Color = 'blue'
            else:
                Color = 'grey'
            SNP_confident = 'False'
            if prevalence >= SNP_prevalence_cutoff:
                SNP_confident = 'True'
            elif prevalence_strict >= SNP_prevalence_cutoff_strict_ALT_freq:
                SNP_confident = 'True'
                pass_num += 1
            Avg_ALT_freq = '%.2f' % (100*statistics.mean([i for i in freqsum if i > 0]))
            prevalence = '%.2f' % (100 * prevalence / total_number_samples)
            prevalence_strict = '%.2f' % (100 * prevalence_strict / total_number_samples)
            alloutput.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(Color,INclinical,clinical_time,str(allsampletime_withmut[0]),SNP_confident,POS_REF_ALT,prevalence,prevalence_strict,geneinfor,'\t'.join(['%.0f'%i for i in depthsum])))
            alloutputfreq.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (Color,INclinical,clinical_time,str(allsampletime_withmut[0]),SNP_confident,POS_REF_ALT, prevalence,prevalence_strict,Avg_ALT_freq,geneinfor, '\t'.join(['%.2f'%(i*100) for i in freqsum])))
    print(pass_num)
    if args.state == 'None':
        f1 = open('%s/alldepthsum.txt'%(args.i),'w')
        f1.write(''.join(alloutput))
        f1.close()
        f1 = open('%s/allfreqsum.txt' % (args.i), 'w')
        f1.write(''.join(alloutputfreq))
        f1.close()
    else:
        f1 = open('%s/%s.depthsum.txt' % (args.i,args.state), 'w')
        f1.write(''.join(alloutput))
        f1.close()
        f1 = open('%s/%s.freqsum.txt' % (args.i,args.state), 'w')
        f1.write(''.join(alloutputfreq))
        f1.close()

def load_primer_coverage(covfile):
    allcov = dict()
    max_primer_cover = 0
    # load all cov
    for lines in open(covfile,'r'):
        if not lines.startswith('sample'):
            sample,primer_cov = lines.split('\n')[0].split('\t')
            primer_cov = int(primer_cov)
            allcov.setdefault(sample,primer_cov)
            max_primer_cover = max(max_primer_cover,primer_cov)
    # evaluate molecules
    newoutput = []
    newoutput.append('sample\tprimer_cover_number\tmolecule_number\n')
    p_bp = (max_primer_cover-1)/max_primer_cover
    allfreq = dict()
    for sample in allcov:
        primer_cov = allcov[sample]
        if primer_cov == 0:
            lambda_mol = 0
        else:
            p = primer_cov/max_primer_cover
            lambda_mol = math.log(max(1 - p,1E-2), p_bp)/max_primer_cover
            allfreq.setdefault(sample, 1 / max(lambda_mol,1E-2))
        newoutput.append('%s\t%s\t%s\n'%(sample,primer_cov,lambda_mol))
    f1 = open(covfile + '.molecule.txt', 'w')
    f1.write(''.join(newoutput))
    f1.close()
    return allfreq

def load_clinical_variant(clinical_file):
    clinical_variant = dict()
    for lines in open(clinical_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            try:
                POS,REF,ALT,year,month,day = lines_set[0:6]
            except ValueError:
                print(lines_set)
                POS, REF, ALT = lines_set[0:3]
                year, month, day = ['2021','01','21']
            mut = '%s\t%s\t%s'%(POS,REF,ALT)
            clinical_variant.setdefault(mut,'%s-%s-%s'%(year,month,day))
    return clinical_variant
################################################### Main ########################################################
# load clinical variants
clinical_variant = load_clinical_variant(args.clinical)
# evaluate molecules and ALT frequency by primer coverage
allfreq = load_primer_coverage(args.cov)
# read and store all SNPs
allSNP = dict()
if args.state == 'None':
    allSNPfiles = glob.glob('%s/*snpfreq.txt' % (args.i))
else:
    allSNPfiles = glob.glob('%s/*_%s*snpfreq.txt' % (args.i,args.state))
total_number_samples = len(allSNPfiles)
this_sample_number = 0
allsamples = []
print('process %s samples' % (total_number_samples))
for SNPfile in allSNPfiles:
    samplename = os.path.split(SNPfile)[-1].split('.mapper1.vcf.final.snpfreq.txt')[0]
    allsamples.append(samplename)
    allSNP = load_snp(SNPfile, allSNP, total_number_samples, this_sample_number)
    this_sample_number += 1

# process all SNPs gene information
database_file = args.ref + '.fna'
Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
allSNP = find_SNP_geneinfor(allSNP)

# output all SNPs infor
print('output SNP summary')
outputSNP(allSNP, allsamples,allfreq,clinical_variant)

import glob,os
allsh = glob.glob('*mapper1.vcf.sh*')
for shfile in allsh:
    if '.err' not in shfile and '.out' not in shfile:
        try:
            f1 = open('%s.err'%(shfile),'r')
            os.system('mv %s* finished/'%(shfile))
        except IOError:
            os.system('cat %s %s>%snew'%('addtitle.txt',shfile,shfile))
            os.system('mv %s finished/'%(shfile))
