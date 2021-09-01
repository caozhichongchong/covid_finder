################################################### END ########################################################
################################################### SET PATH ########################################################
# After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from statistics import median
from datetime import datetime
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-vcf",
                      help="file extension of vcfs with only SNPs(.flt.snp.vcf)",
                      type=str, default='.flt.snp.vcf',
                      metavar='.flt.snp.vcf')
required.add_argument("-fq",
                      help="file extension of fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
# requirement for software calling
optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='prodigal', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# set up path
input_script = args.s
vcf_name = '.raw.vcf'
fastq_name = args.fq
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
Allels_order = ['A','T','G','C']

min_qual_for_call = 40 #Remove sample*candidate that has lower than this quality
neighbour_cov_range = 100
min_cov = 0.75 # at least min_cov * median depth of neighbour_cov_range bp
min_cov_for_call_SNP = 10 #Remove sample*candidate

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def all_ALT(Allels_count):
    allALT = set()
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            allALT.add(Allels_order[alleles])
    return allALT

def outputvcf(output_name,vcf_file_list_freq,vcf_file_list_freq_snp,Sample_name):
    vcf_file_filtered = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write('CHR\tPOS\tREF\tALT\tDepth\tDepth_median_cutoff\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\t%s\n'%(
        '\t'.join(Sample_name))
                            + ''.join(vcf_file_list_freq_snp))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.allfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write(
        'CHR\tPOS\tREF\tALT\tDepth\tDepth_median_cutoff\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\t%s\n' % (
            '\t'.join(Sample_name))
        + ''.join(vcf_file_list_freq))
    vcf_file_filtered.close()

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

def SNP_check_fq(lines_set,vcf_file_list_freq, vcf_file_list_freq_snp,min_cov_for_call_SNP, depth_cutoff):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    REF = lines_set[3]
    allels_set = [REF]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    SNP_quality = lines_set[5]
    allels_set += lines_set[4].split(',')
    allels_set = [x for x in allels_set if x in Allels]
    Total_alleles = len(allels_set)
    Subdepth_all = lines_set[9]
    Qual = float(SNP_quality)
    if Qual >= min_qual_for_call:
        # Keep SNP that has higher than this quality
        Allels_frq = [0, 0, 0, 0, 0]
        Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        Subdepth_forward = [int(i) for i in Subdepth_all.split(':')[-3].split(',')]
        Subdepth_reverse = [int(i) for i in Subdepth_all.split(':')[-2].split(',')]
        for num_allels in range(0, min(len(Subdepth), Total_alleles)):
            allels = allels_set[num_allels]
            forward = Subdepth_forward[num_allels]
            reverse = Subdepth_reverse[num_allels]
            if forward + reverse < min_cov_for_call_SNP:
                forward = 0
                reverse = 0
            Allels_frq[Allels[allels]] += forward + reverse
            Allels_frq_sub[Allels[allels] * 2] += forward
            Allels_frq_sub[Allels[allels] * 2 + 1] += reverse
        total_depth = sum(Allels_frq)
        # find major alt and calculate frequency
        if sum(Allels_frq) > 0:
            allALT = all_ALT(Allels_frq)
            if total_depth >= depth_cutoff:
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                # calculate NS
                gene_info = contig_to_gene(CHR, POS)
                if gene_info != []:
                    Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                    if Ref_seq_chr != 'None':
                        #  observed NS ratio calculated
                        temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                        if codon_start <= POS_gene - 1:
                            Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_chr = Ref_seq_chr
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                temp_snp_line_AA += Ref_seq_aa
                                for ALT in allels_set:
                                    if ALT != REF:
                                        SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                        SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                        SNP_seq_aa = translate(SNP_seq_codon)[0]
                                        temp_snp_line_AA += SNP_seq_aa
                                        temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                        temp_snp_line_NS[-1] += temp_NorS
                # output lines and output major alt
                temp_snp_line.append(CHR)
                temp_snp_line.append(str(POS))
                temp_snp_line.append(REF)
                temp_snp_line.append(','.join(list(allALT)))
                vcf_file_list_freq.append(
                    '\t'.join(temp_snp_line) + '\t%s\t%.1f\t%.1f\t%s\t%s\t%s\n' % (
                        total_depth,depth_cutoff,
                        Qual, '\t'.join(temp_snp_line_NS),
                        temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
                if allALT != set(REF):
                    # a SNP
                    vcf_file_list_freq_snp.append(
                        '\t'.join(temp_snp_line) + '\t%s\t%.1f\t%.1f\t%s\t%s\t%s\n' % (
                            total_depth, depth_cutoff,
                            Qual, '\t'.join(temp_snp_line_NS),
                            temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
    return [vcf_file_list_freq,vcf_file_list_freq_snp]

def load_sample(vcf_file):
    Sample_name = []
    for lines in open(vcf_file.split(args.vcf)[0] + vcf_name, 'r'):
        if lines.startswith('##reference=file:'):
            # set database
            database_file = lines.split('##reference=file:')[1].split('\n')[0]
        if lines.startswith('#CHROM'):
            Sample_name.append(os.path.split(lines.split('\t')[9].split('\n')[0])[-1].split('.sorted.bam')[0])
            break
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    return [database_file,Sample_name]

def SNP_filter(vcf_file,Sample_name,output_name,min_cov_for_call_SNP):
    vcf_file_list_freq = []
    vcf_file_list_freq_snp = []
    alldepth = []
    m = 0
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#") and not lines.startswith("CHR"):
            lines_set = lines.split('\n')[0].split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            alldepth.append(Depth)
            if len(alldepth) > 10:
                depth_cutoff = min_cov * median(alldepth[max(len(alldepth)-neighbour_cov_range,0):])
            else:
                depth_cutoff = min_cov_for_call_SNP
            if '.' not in lines_set[4]:
                # potential SNP
                if Depth >= depth_cutoff:
                    # Remove candidate locations have lower than this coverage
                    m += 1
                    if m % 1000 == 0:
                        print('%s processed %s SNPs' % (datetime.now(), m))
                    vcf_file_list_freq,vcf_file_list_freq_snp = \
                        SNP_check_fq(lines_set,
                                     vcf_file_list_freq,vcf_file_list_freq_snp,
                                     min_cov_for_call_SNP, depth_cutoff)
            else:# NEED CHANGE
                if Depth >= depth_cutoff:
                    # Remove candidate locations have lower than this coverage
                    m += 1
                    if m % 1000 == 0:
                        print('%s processed %s SNPs' % (datetime.now(), m))
                    vcf_file_list_freq,vcf_file_list_freq_snp = \
                        SNP_check_fq(lines_set,
                                     vcf_file_list_freq,vcf_file_list_freq_snp,
                                     min_cov_for_call_SNP, depth_cutoff)
    outputvcf(output_name,vcf_file_list_freq,vcf_file_list_freq_snp,Sample_name)

################################################### Main ########################################################
# run vcf filtering for bowtie
output_name = 'final'
allvcf_file = glob.glob(os.path.join(args.i, '*%s' % (args.vcf)))
print(allvcf_file)
for vcf_file in allvcf_file:
    try:
        f1 = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'r')
    except IOError:
        sample = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        print('%s load database for %s' % (datetime.now(), sample))
        Ref_seq = dict()
        database_file, Sample_name = load_sample(vcf_file)
        if Ref_seq == dict():
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        # SNP filtering
        print('%s start filtering SNPs %s' % (datetime.now(), sample))
        SNP_filter(vcf_file, Sample_name, output_name,
                   min_cov_for_call_SNP)
        print('%s finished output %s' % (datetime.now(), sample))


