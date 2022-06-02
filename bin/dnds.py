import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                      help="input snp sum",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/merge/summary/alldepthsum.txt',
                      metavar='alldepthsum.txt')
parser.add_argument("-ref",
                      help="input ref fasta",
                      type=str, default='/scratch/users/anniz44/scripts/covid/trial/reference.fasta',
                      metavar='reference.fasta')
parser.add_argument("-genelist",
                      help="input ref gene mapping",
                      type=str, default='/scratch/users/anniz44/genomes/covid/SNPcalling/ref_gene.txt',
                      metavar='ref_gene.txt')

################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
################################################### new class #########################################################
__metaclass__ = type

class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'
    def init(self, gene):
        self.gene = gene
        # [[N,S],freq]
        # not observed but predicted SNP pair by codon freq
        self.SNP_pair = {'A-T': [0,0],
                         'A-C': [0,0],
                         'G-C': [0,0],
                         'G-T': [0,0],
                         'A-G': [0,0],
                         'G-A': [0,0]}
        self.SNP_pair_freq = {'A-T': 0,
                         'A-C': 0,
                         'G-C': 0,
                         'G-T': 0,
                         'A-G': 0,
                         'G-A': 0}
        self.protein = ''
        self.NSratio = [0, 0]
    def addprotein(self, aa):
        self.protein += aa
    def addSNP_pair(self, pair, position, unique_snp_count,depth = 0):
        self.SNP_pair_freq[pair] += unique_snp_count
        self.NSratio[position] += unique_snp_count
        if position < 2 and depth == 0:
            # add to NS for each SNP pair of reference genes
            self.SNP_pair[pair][position] += unique_snp_count
    def addpredictSNP_pair(self,refSNP_pair_sum):
        for pair in refSNP_pair_sum:
            self.SNP_pair[pair][0] += refSNP_pair_sum[pair][0]
            self.SNP_pair[pair][1] += refSNP_pair_sum[pair][1]
    def sum_SNP_pair(self):
        self.SNP_pair_sum = {'A-T': [0, 0],
                         'A-C': [0, 0],
                         'G-C': [0, 0],
                         'G-T': [0, 0],
                         'A-G': [0, 0],
                         'G-A': [0, 0]}
        for pair in self.SNP_pair:
            self.SNP_pair_sum[pair][0] += self.SNP_pair[pair][0]
            self.SNP_pair_sum[pair][1] += self.SNP_pair[pair][1]
    def compute_expect(self):
        print(self.SNP_pair)
        print(self.SNP_pair_freq)
        self.expectNSratio = 'No_expect'
        expectNSratio = [0, 0]
        for pair in self.SNP_pair_freq:
            # use selected codon NS ratio (SNP pair) * all genes freq
            expectNSratio[0] += self.SNP_pair[pair][0] * self.SNP_pair_freq[pair]
            expectNSratio[1] += self.SNP_pair[pair][1] * self.SNP_pair_freq[pair]
        if expectNSratio[1] > 0:
            self.expectNSratio = expectNSratio[0] / expectNSratio[1]
        elif self.NSratio[0] == 0:
            # only S observed
            self.expectNSratio = 'expect_None'
        else:
            # only N observed
            self.expectNSratio = 'expect_N_only'

################################################### Function ########################################################


def translate(seq):
    seq = Seq(''.join(seq))
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

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
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
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def expectNSsub(record_name,record_seq,position=0):
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = ''.join(record_seq[(i * 3 + position):((i + 1) * 3 + position)])
        try:
            codon_NSratio = codontable_NSratio[codon]
            temp_SNP_gene.addprotein(codontable[codon])
            for pair in codon_NSratio.SNP_pair:
                temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0],0)
                temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1],0)
        except KeyError:
            pass
    temp_SNP_gene.sum_SNP_pair()
    return [temp_SNP_gene.SNP_pair_sum,temp_SNP_gene.protein,position]

def expectNS(record_name,record_seq):
    Total = int(len(record_seq) / 3)
    temp_result = expectNSsub(record_name, record_seq)
    if len(temp_result[1]) < 0.8 * Total:
        temp_result = expectNSsub(record_name, record_seq,1)
        if len(temp_result[1]) < 0.8 * Total:
            temp_result = expectNSsub(record_name, record_seq,2)
            if len(temp_result[1]) < 0.8 * Total:
                return 'None'
    return temp_result

def causeSNP(seq,position,ALT,Reverse_chr):
    seq = list(seq)
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq[position - 1]=ALT
    return seq

def refgenelist():
    refgene = dict()
    for lines in open(args.genelist,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        refgene.setdefault(lines_set[0],lines_set[1])
    return refgene

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Ref_NSratio=dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = refgene.get(str(record.id),str(record.id))
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_NSratio.setdefault(record_id,
                               expectNS(record_id, record_seq))
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Ref_NSratio,Mapping,Mapping_loci,Reverse]

def load_snpfile(depthfile):
    all_SNP_gene_temp = dict()
    for lines in open(depthfile, 'r'):
        if not lines.startswith('Color'):
            lines_set = lines.split('\t')[:16]
            if lines_set[4] == 'True':
                # good SNP only
                REF = lines_set[6]
                ALT = lines_set[7]
                if ALT in Allels:
                    #position = int(lines_set[12])
                    Chr = lines_set[11]
                    if Chr != 'None':
                        N_or_S = lines_set[13]
                        Chr = refgene.get(Chr,Chr)
                        # check gene
                        SNP_pair = transitions(REF, ALT)
                        # check new gene
                        New_gene = 0
                        if Chr not in all_SNP_gene_temp:
                            SNP_gene_temp = SNP_gene()
                            SNP_gene_temp.init(Chr)
                            all_SNP_gene_temp.setdefault(Chr, SNP_gene_temp)
                            New_gene = 1
                        SNP_gene_temp = all_SNP_gene_temp[Chr]
                        refSNP_pair_sum_all = Ref_NSratio.get(Chr, 'None')
                        if refSNP_pair_sum_all != 'None':
                            if New_gene == 1:
                                # predicted NS store in SNP_pair
                                refSNP_pair = refSNP_pair_sum_all[0]
                                SNP_gene_temp.addpredictSNP_pair(refSNP_pair)
                                SNP_gene_all.addpredictSNP_pair(refSNP_pair)
                            SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[N_or_S],
                                                      1, 2)
                            SNP_gene_all.addSNP_pair(SNP_pair, N_S_set[N_or_S],
                                                      1, 2)
    # sum result
    SNP_gene_all.compute_expect()
    expectNSratio = SNP_gene_all.expectNSratio
    print(expectNSratio)
    alloutput.append('allgenes\t%s\t%s\t%s\n'%(SNP_gene_all.NSratio[0],SNP_gene_all.NSratio[1],expectNSratio))
    for Chr in all_SNP_gene_temp:
        SNP_gene_temp = all_SNP_gene_temp[Chr]
        alloutput.append('%s\t%s\t%s\t%s\n' % (SNP_gene_temp.gene,SNP_gene_temp.NSratio[0], SNP_gene_temp.NSratio[1],expectNSratio))

# calculate codon freq
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
codontable_NSratio = dict()
for codon in codontable:
    SNP_gene_temp = SNP_gene()
    SNP_gene_temp.init(codon)
    codontable_NSratio.setdefault(codon, SNP_gene_temp)
    for position in range(0, 3):
        REF = codon[position]
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position+1, ALT,0)
                temp_NorS = dnORds(translate(codon)[0], translate(new_codon)[0])
                SNP_pair = transitions(REF, ALT)
                SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1, 0)
    SNP_gene_temp.compute_expect()

################################################### Main ########################################################
refgene = refgenelist()
# load database
Ref_seq, Ref_NSratio, Mapping, Mapping_loci, Reverse = loaddatabase(args.ref + '.fna')
# set up all genes dnds
SNP_gene_all = SNP_gene() # all denovo mutation
SNP_gene_all.init('allgenes')
alloutput = ['genename\tN\tS\texpected_ratio\n']
# load snps
load_snpfile(args.i)
f1 = open('%s.NS.txt'%(args.i),'w')
f1.write(''.join(alloutput))
f1.close()