################################################### END ########################################################
################################################### SET PATH ########################################################
import glob
import os
from Bio import SeqIO

allfasta = glob.glob('/scratch/users/amyxiao/projects/COVID19/ALL_GISAID_SEQUENCES_2021-10-13/final_fasta/*.fasta')
reffasta = '/scratch/users/anniz44/genomes/covid/SNPcalling/reference.fa'
outputfolder = '/scratch/users/anniz44/genomes/covid/clinical/'
allfasta.sort()

################################################### Function ########################################################
def read_fasta(fasta,knownseq):
    newoutput = []
    knownseqnew = set()
    for record in SeqIO.parse(fasta, 'fasta'):
        record_seq = str(record.seq)
        if record_seq not in knownseq:
            knownseqnew.add(record_seq)
            newoutput.append('>%s\n%s\n'%(str(record.id),record_seq))
    f1 = open(os.path.join(outputfolder,os.path.split(fasta)[-1]),'w')
    f1.write(''.join(newoutput))
    f1.close()
    return set(list(knownseqnew)+list(knownseq))

def read_ref(fasta,knownseq): # not used
    ref_seq = ''
    for record in SeqIO.parse(fasta, 'fasta'):
        record_seq = str(record.seq)
        if record_seq not in knownseq:
            knownseq.add(record_seq)
            ref_seq = '>reference\n%s\n'%(record_seq)
    return [ref_seq,knownseq]

################################################### Main ########################################################
if False:
    knownseq = set()
    # load ref
    #ref_seq,knownseq = read_ref(reffasta,knownseq)
    # load fasta
    for fasta in allfasta:
        print(fasta)
        knownseq = read_fasta(fasta,knownseq)

# sequence filtering
length_cutoff = 29000# min total length
max_N = 0.05# max total Ns
for fasta in glob.glob('%s/*.fasta'%(outputfolder)):
    print('process %s'%(fasta))
    bad_sample = set()
    for record in SeqIO.parse(fasta, 'fasta'):
        record_seq = str(record.seq)
        record_id = str(record.id)
        record_seq = record_seq.replace('-','')
        record_seq_length = len(record_seq)
        if record_seq_length < length_cutoff or record_seq.count('N') + record_seq.count('n')>max_N*record_seq_length:
            bad_sample.add(record_id)
    print(len(bad_sample))
    newoutput = []
    for record in SeqIO.parse(fasta + '.align', 'fasta'):
        record_seq = str(record.seq)
        record_id = str(record.id)
        if record_id not in bad_sample:
            if record_id.startswith('EPI_ISL_1008203'):
                record_id = 'EPI_ISL_1008203:Globe_2021-01-21'
            newoutput.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(fasta + '.filtered.align', 'w')
    f1.write(''.join(newoutput))
    f1.close()
