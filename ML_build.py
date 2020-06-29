import sys, getopt
import os
from random import random, randint

import Bio
import linecache
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def aln_cds(Parameters):
    seq_file = Parameters[0].split('/')[-1][:Parameters[0].split('/')[-1].index('_nc.fasta')]
    pwd = '/'.join(Parameters[0].split("/")[:-1]) + '/'
    pro_file = seq_file + '_aa.fasta'
    thread_number_iqtree = Parameters[2]
    print seq_file, pro_file
    ####Filter bad sequence###
    outfile_del = open(pwd + seq_file + '.del', 'w')
    dic = {}
    for seq_record in SeqIO.parse(pwd + pro_file, 'fasta'):
        seq_seq = str(seq_record.seq)
        if seq_seq.find('X') != -1:
            seq_seq = seq_seq.replace('X', '*')
            errfile.write('stop codon in string:' + ' ' + seq_file + '\n')
        dic[seq_record.id] = str(seq_seq)
    for seq_record in SeqIO.parse(pwd + seq_file + '_nc.fasta', 'fasta'):
        seq_id = str(seq_record.id)
        seq_seq = str(seq_record.seq)
        if seq_seq.find('n') != -1 or seq_seq.find('K') != -1 or seq_seq.find('W') != -1 or seq_seq.find(
                'N') != -1 or seq_seq.find('M') != -1 or seq_seq.find('Y') != -1:
            outfile_del.close()
            errfile.write('nnnn: found in' + seq_file + '\n')
            break
        outfile_del.write('>' + seq_id + '\n')
        Num = len(seq_seq)
        for start in range(Num):
            if dic[seq_id] + '*' == str(Seq(seq_seq, IUPAC.unambiguous_dna).translate()) or dic[seq_id] == str(
                    Seq(seq_seq, IUPAC.unambiguous_dna).translate()):
                break
            else:
                seq_seq = seq_seq[1:]
        mu = len(seq_seq) % 3
        if mu != 0:
            seq_seq = seq_seq[:-mu]
        if str(seq_seq)[-3:] in ['TAA', 'TAG', 'TGA', 'taa', 'tag', 'tga']:
            outfile_del.write(str(seq_seq)[:-3] + '\n')
        else:
            outfile_del.write(str(seq_seq) + '\n')
    outfile_del.close()

    ###alignment sequence###

    os.system("mafft --localpair --maxiterate 1000 %s_aa.fasta > %s_aa.aln 2> %s_mafft_aln.log" % (
    pwd + seq_file, pwd + seq_file, pwd + seq_file))
    os.system("perl translatorx_vLocal.pl -i %s.del -a %s_aa.aln -o %s 2> tranlatorx.log" % (
    pwd + seq_file, pwd + seq_file, pwd + seq_file))
    len_nc_aln = len((linecache.getline(pwd + seq_file + ".nt_ali.fasta", 2)).strip())
    print len_nc_aln, pwd + seq_file + ".nt.aln"
    ####running iqtree###
    if len_nc_aln != 0:
        os.system('mv %s.nt_* %saln_seqs/%s.nt.aln' % (pwd + seq_file, pwd, seq_file))
        pwd = pwd + '/aln_seqs/'
        outfile = open(pwd + seq_file + '.nt.aln.nex', 'w')
        outfile.write(
            "#nexus\nbegin sets;\n\tcharset p1 = 1-%d \ 3;\n\tcharset p2 = 2-%d \ 3;\n\tcharset p3 = 3-%d \ 3;\nend;" % (
            len_nc_aln, len_nc_aln, len_nc_aln))
        outfile.close()
        os.system("iqtree -s %s -bb 1000 -redo -nt %d -spp %s.nex -m MFPMERGE 1>%s_ML_iqtree.log" % (
        pwd + seq_file + '.nt.aln', thread_number_iqtree, pwd + seq_file + ".nt.aln", pwd + seq_file))


###for non-cds seqs####
def aln_noncds(Parameters):
    seq_file = Parameters[0].split('/')[-1][:Parameters[0].split('/')[-1].index('.fasta')]
    pwd = '/'.join(Parameters[0].split("/")[:-1]) + '/'
    thread_number_iqtree = Parameters[1]
    os.system("mafft --localpair --maxiterate 1000 %s.fasta > %s.aln 2> %s_mafft_aln.log" % (
    pwd + seq_file, pwd + seq_file, pwd + seq_file))
    ####running iqtree###
    os.system('mv %s.aln %saln_seqs/' % (pwd + seq_file, pwd))
    pwd = pwd + '/aln_seqs/'
    random_seed_number = randint(0, 2**32)
    os.system("iqtree -s %s -bb 1000 -redo -nt %d -m MFP -seed %d 1>%s_ML_iqtree.log" % (
    pwd + seq_file + '.aln', thread_number_iqtree, random_seed_number, pwd + seq_file))


###Main Process###
errfile = open('errfile.txt', 'a')
