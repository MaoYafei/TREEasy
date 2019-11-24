import sys, getopt
import threading
import time
import os
import glob
from Spetree import run_Spetree

from ML_build import *
from TREEasy_Help import *

def dependency_check():
    ####Dependence Check####
    try:
        import multiprocessing, linecache, sys, getopt, os, Bio, threading, glob, subprocess
        from ete3 import Tree
        from random import random
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        subprocess.check_call(['which', 'mafft', 'julia', 'java', 'stells-v2', 'mpest', 'iqtree', 'raxml-ng'])
        subprocess.check_call(
            ['ls', 'translatorx_vLocal.pl', 'ML_build.py', 'Spetree.py', 'PhyloNet_3.6.8.jar', 'astral.5.6.3.jar',
             'AMAS.py'])
        print "Passing checking"
    except Exception:
        print """Please run 'python TREEasy_Help.py' to see which dependencies you need to install"""
        return False

    if os.system("mkdir test_test_test_test") != 0 or os.system("rm -r test_test_test_test") != 0:
        print 'You need to have authority to make/delete directory'
        return False
    return True


def run_TREEasy(seq_file, roottaxon, type, species_namefile, gene_namefile, boot_value, Net_num,
                cross_value, thread_number, message_queue=None):
    ####Main Running####
    if message_queue:
        message_queue.put('ML tree running...')

    if type == "CDS":
        L = glob.glob(seq_file + "/*_nc.fasta")
    else:
        L = glob.glob(seq_file + "/*.fasta")
    os.system('rm -r %s/aln_seqs' % (seq_file))
    os.system('mkdir %s/aln_seqs' % (seq_file))
    ThreadList = []
    if type == "CDS":  # CDS
        for i in L:
            Parameters = [i, i.replace("_nc.fasta", "_aa.fasta"), thread_number]
            t = threading.Thread(target=aln_cds, args=(Parameters,))
            t.setDaemon(True)
            ThreadList.append(t)
    else:  # NONCDS
        for i in L:
            Parameters = [i, thread_number]
            t = threading.Thread(target=aln_noncds, args=(Parameters,))
            t.setDaemon(True)
            ThreadList.append(t)

    for t in ThreadList:
        t.start()
    for t in ThreadList:
        t.join()
    print 'ML tree done'
    if message_queue:
        message_queue.put('ML tree done')

    ######Species tree inferrring#######
    run_Spetree(seq_file + '/', species_namefile, gene_namefile, boot_value, roottaxon, Net_num, cross_value, thread_number,
                message_queue=message_queue)


if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "hd:n:r:c:s:g:b:k:t:")
    type = ""
    Net_num = 5
    boot_value = 50


    for op, value in opts:
        if op == "-d":
            seq_file = value
        elif op == "-r":
            roottaxon = value
        elif op == "-c":
            type = value
        elif op == "-s":
            species_namefile = value
        elif op == "-g":
            gene_namefile = value
        elif op == "-b":
            boot_value = int(value)
        elif op == "-n":
            Net_num = int(value)
        elif op == "-k":
            cross_value = int(value)
        elif op == "-t":
            thread_number = int(value)
            if int(value) >= 4:
                thread_number_iqtree = 3
            else:
                thread_number_iqtree = int(value)
        elif op == "-h":
            Help()
            sys.exit()

    ###Parameter Check####
    print len(opts), 'Parameters'

    if len(opts) != 9:
        if ('-h', '') in opts:
            Help()
            sys.exit()
        else:
            print 'can not work, you need to put proper parameters(total 9 parameters), see help(-h)'
            sys.exit()

    if not dependency_check():
        sys.exit()

    run_TREEasy(seq_file, roottaxon, type, species_namefile, gene_namefile, boot_value, Net_num, cross_value, thread_number_iqtree)



