def main():
    ###Checking dependece####
    try:
        import multiprocessing, linecache, sys, getopt, os, Bio, threading, glob, subprocess
        from ete3 import Tree
        from random import random
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        print "ALL good for python modules"
    except Exception:
        print 'You have to install related python modules'
        sys.exit()

    if os.system("mkdir test_test_test_test") != 0 or os.system("rm -r test_test_test_test") != 0:
        print 'You need to have authority to make/delete directory'
        sys.exit()

    try:
        print 'check mafft'
        subprocess.check_call(['which', 'mafft'])
        print "mafft pass"
    except:
        print 'You have to install mafft: https://mafft.cbrc.jp/alignment/software/'
        sys.exit()
    try:
        print 'check IQ-TREE'
        subprocess.check_call(['which', 'iqtree'])
        print "IQ-TREE pass"
    except:
        print 'You have to install IQ-TREE: http://www.iqtree.org/#download'
    try:
        print 'check julia'
        subprocess.check_call(['which', 'julia'])
        print "julia pass"
    except:
        print 'You have to install julia: https://julialang.org/downloads/'
        sys.exit()
    try:
        print 'check java'
        subprocess.check_call(['which', 'java'])
        print "java pass"
    except:
        print 'You have to install java: https://www.java.com/en/download/manual.jsp'
        sys.exit()
    try:
        print 'check stells-v2'
        subprocess.check_call(['which', 'stells-v2'])
        print "stells-v2 pass"
    except:
        print 'You have to install stells-v2: https://github.com/yufengwudcs/STELLS2'
        sys.exit()
    try:
        print 'check mpest'
        subprocess.check_call(['which', 'mpest'])
        print "mpest pass"
    except:
        print 'You have to install mpest: http://faculty.franklin.uga.edu/lliu/mp-est'
        sys.exit()
    try:
        print 'Translatorx/ML_build/Spetree/Phylonet/astral/AMAS.py:'
        subprocess.check_call(
            ['ls', 'translatorx_vLocal.pl', 'ML_build.py', 'Spetree.py', 'PhyloNet_3.6.8.jar', 'astral.5.6.3.jar',
             'AMAS.py'])
        print "All good for scripts"
    except:
        print 'You have to install TRANSLATORX or ML_build.py or Spetree.py or PhyloNet_3.6.8.jar or astral.5.6.3.jar or AMAS.py in current directory'
        sys.exit()
    return;


####Help Page######

def Help():
    print """This python script will help you to infer species tree or Phylonetwork from gene trees, More information in README.file\n 
        Parameters: \n 	
	-d: diretory path, it contains all the fasta files\n 
	-c: type name, it has to be "CDS" or "nonCDS". "CDS" means you will provide nuclear seqs and corresponding protein seqs\n 
        -r: root taxon(s), it should be species name(s), if more than two species, please use "," to link two species names\n 
	-s: speices name file \n 
	-g: gene name file \n 
	-b: Bootstrap value, default 50 \n 
	-k: Cross value, it should be divisible by the numerber of gene trees (wc all_iqtree_rooted.txt)\n 
	-t: Threads number\n 
	-n: Maximum netwrok numbers\n """
    return;


if __name__ == '__main__':
    main()
