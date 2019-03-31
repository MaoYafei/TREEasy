
def main():
###Checking dependece####
	try:
		import multiprocessing,linecache,sys,getopt,os,Bio,threading,glob,subprocess
		from ete3 import Tree
		from random import random
		from Bio import SeqIO
		from Bio.Seq import Seq
		from Bio.Alphabet import IUPAC
	except Exception:
		print 'You have to install related python modules'
		sys.exit()

	if os.system("mkdir test_test_test_test")!=0 or os.system("rm -r test_test_test_test")!=0:
        	print 'You need to have authority to make/delete directory'
        	sys.exit()

	try:
		print 'check mafft/julia/java/stells-v2/mpest:'
		subprocess.check_call(['which', 'mafft','julia','java','stells-v2','mpest'])
	except:
		print 'You have to install mafft/julia/java/stells-v2/mpest'
		sys.exit()

	try:
		print 'Translatorx/ML_build/Spetree/Phylonet/astral/AMAS.py:'
        	subprocess.check_call(['ls', 'translatorx_vLocal.pl','ML_build.py','Spetree.py','PhyloNet_3.6.8.jar','astral.5.6.3.jar','AMAS.py'])
	except:
        	print 'You have to install TRANSLATORX/ML_build.py/Spetree.py/PhyloNet_3.6.8.jar/astral.5.6.3.jar/AMAS.py in current directory'
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
