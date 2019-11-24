import os
from Bio import Phylo
from ete3 import Tree
from random import random, randint
import sys, getopt
import glob
import numpy as np


def bootstrap_check(tree, value):
    for node in tree.iter_search_nodes():
        if node.support > 1 and node.support < value:
            return False
    return True


def spe_namechange(namefile, tree):
    infile = open(namefile, 'r')
    n_spe = 0
    control_mpest[0] = ''
    taxonmap = ""
    for line in infile:
        dic_name = {}
        namelist = line.strip().split()
        control_mpest[0] += ' '.join(namelist) + '\n'
        taxonmap += namelist[0] + ':' + ','.join(namelist[2:]) + ';'
        n_spe += 1
        dic_name[namelist[0]] = namelist[2:]
        if int(namelist[1]) == 1:
            tree = tree.replace(namelist[2], namelist[0])
        else:
            for name in dic_name[namelist[0]]:
                tree = tree.replace(name, namelist[0])
    infile.close()
    taxonmap_phylonet[0] = taxonmap
    control_mpest[1] = n_spe
    return tree


def gene_namechange(pwd, gene_namefile):
    os.system('rm %sall_iqtree_spename.contree' % (pwd))
    infile = open(gene_namefile, 'r')
    dic_name = {}
    for line in infile:
        namelist = line.strip().split()
        dic_name[namelist[1]] = namelist[0]
    infile.close()
    infile = open(pwd + "all_iqtree.contree", 'r')
    outfile = open(pwd + 'all_iqtree_spename.contree', 'a')
    for line in infile:
        tree = line.strip()
        for key in dic_name:
            if key in tree:
                tree = tree.replace(key, dic_name[key])
        outfile.write(tree + '\n')
    outfile.close()
    infile.close()
    return dic_name


def run_Spetree(pwd, species_namefile, gene_namefile, boot_value, roottaxon, Net_num, cross_value, thread_number,
                message_queue=None):
    os.system('cat %saln_seqs/*.contree > %sall_iqtree.contree' % (pwd, pwd))

    infile = open(pwd + "all_iqtree.contree", 'r')

    ###change gene name#####
    global control_mpest
    global taxonmap_phylonet
    control_mpest = ['', 0, 0]
    n_gene = 0
    taxonmap_phylonet = ['']
    dic_name = {}
    dic_name = gene_namechange(pwd, gene_namefile)
    infile = open(pwd + 'all_iqtree_spename.contree', 'r')

    outfile_btstraped = open(pwd + 'all_iqtree_btstraped.txt', 'w')  # Input For ASTRAL, SNAQ
    outfile_rooted = open(pwd + 'all_iqtree_rooted.txt', 'w')  # Input For MPEST,PhyloNet
    outfile_namechange_non_branch = open(pwd + 'all_iqtree_namechange_nonbranch.txt', 'w')  # Input For STELLS2

    ###check bootstrap and rootted tree###
    for line in infile:
        t = Tree(line)
        if bootstrap_check(t, boot_value) == False:
            continue
        n_gene += 1
        s = t.write(format=3)
        s = s.replace('NoName', '')
        outfile_btstraped.write(s + '\n')
        print line
        if "," in roottaxon:
            root_taxon = roottaxon.strip().split(',')
            root_taxon = t.get_common_ancestor(root_taxon)
            try:
                t.set_outgroup(root_taxon)
            except:
                continue
        else:
            t.set_outgroup(t & roottaxon)
        s = t.write(format=3)
        s = s.replace('NoName', '')
        outfile_rooted.write(s + '\n')
        s = t.write(format=8)
        s = s.replace('NoName', '')
        outfile_namechange_non_branch.write(spe_namechange(species_namefile, s) + '\n')
    control_mpest[2] = n_gene
    outfile_btstraped.close()
    outfile_rooted.close()
    outfile_namechange_non_branch.close()

    ######inputs for species tree inferring######

    ####Concat RUNNING####
    random_seed_number = randint(0, 2**32)
    for key in dic_name:
        os.system(
            """ sed -i "s/%s/%s/g" `grep "%s" -rl %saln_seqs/*.aln` """ % (">" + key, ">" + dic_name[key], ">" + key, pwd))

    os.system("rm -rf %sCONCAT" % (pwd))
    os.system("mkdir %sCONCAT" % (pwd))
    os.system("python AMAS.py concat -i %saln_seqs/*.aln -f fasta -d dna" % (pwd))
    os.system("mv concatenated.out %sCONCAT/" % (pwd))
    os.system("iqtree -s %s -bb 1000 -redo -nt %d -m MFP -seed %d 1>%s_ML_iqtree.log" % (
    pwd + "CONCAT/concatenated.out", thread_number, random_seed_number, pwd + "CONCAT/concatenated.out"))

    if message_queue:
        message_queue.put("Raxml running...")

    ####Raxml RUNNING####
    os.system("raxml-ng --all --msa %s  --model GTR+G+FO --tree pars{10} --bs-trees autoMRE{1000} --seed %d --threads %d --prefix %s_raxml_ng --bs-cutoff 0.01" %(pwd+"CONCAT/concatenated.out", random_seed_number,thread_number,pwd + "CONCAT/concatenated.out"))

    if message_queue:
        message_queue.put("ASTRAL running...")

    ###ASTRAL RUNNING###
    os.system('rm -rf %sASTRAL' % (pwd))
    os.system('mkdir %sASTRAL' % (pwd))
    outfile = open(pwd + 'ASTRAL/species_name_ASTRAL.txt', 'w')
    outfile.write(control_mpest[0])
    outfile.close()
    os.system(
        'java -jar astral.5.6.3.jar -i %s -o %sASTRAL/%s -a %sASTRAL/species_name_ASTRAL.txt 2> %sASTRAL/run_ASTRAL.log' \
        % (pwd + "all_iqtree_btstraped.txt", pwd, "ASTRAL_output.txt", pwd, pwd))
    if message_queue:
        message_queue.put("ASTRAL done")
        message_queue.put("MP_EST running...")

    ###MP_EST RUNNING###

    os.system('rm -rf %sMP_EST' % (pwd))
    os.system('mkdir %sMP_EST' % (pwd))
    outfile = open(pwd + 'MP_EST/control.file', 'w')
    outfile.write('%s\n0\n%s\n5\n%s %s\n%s0' % (
    pwd + 'all_iqtree_rooted.txt', str(int(random() * 10000000)), str(control_mpest[2]), str(control_mpest[1]),
    control_mpest[0]))
    outfile.close()
    os.system('mpest %sMP_EST/control.file 1> %sMP_EST/run_MPEST.log' % (pwd, pwd))
    os.system('mv %sall_iqtree_rooted.txt_* %sMP_EST/' % (pwd, pwd))

    Phylo.convert(pwd+'MP_EST/all_iqtree_rooted.txt_besttree.tre', 'nexus', pwd+'MP_EST/all_iqtree_rooted.txt_besttree.tmp.nex', 'newick')
    os.system('head -n1 %sMP_EST/all_iqtree_rooted.txt_besttree.tmp.nex >%sMP_EST/all_iqtree_rooted.txt_besttree.nex' % (pwd, pwd)) 
    if message_queue:
        message_queue.put("MP_EST done")
        message_queue.put("STELLS2 running...")

    ###STELLS2 RUNNING###

    os.system('rm -rf %sSTELLS2' % (pwd))
    os.system('mkdir %sSTELLS2' % (pwd))
    os.system('stells-v2 -t %d -g %sall_iqtree_namechange_nonbranch.txt > %sSTELLS2/STELLS2_output.txt' % (
    thread_number, pwd, pwd))
    os.system('mv %sall_iqtree_namechange_nonbranch.txt-nearopt.trees %sSTELLS2/' % (pwd, pwd))
    os.system('grep "the inferred MLE species tree" %sSTELLS2/STELLS2_output.txt | cut -f 2- -d ":" | sed "s/$/;/" > %sSTELLS2/STELLS2_output_tree.txt' %(pwd, pwd))

    if message_queue:
        message_queue.put("STELLS2 done")
        message_queue.put("SNAQ running...")

    ###TREE COMPARE###
    tree_files = [pwd+'CONCAT/concatenated.out.contree', pwd+'CONCAT/concatenated.out_raxml_ng.raxml.bestTree', pwd+'ASTRAL/ASTRAL_output.txt',pwd+'MP_EST/all_iqtree_rooted.txt_besttree.nex', pwd+'STELLS2/STELLS2_output_tree.txt']
    trees = map(Tree, tree_files)
    n_trees = len(tree_files)
    mat = np.zeros((n_trees, n_trees), dtype=np.float)
    for i in range(n_trees):
        for j in range(n_trees):
            mat[i][j] = (trees[i].robinson_foulds(trees[j], unrooted_trees=True))[0]
    np.savetxt(pwd+'tree_compare.csv', mat, delimiter=",", fmt='%.3e')

    ###SNAQ RUNNING###

    os.system('rm -rf %sSNAQ' % (pwd))
    os.system('mkdir %sSNAQ' % (pwd))
    outfile = open(pwd + 'SNAQ/snaq_con.jl', 'w')
    outfile.write("""using Pkg\nPkg.add("PhyloNetworks")\n#Pkg.update()\nusing PhyloNetworks\nd=readTrees2CF("%sall_iqtree_btstraped.txt");\n 
    T=readTopology("%sASTRAL/ASTRAL_output.txt");\n 
    net2=snaq!(T,d,hmax=%d, filename="net2_snaq");""" % (pwd, pwd, Net_num))
    outfile.close()
    os.system('julia %sSNAQ/snaq_con.jl 1> %sSNAQ/run_SNAQ.log' % (pwd, pwd))
    os.system('mv net2* %sSNAQ/' % (pwd))
    os.system('mv summaryTreesQuartets.txt %sSNAQ/' % (pwd))
    os.system('mv tableCF.txt %sSNAQ/' % (pwd))
    if message_queue:
        message_queue.put("SNAQ done")
        message_queue.put("PHYLONET running...")


    ###PHYLONET RUNNING###

    os.system('rm -rf %sPHYLONET' % (pwd))
    os.system('mkdir %sPHYLONET' % (pwd))
    outfile = open(pwd + 'PHYLONET/phylonet_con.txt', 'w')
    outfile.write('#NEXUS\n\nBEGIN TREES;\n\n')
    infile = open(pwd + 'all_iqtree_rooted.txt', 'r')
    n = 0
    for line in infile:
        outfile.write('TREE gt%d = %s' % (n, line))
        n += 1
    infile.close()
    print n
    outfile.write('\nEND;\n\nBEGIN PHYLONET;\nInferNetwork_ML_CV (all) %d -pl %d -cv %d -a <%s>; \n\nEND;' % (
    Net_num, thread_number, cross_value, taxonmap_phylonet[0][:-1]))
    outfile.close()
    os.system('java -jar PhyloNet_3.6.8.jar %sPHYLONET/phylonet_con.txt 1>%sPHYLONET/PHYLONET_output.txt' % (pwd, pwd))
    if message_queue:
        message_queue.put("PHYLONET done")
        message_queue.put("ALL DONE")


if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "hs:g:b:r:n:k:t:d:")
    boot_value = 50
    Net_num = 5

    for op, value in opts:
        if op == "-b":
            boot_value = int(value)
        elif op == "-d":
            pwd = value + '/'
        elif op == "-r":
            roottaxon = value
        elif op == "-s":
            species_namefile = value
        elif op == "-g":
            gene_namefile = value
        elif op == "-n":
            Net_num = int(value)
        elif op == "-k":
            cross_value = int(value)
        elif op == "-t":
            thread_number = int(value)

    run_Spetree(pwd, species_namefile, gene_namefile, boot_value, roottaxon, Net_num, cross_value, thread_number)

