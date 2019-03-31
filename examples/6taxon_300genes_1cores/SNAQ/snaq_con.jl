Pkg.add("PhyloNetworks")
#Pkg.update()
using PhyloNetworks
d=readTrees2CF("examples/6taxon_300genes_1cores/all_iqtree_btstraped.txt");
 
T=readTopology("examples/6taxon_300genes_1cores/ASTRAC/ASTRAC_output.txt");
 
net2=snaq!(T,d,hmax=3, filename="net2_snaq");