Pkg.add("PhyloNetworks")
#Pkg.update()
using PhyloNetworks
d=readTrees2CF("Yeast_2408_clade2/all_iqtree_btstraped.txt");
 
T=readTopology("Yeast_2408_clade2/ASTRAL/ASTRAL_output.txt");
 
net2=snaq!(T,d,hmax=3, filename="net2_snaq");