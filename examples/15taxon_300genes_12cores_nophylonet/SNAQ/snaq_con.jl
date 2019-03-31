Pkg.add("PhyloNetworks")
#Pkg.update()
using PhyloNetworks
d=readTrees2CF("examples/15taxon_300genes_4cores_nophylonet/all_iqtree_btstraped.txt");
 
T=readTopology("examples/15taxon_300genes_4cores_nophylonet/ASTRAL/ASTRAL_output.txt");
 
net2=snaq!(T,d,hmax=3, filename="net2_snaq");