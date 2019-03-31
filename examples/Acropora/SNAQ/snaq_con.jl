Pkg.add("PhyloNetworks")
#Pkg.update()
using PhyloNetworks
d=readTrees2CF("examples/Acropora/all_iqtree_btstraped.txt");
 
T=readTopology("examples/Acropora/ASTRAC/ASTRAC_output.txt");
 
net2=snaq!(T,d,hmax=3, filename="net2_snaq");