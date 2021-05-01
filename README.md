# Multi-species Regulatory neTwork LEarning (MRTLE) 

MRTLE takes expression data for multiple species, a phyloegentic tree relating the species, and other optional species-specific regulatory information, and learns a regulatory network for each species.

Koch, Christopher, et al. "Inference and evolutionary analysis of genome-scale regulatory networks in large phylogenies." Cell systems 4.5 (2017): 543-558.
https://doi.org/10.1016/j.cels.2017.04.010

## INSTALLATION
1) git clone https://github.com/Roy-lab/mrtle.git 
2) cd mrtle/code 
3) make


## EXAMPLE USAGE OF MRTLE
```./mrtle -f data/speciesconf.txt -x300 -v1 -m data/OGid_members.txt -l data/TFs/TFs_OGs.txt -n data/genelists/AllGenes.txt -d data/yeast_tree_rates.txt -s data/specorder_allclade.txt -b -0.9 -q 4.0```


The above example will run MRTLE using all regulators and targets. Since MRTLE learns regulators on a per-target basis, the algorithm can easily be parallelized by running the algorithm for each target orthogroup (or sets of target orthogroups) separately. For example, to run MRTLE using only OG59 (includes HOG1), we can replace the -n parameter with a file that contains only OG59 as such:

```./mrtle -f data/speciesconf.txt -x300 -v1 -l data/TFs/TFs_OGs.txt -n data/genelists/hog1.txt -d data/yeast_tree_rates.txt -m data/OGid_members.txt -s data/specorder_allclade.txt -b -0.9 -q 4.0```


## PARAMETER EXPLANATIONS
f : config file with one, six column row for each species. Each species' row should have the following species-specific entries:
1. Species Name

2. Filename of the expression data for each species

3. Location to place outputs

4. List of regulators to be used

5. List of target genes to be used

6. List of motifs to be used. This file should have three tab-separated columns, listing the regulator, target, and motif score

x : Maximum # of regulators to be used for a given target.

v : Used for cross-validation. Can be left at 1.

l : List of the orthogroups (id #s) to be considered as regulators. Note: a regulator must also be present in the species-specific list of regulators given in the species-specific config file (parameter f)

n : List of the orthogroups (id #s) to be considered as targets. Note: a target must also be present in the species-specific list of targets given in the species-specific config file (parameter f)

d : The species tree to be used. This file should have 5 columns describing the tree:

1. Child Species

2. Child relation to parent (left or right child)

3. Parent Species

4. Branch-specific gain rate

5. Branch-specific loss rate

m : A file describing the orthology relationships. The first column of this file is of the format OGID{NUMBER}_{DUP}. Each NUMBER represents an orthogroup. For orthogroups with duplications, DUP is the duplication count/id. If there are no duplications in the dataset being used, DUP will always be 1. Gene names must be unique across species	

s : A list of the species present in the OGIDS file (parameter m), in the order they exist in the OGIDS file

b : Beta1. Controls the sparsity of the network 

q : Beta2. Controls how strongly motifs are incorporated as prior. A higher value will result in motifs being valued more strongly
