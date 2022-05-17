There are several different types of calculators used by mothur. All are stored in this grouping, and broken down into smaller groups for each type. Below is a breif description of the types and the commands that use them.


/********************************************************************/
Calculator class is parent for to all the otucalcs. 

OtuCalcs are used by:

collect.single
collect.shared
rarefaction.single
rarefaction.shared
summary.single
summary.shared
dist.shared
tree.shared
get.communitytype
heatmap.sim
venn

/********************************************************************/
ClusterMetric class is parent to all the clustercalcs. 
The clustermetrics are used in the opti method of clustering.

ClusterCalcs are used by:

cluster
cluster.split
cluster.fit
mgcluster

/********************************************************************/
DistCalc class is parent to all distcalcs. The distcalcs are used for finding the distance between sequences.

DistCalcs are used by:

dist.seqs
pairwise.seqs


/********************************************************************/
TreeCalculator is parent to the unifraccalcs.

unifraccalcs are used by:

parsimony
unifrac.weighted
unifrac.unweighted


/********************************************************************/

DiversityCalculator class is parent to all the diversity calculators.

The diversity calculators are used by the estimator.single command. 

https://github.com/chrisquince/DiversityEstimates 

https://www.ncbi.nlm.nih.gov/pubmed/18650928

/********************************************************************/
