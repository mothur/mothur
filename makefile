###################################################
#
# Makefile for mothur
# Creator [Xcode -> Makefile Ver: Feb 14 2007 09:18:41]
# Created: April 16, 2010
#
###################################################

#
# Macros
#

CC = mpic++
CC_OPTIONS = -O3

# if you do not want to use the readline library set to no, default yes.
# make sure you have the library installed

USEREADLINE ?= yes

ifeq  ($(strip $(USEREADLINE)),yes)
    CC_OPTIONS += -DUSE_READLINE
	LNK_OPTIONS += \
      -lreadline\
      -lncurses\
      -L../readline-6.0
endif

USEMPI ?= yes


ifeq  ($(strip $(USEMPI)),yes)
    CC_OPTIONS += -DUSE_MPI
endif

#
# INCLUDE directories for mothur
#

     INCLUDE = -I.

#
# Build mothur
#

mothur : \
		./sharedutilities.o\
		./treegroupscommand.o\
		./bootstrapsharedcommand.o\
		./matrixoutputcommand.o\
		./getoturepcommand.o\
		./screenseqscommand.o\
		./chimera.o\
		./decalc.o\
		./readotucommand.o\
		./readdistcommand.o\
		./commandfactory.o\
		./alignment.o\
		./alignmentcell.o\
		./gotohoverlap.o\
		./overlap.o\
		./needlemanoverlap.o\
		./blastalign.o\
		./noalign.o\
		./suffixdb.o\
		./suffixnodes.o\
		./suffixtree.o\
		./blastdb.o\
		./nast.o\
		./nastreport.o\
		./boneh.o\
		./efron.o\
		./solow.o\
		./unifracweightedcommand.o\
		./weighted.o\
		./unweighted.o\
		./unifracunweightedcommand.o\
		./getsabundcommand.o\
		./getrabundcommand.o\
		./bellerophon.o\
		./pintail.o\
		./sharedanderbergs.o\
		./venncommand.o\
		./venn.o\
		./fullmatrix.o\
		./heatmap.o\
		./heatmapcommand.o\
		./libshuffcommand.o\
		./nocommands.o\
		./sharedbraycurtis.o\
		./sharedkulczynski.o\
		./sharedlennon.o\
		./sharedkulczynskicody.o\
		./sharedmorisitahorn.o\
		./sharedochiai.o\
		./readcolumn.o\
		./readotu.o\
		./readphylip.o\
		./consensuscommand.o\
		./heatmapsimcommand.o\
		./heatmapsim.o\
		./optionparser.o\
		./filterseqscommand.o\
		./goodscoverage.o\
		./sequencedb.o\
		./sharedjackknife.o\
		./sharedmarczewski.o\
		./aligncommand.o\
		./treemap.o\
		./parsimonycommand.o\
		./parsimony.o\
		./seqsummarycommand.o\
		./chimeraccodecommand.o\
		./chimerabellerophoncommand.o\
		./chimeracheckcommand.o\
		./chimeraslayercommand.o\
		./chimerapintailcommand.o\
		./chimeraseqscommand.o\
		./sharedlistvector.o\
		./tree.o\
		./readtree.o\
		./sharedsobscollectsummary.o\
		./deconvolutecommand.o\
		./listseqscommand.o\
		./getseqscommand.o\
		./removeseqscommand.o\
		./systemcommand.o\
		./binsequencecommand.o\
		./distancecommand.o\
		./ace.o\
		./averagelinkage.o\
		./bootstrap.o\
		./calculator.o\
		./chao1.o\
		./cluster.o\
		./clustercommand.o\
		./collect.o\
		./collectcommand.o\
		./collectsharedcommand.o\
		./commandoptionparser.o\
		./completelinkage.o\
		./database.o\
		./engine.o\
		./fastamap.o\
		./fileoutput.o\
		./globaldata.o\
		./groupmap.o\
		./helpcommand.o\
		./inputdata.o\
		./jackknife.o\
		./kmer.o\
		./kmerdb.o\
		./listvector.o\
		./mothur.o\
		./nameassignment.o\
		./npshannon.o\
		./ordervector.o\
		./progress.o\
		./quitcommand.o\
		./rabundvector.o\
		./rarecalc.o\
		./raredisplay.o\
		./rarefact.o\
		./rarefactcommand.o\
		./rarefactsharedcommand.o\
		./sabundvector.o\
		./sequence.o\
		./shannon.o\
		./sharedace.o\
		./sharedchao1.o\
		./sharedcommand.o\
		./sharedjabund.o\
		./sharedjclass.o\
		./sharedjest.o\
		./sharedordervector.o\
		./sharedrabundvector.o\
		./sharedsabundvector.o\
		./sharedsobs.o\
		./sharedsorabund.o\
		./sharedsorclass.o\
		./sharedsorest.o\
		./sharedthetan.o\
		./sharedthetayc.o\
		./simpson.o\
		./singlelinkage.o\
		./sparsematrix.o\
		./summarycommand.o\
		./summarysharedcommand.o\
		./uvest.o\
		./validcalculator.o\
		./validparameter.o\
		./treenode.o\
		./readtreecommand.o\
		./reversecommand.o\
		./trimseqscommand.o\
		./slibshuff.o\
		./libshuff.o\
		./dlibshuff.o\
		./mergefilecommand.o\
		./coverage.o\
		./whittaker.o\
		./preclustercommand.o\
		./otuhierarchycommand.o\
		./setdircommand.o\
		./getgroupcommand.o\
		./getlabelcommand.o\
		./secondarystructurecommand.o\
		./mothurout.o\
		./parselistscommand.o\
		./readblast.o\
		./chimeracheckrdp.o\
		./hclustercommand.o\
		./hcluster.o\
		./getlistcountcommand.o\
		./readcluster.o\
		./ccode.o\
		./taxonomyequalizer.o\
		./phylotypecommand.o\
		./classifyseqscommand.o\
		./parsesffcommand.o\
		./classify.o\
		./phylotree.o\
		./bayesian.o\
		./phylosummary.o\
		./alignmentdb.o\
		./knn.o\
		./distancedb.o\
		./chimeraslayer.o\
		./slayer.o\
		./pcacommand.o\
		./formatcolumn.o\
		./formatphylip.o\
		./mgclustercommand.o\
		./getsharedotucommand.o\
		./maligner.o\
		./chimerarealigner.o\
		./bergerparker.o\
		./bstick.o\
		./sharedkstest.o\
		./qstat.o\
		./shen.o\
		./logsd.o\
		./geom.o\
		./setlogfilecommand.o
	$(CC) $(LNK_OPTIONS) \
		./sharedutilities.o\
		./treegroupscommand.o\
		./bootstrapsharedcommand.o\
		./matrixoutputcommand.o\
		./getoturepcommand.o\
		./screenseqscommand.o\
		./chimera.o\
		./decalc.o\
		./readotucommand.o\
		./readdistcommand.o\
		./commandfactory.o\
		./alignment.o\
		./alignmentcell.o\
		./gotohoverlap.o\
		./overlap.o\
		./needlemanoverlap.o\
		./blastalign.o\
		./noalign.o\
		./suffixdb.o\
		./suffixnodes.o\
		./suffixtree.o\
		./blastdb.o\
		./nast.o\
		./nastreport.o\
		./boneh.o\
		./efron.o\
		./solow.o\
		./unifracweightedcommand.o\
		./weighted.o\
		./unweighted.o\
		./unifracunweightedcommand.o\
		./getsabundcommand.o\
		./getrabundcommand.o\
		./bellerophon.o\
		./pintail.o\
		./sharedanderbergs.o\
		./venncommand.o\
		./venn.o\
		./fullmatrix.o\
		./heatmap.o\
		./heatmapcommand.o\
		./libshuffcommand.o\
		./nocommands.o\
		./sharedbraycurtis.o\
		./sharedkulczynski.o\
		./sharedlennon.o\
		./sharedkulczynskicody.o\
		./sharedmorisitahorn.o\
		./sharedochiai.o\
		./readcolumn.o\
		./readotu.o\
		./readphylip.o\
		./consensuscommand.o\
		./heatmapsimcommand.o\
		./heatmapsim.o\
		./optionparser.o\
		./filterseqscommand.o\
		./goodscoverage.o\
		./sequencedb.o\
		./sharedjackknife.o\
		./sharedmarczewski.o\
		./aligncommand.o\
		./treemap.o\
		./parsimonycommand.o\
		./parsimony.o\
		./seqsummarycommand.o\
		./chimeraccodecommand.o\
		./chimerabellerophoncommand.o\
		./chimeracheckcommand.o\
		./chimeraslayercommand.o\
		./chimerapintailcommand.o\
		./chimeraseqscommand.o\
		./sharedlistvector.o\
		./tree.o\
		./readtree.o\
		./sharedsobscollectsummary.o\
		./deconvolutecommand.o\
		./listseqscommand.o\
		./getseqscommand.o\
		./removeseqscommand.o\
		./systemcommand.o\
		./binsequencecommand.o\
		./distancecommand.o\
		./ace.o\
		./averagelinkage.o\
		./bootstrap.o\
		./calculator.o\
		./chao1.o\
		./cluster.o\
		./clustercommand.o\
		./collect.o\
		./collectcommand.o\
		./collectsharedcommand.o\
		./commandoptionparser.o\
		./completelinkage.o\
		./database.o\
		./engine.o\
		./fastamap.o\
		./fileoutput.o\
		./globaldata.o\
		./groupmap.o\
		./helpcommand.o\
		./inputdata.o\
		./jackknife.o\
		./kmer.o\
		./kmerdb.o\
		./listvector.o\
		./mothur.o\
		./nameassignment.o\
		./npshannon.o\
		./ordervector.o\
		./progress.o\
		./quitcommand.o\
		./rabundvector.o\
		./rarecalc.o\
		./raredisplay.o\
		./rarefact.o\
		./rarefactcommand.o\
		./rarefactsharedcommand.o\
		./sabundvector.o\
		./sequence.o\
		./shannon.o\
		./sharedace.o\
		./sharedchao1.o\
		./sharedcommand.o\
		./sharedjabund.o\
		./sharedjclass.o\
		./sharedjest.o\
		./sharedordervector.o\
		./sharedrabundvector.o\
		./sharedsabundvector.o\
		./sharedsobs.o\
		./sharedsorabund.o\
		./sharedsorclass.o\
		./sharedsorest.o\
		./sharedthetan.o\
		./sharedthetayc.o\
		./simpson.o\
		./singlelinkage.o\
		./sparsematrix.o\
		./summarycommand.o\
		./summarysharedcommand.o\
		./uvest.o\
		./validcalculator.o\
		./validparameter.o\
		./treenode.o\
		./readtreecommand.o\
		./reversecommand.o\
		./trimseqscommand.o\
		./slibshuff.o\
		./libshuff.o\
		./dlibshuff.o\
		./mergefilecommand.o\
		./coverage.o\
		./whittaker.o\
		./preclustercommand.o\
		./otuhierarchycommand.o\
		./setdircommand.o\
		./getgroupcommand.o\
		./getlabelcommand.o\
		./secondarystructurecommand.o\
		./mothurout.o\
		./parselistscommand.o\
		./readblast.o\
		./chimeracheckrdp.o\
		./hclustercommand.o\
		./hcluster.o\
		./getlistcountcommand.o\
		./readcluster.o\
		./ccode.o\
		./taxonomyequalizer.o\
		./phylotypecommand.o\
		./classifyseqscommand.o\
		./parsesffcommand.o\
		./classify.o\
		./phylotree.o\
		./bayesian.o\
		./phylosummary.o\
		./alignmentdb.o\
		./knn.o\
		./distancedb.o\
		./chimeraslayer.o\
		./slayer.o\
		./pcacommand.o\
		./formatcolumn.o\
		./formatphylip.o\
		./mgclustercommand.o\
		./getsharedotucommand.o\
		./maligner.o\
		./chimerarealigner.o\
		./bergerparker.o\
		./bstick.o\
		./sharedkstest.o\
		./qstat.o\
		./shen.o\
		./logsd.o\
		./geom.o\
		./setlogfilecommand.o\
		-o ../Release/mothur

clean : 
		rm \
		./sharedutilities.o\
		./treegroupscommand.o\
		./bootstrapsharedcommand.o\
		./matrixoutputcommand.o\
		./getoturepcommand.o\
		./screenseqscommand.o\
		./chimera.o\
		./decalc.o\
		./readotucommand.o\
		./readdistcommand.o\
		./commandfactory.o\
		./alignment.o\
		./alignmentcell.o\
		./gotohoverlap.o\
		./overlap.o\
		./needlemanoverlap.o\
		./blastalign.o\
		./noalign.o\
		./suffixdb.o\
		./suffixnodes.o\
		./suffixtree.o\
		./blastdb.o\
		./nast.o\
		./nastreport.o\
		./boneh.o\
		./efron.o\
		./solow.o\
		./unifracweightedcommand.o\
		./weighted.o\
		./unweighted.o\
		./unifracunweightedcommand.o\
		./getsabundcommand.o\
		./getrabundcommand.o\
		./bellerophon.o\
		./pintail.o\
		./sharedanderbergs.o\
		./venncommand.o\
		./venn.o\
		./fullmatrix.o\
		./heatmap.o\
		./heatmapcommand.o\
		./libshuffcommand.o\
		./nocommands.o\
		./sharedbraycurtis.o\
		./sharedkulczynski.o\
		./sharedlennon.o\
		./sharedkulczynskicody.o\
		./sharedmorisitahorn.o\
		./sharedochiai.o\
		./readcolumn.o\
		./readotu.o\
		./readphylip.o\
		./consensuscommand.o\
		./heatmapsimcommand.o\
		./heatmapsim.o\
		./optionparser.o\
		./filterseqscommand.o\
		./goodscoverage.o\
		./sequencedb.o\
		./sharedjackknife.o\
		./sharedmarczewski.o\
		./aligncommand.o\
		./treemap.o\
		./parsimonycommand.o\
		./parsimony.o\
		./seqsummarycommand.o\
		./chimeraccodecommand.o\
		./chimerabellerophoncommand.o\
		./chimeracheckcommand.o\
		./chimeraslayercommand.o\
		./chimerapintailcommand.o\
		./chimeraseqscommand.o\
		./sharedlistvector.o\
		./tree.o\
		./readtree.o\
		./sharedsobscollectsummary.o\
		./deconvolutecommand.o\
		./listseqscommand.o\
		./getseqscommand.o\
		./removeseqscommand.o\
		./systemcommand.o\
		./binsequencecommand.o\
		./distancecommand.o\
		./ace.o\
		./averagelinkage.o\
		./bootstrap.o\
		./calculator.o\
		./chao1.o\
		./cluster.o\
		./clustercommand.o\
		./collect.o\
		./collectcommand.o\
		./collectsharedcommand.o\
		./commandoptionparser.o\
		./completelinkage.o\
		./database.o\
		./engine.o\
		./fastamap.o\
		./fileoutput.o\
		./globaldata.o\
		./groupmap.o\
		./helpcommand.o\
		./inputdata.o\
		./jackknife.o\
		./kmer.o\
		./kmerdb.o\
		./listvector.o\
		./mothur.o\
		./nameassignment.o\
		./npshannon.o\
		./ordervector.o\
		./progress.o\
		./quitcommand.o\
		./rabundvector.o\
		./rarecalc.o\
		./raredisplay.o\
		./rarefact.o\
		./rarefactcommand.o\
		./rarefactsharedcommand.o\
		./sabundvector.o\
		./sequence.o\
		./shannon.o\
		./sharedace.o\
		./sharedchao1.o\
		./sharedcommand.o\
		./sharedjabund.o\
		./sharedjclass.o\
		./sharedjest.o\
		./sharedordervector.o\
		./sharedrabundvector.o\
		./sharedsabundvector.o\
		./sharedsobs.o\
		./sharedsorabund.o\
		./sharedsorclass.o\
		./sharedsorest.o\
		./sharedthetan.o\
		./sharedthetayc.o\
		./simpson.o\
		./singlelinkage.o\
		./sparsematrix.o\
		./summarycommand.o\
		./summarysharedcommand.o\
		./uvest.o\
		./validcalculator.o\
		./validparameter.o\
		./treenode.o\
		./readtreecommand.o\
		./reversecommand.o\
		./trimseqscommand.o\
		./slibshuff.o\
		./libshuff.o\
		./dlibshuff.o\
		./mergefilecommand.o\
		./coverage.o\
		./whittaker.o\
		./preclustercommand.o\
		./otuhierarchycommand.o\
		./setdircommand.o\
		./getgroupcommand.o\
		./getlabelcommand.o\
		./secondarystructurecommand.o\
		./mothurout.o\
		./parselistscommand.o\
		./readblast.o\
		./chimeracheckrdp.o\
		./hclustercommand.o\
		./hcluster.o\
		./getlistcountcommand.o\
		./readcluster.o\
		./ccode.o\
		./taxonomyequalizer.o\
		./phylotypecommand.o\
		./classifyseqscommand.o\
		./parsesffcommand.o\
		./classify.o\
		./phylotree.o\
		./bayesian.o\
		./phylosummary.o\
		./alignmentdb.o\
		./knn.o\
		./distancedb.o\
		./chimeraslayer.o\
		./slayer.o\
		./pcacommand.o\
		./formatcolumn.o\
		./formatphylip.o\
		./mgclustercommand.o\
		./getsharedotucommand.o\
		./maligner.o\
		./chimerarealigner.o\
		./bergerparker.o\
		./bstick.o\
		./sharedkstest.o\
		./qstat.o\
		./shen.o\
		./logsd.o\
		./geom.o\
		./setlogfilecommand.o\
		mothur

install : mothur
		#cp mothur ../Release/mothur

#
# Build the parts of mothur
#


# Item # 1 -- sharedutilities --
./sharedutilities.o : sharedutilities.cpp
	$(CC) $(CC_OPTIONS) sharedutilities.cpp -c $(INCLUDE) -o ./sharedutilities.o


# Item # 2 -- treegroupscommand --
./treegroupscommand.o : treegroupscommand.cpp
	$(CC) $(CC_OPTIONS) treegroupscommand.cpp -c $(INCLUDE) -o ./treegroupscommand.o


# Item # 3 -- bootstrapsharedcommand --
./bootstrapsharedcommand.o : bootstrapsharedcommand.cpp
	$(CC) $(CC_OPTIONS) bootstrapsharedcommand.cpp -c $(INCLUDE) -o ./bootstrapsharedcommand.o


# Item # 4 -- matrixoutputcommand --
./matrixoutputcommand.o : matrixoutputcommand.cpp
	$(CC) $(CC_OPTIONS) matrixoutputcommand.cpp -c $(INCLUDE) -o ./matrixoutputcommand.o


# Item # 5 -- getoturepcommand --
./getoturepcommand.o : getoturepcommand.cpp
	$(CC) $(CC_OPTIONS) getoturepcommand.cpp -c $(INCLUDE) -o ./getoturepcommand.o


# Item # 6 -- screenseqscommand --
./screenseqscommand.o : screenseqscommand.cpp
	$(CC) $(CC_OPTIONS) screenseqscommand.cpp -c $(INCLUDE) -o ./screenseqscommand.o


# Item # 7 -- chimera --
./chimera.o : chimera.cpp
	$(CC) $(CC_OPTIONS) chimera.cpp -c $(INCLUDE) -o ./chimera.o


# Item # 8 -- decalc --
./decalc.o : decalc.cpp
	$(CC) $(CC_OPTIONS) decalc.cpp -c $(INCLUDE) -o ./decalc.o


# Item # 9 -- readotucommand --
./readotucommand.o : readotucommand.cpp
	$(CC) $(CC_OPTIONS) readotucommand.cpp -c $(INCLUDE) -o ./readotucommand.o


# Item # 10 -- readdistcommand --
./readdistcommand.o : readdistcommand.cpp
	$(CC) $(CC_OPTIONS) readdistcommand.cpp -c $(INCLUDE) -o ./readdistcommand.o


# Item # 11 -- commandfactory --
./commandfactory.o : commandfactory.cpp
	$(CC) $(CC_OPTIONS) commandfactory.cpp -c $(INCLUDE) -o ./commandfactory.o


# Item # 12 -- alignment --
./alignment.o : alignment.cpp
	$(CC) $(CC_OPTIONS) alignment.cpp -c $(INCLUDE) -o ./alignment.o


# Item # 13 -- alignmentcell --
./alignmentcell.o : alignmentcell.cpp
	$(CC) $(CC_OPTIONS) alignmentcell.cpp -c $(INCLUDE) -o ./alignmentcell.o


# Item # 14 -- gotohoverlap --
./gotohoverlap.o : gotohoverlap.cpp
	$(CC) $(CC_OPTIONS) gotohoverlap.cpp -c $(INCLUDE) -o ./gotohoverlap.o


# Item # 15 -- overlap --
./overlap.o : overlap.cpp
	$(CC) $(CC_OPTIONS) overlap.cpp -c $(INCLUDE) -o ./overlap.o


# Item # 16 -- needlemanoverlap --
./needlemanoverlap.o : needlemanoverlap.cpp
	$(CC) $(CC_OPTIONS) needlemanoverlap.cpp -c $(INCLUDE) -o ./needlemanoverlap.o


# Item # 17 -- blastalign --
./blastalign.o : blastalign.cpp
	$(CC) $(CC_OPTIONS) blastalign.cpp -c $(INCLUDE) -o ./blastalign.o


# Item # 18 -- noalign --
./noalign.o : noalign.cpp
	$(CC) $(CC_OPTIONS) noalign.cpp -c $(INCLUDE) -o ./noalign.o


# Item # 19 -- suffixdb --
./suffixdb.o : suffixdb.cpp
	$(CC) $(CC_OPTIONS) suffixdb.cpp -c $(INCLUDE) -o ./suffixdb.o


# Item # 20 -- suffixnodes --
./suffixnodes.o : suffixnodes.cpp
	$(CC) $(CC_OPTIONS) suffixnodes.cpp -c $(INCLUDE) -o ./suffixnodes.o


# Item # 21 -- suffixtree --
./suffixtree.o : suffixtree.cpp
	$(CC) $(CC_OPTIONS) suffixtree.cpp -c $(INCLUDE) -o ./suffixtree.o


# Item # 22 -- blastdb --
./blastdb.o : blastdb.cpp
	$(CC) $(CC_OPTIONS) blastdb.cpp -c $(INCLUDE) -o ./blastdb.o


# Item # 23 -- nast --
./nast.o : nast.cpp
	$(CC) $(CC_OPTIONS) nast.cpp -c $(INCLUDE) -o ./nast.o


# Item # 24 -- nastreport --
./nastreport.o : nastreport.cpp
	$(CC) $(CC_OPTIONS) nastreport.cpp -c $(INCLUDE) -o ./nastreport.o


# Item # 25 -- boneh --
./boneh.o : boneh.cpp
	$(CC) $(CC_OPTIONS) boneh.cpp -c $(INCLUDE) -o ./boneh.o


# Item # 26 -- efron --
./efron.o : efron.cpp
	$(CC) $(CC_OPTIONS) efron.cpp -c $(INCLUDE) -o ./efron.o


# Item # 27 -- solow --
./solow.o : solow.cpp
	$(CC) $(CC_OPTIONS) solow.cpp -c $(INCLUDE) -o ./solow.o


# Item # 28 -- unifracweightedcommand --
./unifracweightedcommand.o : unifracweightedcommand.cpp
	$(CC) $(CC_OPTIONS) unifracweightedcommand.cpp -c $(INCLUDE) -o ./unifracweightedcommand.o


# Item # 29 -- weighted --
./weighted.o : weighted.cpp
	$(CC) $(CC_OPTIONS) weighted.cpp -c $(INCLUDE) -o ./weighted.o


# Item # 30 -- unweighted --
./unweighted.o : unweighted.cpp
	$(CC) $(CC_OPTIONS) unweighted.cpp -c $(INCLUDE) -o ./unweighted.o


# Item # 31 -- unifracunweightedcommand --
./unifracunweightedcommand.o : unifracunweightedcommand.cpp
	$(CC) $(CC_OPTIONS) unifracunweightedcommand.cpp -c $(INCLUDE) -o ./unifracunweightedcommand.o


# Item # 32 -- getsabundcommand --
./getsabundcommand.o : getsabundcommand.cpp
	$(CC) $(CC_OPTIONS) getsabundcommand.cpp -c $(INCLUDE) -o ./getsabundcommand.o


# Item # 33 -- getrabundcommand --
./getrabundcommand.o : getrabundcommand.cpp
	$(CC) $(CC_OPTIONS) getrabundcommand.cpp -c $(INCLUDE) -o ./getrabundcommand.o


# Item # 34 -- bellerophon --
./bellerophon.o : bellerophon.cpp
	$(CC) $(CC_OPTIONS) bellerophon.cpp -c $(INCLUDE) -o ./bellerophon.o


# Item # 35 -- pintail --
./pintail.o : pintail.cpp
	$(CC) $(CC_OPTIONS) pintail.cpp -c $(INCLUDE) -o ./pintail.o


# Item # 36 -- sharedanderbergs --
./sharedanderbergs.o : sharedanderbergs.cpp
	$(CC) $(CC_OPTIONS) sharedanderbergs.cpp -c $(INCLUDE) -o ./sharedanderbergs.o


# Item # 37 -- venncommand --
./venncommand.o : venncommand.cpp
	$(CC) $(CC_OPTIONS) venncommand.cpp -c $(INCLUDE) -o ./venncommand.o


# Item # 38 -- venn --
./venn.o : venn.cpp
	$(CC) $(CC_OPTIONS) venn.cpp -c $(INCLUDE) -o ./venn.o


# Item # 39 -- fullmatrix --
./fullmatrix.o : fullmatrix.cpp
	$(CC) $(CC_OPTIONS) fullmatrix.cpp -c $(INCLUDE) -o ./fullmatrix.o


# Item # 40 -- heatmap --
./heatmap.o : heatmap.cpp
	$(CC) $(CC_OPTIONS) heatmap.cpp -c $(INCLUDE) -o ./heatmap.o


# Item # 41 -- heatmapcommand --
./heatmapcommand.o : heatmapcommand.cpp
	$(CC) $(CC_OPTIONS) heatmapcommand.cpp -c $(INCLUDE) -o ./heatmapcommand.o


# Item # 42 -- libshuffcommand --
./libshuffcommand.o : libshuffcommand.cpp
	$(CC) $(CC_OPTIONS) libshuffcommand.cpp -c $(INCLUDE) -o ./libshuffcommand.o


# Item # 43 -- nocommands --
./nocommands.o : nocommands.cpp
	$(CC) $(CC_OPTIONS) nocommands.cpp -c $(INCLUDE) -o ./nocommands.o


# Item # 44 -- sharedbraycurtis --
./sharedbraycurtis.o : sharedbraycurtis.cpp
	$(CC) $(CC_OPTIONS) sharedbraycurtis.cpp -c $(INCLUDE) -o ./sharedbraycurtis.o


# Item # 45 -- sharedkulczynski --
./sharedkulczynski.o : sharedkulczynski.cpp
	$(CC) $(CC_OPTIONS) sharedkulczynski.cpp -c $(INCLUDE) -o ./sharedkulczynski.o


# Item # 46 -- sharedlennon --
./sharedlennon.o : sharedlennon.cpp
	$(CC) $(CC_OPTIONS) sharedlennon.cpp -c $(INCLUDE) -o ./sharedlennon.o


# Item # 47 -- sharedkulczynskicody --
./sharedkulczynskicody.o : sharedkulczynskicody.cpp
	$(CC) $(CC_OPTIONS) sharedkulczynskicody.cpp -c $(INCLUDE) -o ./sharedkulczynskicody.o


# Item # 48 -- sharedmorisitahorn --
./sharedmorisitahorn.o : sharedmorisitahorn.cpp
	$(CC) $(CC_OPTIONS) sharedmorisitahorn.cpp -c $(INCLUDE) -o ./sharedmorisitahorn.o


# Item # 49 -- sharedochiai --
./sharedochiai.o : sharedochiai.cpp
	$(CC) $(CC_OPTIONS) sharedochiai.cpp -c $(INCLUDE) -o ./sharedochiai.o


# Item # 50 -- readcolumn --
./readcolumn.o : readcolumn.cpp
	$(CC) $(CC_OPTIONS) readcolumn.cpp -c $(INCLUDE) -o ./readcolumn.o


# Item # 51 -- readotu --
./readotu.o : readotu.cpp
	$(CC) $(CC_OPTIONS) readotu.cpp -c $(INCLUDE) -o ./readotu.o


# Item # 52 -- readphylip --
./readphylip.o : readphylip.cpp
	$(CC) $(CC_OPTIONS) readphylip.cpp -c $(INCLUDE) -o ./readphylip.o


# Item # 53 -- consensuscommand --
./consensuscommand.o : consensuscommand.cpp
	$(CC) $(CC_OPTIONS) consensuscommand.cpp -c $(INCLUDE) -o ./consensuscommand.o


# Item # 54 -- heatmapsimcommand --
./heatmapsimcommand.o : heatmapsimcommand.cpp
	$(CC) $(CC_OPTIONS) heatmapsimcommand.cpp -c $(INCLUDE) -o ./heatmapsimcommand.o


# Item # 55 -- heatmapsim --
./heatmapsim.o : heatmapsim.cpp
	$(CC) $(CC_OPTIONS) heatmapsim.cpp -c $(INCLUDE) -o ./heatmapsim.o


# Item # 56 -- optionparser --
./optionparser.o : optionparser.cpp
	$(CC) $(CC_OPTIONS) optionparser.cpp -c $(INCLUDE) -o ./optionparser.o


# Item # 57 -- filterseqscommand --
./filterseqscommand.o : filterseqscommand.cpp
	$(CC) $(CC_OPTIONS) filterseqscommand.cpp -c $(INCLUDE) -o ./filterseqscommand.o


# Item # 58 -- goodscoverage --
./goodscoverage.o : goodscoverage.cpp
	$(CC) $(CC_OPTIONS) goodscoverage.cpp -c $(INCLUDE) -o ./goodscoverage.o


# Item # 59 -- sequencedb --
./sequencedb.o : sequencedb.cpp
	$(CC) $(CC_OPTIONS) sequencedb.cpp -c $(INCLUDE) -o ./sequencedb.o


# Item # 60 -- sharedjackknife --
./sharedjackknife.o : sharedjackknife.cpp
	$(CC) $(CC_OPTIONS) sharedjackknife.cpp -c $(INCLUDE) -o ./sharedjackknife.o


# Item # 61 -- sharedmarczewski --
./sharedmarczewski.o : sharedmarczewski.cpp
	$(CC) $(CC_OPTIONS) sharedmarczewski.cpp -c $(INCLUDE) -o ./sharedmarczewski.o


# Item # 62 -- aligncommand --
./aligncommand.o : aligncommand.cpp
	$(CC) $(CC_OPTIONS) aligncommand.cpp -c $(INCLUDE) -o ./aligncommand.o


# Item # 63 -- treemap --
./treemap.o : treemap.cpp
	$(CC) $(CC_OPTIONS) treemap.cpp -c $(INCLUDE) -o ./treemap.o


# Item # 64 -- parsimonycommand --
./parsimonycommand.o : parsimonycommand.cpp
	$(CC) $(CC_OPTIONS) parsimonycommand.cpp -c $(INCLUDE) -o ./parsimonycommand.o


# Item # 65 -- parsimony --
./parsimony.o : parsimony.cpp
	$(CC) $(CC_OPTIONS) parsimony.cpp -c $(INCLUDE) -o ./parsimony.o


# Item # 66 -- seqsummarycommand --
./seqsummarycommand.o : seqsummarycommand.cpp
	$(CC) $(CC_OPTIONS) seqsummarycommand.cpp -c $(INCLUDE) -o ./seqsummarycommand.o


# Item # 67 -- chimeraseqscommand --
./chimeraseqscommand.o : chimeraseqscommand.cpp
	$(CC) $(CC_OPTIONS) chimeraseqscommand.cpp -c $(INCLUDE) -o ./chimeraseqscommand.o


# Item # 68 -- sharedlistvector --
./sharedlistvector.o : sharedlistvector.cpp
	$(CC) $(CC_OPTIONS) sharedlistvector.cpp -c $(INCLUDE) -o ./sharedlistvector.o


# Item # 69 -- tree --
./tree.o : tree.cpp
	$(CC) $(CC_OPTIONS) tree.cpp -c $(INCLUDE) -o ./tree.o


# Item # 70 -- readtree --
./readtree.o : readtree.cpp
	$(CC) $(CC_OPTIONS) readtree.cpp -c $(INCLUDE) -o ./readtree.o


# Item # 71 -- sharedsobscollectsummary --
./sharedsobscollectsummary.o : sharedsobscollectsummary.cpp
	$(CC) $(CC_OPTIONS) sharedsobscollectsummary.cpp -c $(INCLUDE) -o ./sharedsobscollectsummary.o


# Item # 72 -- deconvolutecommand --
./deconvolutecommand.o : deconvolutecommand.cpp
	$(CC) $(CC_OPTIONS) deconvolutecommand.cpp -c $(INCLUDE) -o ./deconvolutecommand.o


# Item # 73 -- listseqscommand --
./listseqscommand.o : listseqscommand.cpp
	$(CC) $(CC_OPTIONS) listseqscommand.cpp -c $(INCLUDE) -o ./listseqscommand.o


# Item # 74 -- getseqscommand --
./getseqscommand.o : getseqscommand.cpp
	$(CC) $(CC_OPTIONS) getseqscommand.cpp -c $(INCLUDE) -o ./getseqscommand.o


# Item # 75 -- removeseqscommand --
./removeseqscommand.o : removeseqscommand.cpp
	$(CC) $(CC_OPTIONS) removeseqscommand.cpp -c $(INCLUDE) -o ./removeseqscommand.o


# Item # 76 -- systemcommand --
./systemcommand.o : systemcommand.cpp
	$(CC) $(CC_OPTIONS) systemcommand.cpp -c $(INCLUDE) -o ./systemcommand.o


# Item # 77 -- binsequencecommand --
./binsequencecommand.o : binsequencecommand.cpp
	$(CC) $(CC_OPTIONS) binsequencecommand.cpp -c $(INCLUDE) -o ./binsequencecommand.o


# Item # 78 -- distancecommand --
./distancecommand.o : distancecommand.cpp
	$(CC) $(CC_OPTIONS) distancecommand.cpp -c $(INCLUDE) -o ./distancecommand.o


# Item # 79 -- ace --
./ace.o : ace.cpp
	$(CC) $(CC_OPTIONS) ace.cpp -c $(INCLUDE) -o ./ace.o


# Item # 80 -- averagelinkage --
./averagelinkage.o : averagelinkage.cpp
	$(CC) $(CC_OPTIONS) averagelinkage.cpp -c $(INCLUDE) -o ./averagelinkage.o


# Item # 81 -- bootstrap --
./bootstrap.o : bootstrap.cpp
	$(CC) $(CC_OPTIONS) bootstrap.cpp -c $(INCLUDE) -o ./bootstrap.o


# Item # 82 -- calculator --
./calculator.o : calculator.cpp
	$(CC) $(CC_OPTIONS) calculator.cpp -c $(INCLUDE) -o ./calculator.o


# Item # 83 -- chao1 --
./chao1.o : chao1.cpp
	$(CC) $(CC_OPTIONS) chao1.cpp -c $(INCLUDE) -o ./chao1.o


# Item # 84 -- cluster --
./cluster.o : cluster.cpp
	$(CC) $(CC_OPTIONS) cluster.cpp -c $(INCLUDE) -o ./cluster.o


# Item # 85 -- clustercommand --
./clustercommand.o : clustercommand.cpp
	$(CC) $(CC_OPTIONS) clustercommand.cpp -c $(INCLUDE) -o ./clustercommand.o


# Item # 86 -- collect --
./collect.o : collect.cpp
	$(CC) $(CC_OPTIONS) collect.cpp -c $(INCLUDE) -o ./collect.o


# Item # 87 -- collectcommand --
./collectcommand.o : collectcommand.cpp
	$(CC) $(CC_OPTIONS) collectcommand.cpp -c $(INCLUDE) -o ./collectcommand.o


# Item # 88 -- collectsharedcommand --
./collectsharedcommand.o : collectsharedcommand.cpp
	$(CC) $(CC_OPTIONS) collectsharedcommand.cpp -c $(INCLUDE) -o ./collectsharedcommand.o


# Item # 89 -- commandoptionparser --
./commandoptionparser.o : commandoptionparser.cpp
	$(CC) $(CC_OPTIONS) commandoptionparser.cpp -c $(INCLUDE) -o ./commandoptionparser.o


# Item # 90 -- completelinkage --
./completelinkage.o : completelinkage.cpp
	$(CC) $(CC_OPTIONS) completelinkage.cpp -c $(INCLUDE) -o ./completelinkage.o


# Item # 91 -- database --
./database.o : database.cpp
	$(CC) $(CC_OPTIONS) database.cpp -c $(INCLUDE) -o ./database.o


# Item # 92 -- engine --
./engine.o : engine.cpp
	$(CC) $(CC_OPTIONS) engine.cpp -c $(INCLUDE) -o ./engine.o


# Item # 93 -- fastamap --
./fastamap.o : fastamap.cpp
	$(CC) $(CC_OPTIONS) fastamap.cpp -c $(INCLUDE) -o ./fastamap.o


# Item # 94 -- fileoutput --
./fileoutput.o : fileoutput.cpp
	$(CC) $(CC_OPTIONS) fileoutput.cpp -c $(INCLUDE) -o ./fileoutput.o


# Item # 95 -- globaldata --
./globaldata.o : globaldata.cpp
	$(CC) $(CC_OPTIONS) globaldata.cpp -c $(INCLUDE) -o ./globaldata.o


# Item # 96 -- groupmap --
./groupmap.o : groupmap.cpp
	$(CC) $(CC_OPTIONS) groupmap.cpp -c $(INCLUDE) -o ./groupmap.o


# Item # 97 -- helpcommand --
./helpcommand.o : helpcommand.cpp
	$(CC) $(CC_OPTIONS) helpcommand.cpp -c $(INCLUDE) -o ./helpcommand.o


# Item # 98 -- inputdata --
./inputdata.o : inputdata.cpp
	$(CC) $(CC_OPTIONS) inputdata.cpp -c $(INCLUDE) -o ./inputdata.o


# Item # 99 -- jackknife --
./jackknife.o : jackknife.cpp
	$(CC) $(CC_OPTIONS) jackknife.cpp -c $(INCLUDE) -o ./jackknife.o


# Item # 100 -- kmer --
./kmer.o : kmer.cpp
	$(CC) $(CC_OPTIONS) kmer.cpp -c $(INCLUDE) -o ./kmer.o


# Item # 101 -- kmerdb --
./kmerdb.o : kmerdb.cpp
	$(CC) $(CC_OPTIONS) kmerdb.cpp -c $(INCLUDE) -o ./kmerdb.o


# Item # 102 -- listvector --
./listvector.o : listvector.cpp
	$(CC) $(CC_OPTIONS) listvector.cpp -c $(INCLUDE) -o ./listvector.o


# Item # 103 -- mothur --
./mothur.o : mothur.cpp
	$(CC) $(CC_OPTIONS) mothur.cpp -c $(INCLUDE) -o ./mothur.o


# Item # 104 -- nameassignment --
./nameassignment.o : nameassignment.cpp
	$(CC) $(CC_OPTIONS) nameassignment.cpp -c $(INCLUDE) -o ./nameassignment.o


# Item # 105 -- npshannon --
./npshannon.o : npshannon.cpp
	$(CC) $(CC_OPTIONS) npshannon.cpp -c $(INCLUDE) -o ./npshannon.o


# Item # 106 -- ordervector --
./ordervector.o : ordervector.cpp
	$(CC) $(CC_OPTIONS) ordervector.cpp -c $(INCLUDE) -o ./ordervector.o


# Item # 107 -- progress --
./progress.o : progress.cpp
	$(CC) $(CC_OPTIONS) progress.cpp -c $(INCLUDE) -o ./progress.o


# Item # 108 -- quitcommand --
./quitcommand.o : quitcommand.cpp
	$(CC) $(CC_OPTIONS) quitcommand.cpp -c $(INCLUDE) -o ./quitcommand.o


# Item # 109 -- rabundvector --
./rabundvector.o : rabundvector.cpp
	$(CC) $(CC_OPTIONS) rabundvector.cpp -c $(INCLUDE) -o ./rabundvector.o


# Item # 110 -- rarecalc --
./rarecalc.o : rarecalc.cpp
	$(CC) $(CC_OPTIONS) rarecalc.cpp -c $(INCLUDE) -o ./rarecalc.o


# Item # 111 -- raredisplay --
./raredisplay.o : raredisplay.cpp
	$(CC) $(CC_OPTIONS) raredisplay.cpp -c $(INCLUDE) -o ./raredisplay.o


# Item # 112 -- rarefact --
./rarefact.o : rarefact.cpp
	$(CC) $(CC_OPTIONS) rarefact.cpp -c $(INCLUDE) -o ./rarefact.o


# Item # 113 -- rarefactcommand --
./rarefactcommand.o : rarefactcommand.cpp
	$(CC) $(CC_OPTIONS) rarefactcommand.cpp -c $(INCLUDE) -o ./rarefactcommand.o


# Item # 114 -- rarefactsharedcommand --
./rarefactsharedcommand.o : rarefactsharedcommand.cpp
	$(CC) $(CC_OPTIONS) rarefactsharedcommand.cpp -c $(INCLUDE) -o ./rarefactsharedcommand.o


# Item # 115 -- sabundvector --
./sabundvector.o : sabundvector.cpp
	$(CC) $(CC_OPTIONS) sabundvector.cpp -c $(INCLUDE) -o ./sabundvector.o


# Item # 116 -- sequence --
./sequence.o : sequence.cpp
	$(CC) $(CC_OPTIONS) sequence.cpp -c $(INCLUDE) -o ./sequence.o


# Item # 117 -- shannon --
./shannon.o : shannon.cpp
	$(CC) $(CC_OPTIONS) shannon.cpp -c $(INCLUDE) -o ./shannon.o


# Item # 118 -- sharedace --
./sharedace.o : sharedace.cpp
	$(CC) $(CC_OPTIONS) sharedace.cpp -c $(INCLUDE) -o ./sharedace.o


# Item # 119 -- sharedchao1 --
./sharedchao1.o : sharedchao1.cpp
	$(CC) $(CC_OPTIONS) sharedchao1.cpp -c $(INCLUDE) -o ./sharedchao1.o


# Item # 120 -- sharedcommand --
./sharedcommand.o : sharedcommand.cpp
	$(CC) $(CC_OPTIONS) sharedcommand.cpp -c $(INCLUDE) -o ./sharedcommand.o


# Item # 121 -- sharedjabund --
./sharedjabund.o : sharedjabund.cpp
	$(CC) $(CC_OPTIONS) sharedjabund.cpp -c $(INCLUDE) -o ./sharedjabund.o


# Item # 122 -- sharedjclass --
./sharedjclass.o : sharedjclass.cpp
	$(CC) $(CC_OPTIONS) sharedjclass.cpp -c $(INCLUDE) -o ./sharedjclass.o


# Item # 123 -- sharedjest --
./sharedjest.o : sharedjest.cpp
	$(CC) $(CC_OPTIONS) sharedjest.cpp -c $(INCLUDE) -o ./sharedjest.o


# Item # 124 -- sharedordervector --
./sharedordervector.o : sharedordervector.cpp
	$(CC) $(CC_OPTIONS) sharedordervector.cpp -c $(INCLUDE) -o ./sharedordervector.o


# Item # 125 -- sharedrabundvector --
./sharedrabundvector.o : sharedrabundvector.cpp
	$(CC) $(CC_OPTIONS) sharedrabundvector.cpp -c $(INCLUDE) -o ./sharedrabundvector.o


# Item # 126 -- sharedsabundvector --
./sharedsabundvector.o : sharedsabundvector.cpp
	$(CC) $(CC_OPTIONS) sharedsabundvector.cpp -c $(INCLUDE) -o ./sharedsabundvector.o


# Item # 127 -- sharedsobs --
./sharedsobs.o : sharedsobs.cpp
	$(CC) $(CC_OPTIONS) sharedsobs.cpp -c $(INCLUDE) -o ./sharedsobs.o


# Item # 128 -- sharedsorabund --
./sharedsorabund.o : sharedsorabund.cpp
	$(CC) $(CC_OPTIONS) sharedsorabund.cpp -c $(INCLUDE) -o ./sharedsorabund.o


# Item # 129 -- sharedsorclass --
./sharedsorclass.o : sharedsorclass.cpp
	$(CC) $(CC_OPTIONS) sharedsorclass.cpp -c $(INCLUDE) -o ./sharedsorclass.o


# Item # 130 -- sharedsorest --
./sharedsorest.o : sharedsorest.cpp
	$(CC) $(CC_OPTIONS) sharedsorest.cpp -c $(INCLUDE) -o ./sharedsorest.o


# Item # 131 -- sharedthetan --
./sharedthetan.o : sharedthetan.cpp
	$(CC) $(CC_OPTIONS) sharedthetan.cpp -c $(INCLUDE) -o ./sharedthetan.o


# Item # 132 -- sharedthetayc --
./sharedthetayc.o : sharedthetayc.cpp
	$(CC) $(CC_OPTIONS) sharedthetayc.cpp -c $(INCLUDE) -o ./sharedthetayc.o


# Item # 133 -- simpson --
./simpson.o : simpson.cpp
	$(CC) $(CC_OPTIONS) simpson.cpp -c $(INCLUDE) -o ./simpson.o


# Item # 134 -- singlelinkage --
./singlelinkage.o : singlelinkage.cpp
	$(CC) $(CC_OPTIONS) singlelinkage.cpp -c $(INCLUDE) -o ./singlelinkage.o


# Item # 135 -- sparsematrix --
./sparsematrix.o : sparsematrix.cpp
	$(CC) $(CC_OPTIONS) sparsematrix.cpp -c $(INCLUDE) -o ./sparsematrix.o


# Item # 136 -- summarycommand --
./summarycommand.o : summarycommand.cpp
	$(CC) $(CC_OPTIONS) summarycommand.cpp -c $(INCLUDE) -o ./summarycommand.o


# Item # 137 -- summarysharedcommand --
./summarysharedcommand.o : summarysharedcommand.cpp
	$(CC) $(CC_OPTIONS) summarysharedcommand.cpp -c $(INCLUDE) -o ./summarysharedcommand.o


# Item # 138 -- uvest --
./uvest.o : uvest.cpp
	$(CC) $(CC_OPTIONS) uvest.cpp -c $(INCLUDE) -o ./uvest.o


# Item # 139 -- validcalculator --
./validcalculator.o : validcalculator.cpp
	$(CC) $(CC_OPTIONS) validcalculator.cpp -c $(INCLUDE) -o ./validcalculator.o


# Item # 140 -- validparameter --
./validparameter.o : validparameter.cpp
	$(CC) $(CC_OPTIONS) validparameter.cpp -c $(INCLUDE) -o ./validparameter.o


# Item # 141 -- treenode --
./treenode.o : treenode.cpp
	$(CC) $(CC_OPTIONS) treenode.cpp -c $(INCLUDE) -o ./treenode.o


# Item # 142 -- readtreecommand --
./readtreecommand.o : readtreecommand.cpp
	$(CC) $(CC_OPTIONS) readtreecommand.cpp -c $(INCLUDE) -o ./readtreecommand.o


# Item # 143 -- reversecommand --
./reversecommand.o : reversecommand.cpp
	$(CC) $(CC_OPTIONS) reversecommand.cpp -c $(INCLUDE) -o ./reversecommand.o


# Item # 144 -- trimseqscommand --
./trimseqscommand.o : trimseqscommand.cpp
	$(CC) $(CC_OPTIONS) trimseqscommand.cpp -c $(INCLUDE) -o ./trimseqscommand.o


# Item # 145 -- slibshuff --
./slibshuff.o : slibshuff.cpp
	$(CC) $(CC_OPTIONS) slibshuff.cpp -c $(INCLUDE) -o ./slibshuff.o


# Item # 146 -- libshuff --
./libshuff.o : libshuff.cpp
	$(CC) $(CC_OPTIONS) libshuff.cpp -c $(INCLUDE) -o ./libshuff.o


# Item # 147 -- dlibshuff --
./dlibshuff.o : dlibshuff.cpp
	$(CC) $(CC_OPTIONS) dlibshuff.cpp -c $(INCLUDE) -o ./dlibshuff.o


# Item # 148 -- mergefilecommand --
./mergefilecommand.o : mergefilecommand.cpp
	$(CC) $(CC_OPTIONS) mergefilecommand.cpp -c $(INCLUDE) -o ./mergefilecommand.o


# Item # 149 -- coverage --
./coverage.o : coverage.cpp
	$(CC) $(CC_OPTIONS) coverage.cpp -c $(INCLUDE) -o ./coverage.o


# Item # 150 -- whittaker --
./whittaker.o : whittaker.cpp
	$(CC) $(CC_OPTIONS) whittaker.cpp -c $(INCLUDE) -o ./whittaker.o


# Item # 151 -- preclustercommand --
./preclustercommand.o : preclustercommand.cpp
	$(CC) $(CC_OPTIONS) preclustercommand.cpp -c $(INCLUDE) -o ./preclustercommand.o


# Item # 152 -- otuhierarchycommand --
./otuhierarchycommand.o : otuhierarchycommand.cpp
	$(CC) $(CC_OPTIONS) otuhierarchycommand.cpp -c $(INCLUDE) -o ./otuhierarchycommand.o


# Item # 153 -- setdircommand --
./setdircommand.o : setdircommand.cpp
	$(CC) $(CC_OPTIONS) setdircommand.cpp -c $(INCLUDE) -o ./setdircommand.o


# Item # 154 -- getgroupcommand --
./getgroupcommand.o : getgroupcommand.cpp
	$(CC) $(CC_OPTIONS) getgroupcommand.cpp -c $(INCLUDE) -o ./getgroupcommand.o


# Item # 155 -- getlabelcommand --
./getlabelcommand.o : getlabelcommand.cpp
	$(CC) $(CC_OPTIONS) getlabelcommand.cpp -c $(INCLUDE) -o ./getlabelcommand.o


# Item # 156 -- secondarystructurecommand --
./secondarystructurecommand.o : secondarystructurecommand.cpp
	$(CC) $(CC_OPTIONS) secondarystructurecommand.cpp -c $(INCLUDE) -o ./secondarystructurecommand.o


# Item # 157 -- mothurout --
./mothurout.o : mothurout.cpp
	$(CC) $(CC_OPTIONS) mothurout.cpp -c $(INCLUDE) -o ./mothurout.o


# Item # 158 -- parselistscommand --
./parselistscommand.o : parselistscommand.cpp
	$(CC) $(CC_OPTIONS) parselistscommand.cpp -c $(INCLUDE) -o ./parselistscommand.o


# Item # 159 -- readblast --
./readblast.o : readblast.cpp
	$(CC) $(CC_OPTIONS) readblast.cpp -c $(INCLUDE) -o ./readblast.o


# Item # 160 -- chimeracheckrdp --
./chimeracheckrdp.o : chimeracheckrdp.cpp
	$(CC) $(CC_OPTIONS) chimeracheckrdp.cpp -c $(INCLUDE) -o ./chimeracheckrdp.o


# Item # 161 -- hclustercommand --
./hclustercommand.o : hclustercommand.cpp
	$(CC) $(CC_OPTIONS) hclustercommand.cpp -c $(INCLUDE) -o ./hclustercommand.o


# Item # 162 -- hcluster --
./hcluster.o : hcluster.cpp
	$(CC) $(CC_OPTIONS) hcluster.cpp -c $(INCLUDE) -o ./hcluster.o


# Item # 163 -- getlistcountcommand --
./getlistcountcommand.o : getlistcountcommand.cpp
	$(CC) $(CC_OPTIONS) getlistcountcommand.cpp -c $(INCLUDE) -o ./getlistcountcommand.o


# Item # 164 -- readcluster --
./readcluster.o : readcluster.cpp
	$(CC) $(CC_OPTIONS) readcluster.cpp -c $(INCLUDE) -o ./readcluster.o


# Item # 165 -- ccode --
./ccode.o : ccode.cpp
	$(CC) $(CC_OPTIONS) ccode.cpp -c $(INCLUDE) -o ./ccode.o


# Item # 166 -- taxonomyequalizer --
./taxonomyequalizer.o : taxonomyequalizer.cpp
	$(CC) $(CC_OPTIONS) taxonomyequalizer.cpp -c $(INCLUDE) -o ./taxonomyequalizer.o


# Item # 167 -- phylotypecommand --
./phylotypecommand.o : phylotypecommand.cpp
	$(CC) $(CC_OPTIONS) phylotypecommand.cpp -c $(INCLUDE) -o ./phylotypecommand.o


# Item # 168 -- classifyseqscommand --
./classifyseqscommand.o : classifyseqscommand.cpp
	$(CC) $(CC_OPTIONS) classifyseqscommand.cpp -c $(INCLUDE) -o ./classifyseqscommand.o


# Item # 169 -- classify --
./classify.o : classify.cpp
	$(CC) $(CC_OPTIONS) classify.cpp -c $(INCLUDE) -o ./classify.o


# Item # 170 -- phylotree --
./phylotree.o : phylotree.cpp
	$(CC) $(CC_OPTIONS) phylotree.cpp -c $(INCLUDE) -o ./phylotree.o


# Item # 171 -- bayesian --
./bayesian.o : bayesian.cpp
	$(CC) $(CC_OPTIONS) bayesian.cpp -c $(INCLUDE) -o ./bayesian.o


# Item # 172 -- alignmentdb --
./alignmentdb.o : alignmentdb.cpp
	$(CC) $(CC_OPTIONS) alignmentdb.cpp -c $(INCLUDE) -o ./alignmentdb.o


# Item # 173 -- knn --
./knn.o : knn.cpp
	$(CC) $(CC_OPTIONS) knn.cpp -c $(INCLUDE) -o ./knn.o


# Item # 174 -- distancedb --
./distancedb.o : distancedb.cpp
	$(CC) $(CC_OPTIONS) distancedb.cpp -c $(INCLUDE) -o ./distancedb.o


# Item # 175 -- chimeraslayer --
./chimeraslayer.o : chimeraslayer.cpp
	$(CC) $(CC_OPTIONS) chimeraslayer.cpp -c $(INCLUDE) -o ./chimeraslayer.o


# Item # 176 -- slayer --
./slayer.o : slayer.cpp
	$(CC) $(CC_OPTIONS) slayer.cpp -c $(INCLUDE) -o ./slayer.o


# Item # 177 -- pcacommand --
./pcacommand.o : pcacommand.cpp
	$(CC) $(CC_OPTIONS) pcacommand.cpp -c $(INCLUDE) -o ./pcacommand.o


# Item # 178 -- formatcolumn --
./formatcolumn.o : formatcolumn.cpp
	$(CC) $(CC_OPTIONS) formatcolumn.cpp -c $(INCLUDE) -o ./formatcolumn.o


# Item # 179 -- formatphylip --
./formatphylip.o : formatphylip.cpp
	$(CC) $(CC_OPTIONS) formatphylip.cpp -c $(INCLUDE) -o ./formatphylip.o


# Item # 180 -- mgclustercommand --
./mgclustercommand.o : mgclustercommand.cpp
	$(CC) $(CC_OPTIONS) mgclustercommand.cpp -c $(INCLUDE) -o ./mgclustercommand.o


# Item # 181 -- getsharedotucommand --
./getsharedotucommand.o : getsharedotucommand.cpp
	$(CC) $(CC_OPTIONS) getsharedotucommand.cpp -c $(INCLUDE) -o ./getsharedotucommand.o


# Item # 182 -- maligner --
./maligner.o : maligner.cpp
	$(CC) $(CC_OPTIONS) maligner.cpp -c $(INCLUDE) -o ./maligner.o


# Item # 183 -- chimerarealigner --
./chimerarealigner.o : chimerarealigner.cpp
	$(CC) $(CC_OPTIONS) chimerarealigner.cpp -c $(INCLUDE) -o ./chimerarealigner.o


# Item # 184 -- bergerparker --
./bergerparker.o : bergerparker.cpp
	$(CC) $(CC_OPTIONS) bergerparker.cpp -c $(INCLUDE) -o ./bergerparker.o


# Item # 185 -- bstick --
./bstick.o : bstick.cpp
	$(CC) $(CC_OPTIONS) bstick.cpp -c $(INCLUDE) -o ./bstick.o


# Item # 186 -- sharedkstest --
./sharedkstest.o : sharedkstest.cpp
	$(CC) $(CC_OPTIONS) sharedkstest.cpp -c $(INCLUDE) -o ./sharedkstest.o


# Item # 187 -- qstat --
./qstat.o : qstat.cpp
	$(CC) $(CC_OPTIONS) qstat.cpp -c $(INCLUDE) -o ./qstat.o


# Item # 188 -- shen --
./shen.o : shen.cpp
	$(CC) $(CC_OPTIONS) shen.cpp -c $(INCLUDE) -o ./shen.o


# Item # 189 -- logsd --
./logsd.o : logsd.cpp
	$(CC) $(CC_OPTIONS) logsd.cpp -c $(INCLUDE) -o ./logsd.o


# Item # 190 -- geom --
./geom.o : geom.cpp
	$(CC) $(CC_OPTIONS) geom.cpp -c $(INCLUDE) -o ./geom.o

# Item # 191 -- parsesffcommand --
./parsesffcommand.o : parsesffcommand.cpp
	$(CC) $(CC_OPTIONS) parsesffcommand.cpp -c $(INCLUDE) -o ./parsesffcommand.o

# Item # 192 -- chimeraccodecommand --
./chimeraccodecommand.o : chimeraccodecommand.cpp
	$(CC) $(CC_OPTIONS) chimeraccodecommand.cpp -c $(INCLUDE) -o ./chimeraccodecommand.o

# Item # 193 -- chimeracheckcommand --
./chimeracheckcommand.o : chimeracheckcommand.cpp
	$(CC) $(CC_OPTIONS) chimeracheckcommand.cpp -c $(INCLUDE) -o ./chimeracheckcommand.o


# Item # 194 -- chimeraslayercommand --
./chimeraslayercommand.o : chimeraslayercommand.cpp
	$(CC) $(CC_OPTIONS) chimeraslayercommand.cpp -c $(INCLUDE) -o ./chimeraslayercommand.o

# Item # 195 -- chimerapintailcommand --
./chimerapintailcommand.o : chimerapintailcommand.cpp
	$(CC) $(CC_OPTIONS) chimerapintailcommand.cpp -c $(INCLUDE) -o ./chimerapintailcommand.o

# Item # 196 -- chimerabellerophoncommand --
./chimerabellerophoncommand.o : chimerabellerophoncommand.cpp
	$(CC) $(CC_OPTIONS) chimerabellerophoncommand.cpp -c $(INCLUDE) -o ./chimerabellerophoncommand.o

# Item # 197 -- phylosummary --
./phylosummary.o : phylosummary.cpp
	$(CC) $(CC_OPTIONS) phylosummary.cpp -c $(INCLUDE) -o ./phylosummary.o

# Item # 198 -- setlogfilecommand --
./setlogfilecommand.o : setlogfilecommand.cpp
	$(CC) $(CC_OPTIONS) setlogfilecommand.cpp -c $(INCLUDE) -o ./setlogfilecommand.o



##### END RUN ####
