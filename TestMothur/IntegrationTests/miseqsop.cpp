//
//  miseqsop.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/6/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "miseqsop.hpp"
#include "makecontigscommand.h"
#include "screenseqscommand.h"
#include "deconvolutecommand.h"
#include "countseqscommand.h"
#include "aligncommand.h"
#include "filterseqscommand.h"
#include "preclustercommand.h"
#include "setseedcommand.h"
#include "seqsummarycommand.h"
#include "chimeravsearchcommand.h"
#include "removeseqscommand.h"
#include "classifyseqscommand.h"
#include "removelineagecommand.h"
#include "getgroupscommand.h"
#include "seqerrorcommand.h"
#include "distancecommand.h"
#include "clustercommand.h"
#include "sharedcommand.h"
#include "rarefactcommand.h"
#include "removegroupscommand.h"
#include "classifyotucommand.h"
#include "phylotypecommand.h"
#include "clearcutcommand.h"
#include "renamefilecommand.h"
#include "countgroupscommand.h"
#include "subsamplecommand.h"
#include "summarycommand.h"
#include "heatmapcommand.h"
#include "matrixoutputcommand.h"
#include "heatmapsimcommand.h"
#include "venncommand.h"
#include "treegroupscommand.h"
#include "parsimonycommand.h"
#include "pcoacommand.h"
#include "nmdscommand.h"
#include "amovacommand.h"
#include "homovacommand.h"
#include "corraxescommand.h"
#include "getmetacommunitycommand.h"
#include "metastatscommand.h"
#include "lefsecommand.h"
#include "phylodiversitycommand.h"
#include "unifracunweightedcommand.h"
#include "unifracweightedcommand.h"
#include "setlogfilecommand.h"

/**************************************************************************************************/
TestMiSeqSOP::TestMiSeqSOP() {  //setup
    inputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    outputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    setDirInputs = "inputdir=" + inputDir + ", outputdir=" + outputDir;
    current = CurrentFile::getInstance();
    m = MothurOut::getInstance();
    
    //remove shortcut files
    Utils util;
    util.mothurRemove(inputDir+"silva.v4.8mer");
    util.mothurRemove(inputDir+"trainset9_032012.pds.8mer");
    util.mothurRemove(inputDir+"trainset9_032012.pds.trainset9_032012.pds.8mer.numNonZero");
    util.mothurRemove(inputDir+"trainset9_032012.pds.trainset9_032012.pds.8mer.prob");
    util.mothurRemove(inputDir+"trainset9_032012.pds.tree.sum");
    util.mothurRemove(inputDir+"trainset9_032012.pds.tree.train");

}
/**************************************************************************************************/
TestMiSeqSOP::~TestMiSeqSOP() {
}
/**************************************************************************************************/

TEST(Test_Integration_MISEQSOP, MISEQMakeContigs2Precluster) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.MISEQMakeContigs2Precluster.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //make.contigs
    inputs = test.setDirInputs + ", file=stability.files";
    test.m->mothurOut("Running command: make.contigs(" + inputs + ")\n");
    
    Command* makeContigsCommand = new MakeContigsCommand(inputs);
    makeContigsCommand->execute();
    delete makeContigsCommand;
    
    //screen.seqs
    inputs = test.setDirInputs + ", fasta=current, group=current, maxambig=0, maxlength=275";
    test.m->mothurOut("Running command: screen.seqs(" + inputs + ")\n");
    
    Command* screenSeqsCommand = new ScreenSeqsCommand(inputs);
    screenSeqsCommand->execute();
    delete screenSeqsCommand;
    
    //unique.seqs
    inputs = test.setDirInputs + ", fasta=current";
    test.m->mothurOut("Running command: unique.seqs(" + inputs + ")\n");
    
    Command* deconvoluteCommand = new DeconvoluteCommand(inputs);
    deconvoluteCommand->execute();
    delete deconvoluteCommand;

    //count.seqs
    inputs = test.setDirInputs + ", name=current, group=current";
    test.m->mothurOut("Running command: count.seqs(" + inputs + ")\n");
    
    Command* countSeqsCommand = new CountSeqsCommand(inputs);
    countSeqsCommand->execute();
    delete countSeqsCommand;
    
    //align.seqs
    inputs = test.setDirInputs + ", fasta=current, reference=silva.v4.fasta";
    test.m->mothurOut("Running command: align.seqs(" + inputs + ")\n");
    
    Command* alignSeqsCommand = new AlignCommand(inputs);
    alignSeqsCommand->execute();
    delete alignSeqsCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, count=current";
    test.m->mothurOut("Running command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;

    //screen.seqs
    inputs = test.setDirInputs + ", fasta=current, count=current, start=1968, end=11546, maxhomop=8";
    test.m->mothurOut("Running command: screen.seqs(" + inputs + ")\n");
    
    Command* screenSeqsCommand2 = new ScreenSeqsCommand(inputs);
    screenSeqsCommand2->execute();
    delete screenSeqsCommand2;
    
    //filter.seqs
    inputs = test.setDirInputs + ", fasta=current, vertical=T, trump=.";
    test.m->mothurOut("Running command: filter.seqs(" + inputs + ")\n");
    
    Command* filterSeqsCommand = new FilterSeqsCommand(inputs);
    filterSeqsCommand->execute();
    delete filterSeqsCommand;
    
    //unique.seqs
    inputs = test.setDirInputs + ", fasta=current, count=current";
    test.m->mothurOut("Running command: unique.seqs(" + inputs + ")\n");
    
    Command* deconvoluteCommand2 = new DeconvoluteCommand(inputs);
    deconvoluteCommand2->execute();
    delete deconvoluteCommand2;

    //pre.cluster
    inputs = test.setDirInputs + ", fasta=current, count=current, diffs=2";
    test.m->mothurOut("Running command: pre.cluster(" + inputs + ")\n");
    
    Command* preClusterCommand = new PreClusterCommand(inputs);
    preClusterCommand->execute();
    delete preClusterCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");

}

TEST(Test_Integration_MISEQSOP, RemoveChimerasContaminants) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.RemoveChimerasContaminants.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //chimera.vsearch
    inputs = test.setDirInputs + ", fasta=/Users/sarahwestcott/desktop/Miseq_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=/Users/sarahwestcott/desktop/Miseq_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t";
    test.m->mothurOut("Running command: chimera.vsearch(" + inputs + ")\n");
    
    Command* chimeraCommand = new ChimeraVsearchCommand(inputs);
    chimeraCommand->execute();
    delete chimeraCommand;

    //remove.seqs
    inputs = test.setDirInputs + ", fasta=current";
    test.m->mothurOut("Running command: remove.seqs(" + inputs + ")\n");
    
    Command* removeSeqsCommand = new RemoveSeqsCommand(inputs);
    removeSeqsCommand->execute();
    delete removeSeqsCommand;
    
    //classify.seqs(fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax
    inputs = test.setDirInputs + ", fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax";
    test.m->mothurOut("Running command: classify.seqs(" + inputs + ")\n");
    
    Command* classifySeqsCommand = new ClassifySeqsCommand(inputs);
    classifySeqsCommand->execute();
    delete classifySeqsCommand;
    
    //remove.lineage
    inputs = test.setDirInputs + ", fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota";
    test.m->mothurOut("Running command: remove.lineage(" + inputs + ")\n");
    
    Command* removeLineageCommand = new RemoveLineageCommand(inputs);
    removeLineageCommand->execute();
    delete removeLineageCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, count=current";
    test.m->mothurOut("Running command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_MISEQSOP, AssessingErrorRates) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.AssessingErrorRates.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //get.groups
    inputs = test.setDirInputs + ", fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, groups=Mock";
    test.m->mothurOut("Running command: get.groups(" + inputs + ")\n");
    
    Command* getGroupsCommand = new GetGroupsCommand(inputs);
    getGroupsCommand->execute();
    delete getGroupsCommand;
    
    //seq.error
    inputs = test.setDirInputs + ", fasta=current, count=current, reference=HMP_MOCK.v35.fasta, aligned=F";
    test.m->mothurOut("Running command: seq.error(" + inputs + ")\n");
    
    Command* seqErrorCommand = new SeqErrorCommand(inputs);
    seqErrorCommand->execute();
    delete seqErrorCommand;
    
    //dist.seqs
    inputs = test.setDirInputs + ", fasta=current, cutoff=0.03";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand = new DistanceCommand(inputs);
    distSeqsCommand->execute();
    delete distSeqsCommand;

    //cluster
    inputs = test.setDirInputs + ", column=current, count=current";
    test.m->mothurOut("Running command: cluster(" + inputs + ")\n");
    
    Command* clusterCommand = new ClusterCommand(inputs);
    clusterCommand->execute();
    delete clusterCommand;
    
    //make.shared
    inputs = test.setDirInputs + ", list=current, count=current, label=0.03";
    test.m->mothurOut("Running command: make.shared(" + inputs + ")\n");
    
    Command* makeSharedCommand = new SharedCommand(inputs);
    makeSharedCommand->execute();
    delete makeSharedCommand;

    //rarefaction.single
    inputs = test.setDirInputs + ", shared=current";
    test.m->mothurOut("Running command: rarefaction.single(" + inputs + ")\n");
    
    Command* rarefactionSingleCommand = new RareFactCommand(inputs);
    rarefactionSingleCommand->execute();
    delete rarefactionSingleCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_MISEQSOP, PreparingForAnalysis) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.PreparingForAnalysis.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //remove.groups
    inputs = test.setDirInputs + ", fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, groups=Mock";
    test.m->mothurOut("Running command: remove.groups(" + inputs + ")\n");
    
    Command* removeGroupsCommand = new RemoveGroupsCommand(inputs);
    removeGroupsCommand->execute();
    delete removeGroupsCommand;
    
    //dist.seqs
    inputs = test.setDirInputs + ", fasta=current, cutoff=0.03";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand = new DistanceCommand(inputs);
    distSeqsCommand->execute();
    delete distSeqsCommand;
    
    //cluster
    inputs = test.setDirInputs + ", column=current, count=current";
    test.m->mothurOut("Running command: cluster(" + inputs + ")\n");
    
    Command* clusterCommand = new ClusterCommand(inputs);
    clusterCommand->execute();
    delete clusterCommand;
    
    //make.shared
    inputs = test.setDirInputs + ", list=current, count=current, label=0.03";
    test.m->mothurOut("Running command: make.shared(" + inputs + ")\n");
    
    Command* makeSharedCommand = new SharedCommand(inputs);
    makeSharedCommand->execute();
    delete makeSharedCommand;
    
    //rename.file
    inputs = test.setDirInputs + ", column=current, fasta=current, count=current, taxonomy=current, list=current, shared=current, tree=current, prefix=final, deleteold=f";
    test.m->mothurOut("Running command: rename.file(" + inputs + ")\n");
    
    Command* renameFileCommand = new RenameFileCommand(inputs);
    renameFileCommand->execute();
    delete renameFileCommand;

    //classify.otu
    inputs = test.setDirInputs + ", list=current, count=current, taxonomy=current, label=0.03";
    test.m->mothurOut("Running command: classify.otu(" + inputs + ")\n");
    
    Command* classifyOtuCommand = new ClassifyOtuCommand(inputs);
    classifyOtuCommand->execute();
    delete classifyOtuCommand;
    
    //phylotype
    inputs = test.setDirInputs + ", taxonomy=current, label=1";
    test.m->mothurOut("Running command: classify.otu(" + inputs + ")\n");
    
    Command* phylotypeCommand = new PhylotypeCommand(inputs);
    phylotypeCommand->execute();
    delete phylotypeCommand;
    
    //make.shared
    inputs = test.setDirInputs + ", list=current, count=current";
    test.m->mothurOut("Running command: make.shared(" + inputs + ")\n");
    
    Command* makeSharedCommand2 = new SharedCommand(inputs);
    makeSharedCommand2->execute();
    delete makeSharedCommand2;
    
    //classify.otu
    inputs = test.setDirInputs + ", list=current, count=current, taxonomy=current";
    test.m->mothurOut("Running command: classify.otu(" + inputs + ")\n");
    
    Command* classifyOtuCommand2 = new ClassifyOtuCommand(inputs);
    classifyOtuCommand2->execute();
    delete classifyOtuCommand2;
    
    //dist.seqs
    inputs = test.setDirInputs + ", fasta=current, output=lt";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand2 = new DistanceCommand(inputs);
    distSeqsCommand2->execute();
    delete distSeqsCommand2;

    //clearcut
    inputs = test.setDirInputs + ", phylip=current";
    test.m->mothurOut("Running command: clearcut(" + inputs + ")\n");
    
    Command* clearcutCommand = new ClearcutCommand(inputs);
    clearcutCommand->execute();
    delete clearcutCommand;
    
    //count.groups
    inputs = test.setDirInputs + ", shared=current";
    test.m->mothurOut("Running command: count.groups(" + inputs + ")\n");
    
    Command* countGroupsCommand = new CountGroupsCommand(inputs);
    countGroupsCommand->execute();
    delete countGroupsCommand;
    
    //sub.sample
    inputs = test.setDirInputs + ", shared=current";
    test.m->mothurOut("Running command: sub.sample(" + inputs + ")\n");
    
    Command* subsampleCommand = new SubSampleCommand(inputs);
    subsampleCommand->execute();
    delete subsampleCommand;

    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_MISEQSOP, OTUBasedAnalysis) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.OTUBasedAnalysis.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //rarefaction.single
    inputs = test.setDirInputs + ", shared=final.opti_mcc.shared, calc=sobs, freq=100";
    test.m->mothurOut("Running command: rarefaction.single(" + inputs + ")\n");
    
    Command* rarefactionSingleCommand = new RareFactCommand(inputs);
    rarefactionSingleCommand->execute();
    delete rarefactionSingleCommand;
    
    //summary.single
    inputs = test.setDirInputs + ", shared=final.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T";
    test.m->mothurOut("Running command: summary.single(" + inputs + ")\n");
    
    Command* summarySingleCommand = new SummaryCommand(inputs);
    summarySingleCommand->execute();
    delete summarySingleCommand;
    
    //heatmap.bin
    inputs = test.setDirInputs + ", shared=final.opti_mcc.0.03.subsample.shared, scale=log2, numotu=50";
    test.m->mothurOut("Running command: heatmap.bin(" + inputs + ")\n");
    
    Command* heatmapCommand = new HeatMapCommand(inputs);
    heatmapCommand->execute();
    delete heatmapCommand;
    
    //venn
    inputs = test.setDirInputs + ", shared=final.opti_mcc.0.03.subsample.shared, groups=F3D0-F3D1-F3D2-F3D3";
    test.m->mothurOut("Running command: venn(" + inputs + ")\n");
    
    Command* vennCommand = new VennCommand(inputs);
    vennCommand->execute();
    delete vennCommand;
    
    //dist.shared
    inputs = test.setDirInputs + ", shared=final.opti_mcc.shared, calc=thetayc-jclass, subsample=T";
    test.m->mothurOut("Running command: dist.shared(" + inputs + ")\n");
    
    Command* distSharedCommand = new MatrixOutputCommand(inputs);
    distSharedCommand->execute();
    delete distSharedCommand;
    
    //heatmap.sim
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.jclass.0.03.lt.ave.dist";
    test.m->mothurOut("Running command: heatmap.sim(" + inputs + ")\n");
    
    Command* heatmapSimCommand = new HeatMapSimCommand(inputs);
    heatmapSimCommand->execute();
    delete heatmapSimCommand;
    
    //heatmap.sim
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist";
    test.m->mothurOut("Running command: heatmap.sim(" + inputs + ")\n");
    
    Command* heatmapSimCommand2 = new HeatMapSimCommand(inputs);
    heatmapSimCommand2->execute();
    delete heatmapSimCommand2;
    
    //tree.shared - phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current, calc=thetayc";
    test.m->mothurOut("Running command: tree.shared(" + inputs + ")\n");
    
    Command* treeSharedCommand = new TreeGroupCommand(inputs);
    treeSharedCommand->execute();
    delete treeSharedCommand;
    
    //parsimony - tree=stability.opti_mcc.thetayc.0.03.lt.ave.tre
    inputs = test.setDirInputs + ", tree=current, group=mouse.time.design";
    test.m->mothurOut("Running command: parsimony(" + inputs + ")\n");
    
    Command* parsimonyCommand = new ParsimonyCommand(inputs);
    parsimonyCommand->execute();
    delete parsimonyCommand;
    
    //pcoa - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current";
    test.m->mothurOut("Running command: pcoa(" + inputs + ")\n");
    
    Command* pcoaCommand = new PCOACommand(inputs);
    pcoaCommand->execute();
    delete pcoaCommand;
    
    //nmds - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current";
    test.m->mothurOut("Running command: nmds(" + inputs + ")\n");
    
    Command* nmdsCommand = new NMDSCommand(inputs);
    nmdsCommand->execute();
    delete nmdsCommand;
    
    //nmds - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current, mindim=3, maxdim=3";
    test.m->mothurOut("Running command: nmds(" + inputs + ")\n");
    
    Command* nmdsCommand2 = new NMDSCommand(inputs);
    nmdsCommand2->execute();
    delete nmdsCommand2;
    
    //amova - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current, design=mouse.time.design";
    test.m->mothurOut("Running command: amova(" + inputs + ")\n");
    
    Command* amovaCommand = new AmovaCommand(inputs);
    amovaCommand->execute();
    delete amovaCommand;
    
    //homova - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current, design=mouse.time.design";
    test.m->mothurOut("Running command: homova(" + inputs + ")\n");
    
    Command* homovaCommand = new HomovaCommand(inputs);
    homovaCommand->execute();
    delete homovaCommand;
    
    //corr.axes
    inputs = test.setDirInputs + ", axes=final.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, shared=final.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3";
    test.m->mothurOut("Running command: corr.axes(" + inputs + ")\n");
    
    Command* corrAxesCommand = new CorrAxesCommand(inputs);
    corrAxesCommand->execute();
    delete corrAxesCommand;
    
    //corr.axes(axes=stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3)
    inputs = test.setDirInputs + ", axes=final.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3";
    test.m->mothurOut("Running command: corr.axes(" + inputs + ")\n");
    
    Command* corrAxesCommand2 = new CorrAxesCommand(inputs);
    corrAxesCommand2->execute();
    delete corrAxesCommand2;
    
    //get.communitytype
    inputs = test.setDirInputs + ", shared=final.opti_mcc.0.03.subsample.shared";
    test.m->mothurOut("Running command: get.communitytype(" + inputs + ")\n");
    
    Command* getCommunityTypeCommand = new GetMetaCommunityCommand(inputs);
    getCommunityTypeCommand->execute();
    delete getCommunityTypeCommand;

    //metastats
    inputs = test.setDirInputs + ", shared=current, design=current";
    test.m->mothurOut("Running command: metastats(" + inputs + ")\n");
    
    Command* metastatsCommand = new MetaStatsCommand(inputs);
    metastatsCommand->execute();
    delete metastatsCommand;
    
    //lefse
    inputs = test.setDirInputs + ", shared=current, design=current";
    test.m->mothurOut("Running command: lefse(" + inputs + ")\n");
    
    Command* lefseCommand = new LefseCommand(inputs);
    lefseCommand->execute();
    delete lefseCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_MISEQSOP, PhylogenyBasedAnalysis) {
    TestMiSeqSOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration.PhylogenyBasedAnalysis.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;

    
    //phylo.diversity
    inputs = test.setDirInputs + ", tree=final.phylip.tre, count=final.count_table, rarefy=T";
    test.m->mothurOut("Running command: phylo.diversity(" + inputs + ")\n");
    
    Command* phyloDiversityCommand = new PhyloDiversityCommand(inputs);
    phyloDiversityCommand->execute();
    delete phyloDiversityCommand;
    
    //unifrac.unweighted
    inputs = test.setDirInputs + ", tree=current, count=current, distance=lt, random=F, subsample=T";
    test.m->mothurOut("Running command: unifrac.unweighted(" + inputs + ")\n");
    
    Command* unifracUnweightedCommand = new UnifracUnweightedCommand(inputs);
    unifracUnweightedCommand->execute();
    delete unifracUnweightedCommand;
    
    //unifrac.weighted
    inputs = test.setDirInputs + ", tree=current, count=current, distance=lt, random=F, subsample=T";
    test.m->mothurOut("Running command: unifrac.weighted(" + inputs + ")\n");
    
    Command* unifracWeightedCommand = new UnifracWeightedCommand(inputs);
    unifracWeightedCommand->execute();
    delete unifracWeightedCommand;

    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

/**************************************************************************************************/


