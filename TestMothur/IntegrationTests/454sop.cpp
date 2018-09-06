//
//  454sop.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/17/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "454sop.hpp"
#include "sffinfocommand.h"
#include "trimflowscommand.h"
#include "trimseqscommand.h"
#include "shhhercommand.h"
#include "screenseqscommand.h"
#include "deconvolutecommand.h"
#include "countseqscommand.h"
#include "aligncommand.h"
#include "filterseqscommand.h"
#include "preclustercommand.h"
#include "setseedcommand.h"
#include "seqsummarycommand.h"
#include "chimerauchimecommand.h"
#include "removeseqscommand.h"
#include "classifyseqscommand.h"
#include "removelineagecommand.h"
#include "getgroupscommand.h"
#include "getseqscommand.h"
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
#include "systemcommand.h"
#include "summarycommand.h"
#include "collectcommand.h"
#include "summarysharedcommand.h"

/**************************************************************************************************/
Test454SOP::Test454SOP() {  //setup
    inputDir = "/Users/sarahwestcott/desktop/454_sop/";
    outputDir = "/Users/sarahwestcott/desktop/454_sop/";
    setDirInputs = "inputdir=" + inputDir + ", outputdir=" + outputDir;
    current = CurrentFile::getInstance();
    m = MothurOut::getInstance();
    
    //remove shortcut files
    Utils util;
    util.mothurRemove(inputDir+"silva.bacteria.8mer");
    util.mothurRemove(inputDir+"trainset9_032012.pds.8mer");
    util.mothurRemove(inputDir+"trainset9_032012.pds.trainset9_032012.pds.8mer.numNonZero");
    util.mothurRemove(inputDir+"trainset9_032012.pds.trainset9_032012.pds.8mer.prob");
    util.mothurRemove(inputDir+"trainset9_032012.pds.tree.sum");
    util.mothurRemove(inputDir+"trainset9_032012.pds.tree.train");
    
}
/**************************************************************************************************/
Test454SOP::~Test454SOP() {}
/**************************************************************************************************/

TEST(Test_Integration_454SOP, SFF) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.SFF.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //sffinfo
    inputs = test.setDirInputs + ", sff=GQY1XT001.sff, flow=t";
    test.m->mothurOut("Running command: sff.info(" + inputs + ")\n");
    
    Command* sffCommand = new SffInfoCommand(inputs);
    sffCommand->execute();
    delete sffCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current";
    test.m->mothurOut("\nRunning command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, ReducingSequencingErrors) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.ReducingSequencingErrors.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //trim.flows
    inputs = test.setDirInputs + ", flow=GQY1XT001.flow, oligos=GQY1XT001.oligos, pdiffs=2, bdiffs=1";
    test.m->mothurOut("Running command: trim.flows(" + inputs + ")\n");
    
    Command* trimFlowsCommand = new TrimFlowsCommand(inputs);
    trimFlowsCommand->execute();
    delete trimFlowsCommand;
    
    //shhh.flows
    inputs = test.setDirInputs + ", file=GQY1XT001.flow.files";
    test.m->mothurOut("\nRunning command: shhh.flows(" + inputs + ")\n");
    
    Command* shhhFlowsCommand = new ShhherCommand(inputs);
    shhhFlowsCommand->execute();
    delete shhhFlowsCommand;
    
    //trim.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current, oligos=current, pdiffs=2, bdiffs=1, maxhomop=8, minlength=200, flip=T";
    test.m->mothurOut("Running command: trim.seqs(" + inputs + ")\n");
    
    Command* trimSeqsCommand = new TrimSeqsCommand(inputs);
    trimSeqsCommand->execute();
    delete trimSeqsCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current";
    test.m->mothurOut("\nRunning command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    //trim.seqs
    inputs = test.setDirInputs + ", fasta=current, qfile=current, oligos=current, maxambig=0, maxhomop=8, flip=T, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50";
    test.m->mothurOut("Running command: trim.seqs(" + inputs + ")\n");
    
    Command* trimSeqsCommand2 = new TrimSeqsCommand(inputs);
    trimSeqsCommand2->execute();
    delete trimSeqsCommand2;
    
    //trim.seqs
    inputs = test.setDirInputs + ", fasta=current, oligos=current, name=current, maxambig=0, maxhomop=8, flip=T, bdiffs=1, pdiffs=2, keepfirst=200";
    test.m->mothurOut("Running command: trim.seqs(" + inputs + ")\n");
    
    Command* trimSeqsCommand3 = new TrimSeqsCommand(inputs);
    trimSeqsCommand3->execute();
    delete trimSeqsCommand3;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current";
    test.m->mothurOut("\nRunning command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand2 = new SeqSummaryCommand(inputs);
    summarySeqsCommand2->execute();
    delete summarySeqsCommand2;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, ProcessingImprovedSequences) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.ProcessingImprovedSequences.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //unique.seqs
    inputs = test.setDirInputs + ", fasta=GQY1XT001.shhh.trim.fasta, name=GQY1XT001.shhh.trim.names";
    test.m->mothurOut("Running command: unique.seqs(" + inputs + ")\n");
    
    Command* deconvoluteCommand = new DeconvoluteCommand(inputs);
    deconvoluteCommand->execute();
    delete deconvoluteCommand;

    //align.seqs
    inputs = test.setDirInputs + ", fasta=current, reference=silva.bacteria.fasta";
    test.m->mothurOut("Running command: align.seqs(" + inputs + ")\n");
    
    Command* alignSeqsCommand = new AlignCommand(inputs);
    alignSeqsCommand->execute();
    delete alignSeqsCommand;

    //screen.seqs
    inputs = test.setDirInputs + ", fasta=current, group=current, name=current, end=27659, optimize=start, criteria=95";
    test.m->mothurOut("Running command: screen.seqs(" + inputs + ")\n");
    
    Command* screenSeqsCommand = new ScreenSeqsCommand(inputs);
    screenSeqsCommand->execute();
    delete screenSeqsCommand;
    
    //filter.seqs
    inputs = test.setDirInputs + ", fasta=current, vertical=T, trump=.";
    test.m->mothurOut("Running command: filter.seqs(" + inputs + ")\n");
    
    Command* filterSeqsCommand = new FilterSeqsCommand(inputs);
    filterSeqsCommand->execute();
    delete filterSeqsCommand;
    
    //unique.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current";
    test.m->mothurOut("Running command: unique.seqs(" + inputs + ")\n");
    
    Command* deconvoluteCommand2 = new DeconvoluteCommand(inputs);
    deconvoluteCommand2->execute();
    delete deconvoluteCommand2;
    
    //pre.cluster
    inputs = test.setDirInputs + ", fasta=current, name=current, group=current, diffs=2";
    test.m->mothurOut("Running command: pre.cluster(" + inputs + ")\n");
    
    Command* preClusterCommand = new PreClusterCommand(inputs);
    preClusterCommand->execute();
    delete preClusterCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current";
    test.m->mothurOut("\nRunning command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, RemoveChimerasContaminants) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.RemoveChimerasContaminants.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //chimera.uchime
    inputs = test.setDirInputs + ", fasta=/Users/sarahwestcott/desktop/454_sop/GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.fasta, group=/Users/sarahwestcott/desktop/454_sop/GQY1XT001.shhh.good.groups, name=/Users/sarahwestcott/desktop/454_sop/GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.names";
    test.m->mothurOut("Running command: chimera.uchime(" + inputs + ")\n");
    
    Command* chimeraCommand = new ChimeraUchimeCommand(inputs);
    chimeraCommand->execute();
    delete chimeraCommand;
    
    //remove.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current, group=current";
    test.m->mothurOut("Running command: remove.seqs(" + inputs + ")\n");
    
    Command* removeSeqsCommand = new RemoveSeqsCommand(inputs);
    removeSeqsCommand->execute();
    delete removeSeqsCommand;
    
    //classify.seqs(fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax
    inputs = test.setDirInputs + ", fasta=current, name=current, group=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax";
    test.m->mothurOut("Running command: classify.seqs(" + inputs + ")\n");
    
    Command* classifySeqsCommand = new ClassifySeqsCommand(inputs);
    classifySeqsCommand->execute();
    delete classifySeqsCommand;
    
    //remove.lineage
    inputs = test.setDirInputs + ", fasta=current, name=current, group=current, taxonomy=current, taxon=Mitochondria-Chloroplast-Archaea-Eukaryota-unknown";
    test.m->mothurOut("Running command: remove.lineage(" + inputs + ")\n");
    
    Command* removeLineageCommand = new RemoveLineageCommand(inputs);
    removeLineageCommand->execute();
    delete removeLineageCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current, name=current";
    test.m->mothurOut("Running command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, AssessingErrorRates) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.AssessingErrorRates.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //filter.seqs
    inputs = test.setDirInputs + ", fasta=HMP_MOCK.v35.align, hard=GQY1XT001.filter, vertical=T, trump=.";
    test.m->mothurOut("Running command: filter.seqs(" + inputs + ")\n");
    
    Command* filterSeqsCommand = new FilterSeqsCommand(inputs);
    filterSeqsCommand->execute();
    delete filterSeqsCommand;
    
    //get.groups
    inputs = test.setDirInputs + ", fasta=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta, group=GQY1XT001.shhh.good.pick.pick.groups, name=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.names, groups=MOCK.GQY1XT001";
    test.m->mothurOut("Running command: get.groups(" + inputs + ")\n");
    
    Command* getGroupsCommand = new GetGroupsCommand(inputs);
    getGroupsCommand->execute();
    delete getGroupsCommand;
    
    //seq.error
    inputs = test.setDirInputs + ", fasta=current, name=current, reference=HMP_MOCK.v35.filter.fasta";
    test.m->mothurOut("Running command: seq.error(" + inputs + ")\n");
    
    Command* seqErrorCommand = new SeqErrorCommand(inputs);
    seqErrorCommand->execute();
    delete seqErrorCommand;
    
    //system
    inputs = "grep -c /""2$/"" /Users/sarahwestcott/desktop/454_sop/GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.pick.error.chimera";
    test.m->mothurOut("Running command: system(" + inputs + ")\n");
    
    Command* systemCommand = new SystemCommand(inputs);
    systemCommand->execute();
    delete systemCommand;
    
    //system
    inputs = "cut -f 2 *error.summary | sort | uniq > /Users/sarahwestcott/desktop/454_sop/mock.accnos";
    test.m->mothurOut("Running command: system(" + inputs + ")\n");
    
    Command* systemCommand2 = new SystemCommand(inputs);
    systemCommand2->execute();
    delete systemCommand2;
    
    //get.seqs
    inputs = test.setDirInputs + ", accnos=current, fasta=HMP_MOCK.v35.filter.fasta";
    test.m->mothurOut("Running command: remove.seqs(" + inputs + ")\n");
    
    Command* getSeqsCommand = new GetSeqsCommand(inputs);
    getSeqsCommand->execute();
    delete getSeqsCommand;
    
    //dist.seqs
    inputs = test.setDirInputs + ", fasta=current, output=lt";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand = new DistanceCommand(inputs);
    distSeqsCommand->execute();
    delete distSeqsCommand;

    //cluster
    inputs = test.setDirInputs + ", phylip=current";
    test.m->mothurOut("Running command: cluster(" + inputs + ")\n");
    
    Command* clusterCommand = new ClusterCommand(inputs);
    clusterCommand->execute();
    delete clusterCommand;

    //dist.seqs
    inputs = test.setDirInputs + ", fasta=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.pick.fasta, cutoff=0.03";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand2 = new DistanceCommand(inputs);
    distSeqsCommand2->execute();
    delete distSeqsCommand2;
    
    //cluster
    inputs = test.setDirInputs + ", column=current, name=current";
    test.m->mothurOut("Running command: cluster(" + inputs + ")\n");
    
    Command* clusterCommand2 = new ClusterCommand(inputs);
    clusterCommand2->execute();
    delete clusterCommand2;

    //summary.single
    inputs = test.setDirInputs + ", calc=sobs, label=0.03";
    test.m->mothurOut("Running command: summary.single(" + inputs + ")\n");
    
    Command* summaryCommand = new SummaryCommand(inputs);
    summaryCommand->execute();
    delete summaryCommand;
    
    //summary.single
    inputs = test.setDirInputs + ", calc=sobs, label=0.03, subsample=4419";
    test.m->mothurOut("Running command: summary.single(" + inputs + ")\n");
    
    Command* summaryCommand2 = new SummaryCommand(inputs);
    summaryCommand2->execute();
    delete summaryCommand2;
    
    //label	method	sobs
    //0.03	ave	26.227000
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, PreparingForAnalysis) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.PreparingForAnalysis.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //remove.groups
    inputs = test.setDirInputs + ", fasta=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta, group=GQY1XT001.shhh.good.pick.pick.groups, name=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pick.names, taxonomy=GQY1XT001.shhh.trim.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, groups=MOCK.GQY1XT001";
    test.m->mothurOut("Running command: remove.groups(" + inputs + ")\n");
    
    Command* removeGroupsCommand = new RemoveGroupsCommand(inputs);
    removeGroupsCommand->execute();
    delete removeGroupsCommand;
    
    //rename.file(fasta=current, name=current, group=current, prefix=final, deleteold=false)
    inputs = test.setDirInputs + ", fasta=current, name=current, group=current, prefix=final, deleteold=false";
    test.m->mothurOut("Running command: rename.file(" + inputs + ")\n");
    
    Command* renameFileCommand = new RenameFileCommand(inputs);
    renameFileCommand->execute();
    delete renameFileCommand;
    
    //dist.seqs
    inputs = test.setDirInputs + ", fasta=current, cutoff=0.03";
    test.m->mothurOut("Running command: dist.seqs(" + inputs + ")\n");
    
    Command* distSeqsCommand = new DistanceCommand(inputs);
    distSeqsCommand->execute();
    delete distSeqsCommand;
    
    //cluster
    inputs = test.setDirInputs + ", column=current, name=current";
    test.m->mothurOut("Running command: cluster(" + inputs + ")\n");
    
    Command* clusterCommand = new ClusterCommand(inputs);
    clusterCommand->execute();
    delete clusterCommand;
    
    //make.shared
    inputs = test.setDirInputs + ", list=current, group=current, label=0.03";
    test.m->mothurOut("Running command: make.shared(" + inputs + ")\n");
    
    Command* makeSharedCommand = new SharedCommand(inputs);
    makeSharedCommand->execute();
    delete makeSharedCommand;
    
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
    
    //classify.otu
    inputs = test.setDirInputs + ", list=current, name=current, taxonomy=current, label=0.03";
    test.m->mothurOut("Running command: classify.otu(" + inputs + ")\n");
    
    Command* classifyOtuCommand = new ClassifyOtuCommand(inputs);
    classifyOtuCommand->execute();
    delete classifyOtuCommand;
    
    //phylotype
    inputs = test.setDirInputs + ", name=current, taxonomy=current, label=1";
    test.m->mothurOut("Running command: classify.otu(" + inputs + ")\n");
    
    Command* phylotypeCommand = new PhylotypeCommand(inputs);
    phylotypeCommand->execute();
    delete phylotypeCommand;
    
    //make.shared
    inputs = test.setDirInputs + ", list=current, group=current";
    test.m->mothurOut("Running command: make.shared(" + inputs + ")\n");
    
    Command* makeSharedCommand2 = new SharedCommand(inputs);
    makeSharedCommand2->execute();
    delete makeSharedCommand2;
    
    //classify.otu
    inputs = test.setDirInputs + ", list=current, name=current, taxonomy=current";
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
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, OTUBasedAnalysis) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.OTUBasedAnalysis.log";
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
    inputs = test.setDirInputs + ", shared=current, calc=nseqs-coverage-sobs-invsimpson, subsample=T";
    test.m->mothurOut("Running command: summary.single(" + inputs + ")\n");
    
    Command* summarySingleCommand = new SummaryCommand(inputs);
    summarySingleCommand->execute();
    delete summarySingleCommand;
    
    //collect.single
    inputs = test.setDirInputs + ", shared=current, calc=chao-invsimpson, freq=100";
    test.m->mothurOut("Running command: collect.single(" + inputs + ")\n");
    
    Command* collectSingleCommand = new CollectCommand(inputs);
    collectSingleCommand->execute();
    delete collectSingleCommand;
    
    //heatmap.bin
    inputs = test.setDirInputs + ", shared=current, scale=log2, numotu=50";
    test.m->mothurOut("Running command: heatmap.bin(" + inputs + ")\n");
    
    Command* heatmapCommand = new HeatMapCommand(inputs);
    heatmapCommand->execute();
    delete heatmapCommand;
    
    //heatmap.sim
    inputs = test.setDirInputs + ", calc=jclass-thetayc";
    test.m->mothurOut("Running command: heatmap.sim(" + inputs + ")\n");
    
    Command* heatmapSimCommand = new HeatMapSimCommand(inputs);
    heatmapSimCommand->execute();
    delete heatmapSimCommand;
    
    //venn
    inputs = test.setDirInputs + ", groups=F003D000-F003D002-F003D004-F003D006";
    test.m->mothurOut("Running command: venn(" + inputs + ")\n");
    
    Command* vennCommand = new VennCommand(inputs);
    vennCommand->execute();
    delete vennCommand;
    
    //summary.shared
    inputs = test.setDirInputs + ", shared=current, calc=sharedchao, groups=F003D000-F003D002-F003D004-F003D006, all=T";
    test.m->mothurOut("Running command: summary.shared(" + inputs + ")\n");
    
    Command* summarySharedCommand = new SummarySharedCommand(inputs);
    summarySharedCommand->execute();
    delete summarySharedCommand;
    
    //tree.shared -
    inputs = test.setDirInputs + ", calc=thetayc-jclass, subsample=T";
    test.m->mothurOut("Running command: tree.shared(" + inputs + ")\n");
    
    Command* treeSharedCommand = new TreeGroupCommand(inputs);
    treeSharedCommand->execute();
    delete treeSharedCommand;

    //parsimony - tree=final.opti_mcc.thetayc.0.03.lt.ave.tre
    inputs = test.setDirInputs + ", tree=final.opti_mcc.thetayc.0.03.lt.ave.tre, group=mouse.sex_time.design, groups=all";
    test.m->mothurOut("Running command: parsimony(" + inputs + ")\n");
    
    Command* parsimonyCommand = new ParsimonyCommand(inputs);
    parsimonyCommand->execute();
    delete parsimonyCommand;
    
    //unifrac.unweighted
    inputs = test.setDirInputs + ", tree=current, group=current, random=T, groups=all";
    test.m->mothurOut("Running command: unifrac.unweighted(" + inputs + ")\n");
    
    Command* unifracUnweightedCommand = new UnifracUnweightedCommand(inputs);
    unifracUnweightedCommand->execute();
    delete unifracUnweightedCommand;
    
    //unifrac.weighted
    inputs = test.setDirInputs + ", tree=current, group=current, random=T";
    test.m->mothurOut("Running command: unifrac.weighted(" + inputs + ")\n");
    
    Command* unifracWeightedCommand = new UnifracWeightedCommand(inputs);
    unifracWeightedCommand->execute();
    delete unifracWeightedCommand;
    
    //dist.shared
    inputs = test.setDirInputs + ", shared=final.opti_mcc.shared, calc=thetayc-jclass, subsample=T";
    test.m->mothurOut("Running command: dist.shared(" + inputs + ")\n");
    
    Command* distSharedCommand = new MatrixOutputCommand(inputs);
    distSharedCommand->execute();
    delete distSharedCommand;
    
    //pcoa - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist";
    test.m->mothurOut("Running command: pcoa(" + inputs + ")\n");
    
    Command* pcoaCommand = new PCOACommand(inputs);
    pcoaCommand->execute();
    delete pcoaCommand;
    
    //pcoa - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.jclassc.0.03.lt.ave.dist";
    test.m->mothurOut("Running command: pcoa(" + inputs + ")\n");
    
    Command* pcoaCommand2 = new PCOACommand(inputs);
    pcoaCommand2->execute();
    delete pcoaCommand2;
    
    //nmds - phylip=stability.opti_mcc.jclass.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.jclass.0.03.lt.ave.dist";
    test.m->mothurOut("Running command: nmds(" + inputs + ")\n");
    
    Command* nmdsCommand = new NMDSCommand(inputs);
    nmdsCommand->execute();
    delete nmdsCommand;
    
    //nmds - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3";
    test.m->mothurOut("Running command: nmds(" + inputs + ")\n");
    
    Command* nmdsCommand2 = new NMDSCommand(inputs);
    nmdsCommand2->execute();
    delete nmdsCommand2;
    
    //amova - phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist
    inputs = test.setDirInputs + ", phylip=current, design=mouse.sex_time.design";
    test.m->mothurOut("Running command: amova(" + inputs + ")\n");
    
    Command* amovaCommand = new AmovaCommand(inputs);
    amovaCommand->execute();
    delete amovaCommand;
    
    //corr.axes
    inputs = test.setDirInputs + ", axes=final.opti_mcc.thetayc.0.03.lt.ave.nmds.axes, shared=final.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3";
    test.m->mothurOut("Running command: corr.axes(" + inputs + ")\n");
    
    Command* corrAxesCommand = new CorrAxesCommand(inputs);
    corrAxesCommand->execute();
    delete corrAxesCommand;
    
    //corr.axes(axes=stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3)
    inputs = test.setDirInputs + ", axes=final.opti_mcc.thetayc.0.03.lt.ave.nmds.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3";
    test.m->mothurOut("Running command: corr.axes(" + inputs + ")\n");
    
    Command* corrAxesCommand2 = new CorrAxesCommand(inputs);
    corrAxesCommand2->execute();
    delete corrAxesCommand2;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

TEST(Test_Integration_454SOP, PhylogenyBasedAnalysis) {
    Test454SOP test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_454SOP.PhylogenyBasedAnalysis.log";
    test.m->mothurOut("Running command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //metastats
    inputs = test.setDirInputs + ", shared=final.opti_mcc.0.03.subsample.shared, design=mouse.sex_time.design";
    test.m->mothurOut("Running command: metastats(" + inputs + ")\n");
    
    Command* metastatsCommand = new MetaStatsCommand(inputs);
    metastatsCommand->execute();
    delete metastatsCommand;
    
    //phylo.diversity
    inputs = test.setDirInputs + ", tree=final.phylip.tre, name=final.names, group=final.groups, rarefy=T";
    test.m->mothurOut("Running command: phylo.diversity(" + inputs + ")\n");
    
    Command* phyloDiversityCommand = new PhyloDiversityCommand(inputs);
    phyloDiversityCommand->execute();
    delete phyloDiversityCommand;
    
    //unifrac.unweighted
    inputs = test.setDirInputs + ", tree=current, name=current, group=current, distance=lt, random=F, subsample=T";
    test.m->mothurOut("Running command: unifrac.unweighted(" + inputs + ")\n");
    
    Command* unifracUnweightedCommand = new UnifracUnweightedCommand(inputs);
    unifracUnweightedCommand->execute();
    delete unifracUnweightedCommand;
    
    //unifrac.weighted
    inputs = test.setDirInputs + ", tree=current, name=current, group=current, distance=lt, random=F, subsample=T";
    test.m->mothurOut("Running command: unifrac.weighted(" + inputs + ")\n");
    
    Command* unifracWeightedCommand = new UnifracWeightedCommand(inputs);
    unifracWeightedCommand->execute();
    delete unifracWeightedCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
}

/**************************************************************************************************/


