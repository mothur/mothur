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

/**************************************************************************************************/
TestMiSeqSOP::TestMiSeqSOP() {  //setup
    inputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    outputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    setDirInputs = "inputdir=" + inputDir + ", outputdir=" + outputDir;
    current = CurrentFile::getInstance();
    m = MothurOut::getInstance();
}
/**************************************************************************************************/
TestMiSeqSOP::~TestMiSeqSOP() {
}
/**************************************************************************************************/

TEST(Test_Integration, MISEQMakeContigs2Precluster) {
    TestMiSeqSOP test;

    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("Running command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
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

TEST(Test_Integration, RemoveChimerasContaminants) {
    TestMiSeqSOP test;
    
}

/**************************************************************************************************/


