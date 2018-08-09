//
//  testintegrationalignseqscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testintegrationalignseqs.hpp"
#include "aligncommand.h"
#include "setlogfilecommand.h"
#include "setseedcommand.h"
#include "seqsummarycommand.h"

/**************************************************************************************************/
TestAlignSeqsIntegration::TestAlignSeqsIntegration() {  //setup
    inputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    outputDir = "/Users/sarahwestcott/desktop/Miseq_sop/";
    setDirInputs = "inputdir=" + inputDir + ", outputdir=" + outputDir;
    current = CurrentFile::getInstance();
    m = MothurOut::getInstance();
    
    //remove shortcut files
    Utils util; util.mothurRemove(inputDir+"silva.v4.8mer");
    
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(); //Fasta, name, group returned
}
/**************************************************************************************************/
TestAlignSeqsIntegration::~TestAlignSeqsIntegration() { }
/**************************************************************************************************/

TEST(Test_Integration_AlignSeqs, testOptions) {
    TestAlignSeqsIntegration test;
    
    test.m->mothurOut("/******************************************/\n");
    test.current->setMothurCalling(true);
    
    string inputs = "seed=123456";
    test.m->mothurOut("\nRunning command: set.seed(" + inputs + ")\n");
    
    Command* setSeedCommand = new SetSeedCommand(inputs);
    setSeedCommand->execute();
    delete setSeedCommand;
    
    //set.logfile
    inputs = test.setDirInputs + ", name=Test_Integration_AlignSeqs.testOptions.log";
    test.m->mothurOut("\nRunning command: set.logfile(" + inputs + ")\n");
    
    Command* setLogfileCommand = new SetLogFileCommand(inputs);
    setLogfileCommand->execute();
    delete setLogfileCommand;
    
    //align.seqs
    inputs = test.setDirInputs + ", fasta=" + test.filenames[0] + ", reference=silva.v4.fasta";
    test.m->mothurOut("\nRunning command: align.seqs(" + inputs + ")\n");
    
    Command* alignSeqsCommand = new AlignCommand(inputs);
    alignSeqsCommand->execute();
    delete alignSeqsCommand;
    
    //summary.seqs
    inputs = test.setDirInputs + ", fasta=current";
    test.m->mothurOut("\nRunning command: summary.seqs(" + inputs + ")\n");
    
    Command* summarySeqsCommand = new SeqSummaryCommand(inputs);
    summarySeqsCommand->execute();
    delete summarySeqsCommand;
    
    test.current->setMothurCalling(false);
    test.m->mothurOut("/******************************************/\n");
    
}

/**************************************************************************************************/


