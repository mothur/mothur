//
//  testrenameseqscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 9/29/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testrenameseqscommand__
#define __Mothur__testrenameseqscommand__

#include "renameseqscommand.h"

class TestRenameSeqsCommand : public RenameSeqsCommand {
    
public:
    
    MothurOut* m;
    
    //private functions
    using RenameSeqsCommand::readQual;
    using RenameSeqsCommand::readContigs;
    using RenameSeqsCommand::readFasta;
    using RenameSeqsCommand::processFile;
    using RenameSeqsCommand::readMapFile;
    using RenameSeqsCommand::readFiles;
    
    //private variables
    using RenameSeqsCommand::fastaFile;
    using RenameSeqsCommand::nameFile;
    using RenameSeqsCommand::groupfile;
    using RenameSeqsCommand::countfile;
    using RenameSeqsCommand::qualfile;
    using RenameSeqsCommand::contigsfile;
    using RenameSeqsCommand::fileFile;
    using RenameSeqsCommand::mapFile;
    
};

#endif /* defined(__Mothur__testrenameseqscommand__) */
