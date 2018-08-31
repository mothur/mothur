//
//  testrenamefilecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/4/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testrenamefilecommand__
#define __Mothur__testrenamefilecommand__

#include "renamefilecommand.h"

class TestRenameFileCommand : public RenameFileCommand {
    
public:
    
    TestRenameFileCommand();
    ~TestRenameFileCommand();
    

    MothurOut* m;
    vector<string> filenames;
    
    //private functions
    using RenameFileCommand::getNewName;
    using RenameFileCommand::renameOrCopy;
    
    //private variables
    using RenameFileCommand::prefix;
    using RenameFileCommand::mothurGenerated;
    using RenameFileCommand::outputfile;
    using RenameFileCommand::deleteOld;
    
};


#endif /* defined(__Mothur__testrenamefilecommand__) */
