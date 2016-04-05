//
//  testremovegroupscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 7/30/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testremovegroupscommand__
#define __Mothur__testremovegroupscommand__

#include "removegroupscommand.h"

class TestRemoveGroupsCommand : public RemoveGroupsCommand {
    
public:
    
    using RemoveGroupsCommand::readFasta;
    using RemoveGroupsCommand::readName;
    using RemoveGroupsCommand::readGroup;
    using RemoveGroupsCommand::readCount;
    using RemoveGroupsCommand::readList;
    using RemoveGroupsCommand::readTax;
    using RemoveGroupsCommand::fillNames;
    using RemoveGroupsCommand::readShared;
    using RemoveGroupsCommand::readDesign;
    using RemoveGroupsCommand::readPhylip;
    using RemoveGroupsCommand::readColumn;
    
    
};


#endif /* defined(__Mothur__testremovegroupscommand__) */
