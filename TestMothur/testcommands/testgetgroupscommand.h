//
//  testgetgroupscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 7/30/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testgetgroupscommand__
#define __Mothur__testgetgroupscommand__

#include "getgroupscommand.h"

class TestGetGroupsCommand : public GetGroupsCommand {
    
public:
    
    using GetGroupsCommand::readFasta;
    using GetGroupsCommand::readName;
    using GetGroupsCommand::readGroup;
    using GetGroupsCommand::readCount;
    using GetGroupsCommand::readList;
    using GetGroupsCommand::readTax;
    using GetGroupsCommand::fillNames;
    using GetGroupsCommand::readShared;
    using GetGroupsCommand::readDesign;
    using GetGroupsCommand::readPhylip;
    using GetGroupsCommand::readColumn;

    
};


#endif /* defined(__Mothur__testgetgroupscommand__) */
