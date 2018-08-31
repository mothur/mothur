//
//  testbiominfocommand.h
//  Mothur
//
//  Created by Sarah Westcott on 8/18/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testbiominfocommand__
#define __Mothur__testbiominfocommand__

#include "biominfocommand.h"

class TestBiomInfoCommand : public BiomInfoCommand {
    
public:
    
    using BiomInfoCommand::getTag;
    using BiomInfoCommand::getName;
    using BiomInfoCommand::getTaxonomy;
    using BiomInfoCommand::readRows;
    using BiomInfoCommand::getDims;
    using BiomInfoCommand::readData;
    using BiomInfoCommand::getNamesAndTaxonomies;

};


#endif /* defined(__Mothur__testbiominfocommand__) */
