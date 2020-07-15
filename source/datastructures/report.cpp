//
//  report.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/15/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "report.hpp"

/**************************************************************************************************

vector<string> Report::readHeaders(ifstream& repFile, string repFileName){
    try {
        
        util.openInputFile(repFileName, repFile);
        string headers = util.getline(repFile);
        
        vector<string> sHeaders = util.splitWhiteSpace(headers);
        
        return sHeaders;
    }
    catch(exception& e) {
        m->errorOut(e, "Report", "readHeaders");
        exit(1);
    }
}
/**************************************************************************************************/

