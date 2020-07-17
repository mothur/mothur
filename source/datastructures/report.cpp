//
//  report.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/15/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "report.hpp"

/**************************************************************************************************/

vector<string> Report::readHeaders(ifstream& repFile){
    try {
        
        util.openInputFile(reportFileName, repFile);
        string headers = util.getline(repFile);
        
        reportHeaders = util.splitWhiteSpace(headers);
        
        return reportHeaders;
    }
    catch(exception& e) {
        m->errorOut(e, "Report", "readHeaders");
        exit(1);
    }
}
/**************************************************************************************************/

void Report::printHeaders(ofstream& repFile){
    try {
        for (int i = 0; i < reportHeaders.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            repFile << reportHeaders[i] << '\t';
        }
        repFile << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "Report", "printHeaders");
        exit(1);
    }
}
/**************************************************************************************************/

/**************************************************************************************************/

