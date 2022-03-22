//
//  report.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/15/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef report_hpp
#define report_hpp

#include "utils.hpp"
#include "mothurout.h"

/**************************************************************************************************/

class Report {

public:
    
    Report() { m = MothurOut::getInstance(); }
    virtual ~Report() = default;
    
    virtual void read(ifstream&) = 0;
    
    vector<string> getHeaders() { return reportHeaders; }
    vector<string> readHeaders(ifstream&);
    void printHeaders(ofstream&);

    
protected:
    
    virtual void fillHeaders() = 0;
        
    MothurOut* m;
    Utils util;
    
    vector<string> reportHeaders;
    
};

/**************************************************************************************************/


#endif /* report_hpp */
