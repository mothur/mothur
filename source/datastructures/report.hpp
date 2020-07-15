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

/**************************************************************************************************/

class Report {

public:
    
    Report() { m = MothurOut::getInstance(); }
    virtual ~Report() {}
    
    virtual int read(ifstream&) = 0;
    
    vector<string> readHeaders(ifstream& repFile, string repFileName) {
        util.openInputFile(repFileName, repFile);
        string headers = util.getline(repFile);
        
        return (util.splitWhiteSpace(headers));
    }

    
protected:
        
    MothurOut* m;
    Utils util;
    
};

/**************************************************************************************************/


#endif /* report_hpp */
