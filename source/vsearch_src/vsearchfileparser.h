//
//  vsearchfileparser.h
//  Mothur
//
//  Created by Sarah Westcott on 10/13/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__vsearchfileparser__
#define __Mothur__vsearchfileparser__

#include "mothurout.h"


/**************************************************************************************************/

class VsearchFileParser {
    
    public:
        VsearchFileParser(string f);
        VsearchFileParser(string f, string n, string format);
        ~VsearchFileParser(){}
    
        string getVsearchFile();
    
    private:
        MothurOut* m;
        string fastafile, namefile, countfile;
        string getNamesFile(string& inputFile);
        string createVsearchFasta(string, map<string, int>&);
    
    
};

/**************************************************************************************************/


#endif /* defined(__Mothur__vsearchfileparser__) */
