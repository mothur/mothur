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
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"

/**************************************************************************************************/

class VsearchFileParser {
    
#ifdef UNIT_TEST
    friend class TestVsearchFileParser;
#endif
    
    public:
        VsearchFileParser();
        VsearchFileParser(string f); //fastafile
        VsearchFileParser(string f, string n, string format); //fastafile, name or count file, format("name" or "count")
        ~VsearchFileParser(){}
    
        string getVsearchFile(); //creates vsearch formatted file. degapped, abundance appended, and sorted.
        string getCountFile(); //returns count filename. If no names or count file is given, then a count file is generated. 
        ListVector createListFile(string uc, int numbins, string label, map<string, int>& counts); //takes vsearch output and creates listvector
        int getNumBins(string logfile);
        void setFastaFile(string f) { fastafile = f; }
    
    private:
        MothurOut* m;
        map<string, string> nameMap;
        map<string, int> counts;
        string fastafile, namefile, countfile, format;
        string getNamesFile(string& inputFile);
        string createVsearchFasta(string);
        string removeAbundances(string seqName);
    
    
};

/**************************************************************************************************/


#endif /* defined(__Mothur__vsearchfileparser__) */
