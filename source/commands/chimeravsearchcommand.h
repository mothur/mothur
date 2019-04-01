//
//  chimeravsearchcommand.h
//  Mothur
//
//  Created by Sarah Westcott on 6/16/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__chimeravsearchcommand__
#define __Mothur__chimeravsearchcommand__

#include "command.hpp"
#include "sequenceparser.h"
#include "counttable.h"
#include "sequencecountparser.h"

/***********************************************************/

class ChimeraVsearchCommand : public Command {
public:
    ChimeraVsearchCommand(string);
    ChimeraVsearchCommand();
    ~ChimeraVsearchCommand() {}
    
    vector<string> setParameters();
    string getCommandName()			{ return "chimera.vsearch";		}
    string getCommandCategory()		{ return "Sequence Processing"; }
    
    string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);
    string getCitation() { return "vsearch by https://github.com/torognes/vsearch.\nhttp://www.mothur.org/wiki/Chimera.vsearch\n"; }
    string getDescription()		{ return "detect chimeric sequences"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    vector<int> processIDS;   //processid
    int driver(string, string, string, string, int&);
    
    bool abort, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, ucl, useMindiffs, hasCount, dups;
    string fastafile, templatefile, outputDir, countfile, abskew, minh, mindiv, xn, dn, mindiffs, vsearchLocation;
    int processors;
    
    vector<string> outputNames;
    
    string getCountFile(string&);
    int readFasta(string, map<string, string>&);
    int deconvoluteResults(string, string, string, long long&);
    int driverGroups(map<string, vector<string> >, string, string, string, string, string);
    int getSeqs(map<string, int>& nameMap, string thisGroupsFormattedOutputFilename, string tag, string tag2, long long& numSeqs, string thisGroupsFastaFile);

    int prepFile(string filename, string);
};

/***********************************************************/

#endif
