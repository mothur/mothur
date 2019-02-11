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
    
    bool abort, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, ucl, useMindiffs, hasCount, hasName, dups;
    string fastafile, groupfile, templatefile, outputDir, namefile, countfile, abskew, minh, mindiv, xn, dn, mindiffs, vsearchLocation;
    int processors;
    
    SequenceParser* sparser;
    SequenceCountParser* cparser;
    vector<string> outputNames;
    vector<string> fastaFileNames;
    vector<string> nameFileNames;
    vector<string> groupFileNames;
    
    string getNamesFile(string&);
    int readFasta(string, map<string, string>&);
    int deconvoluteResults(map<string, string>&, string, string, string);
    int driverGroups(string, string, string, string, string, int, int, vector<string>);
    int prepFile(string filename, string);
    
    
};

/***********************************************************/

#endif
