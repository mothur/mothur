//
//  mergesfffilecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 1/31/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#ifndef Mothur_mergesfffilecommand_h
#define Mothur_mergesfffilecommand_h

#include "command.hpp"
#include "sffheader.hpp"
#include "sffread.hpp"

/**********************************************************/
class MergeSfffilesCommand : public Command {
	
public:
	MergeSfffilesCommand(string);
    ~MergeSfffilesCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "merge.sfffiles";					}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();
    string getOutputPattern(string);
	string getCitation() { return "http://www.mothur.org/wiki/merge.sfffiles"; }
	string getDescription()		{ return "merge individual sfffiles into a single .sff file"; }
    
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	string sffFilename,  file, currentFileName;
	vector<string> filenames, outputNames;
	bool abort, keyTrim;
    int numTotalReads, allFilesnumFlowReads, allFileskeyLength;
    string outputFile, outputFileHeader;
    vector<SffCommonHeader*> commonHeaders;
    
	//extract sff file functions
	int mergeSffInfo(string, ofstream&);
    void adjustCommonHeader();
    void readFile();
};

/**********************************************************/

#endif
