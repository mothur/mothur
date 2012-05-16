#ifndef Mothur_makecontigscommand_h
#define Mothur_makecontigscommand_h

//
//  makecontigscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"

struct fastqRead {
	vector<int> scores;
	string name;
	string sequence;
	
	fastqRead() { name = ""; sequence = ""; scores.clear(); };
	fastqRead(string n, string s, vector<int> sc) : name(n), sequence(s), scores(sc) {};
	~fastqRead() {};
};

/**************************************************************************************************/

class MakeContigsCommand : public Command {
public:
    MakeContigsCommand(string);
    MakeContigsCommand();
    ~MakeContigsCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.contigs";			}
    string getCommandCategory()		{ return "Sequence Processing";		} 
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
    string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Make.contigs"; }
    string getDescription()		{ return "description"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    string outputDir, ffastqfile, rfastqfile, align;
	float match, misMatch, gapOpen, gapExtend;
	int processors, longestBase;
    vector<string> outputNames;
    
    fastqRead readFastq(ifstream&);
    vector< vector<string> > readFastqFiles(int&);
    bool checkReads(fastqRead&, fastqRead&);
};

/**************************************************************************************************/







#endif
