#ifndef CHIMERAPERSEUSCOMMAND_H
#define CHIMERAPERSEUSCOMMAND_H


/*
 *  chimeraperseuscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/26/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */



#include "mothur.h"
#include "command.hpp"
#include "sequenceparser.h"
#include "sequencecountparser.h"
#include "myPerseus.h"
#include "counttable.h"

/***********************************************************/
class ChimeraPerseusCommand : public Command {
public:
	ChimeraPerseusCommand(string);
	ChimeraPerseusCommand();
	~ChimeraPerseusCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.perseus";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Quince C, Lanzen A, Davenport RJ, Turnbaugh PJ (2011).  Removing noise from pyrosequenced amplicons.  BMC Bioinformatics  12:38.\nEdgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011), UCHIME improves sensitivity and speed of chimera detection.  Bioinformatics 27:2194.\nhttp://www.mothur.org/wiki/Chimera.perseus\n"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	bool abort, hasName, hasCount, dups;
	string fastafile, groupfile, countfile, outputDir, namefile;
	int processors, alignLength, numChimeras;
	double cutoff, alpha, beta;
    SequenceParser* parser;
    SequenceCountParser* cparser;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> groupFileNames;
	
	string getNamesFile(string&);
	vector<seqData> readFiles(string, string);
    vector<seqData> readFiles(string inputFile, CountTable* ct);
    int deconvoluteResults(map<string, string>&, string, string);
	int createProcessesGroups(string, string, string, vector<string>, string, string, string);

};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct perseusData {
	string outputFName;
	string accnos;
    string countlist;
    string outputFNameGroup;
    string accnosGroup;
    string countlistGroup;
	MothurOut* m;
	int start;
	int end;
    bool hasName, hasCount, dups;
	long long threadID, count, numChimeras, chimerasInGroup;
	double alpha, beta, cutoff;
	vector<string> groups;
    SequenceParser* parser;
    SequenceCountParser* cparser;
    vector<seqData> sequences;
    int alignLength;
	
	perseusData(){}
	perseusData(SequenceParser* sp, SequenceCountParser* cp, string cl, string ac, bool dps, bool hn, bool hc, double a, double b, double c, string o, vector<string> gr, MothurOut* mout, int st, int en) {
		alpha = a;
		beta = b;
		cutoff = c;
        outputFName = o; outputFNameGroup = o;
        countlist = cl; countlistGroup = cl;
        accnos = ac; accnosGroup = ac;
		m = mout;
		start = st;
		end = en;
		groups = gr;
        hasName = hn;
        hasCount = hc;
        dups = dps;
		count = 0;
		numChimeras = 0;
        parser = sp;
        cparser = cp;
        alignLength = 0;
	}
    
    perseusData(SequenceParser* sp, SequenceCountParser* cp, string cl, string ac, bool dps, bool hn, bool hc, double a, double b, double c, string o, MothurOut* mout) {
        alpha = a;
        beta = b;
        cutoff = c;
        outputFName = o; outputFNameGroup = o;
        countlist = cl; countlistGroup = cl;
        accnos = ac; accnosGroup = ac;
        m = mout;
        start = 0;
        end = 0;
        hasName = hn;
        hasCount = hc;
        dups = dps;
        count = 0;
        numChimeras = 0;
        parser = sp;
        cparser = cp;
        alignLength = 0;
    }
};
/**************************************************************************************************/

#endif


