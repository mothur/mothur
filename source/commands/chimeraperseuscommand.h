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
	int processors, alignLength;
	double cutoff, alpha, beta;
    vector<string> outputNames;
	
	string getNamesFile(string&);
	vector<seqData> readFiles(string, string);
    vector<seqData> readFiles(string inputFile, CountTable* ct);
	int deconvoluteResults(map<string, string>&, string, string);
	int createProcessesGroups(string, string, string, vector<string>, string, string, string, int&);
};
/***********************************************************/

#endif


