#ifndef CHIMERAUCHIMECOMMAND_H
#define CHIMERAUCHIMECOMMAND_H


/*
 *  chimerauchimecommand.h
 *  Mothur
 *
 *  Created by westcott on 5/13/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"

/***********************************************************/

class ChimeraUchimeCommand : public Command {
public:
	ChimeraUchimeCommand(string);
	ChimeraUchimeCommand();
	~ChimeraUchimeCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.uchime";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "uchime by Robert C. Edgar\nhttp://drive5.com/uchime\nThis code is donated to the public domain.\nhttp://www.mothur.org/wiki/Chimera.uchime\nEdgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011), UCHIME improves sensitivity and speed of chimera detection, Bioinformatics, in press.\n"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	vector<int> processIDS;   //processid
	int driver(string, string, string, string);
	int createProcesses(string, string, string, string);
		
	bool abort, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, useChunks, useMinchunk, useIdsmoothwindow, useMinsmoothid, useMaxp, skipgaps, skipgaps2, useMinlen, useMaxlen, ucl, useQueryfract;
	string fastafile, templatefile, outputDir, namefile, abskew, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, minlen, maxlen, queryfract;
	int processors;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	
};

/***********************************************************/

#endif


