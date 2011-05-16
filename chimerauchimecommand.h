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
	string getCitation() { return "http://drive5.com/uchime/ \nhttp://www.mothur.org/wiki/Chimera.uchime"; }
	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	vector<int> processIDS;   //processid
	int driver(string, string, string);
	int createProcesses(string, string, string);
	
#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
#endif
	
	bool abort;
	string fastafile, templatefile, outputDir, namefile;
	int processors;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	
};

/***********************************************************/

#endif


