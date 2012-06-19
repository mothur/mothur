#ifndef CHOPSEQSCOMMAND_H
#define CHOPSEQSCOMMAND_H

/*
 *  chopseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 5/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"

class ChopSeqsCommand : public Command {
	
	public:
	
		ChopSeqsCommand(string);
		ChopSeqsCommand();	
		~ChopSeqsCommand(){};
	
		vector<string> setParameters();
		string getCommandName()			{ return "chop.seqs";		}
		string getCommandCategory()		{ return "Sequence Processing"; }
		string getOutputFileNameTag(string, string);
        string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Chops.seqs"; }
		string getDescription()		{ return "trim sequence length"; }
	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }		
	
	private:
		string fastafile, outputDir, keep;
		bool abort, countGaps, Short;
		int numbases;
		vector<string> outputNames;
		
		string getChopped(Sequence);
};

#endif


