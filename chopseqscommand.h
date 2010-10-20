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
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();	
		
	private:
		string fastafile, outputDir, keep;
		bool abort, countGaps, Short;
		int numbases;
		vector<string> outputNames;
		map<string, vector<string> > outputTypes;
		
		string getChopped(Sequence);
		
		
};

#endif


