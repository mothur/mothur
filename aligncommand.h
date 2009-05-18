#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H

/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "globaldata.hpp"

class AlignCommand : public Command {
	
	public:
		AlignCommand();	
		~AlignCommand();
		int execute();	
	
	private:
		GlobalData* globaldata;
		string candidateFileName, templateFileName, distanceFileName;
		int kmerSize;
		float match, misMatch, gapOpen, gapExtend;
		ofstream out;
		ifstream in;

};



#endif