#ifndef SECONDARYSTRUCTURECHECKERCOMMAND_H
#define SECONDARYSTRUCTURECHECKERCOMMAND_H

/*
 *  secondarystructurecommand.h
 *  Mothur
 *
 *  Created by westcott on 9/18/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"

/**************************************************************************************************/

struct statData {
	int pound;
	int tilde;
	int dash;
	int plus;
	int equal;
	int loop;
	int total;
	statData() : pound(0), loop(0), tilde(0), dash(0), plus(0), equal(0), total(0) {};	
};

/**************************************************************************************************/


class AlignCheckCommand : public Command {
	
	public:
	
		AlignCheckCommand(string);	
		~AlignCheckCommand(){};
		int execute();
		void help();	
		
	private:
		vector<int> structMap;
		string mapfile, fastafile, outputDir;
		bool abort;
		int seqLength;
		
		void readMap();
		statData getStats(string sequence);
};

/**************************************************************************************************/
#endif

