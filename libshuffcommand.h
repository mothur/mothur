#ifndef LIBSHUFFCOMMAND_H
#define LIBSHUFFCOMMAND_H

/*
 *  libshuffcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "fullmatrix.h"
#include "libshuff.h"


class GlobalData;

class LibShuffCommand : public Command {
	
	public:
		LibShuffCommand();	
		~LibShuffCommand(){};
		int execute();	
	
	private:
		vector<string> groupNames;
		
		void setGroups();
		void printCoverageFile();
		void printSummaryFile();

		GlobalData* globaldata;
		FullMatrix* matrix;
		Libshuff* form;
		float cutOff, step;
		int numGroups, numComp, iters;
		string coverageFile, summaryFile;
		vector<vector<int> > pValueCounts;
		vector<vector<double> > savedDXYValues;
		vector<vector<vector<double> > > savedMinValues;
};

#endif
