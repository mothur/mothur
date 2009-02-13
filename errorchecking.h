#ifndef ERRORCHECKING_H
#define ERRORCHECKING_H
/*
 *  errorchecking.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <iostream>
#include <map>
#include "globaldata.hpp"
#include "validcommands.h"
#include "validparameter.h"
#include "utilities.hpp"

class ErrorCheck {
	public:
		ErrorCheck();
		~ErrorCheck();
		bool checkInput(string);
	
	private: 
		GlobalData* globaldata;
		ValidCommands* validCommand;
		ValidParameters* validParameter;
		void validateReadFiles();
		void validateReadDist();
		void validateReadPhil();
		void validateParseFiles();
		void validateTreeFiles();
		void clear();
		void refresh();
		string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, fastafile, treefile, sharedfile, cutoff, format; 
		string precision, method, fileroot, label, line, iters, jumble, freq, single, rarefaction, shared, summary, randomtree;
		string commandName, optionText;
		bool errorFree;
		
};
#endif
