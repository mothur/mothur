#ifndef PARSIMONY_H
#define PARSIMONY_H


/*
 *  parsimony.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "treemap.h"
#include "globaldata.hpp"

/***********************************************************************/

class Parsimony : public TreeCalculator  {
	
	public:
		Parsimony(TreeMap* t) : tmap(t) {};
		~Parsimony() {};
		EstOutput getValues(Tree*, int, string);
		//EstOutput getValues(Tree*, string, string) { return data; }
		
	private:
		struct linePair {
			int start;
			int num;
			linePair(int i, int j) : start(i), num(j) {}
		};
		vector<linePair> lines;
	
		GlobalData* globaldata;
		EstOutput data;
		TreeMap* tmap;
		int processors;
		string outputDir;
	
		EstOutput driver(Tree*, vector< vector<string> >, int, int); 
		EstOutput createProcesses(Tree*, vector< vector<string> >);
};

/***********************************************************************/

#endif
