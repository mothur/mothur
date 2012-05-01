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

/***********************************************************************/

class Parsimony : public TreeCalculator  {
	
	public:
		Parsimony() {};
		~Parsimony() {};
		EstOutput getValues(Tree*, int, string);
		
	private:
		struct linePair {
			int start;
			int num;
			linePair(int i, int j) : start(i), num(j) {}
		};
		vector<linePair> lines;
	
		EstOutput data;
		int processors;
		string outputDir;
	
		EstOutput driver(Tree*, vector< vector<string> >, int, int, TreeMap*); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, TreeMap*);
};

/***********************************************************************/

#endif
