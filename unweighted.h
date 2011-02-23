#ifndef UNWEIGHTED_H
#define UNWEIGHTED_H


/*
 *  unweighted.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "treemap.h"

/***********************************************************************/

class Unweighted : public TreeCalculator  {
	
	public:
		Unweighted(TreeMap* t) : tmap(t) {};
		~Unweighted() {};
		EstOutput getValues(Tree*, int, string);
		EstOutput getValues(Tree*, string, string, int, string);
		
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
		map< vector<string>, set<int> > rootForGrouping;  //maps a grouping combo to the roots for that combo
		
		EstOutput driver(Tree*, vector< vector<string> >, int, int); 
		EstOutput createProcesses(Tree*, vector< vector<string> >);
		EstOutput driver(Tree*, vector< vector<string> >, int, int, bool); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, bool);
		int getRoot(Tree*, int, vector<string>);
};

/***********************************************************************/

#endif
