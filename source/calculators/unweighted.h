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
#include "counttable.h"

/***********************************************************************/

class Unweighted : public TreeCalculator  {
	
	public:
        Unweighted(bool r, vector<string> g);
		~Unweighted() {};
		EstOutput getValues(Tree*, int, string);
		EstOutput getValues(Tree*, vector<vector<int> >&, int, string);
		
	private:
		vector< vector<string> > namesOfGroupCombos;
        vector<string> Groups;
		int processors;
		string outputDir;
		bool includeRoot;
		
		EstOutput createProcesses(Tree*, CountTable*);
		EstOutput createProcesses(Tree*, vector<vector<int> >&, CountTable*);
};

/**************************************************************************************************/

#endif
