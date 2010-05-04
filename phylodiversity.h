#ifndef PHYLODIVERSITY_H
#define PHYLODIVERSITY_H


/*
 *  phylodiversity.h
 *  Mothur
 *
 *  Created by westcott on 4/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "treemap.h"
#include "globaldata.hpp"
#include "mothurout.h"

typedef vector<double> EstOutput; 

/***********************************************************************/

class PhyloDiversity  {
	
	public:
		PhyloDiversity(TreeMap* t) : tmap(t) { globaldata = GlobalData::getInstance();  m = MothurOut::getInstance(); }
		~PhyloDiversity() {};
		
		EstOutput getValues(Tree*, vector<int>);
		void setTotalGroupBranchLengths(Tree*);
		
	private:
		GlobalData* globaldata;
		MothurOut* m;
		EstOutput data;
		TreeMap* tmap;
		map<string, float> groupTotals;
};

/***********************************************************************/


#endif

