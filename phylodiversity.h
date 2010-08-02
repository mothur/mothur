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


/***********************************************************************/

class PhyloDiversity  {
	
	public:
		PhyloDiversity(TreeMap* t) : tmap(t) { globaldata = GlobalData::getInstance();  m = MothurOut::getInstance(); }
		~PhyloDiversity() {};
		
		//int getValues(Tree*, vector<int>, vector< vector< float> >&);
		
		
	private:
		GlobalData* globaldata;
		MothurOut* m;
		TreeMap* tmap;
};

/***********************************************************************/


#endif

