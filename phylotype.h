#ifndef PHYLOTYPE_H
#define PHYLOTYPE_H

/*
 *  phylotype.h
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "classify.h"

/**************************************************************************************************/

class PhyloType : public Classify {
	
public:
	PhyloType(string, string, string, int, int, int, int, int);
	~PhyloType() {};
	
	string getTaxonomy(Sequence*);
	
	
private:
	
};

/**************************************************************************************************/

#endif

