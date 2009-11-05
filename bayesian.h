#ifndef BAYESIAN_H
#define BAYESIAN_H

/*
 *  bayesian.h
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "classify.h"

/**************************************************************************************************/

class Bayesian : public Classify {
	
public:
	Bayesian(string, string, string, int, int, int, int, int);
	~Bayesian() {};
	
	string getTaxonomy(Sequence*);
	
	
private:
	
};

/**************************************************************************************************/

#endif

