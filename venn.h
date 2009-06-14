#ifndef VENN_H
#define VENN_H
/*
 *  venn.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sabundvector.hpp"
#include "sharedrabundvector.h"
#include "datavector.hpp"
#include "globaldata.hpp"
#include "calculator.h"


/***********************************************************************/

class Venn {
public:
	Venn();
	~Venn(){};

	void getPic(SAbundVector*, vector<Calculator*>);
	void getPic(vector<SharedRAbundVector*>, vector<Calculator*>);

private:
	GlobalData* globaldata;
	Calculator* singleCalc;
	string groupComb;
	ofstream outsvg;
};

/***********************************************************************/

#endif

