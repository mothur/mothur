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
	Venn(string);
	~Venn(){};

	vector<string> getPic(SAbundVector*, vector<Calculator*>);
	vector<string> getPic(vector<SharedRAbundVector*>, vector<Calculator*>);

private:
	GlobalData* globaldata;
	Calculator* singleCalc;
	string groupComb, outputDir;
	ofstream outsvg;
	MothurOut* m;
};

/***********************************************************************/

#endif

