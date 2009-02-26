#ifndef WEIGHTED_H
#define WEIGHTED_H


/*
 *  weighted.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "treemap.h"

/***********************************************************************/

class Weighted : public TreeCalculator  {
	
	public:
		Weighted(TreeMap* t) : tmap(t) {};
		~Weighted() {};
		EstOutput getValues(Tree*);
		EstOutput getValues(Tree*, string, string);
		
	private:
		GlobalData* globaldata;
		EstOutput data;
		TreeMap* tmap;
		map<string, int>::iterator it;
		map<string, float> WScore; //a score for each group combination i.e. AB, AC, BC.
};

/***********************************************************************/


#endif