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
		Parsimony(TreeMap* t) : tmap(t) {};
		~Parsimony() {};
		EstOutput getValues(Tree*);
		
	private:
		EstOutput data;
		TreeMap* tmap;
};

/***********************************************************************/
#endif