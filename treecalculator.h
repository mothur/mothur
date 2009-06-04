#ifndef TREECALCULATOR_H
#define TREECALCULATOR_H

/*
 *  treecalculator.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "tree.h"

/* The tree calculator class is the parent class for tree calculators in mothur. */ 

typedef vector<double> EstOutput;

/***********************************************************************/

class TreeCalculator {

public:
	TreeCalculator(){};
	TreeCalculator(string n) : name(n) {};
	~TreeCalculator(){};
	virtual EstOutput getValues(Tree*) = 0;	
	virtual EstOutput getValues(Tree*, string, string) = 0;
	
	virtual string getName()		{	return name;	}
		
protected:
	EstOutput data;
	string name;

};

/***********************************************************************/

#endif