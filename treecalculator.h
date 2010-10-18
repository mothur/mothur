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
#include "mothurout.h"

/* The tree calculator class is the parent class for tree calculators in mothur. */ 

typedef vector<double> EstOutput;

/***********************************************************************/

class TreeCalculator {

public:
	TreeCalculator(){ m = MothurOut::getInstance(); }
	TreeCalculator(string n) : name(n) {};
	
	virtual ~TreeCalculator(){};
	virtual EstOutput getValues(Tree*) { return data; }	
	virtual EstOutput getValues(Tree*, int, string) { return data; }	
	virtual EstOutput getValues(Tree*, string, string) { return data; }
	virtual EstOutput getValues(Tree*, string, string, vector<double>&) { return data; }
	
	virtual string getName()		{	return name;	}
		
protected:
	EstOutput data;
	string name;
	MothurOut* m;

};

/***********************************************************************/

#endif
