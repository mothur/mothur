#ifndef JACKKNIFE_H
#define JACKKNIFE_H

/*
 *  jacknife.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/* This class implements the JackKnife estimator on single group. 
It is a child of the calculator class. */

/***********************************************************************/

class Jackknife : public Calculator  {
	
public:
	Jackknife() : Calculator("jackknife", 3, false) {	getAMatrix(); };
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Jackknife"; }

private:
	static const int maxOrder = 30;
	vector<vector<double> > aMat;

	void getAMatrix();
	double CN(double);
};

/***********************************************************************/

#endif
