#ifndef SHAREDJACKKNIFE_H
#define SHAREDJACKKNIFE_H

/*
 *  sharedjackknife.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/*This class implements the SharedJackknife estimator. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class SharedJackknife : public Calculator  {
	
public:
	SharedJackknife() : numGroups(-1), callCount(0), count(0), currentCallDone(true), Calculator("sharedjackknife", 3, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Sharedjackknife"; }
	
private:
	int numGroups, callCount, count;
	bool currentCallDone;
	vector<RAbundVector*> groups;
	double simpson(vector<int>, double, int);
	double* jackknife();
    double getConfLimit(int row, int col);
};

/***********************************************************************/

#endif

