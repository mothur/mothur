#ifndef LOGSD_H
#define LOGSD_H

/*
 *  logsd.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the LogSD estimator on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class LogSD : public Calculator  {
	
public:
	LogSD() : Calculator("logseries", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) { return data; };

private:
	double logS(double);
	RAbundVector rdata;
};

/***********************************************************************/

#endif





