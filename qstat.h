#ifndef QSTAT_H
#define QSTAT_H
/*
 *  qstat.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/4/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the LogSD estimator on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class QStat : public Calculator  {
	
public:
	QStat() : Calculator("qstat", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};

private:
	RAbundVector rdata;
};

/***********************************************************************/

#endif

