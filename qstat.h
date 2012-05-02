#ifndef QSTAT_H
#define QSTAT_H
/*
 *  qstat.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the q statistic on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class QStat : public Calculator  {
	
public:
	QStat() : Calculator("qstat", 1, false) {};

	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Qstat"; }

private:
	RAbundVector rdata;
};

/***********************************************************************/

#endif

