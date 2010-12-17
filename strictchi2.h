#ifndef STRICTCHI2_H
#define STRICTCHI2_H

/*
 *  strictchi2.h
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class StrictChi2 : public Calculator  {
	
public:
	StrictChi2() :  Calculator("strictchi2", 1, false, true) {};  //the true means this calculator needs all groups to calculate the pair value
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif