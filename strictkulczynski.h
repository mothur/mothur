#ifndef STRICTKULCZYNSKI_H
#define STRICTKULCZYNSKI_H

/*
 *  strictkulczynski.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class StrictKulczynski : public Calculator  {
	
public:
	StrictKulczynski() :  Calculator("strictkulczynski", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif



