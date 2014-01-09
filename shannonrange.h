//
//  shannonrange.h
//  Mothur
//
//  Created by SarahsWork on 1/3/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#ifndef Mothur_shannonrange_h
#define Mothur_shannonrange_h

#include "calculator.h"

/***********************************************************************/

class RangeShannon : public Calculator  {
	
public:
	RangeShannon() : Calculator("rangeshannon", 3, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
    EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/rangeshannon"; }
};

/***********************************************************************/



#endif
