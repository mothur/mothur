//
//  sharedrjsd.h
//  Mothur
//
//  Created by Sarah Westcott on 1/21/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#ifndef Mothur_sharedrjsd_h
#define Mothur_sharedrjsd_h

#include "calculator.h"

/***********************************************************************/
//Jensen-Shannon divergence (JSD)
class RJSD : public Calculator  {
	
public:
	RJSD() :  Calculator("rjsd", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/RJSD"; }
private:
	
};

/***********************************************************************/



#endif
