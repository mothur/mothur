//
//  sharedjsd.h
//  Mothur
//
//  Created by SarahsWork on 12/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_sharedjsd_h
#define Mothur_sharedjsd_h

#include "calculator.h"

/***********************************************************************/
//Jensen-Shannon divergence (JSD)
class JSD : public Calculator  {
	
public:
	JSD() :  Calculator("jsd", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/JSD"; }
private:
	
};

/***********************************************************************/


#endif
