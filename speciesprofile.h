#ifndef SPECIESPROFILE_H
#define SPECIESPROFILE_H

/*
 *  speciesprofile.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class SpeciesProfile : public Calculator  {
	
public:
	SpeciesProfile() :  Calculator("speciesprofile", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif


