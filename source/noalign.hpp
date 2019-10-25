#ifndef NOALIGN_HPP
#define NOALIGN_HPP

/*
 *  noalign.hpp
 *  
 *
 *  Created by Pat Schloss on 2/19/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"

/**************************************************************************************************/

class NoAlign : public Alignment {
	
public:
	NoAlign();
	~NoAlign();
	void align(string, string, bool createBaseMap);
	
private:	
};

/**************************************************************************************************/


#endif
