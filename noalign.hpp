#ifndef NOALIGN_HPP
#define NOALIGN_HPP

/*
 *  noalign.hpp
 *  
 *
 *  Created by Pat Schloss on 2/19/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
using namespace std;

#include "mothur.h"

/**************************************************************************************************/

class NoAlign : public Alignment {
	
public:
	NoAlign();
	~NoAlign();
	void align(string, string);
	
private:	
};

/**************************************************************************************************/


#endif
