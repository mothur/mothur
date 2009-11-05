/*
 *  database.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "sequence.hpp"
#include "database.hpp"

/**************************************************************************************************/

Database::Database(){		
	longest = 0;
	numSeqs = 0;
}
/**************************************************************************************************/

Database::~Database(){}

/**************************************************************************************************/

float Database::getSearchScore()	{	return searchScore;		}	//	we're assuming that the search is already done


/**************************************************************************************************/

int Database::getLongestBase()	{	return longest+1;		}	

/**************************************************************************************************/
