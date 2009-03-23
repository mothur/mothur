#ifndef COVERAGE_H
#define COVERAGE_H

/*
 *  coverage.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "fullmatrix.h"
#include "globaldata.hpp"

using namespace std;

/***********************************************************************/

class Coverage  {
	
	public: 
		Coverage();
		~Coverage(){};
		void getValues(FullMatrix*, vector< vector< vector<float> > >&, vector<float>, string);	//matrix, container for results, vector of distances, mode - for random matrices
		void getValues(FullMatrix*, vector< vector< vector<float> > >&, vector<float>); //for user matrix
		
	private:
		GlobalData* globaldata;
		int numGroups, numUserGroups;

};


/***********************************************************************/



#endif