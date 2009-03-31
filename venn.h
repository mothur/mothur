#ifndef VENN_H
#define VENN_H
/*
 *  venn.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include "ordervector.hpp"
#include "rabundvector.hpp"
#include "sharedrabundvector.h"
#include "sharedordervector.h"
#include "datavector.hpp"
#include "globaldata.hpp"

/***********************************************************************/

class Venn {
	
	public:
		Venn();
		~Venn(){};
	
		void getPic(OrderVector*);
		void getPic(SharedOrderVector*);

	private:
		void getSharedVectors(SharedOrderVector*);
		
		RAbundVector rabund;
		GlobalData* globaldata;
		vector<SharedRAbundVector*> lookup;
		string format, groupComb;
		ofstream outsvg;

			
};
/***********************************************************************/

#endif

