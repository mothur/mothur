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
#include "Sabundvector.hpp"
#include "sharedrabundvector.h"
#include "sharedordervector.h"
#include "datavector.hpp"
#include "globaldata.hpp"
#include "calculator.h"

/***********************************************************************/

class Venn {
	
	public:
		Venn();
		~Venn(){};
	
		void getPic(OrderVector*, vector<Calculator*>);
		void getPic(SharedOrderVector*, vector<Calculator*>);

	private:
		void getSharedVectors(SharedOrderVector*);
		
		SAbundVector* sabund;
		GlobalData* globaldata;
		vector<SharedRAbundVector*> lookup;
		string format, groupComb;
		ofstream outsvg;

			
};
/***********************************************************************/

#endif

