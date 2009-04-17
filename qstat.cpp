/*
 *  qstat.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "qstat.h"


/***********************************************************************/

EstOutput QStat::getValues(SAbundVector* rank){
	try {
		
		/*test data VVV
		int dstring[] = {0,0,1,4,2,0,2,1,1,1,1,1,0,1,1,2,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
		vector<int> dvec;
		for(int i = 0; i < 171; i++)
			dvec.push_back(dstring[i]);
		int mr = 170;
		int nb = 29;
		int ns = 884;
		SAbundVector rankw = SAbundVector(dvec, mr,nb,ns);
		SAbundVector *rank = &rankw;*/
		data.resize(1,0);
		int numSpec = rank->getNumBins();
		int r1 = -1;
		int r3 = -1;
		int r1Ind = 0;
		int r3Ind = 0;
		int sumSpec = 0;
		int iqSum = 0;
		for(int i = 1; i < rank->size(); i++)
		{
			if(r1 != -1 && r3 != -1)
				i = rank->size();
				
			sumSpec += rank->get(i);
			
			if(r1 == -1 && sumSpec >= numSpec*.25)
			{
				r1 = rank->get(i);
				r1Ind = i;
			}
			else if(r3 == -1 && sumSpec >= numSpec*.75)
			{
				r3 = rank->get(i);
				r3Ind = i;
			}
			else if(sumSpec >= numSpec*.25 && sumSpec < numSpec*.75)
				iqSum += rank->get(i);
		}
		
		double qstat = (.5*r1 + iqSum + .5*r3)/log((double)r3Ind/r1Ind);
		//cout << "QStat:\nQStatistic = " << qstat << "\n\n";
		
		data[0] = qstat;
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the NPShannon class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the NPShannon class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

