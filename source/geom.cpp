/*
 *  geom.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "geom.h"

/***********************************************************************/

double Geom::kEq(double k, double spec){
	return k/(1-k)*pow(1-k, spec)/(1-pow(1-k, spec));
}

RAbundVector Geom::getRAbundVector(SAbundVector* rank){
		vector <int> rData;
		int mr = 1;
		int nb = 0;
		int ns = 0;
		
		for(int i = rank->size()-1; i > 0; i--) {
			int cur = rank->get(i);
			if(mr == 1 && cur > 0)
				mr = i;
			nb += cur;
			ns += i*cur;
			for(int j = 0; j < cur; j++)
				rData.push_back(i);
		}
		
		RAbundVector rav = RAbundVector(rData, mr, nb, ns);
		rav.setLabel(rank->getLabel());
		return rav;
}

/***********************************************************************************/

/***********************************************************************************/
EstOutput Geom::getValues(SAbundVector* rank){
	try {
		data.resize(3,0);
		
		rdata = getRAbundVector(rank);
		double numInd = rdata.getNumSeqs();
		double numSpec = rdata.getNumBins();
		double min = rdata.get(rdata.size()-1);
		double k = .5;
		double step = .49999;
		
		while(fabs(min - numInd*kEq(k, (double)numSpec)) > .0001) { //This uses a binary search to find the value of k.
			if(numInd*kEq(k, numSpec) > min)
				k += step;
			else
				k -= step;
			step /= 2;
		}
		
		double cK = 1/(1-pow(1-k, numSpec));
		double sumExp = 0;
		double sumObs = 0;
		double maxDiff = 0;
		
		for(int i = 0; i < numSpec; i++)
		{
			sumObs += rdata.get(i);
			sumExp += numInd*cK*k*pow(1-k, i);
			
			double diff = fabs(sumObs-sumExp);
			if(diff > maxDiff)	{	maxDiff = diff;		}
			
		}


		data[0] = maxDiff/numInd;
		data[1] = 0.886/sqrt(numSpec);
		data[2] = 1.031/sqrt(numSpec);

		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Geom", "getValues");
		exit(1);
	}
}

/***********************************************************************/







