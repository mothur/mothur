/*
 *  bstick.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bstick.h"


/***********************************************************************/
double BStick::invSum(int index, double numSpec)
{
	double sum = 0;
	for(int i = index; i <= numSpec; i++)
		sum += 1/(double)i;
	return sum;
}

RAbundVector BStick::getRAbundVector(SAbundVector* rank){
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

/***************************************************************************/

/***************************************************************************/
EstOutput BStick::getValues(SAbundVector* rank){
	try {
		data.resize(3,0);
		rdata = getRAbundVector(rank);
		double numInd = (double)rdata.getNumSeqs();
		double numSpec = (double)rdata.getNumBins();
		
		double sumExp = 0;
		double sumObs = 0;
		double maxDiff = 0;

		for(int i = 0; i < rdata.size(); i++) {
			sumObs += rdata.get(i);
			sumExp += numInd/numSpec*invSum(i+1,numSpec);
			double diff = fabs(sumObs-sumExp);
			if(diff > maxDiff)
				maxDiff = diff;
		}
		

		data[0] = maxDiff/numInd;
		data[1] = 0.886/sqrt(rdata.size());
		data[2] = 1.031/sqrt(rdata.size());

		/*m->mothurOut(critVal); m->mothurOutEndLine();
		m->mothurOut("If D-Statistic is less than the critical value then the data fits the Broken Stick model w/ 95% confidence.\n\n");*/
		

		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "BStick", "getValues");
		exit(1);
	}
}

/***********************************************************************/


