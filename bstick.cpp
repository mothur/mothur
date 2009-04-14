/*
 *  bstick.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
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
		
		for(int i = rank->size()-1; i > 0; i--)
		{
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
		data.resize(2,0);
		rdata = getRAbundVector(rank);
		double numInd = (double)rdata.getNumSeqs();
		double numSpec = (double)rdata.getNumBins();
		
		double sumExp = 0;
		double sumObs = 0;
		double maxDiff = 0;

		for(int i = 0; i < rdata.size(); i++)
		{
			sumObs += rdata.get(i);
			sumExp += numInd/numSpec*invSum(i+1,numSpec);
			double diff = fabs(sumObs-sumExp);
			if(diff > maxDiff)
				maxDiff = diff;
		}
		
		double DStatistic = maxDiff/numInd;
		double critVal = 0;
		/*cout << "BStick:\n";
		cout << "D-Statistic = " << DStatistic << "\n";
		cout << "Critical value for 95% confidence interval = ";*/
		if(rdata.size() > 20)
		{
			critVal = .886/sqrt(rdata.size());
		}
		else
		{
			KOSTable table;
			critVal = table.getConfLimit(numSpec);
		}
		/*cout << critVal << "\n";
		cout << "If D-Statistic is less than the critical value then the data fits the Broken Stick model w/ 95% confidence.\n\n";*/
		
		data[0] = DStatistic;
		data[1] = critVal;
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		
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


