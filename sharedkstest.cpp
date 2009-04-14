/*
 *  kstest.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sharedkstest.h"

/***********************************************************************/

EstOutput SharedKSTest::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	try {
		data.resize(2,0);
		
		//Must return shared1 and shared2 to original order at conclusion of kstest
		vector <individual> initData1 = shared1->getData();
		vector <individual> initData2 = shared2->getData();
		shared1->sortD();
		shared2->sortD();

		int numNZ1 = shared1->numNZ();
		int numNZ2 = shared2->numNZ();
		double numInd1 = (double)shared1->getNumSeqs();
		double numInd2 = (double)shared2->getNumSeqs();
		
		double maxDiff = -1;
		double sum1 = 0;
		double sum2 = 0;
		for(int i = 1; i < shared1->getNumBins(); i++)
		{
			sum1 += shared1->get(i).abundance;
			sum2 += shared2->get(i).abundance;
			double diff = fabs((double)sum1/numInd1 - (double)sum2/numInd2);
			if(diff > maxDiff)
				maxDiff = diff;
		}
		
		double DStatistic = maxDiff*numNZ1*numNZ2;
		double a = pow((double)(numNZ1 + numNZ2)/(numNZ1*numNZ2),.5);
		//double pVal = exp(-2*pow(maxDiff/a,2));
		double critVal = 1.36*a*numNZ1*numNZ2;
		
		/*cout << "Kolmogorov-Smirnov 2-sample test:\n";
		if(numNZ1 > 25 && numNZ2 > 25) //If the sample sizes are both bigger than 25.
			cout << "P-Value = " << pVal << "\nP-Value is the probability that the data sets are significantly different.\n";
		else
		{	
			//cout << "90% Confidence Critical Value = " << 1.22*a*numNZ1*numNZ2 << "\n";
			cout << "D-Statistic = " << DStatistic << "\n";
			cout << "95% Confidence Critical Value = " << critVal << "\n";
			cout << "If D-Statistic is greater than the critical value then the data sets are significantly different at the 95% confidence level.\n\n";
		}*/
		
		shared1->setData(initData1);
		shared2->setData(initData2);
		
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
