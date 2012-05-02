/*
 *  kstest.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedkstest.h"

/***********************************************************************/

EstOutput KSTest::getValues(vector<SharedRAbundVector*> shared){
	try {
		data.resize(3,0);

		//Must return shared1 and shared2 to original order at conclusion of kstest
		vector <individual> initData1 = shared[0]->getData();
		vector <individual> initData2 = shared[1]->getData();
		shared[0]->sortD();
		shared[1]->sortD();

		int numNZ1 = shared[0]->numNZ();
		int numNZ2 = shared[1]->numNZ();
		double numInd1 = (double)shared[0]->getNumSeqs();
		double numInd2 = (double)shared[1]->getNumSeqs();
		
		double maxDiff = -1;
		double sum1 = 0;
		double sum2 = 0;
		for(int i = 1; i < shared[0]->getNumBins(); i++)
		{
			sum1 += shared[0]->get(i).abundance;
			sum2 += shared[1]->get(i).abundance;
			double diff = fabs((double)sum1/numInd1 - (double)sum2/numInd2);
			if(diff > maxDiff)
				maxDiff = diff;
		}
		
		double DStatistic = maxDiff*numNZ1*numNZ2;
		double a = pow((double)(numNZ1 + numNZ2)/(numNZ1*numNZ2),.5);
		//double pVal = exp(-2*pow(maxDiff/a,2));
		double critVal = 1.36*a*numNZ1*numNZ2;
		
		shared[0]->setData(initData1);
		shared[1]->setData(initData2);
		
		data[0] = DStatistic;
		data[1] = critVal;
		data[2] = 0;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }

		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "KSTest", "getValues");
		exit(1);
	}
}

/***********************************************************************/
