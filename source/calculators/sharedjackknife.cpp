/*
 *  sharedjackknife.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedjackknife.h"

/***************************************************************************************

***************************************************************************************/
double SharedJackknife::simpson(vector<int> abunds, double numInd, int numBins){
	double denom = numInd*(numInd-1);
	double sum = 0;
	for(int i = 0; i < numBins; i++)
		sum += (double)abunds[i]*((double)abunds[i]-1)/denom;
	
	return sum;
}

/*****************************************************************************************/

double* SharedJackknife::jackknife(){		
	int numBins = groups.at(0)->getNumBins()-1;
	vector<int> cArray(numBins);
	for(int i = 0; i < numBins; i++)
		cArray[i] = 0;

	double numInd = 0;
	for(int i = 0; i < numGroups; i++)
		for(int j = 0; j < numBins; j++) {
			int curAbund = groups.at(i)->get(j+1);
			cArray[j] += curAbund;
			numInd += (double)curAbund;
		}

	double baseD = 1/simpson(cArray, numInd, numBins);
	
	vector<double> pseudoVals(numBins);
	double jackknifeEstimate = 0;
	for(int i = 0; i < numGroups; i++) {
		for(int j = 0; j < numBins-1; j++) {
			int abundDiff = -groups.at(i)->get(j+1);
			if(i > 0)
				abundDiff += groups.at(i-1)->get(j+1);

			cArray[j] += abundDiff;
			numInd += abundDiff;	
		}
		
		double curD = 1/simpson(cArray, numInd, numBins);
		pseudoVals[i] = (double)numGroups*(baseD - curD) + curD;
		jackknifeEstimate += pseudoVals[i];
	}
	jackknifeEstimate /= (double)numGroups;
		
	double variance = 0;
	for(int i = 0; i < numGroups; i++)
		variance += pow(pseudoVals[i]-jackknifeEstimate, 2);
		
	variance /= (double)numGroups*((double)numGroups-1);
	double stErr = sqrt(variance);
	double confLimit = 0;
	if(numGroups <= 30)
		confLimit = getConfLimit(numGroups-1, 1);
	else
		confLimit = 1.645;
	
	confLimit *= stErr;
	
	double* rdata = new double[3];
	rdata[0] = baseD;
	rdata[1] = jackknifeEstimate - confLimit;
	rdata[2] = jackknifeEstimate + confLimit;
	
	return rdata;
}

/***********************************************************************/
double SharedJackknife::getConfLimit(int row, int col) //Rows are the degrees of freedom
{
    //Found on http://www.vgtu.lt/leidiniai/elektroniniai/Probability.pdf/Table%203.pdf
    
    //Confidence Level        .90    .95     .975     .99    .995     .999    .9995
    double values[30][7] = {{3.078, 6.314,	12.706,	31.821,	63.656,	318.289, 636.578},
        {1.886,	2.920,	4.303,	6.965,	9.925,	22.328,	31.600},
        {1.638,	2.353,	3.182,	4.541,	5.841,	10.214,	12.924},
        {1.533,	2.132,	2.776,	3.747,	4.604,	7.173,	8.610},
								{1.476,	2.015,	2.571,	3.365,	4.032,	5.894,	6.869},
								{1.440,	1.943,	2.447,	3.143,	3.707,	5.208,	5.959},
								{1.415,	1.895,	2.365,	2.998,	3.499,	4.785,	5.408},
								{1.397,	1.860,	2.306,	2.896,	3.355,	4.501,	5.041},
								{1.383,	1.833,	2.262,	2.821,	3.250,	4.297,	4.781},
								{1.372,	1.812,	2.228,	2.764,	3.169,	4.144,	4.587},
								{1.363,	1.796,	2.201,	2.718,	3.106,	4.025,	4.437},
								{1.356,	1.782,	2.179,	2.681,	3.055,	3.930,	4.318},
								{1.350,	1.771,	2.160,	2.650,	3.012,	3.852,	4.221},
								{1.345,	1.761,	2.145,	2.624,	2.977,	3.787,	4.140},
								{1.341,	1.753,	2.131,	2.602,	2.947,	3.733,	4.073},
								{1.337,	1.746,	2.120,	2.583,	2.921,	3.686,	4.015},
								{1.333,	1.740,	2.110,	2.567,	2.898,	3.646,	3.965},
								{1.330,	1.734,	2.101,	2.552,	2.878,	3.610,	3.922},
								{1.328,	1.729,	2.093,	2.539,	2.861,	3.579,	3.883},
								{1.325,	1.725,	2.086,	2.528,	2.845,	3.552,	3.850},
								{1.323,	1.721,	2.080,	2.518,	2.831,	3.527,	3.819},
								{1.321,	1.717,	2.074,	2.508,	2.819,	3.505,	3.792},
								{1.319,	1.714,	2.069,	2.500,	2.807,	3.485,	3.768},
								{1.318,	1.711,	2.064,	2.492,	2.797,	3.467,	3.745},
								{1.316,	1.708,	2.060,	2.485,	2.787,	3.450,	3.725},
								{1.315,	1.706,	2.056,	2.479,	2.779,	3.435,	3.707},
								{1.314,	1.703,	2.052,	2.473,	2.771,	3.421,	3.689},
								{1.313,	1.701,	2.048,	2.467,	2.763,	3.408,	3.674},
								{1.311,	1.699,	2.045,	2.462,	2.756,	3.396,	3.660},
								{1.310,	1.697,	2.042,	2.457,	2.750,	3.385,	3.646}};
    
    return values[row][col];
    
}

/***********************************************************************/

/************************************************************************************************/

EstOutput SharedJackknife::getValues(vector<SharedRAbundVector*> vectorShared){ //Fix this for collect, mistake was that it was made with summary in mind.
	try {
		//if(numGroups == -1) { numGroups = m->getNumGroups(); }

		if(callCount == numGroups*(numGroups-1)/2) {
			currentCallDone = true;
			callCount = 0;
		}
		callCount++;

		if(currentCallDone) {
			groups.clear();	
			currentCallDone = false;
		}
		
		if(groups.size() != numGroups) {	
            if(groups.size() == 0)
				groups.push_back(vectorShared[0]);
            groups.push_back(vectorShared[1]);
            
		}
		
		if(groups.size() == numGroups && callCount < numGroups) {
			data.resize(3,0);

			double* rdata = jackknife();
			data[0] = rdata[0];
			data[1] = rdata[1];
			data[2] = rdata[2];
		
			if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
			if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
			if (isnan(data[2]) || isinf(data[0])) { data[2] = 0; }
			
			return data;
		}
		
		data.resize(3,0);
		data[0] = 0;
		data[1] = 0;
		data[2] = 0;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		return data;	
	}
		
	catch(exception& e) {
		m->errorOut(e, "SharedJackknife", "getValues");
		exit(1);
	}
}

/***********************************************************************/

