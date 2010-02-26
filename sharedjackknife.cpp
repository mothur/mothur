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
			int curAbund = groups.at(i)->get(j+1).abundance;
			cArray[j] += curAbund;
			numInd += (double)curAbund;
		}

	double baseD = 1/simpson(cArray, numInd, numBins);
	
	vector<double> pseudoVals(numBins);
	double jackknifeEstimate = 0;
	for(int i = 0; i < numGroups; i++) {
		for(int j = 0; j < numBins-1; j++) {
			int abundDiff = -groups.at(i)->get(j+1).abundance;
			if(i > 0)
				abundDiff += groups.at(i-1)->get(j+1).abundance;

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
	TDTable table;
	double confLimit = 0;
	if(numGroups <= 30)
		confLimit = table.getConfLimit(numGroups-1, 1);
	else
		confLimit = 1.645;
	
	confLimit *= stErr;
	
	double* rdata = new double[3];
	rdata[0] = baseD;
	rdata[1] = jackknifeEstimate - confLimit;
	rdata[2] = jackknifeEstimate + confLimit;
	
	return rdata;
}

/************************************************************************************************
************************************************************************************************/

EstOutput SharedJackknife::getValues(vector<SharedRAbundVector*> vectorShared){ //Fix this for collect, mistake was that it was made with summary in mind.
	try {
		SharedRAbundVector* shared1 = vectorShared[0];
		SharedRAbundVector* shared2 = vectorShared[1];
		if(numGroups == -1) {
			globaldata = GlobalData::getInstance();
			numGroups = globaldata->Groups.size();
		}

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
				groups.push_back(shared1);
			groups.push_back(shared2);
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

