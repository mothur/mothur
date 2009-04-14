/*
 *  DLibshuff.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "dlibshuff.h"

/***********************************************************************/

DLibshuff::DLibshuff(FullMatrix* D, int it, float step, float co) : Libshuff(D, it, step, co){

	numDXs = int(cutOff / stepSize);

}

/***********************************************************************/

float DLibshuff::evaluatePair(int i, int j){
	return dCalculate(i,j);
}

/***********************************************************************/

vector<vector<double> > DLibshuff::evaluateAll(){
	savedMins.resize(numGroups);
	vector<vector<double> > dCXYValues(numGroups);
	
	for(int i=0;i<numGroups;i++){
		savedMins[i].resize(numGroups);
		dCXYValues[i].resize(numGroups);
		for(int j=0;j<numGroups;j++){
			if(i!=j){	dCXYValues[i][j] = dCalculate(i,j);	}
			savedMins[i][i] = minX;
			savedMins[i][j] = minXY;
		}
	}
	
	return dCXYValues;
}

/***********************************************************************/

double DLibshuff::dCalculate(int x, int y){
	
	minX = getMinX(x);
	minXY = getMinXY(x, y);
	
	vector<int> nx = calcN(minX);
	vector<int> nxy = calcN(minXY);

	double sum = 0;

	for(int i=0;i<numDXs;i++){
		float h = (nx[i] - nxy[i]) / (float) groupSizes[x];
		sum += h * h * stepSize;
	}

	return sum;
}

/***********************************************************************/

vector<int> DLibshuff::calcN(vector<double> minVector){

	vector<int> counts(numDXs,0);
	
	int precision = int(1 / stepSize);
	
	for(int i=0;i<minVector.size();i++){
		int bin = int (precision * minVector[i]);
		if(bin < numDXs){	counts[bin]++;	}
	}

	for(int i=1;i<numDXs;i++){
		counts[i] += counts[i-1];
	}
	
	return counts;
}

/***********************************************************************/
