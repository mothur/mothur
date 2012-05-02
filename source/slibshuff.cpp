/*
 *  slibshuff.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "slibshuff.h"

/***********************************************************************/

SLibshuff::SLibshuff(FullMatrix* D, int it, float co) : Libshuff(D, it, 0, co){}

/***********************************************************************/

float SLibshuff::evaluatePair(int i, int j){
	return sCalculate(i,j);
}

/***********************************************************************/

vector<vector<double> > SLibshuff::evaluateAll(){
	try{
		savedMins.resize(numGroups);
		vector<vector<double> > dCXYValues(numGroups);

		for(int i=0;i<numGroups;i++){
			dCXYValues[i].resize(numGroups);
			savedMins[i].resize(numGroups);
			for(int j=0;j<numGroups;j++){
				if(i!=j){
					dCXYValues[i][j] = sCalculate(i,j);	
					savedMins[i][j] = minXY;
				}

				if(savedMins[i][i].size() == 0){
					savedMins[i][i] = minX;
				}

			}
		}		
		return dCXYValues;
	}
	catch(exception& e) {
		m->errorOut(e, "SLibshuff", "evaluateAll");
		exit(1);
	}
}

/***********************************************************************/

double SLibshuff::sCalculate(int x, int y){
	try{
		double sum = 0.0,t=0.0;
		
		minX = getMinX(x);
		
		if (m->control_pressed) { return sum; }
		
		minXY = getMinXY(x,y);
		
		if (m->control_pressed) { return sum; }

		sort(minX.begin(), minX.end());
		
		if (m->control_pressed) { return sum; }
		
		sort(minXY.begin(), minXY.end());
		
		if (m->control_pressed) { return sum; }

		int ix=0,iy=0;
		while( (ix < groupSizes[x]) && (iy < groupSizes[x]) ) {
			double h = (ix-iy)/double(groupSizes[x]);
			
			if(minX[ix] < minXY[iy]) {
				sum += (minX[ix] - t)*h*h;
				t = minX[ix++];
			}
			else {
				sum += (minXY[iy] - t)*h*h;
				t = minXY[iy++];
			}
			
		}
		
		if(ix < groupSizes[x]) {
			
			while(ix < groupSizes[x]) {
				double h = (ix-iy)/double(groupSizes[x]);
				sum += (minX[ix] - t)*h*h;
				t = minX[ix++];
			}
			
		}
		else {
			
			while(iy < groupSizes[x]) {
				double h = (ix-iy)/double(groupSizes[x]);
				sum += (minXY[iy] - t)*h*h;
				t = minXY[iy++];
			}
			
		}
		
		return sum;
	}
	catch(exception& e) {
		m->errorOut(e, "SLibshuff", "sCalculate");
		exit(1);
	}
}

/***********************************************************************/
