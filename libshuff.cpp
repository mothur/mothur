/*
 *  libshuffform.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "libshuff.h"

/***********************************************************************/

void swap(int& i,int& j){	int t = i;  i = j;  j = t;	}

/***********************************************************************/

Libshuff::Libshuff(FullMatrix* D, int it, float step, float co) : matrix(D), iters(it), stepSize(step), cutOff(co){
	try{
		groupNames = matrix->getGroups();
		groupSizes = matrix->getSizes();
		numGroups = matrix->getNumGroups();
		initializeGroups(matrix);
	}
	catch(exception& e) {
		errorOut(e, "Libshuff", "Libshuff");
		exit(1);
	}
}

/***********************************************************************/

void Libshuff::initializeGroups(FullMatrix* matrix){
	try{
		groups.resize(numGroups);
		savedGroups.resize(numGroups);
		
		savedGroups.resize(numGroups);
		for(int i=0;i<numGroups;i++) {
			groups[i].resize(groupSizes[i]);
			savedGroups[i].resize(groupSizes[i]);
		}
		int index=0;
		for(int i=0;i<numGroups;i++){

			for(int j=0;j<groupSizes[i];j++){
				savedGroups[i][j] = groups[i][j] = index++;
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "Libshuff", "initializeGroups");
		exit(1);
	}
}

/***********************************************************************/

vector<vector<vector<double> > > Libshuff::getSavedMins(){
	return savedMins;
}

/***********************************************************************/

vector<double> Libshuff::getMinX(int x){
	try{

		vector<double> minX(groupSizes[x], 0);
		for(int i=0;i<groupSizes[x];i++){
			minX[i] = (groupSizes[x] > 1 ? (i==0 ? matrix->get(groups[x][0], groups[x][1]) : matrix->get(groups[x][i], groups[x][0])) : 0.0); //get the first value in row i of this block
			//minX[i] = matrix->get(groups[x][i], groups[x][0]);
			for(int j=0;j<groupSizes[x];j++){
				if(i != j)	{
					double dx = matrix->get(groups[x][i], groups[x][j]);
					if(dx < minX[i]){	minX[i] = dx;	}
				}
			}
		}
		return minX;
	}
	catch(exception& e) {
		errorOut(e, "Libshuff", "getMinX");
		exit(1);
	}
}

/***********************************************************************/

vector<double> Libshuff::getMinXY(int x, int y){
	try{
		vector<double> minXY(groupSizes[x], 0);

		for(int i=0;i<groupSizes[x];i++){
			minXY[i] = matrix->get(groups[x][i], groups[y][0]);
			for(int j=0;j<groupSizes[y];j++){
				double dxy = matrix->get(groups[x][i], groups[y][j]);
				if(dxy<minXY[i]){	minXY[i] = dxy;	}
			}
		}
		return minXY;
	}
	catch(exception& e) {
		errorOut(e, "Libshuff", "getMinXY");
		exit(1);
	}
}

/***********************************************************************/

void Libshuff::randomizeGroups(int x, int y){
	try{
		int nv = groupSizes[x]+groupSizes[y];
		vector<int> v(nv);
		
		int index=0;
		for(int k=0;k<groupSizes[x];k++)	{	v[index++] = groups[x][k];	}
		for(int k=0;k<groupSizes[y];k++)	{	v[index++] = groups[y][k];	}
		
		for(int k=nv-1;k>0;k--){
			int z = (int)(rand() % k);
			swap(v[z],v[k]);
		}
		
		index=0;
		for(int k=0;k<groupSizes[x];k++)	{	groups[x][k]=v[index++];	}
		for(int k=0;k<groupSizes[y];k++)	{	groups[y][k]=v[index++];	}
	}
	catch(exception& e) {
		errorOut(e, "Libshuff", "randomizeGroups");
		exit(1);
	}
}

/***********************************************************************/

void Libshuff::resetGroup(int x){
	
	for(int k=0;k<groupSizes[x];k++)	{	groups[x][k] = savedGroups[x][k];	}
	
}

/***********************************************************************/
