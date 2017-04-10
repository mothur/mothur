/*
 *  calculator.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/
int VecCalc::sumElements(vector<int> vec){
		return sumElements(vec,0);
}
/***********************************************************************/
int VecCalc::sumElements(vector<int> vec, int index){
	
		int sum = 0;
		for(int i = index; i < vec.size(); i++)
			sum += vec.at(i);
		return sum;
	
}

/***********************************************************************/
double VecCalc::sumElements(vector<double> vec){
		double sum = 0;
		for(int i = 0; i < vec.size(); i++)
			sum += vec.at(i);
		return sum;
	
}
/***********************************************************************/
double VecCalc::sumElements(vector<double> vec, int index){
	
		double sum = 0;
		for(int i = index; i < vec.size(); i++)
			sum += vec.at(i);
		return sum;
	
}
/***********************************************************************/
int VecCalc::numNZ(vector<int> vec){
	
		int numNZ = 0;
		for(int i = 0; i < vec.size(); i++)
			if(vec.at(i) != 0)
				numNZ++;
		return numNZ;
	
}
/***********************************************************************/
double VecCalc::numNZ(vector<double> vec){
	
		double numNZ = 0;
		for(int i = 0; i < vec.size(); i++)
			if(vec.at(i) != 0)
				numNZ++;
		return numNZ;
	}



