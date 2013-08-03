#ifndef ALIGNNODE
#define ALIGNNODE

/*
 *  alignNode.h
 *  bayesian
 *
 *  Created by Pat Schloss on 10/11/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "taxonomynode.h"

/**************************************************************************************************/

struct thetaAlign {
	thetaAlign() : A(0), T(0), G(0), C(0), gap(0){}
	unsigned int A;
	unsigned int T;
	unsigned int G;
	unsigned int C;
	unsigned int gap;
};

/**************************************************************************************************/

class AlignNode : public TaxonomyNode {
	
public:
	AlignNode(string, int);
	int loadSequence(string&);
	int checkTheta();
    void printTheta();
	double getPxGivenkj_D_j(string& query);	//P(x | k_j, D, j)
	double getSimToConsensus(string& query);
	vector<thetaAlign> getTheta()	{	return theta;	}
	int addThetas(vector<thetaAlign>, int);
	
private:
	vector<thetaAlign> theta;
	vector<unsigned int> columnCounts;
	int alignLength;
};

/**************************************************************************************************/

#endif

