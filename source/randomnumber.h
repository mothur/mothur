#ifndef RANDOMNUMBER
#define RANDOMNUMBER

/*
 *  randomnumber.h
 *  
 *
 *  Created by Pat Schloss on 7/6/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */


/**************************************************************************************************/

#include "mothur.h"

/**************************************************************************************************/

class RandomNumberGenerator {
	
public:
	RandomNumberGenerator();
    float randomUniform();
	float randomExp();
	float randomNorm();
	float randomGamma(float);
	vector<float> randomDirichlet(vector<float> alphas);
	
};

/**************************************************************************************************/

#endif
