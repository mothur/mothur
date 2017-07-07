#ifndef RANDOMNUMBER_H
#define RANDOMNUMBER_H

/*
 *  randomnumber.cpp
 *  
 *
 *  Created by Pat Schloss on 7/6/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "randomnumber.h"
#include <cmath>

/**************************************************************************************************/

RandomNumberGenerator::RandomNumberGenerator(){
//    m->setRandomSeed( (unsigned)time( NULL ) );
}

/**************************************************************************************************/

float RandomNumberGenerator::randomUniform(){
	
	float randUnif = 0.0000;
	
	while(randUnif == 0.0000){
		
		randUnif = rand() / (float)RAND_MAX;
		
	}
	
	return randUnif;
}

/**************************************************************************************************/
//
//Code shamelessly swiped and modified from Numerical Recipes in C++
//
/**************************************************************************************************/

float RandomNumberGenerator::randomExp(){
	
	float randExp = 0.0000;
	
	while(randExp == 0.0000){
		
		randExp = -log(randomUniform());
		
	}
	
	return randExp;
}

/**************************************************************************************************/
//
//Code shamelessly swiped and modified from Numerical Recipes in C++
//
/**************************************************************************************************/

float RandomNumberGenerator::randomNorm(){
	
	float x, y, rsquare;
	
	do{
		x = 2.0 * randomUniform() - 1.0;
		y = 2.0 * randomUniform() - 1.0;
	
		rsquare = x * x + y * y;
	} while(rsquare >= 1.0 || rsquare == 0.0);
	
	float fac = sqrt(-2.0 * log(rsquare)/rsquare);

	return x * fac;
}


/**************************************************************************************************/
/*
 *	Slightly modified version of:
 *
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2005 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    float rgamma(float a, float scale);
 *
 *  DESCRIPTION
 *
 *    Random variates from the gamma distribution.
 *
 *  REFERENCES
 *
 *    [1] Shape parameter a >= 1.  Algorithm GD in:
 *
 *        Ahrens, J.H. and Dieter, U. (1982).
 *        Generating gamma variates by a modified
 *        rejection technique.
 *        Comm. ACM, 25, 47-54.
 *
 *
 *    [2] Shape parameter 0 < a < 1. Algorithm GS in:
 *
 *        Ahrens, J.H. and Dieter, U. (1974).
 *        Computer methods for sampling from gamma, beta,
 *        poisson and binomial distributions.
 *        Computing, 12, 223-246.
 *
 *    Input: a = parameter (mean) of the standard gamma distribution.
 *    Output: a variate from the gamma(a)-distribution
 */

float RandomNumberGenerator::randomGamma(float a)
{
	/* Constants : */
    const static float sqrt32 = 5.656854;
    const static float exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
	float scale = 1.0;
	
    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
     */
    const static float q1 = 0.04166669;
    const static float q2 = 0.02083148;
    const static float q3 = 0.00801191;
    const static float q4 = 0.00144121;
    const static float q5 = -7.388e-5;
    const static float q6 = 2.4511e-4;
    const static float q7 = 2.424e-4;
	
    const static float a1 = 0.3333333;
    const static float a2 = -0.250003;
    const static float a3 = 0.2000062;
    const static float a4 = -0.1662921;
    const static float a5 = 0.1423657;
    const static float a6 = -0.1367177;
    const static float a7 = 0.1233795;
	
    /* State variables [FIXME for threading!] :*/
    static float aa = 0.;
    static float aaa = 0.;
    static float s, s2, d;    /* no. 1 (step 1) */
    static float q0, b, si, c;/* no. 2 (step 4) */
	
    float e, p, q, r, t, u, v, w, x, ret_val;
	
    if (a <= 0.0 || scale <= 0.0){	cout << "error alpha or scale parameter are less than zero." << endl; exit(1);	}
		
    if (a < 1.) { /* GS algorithm for parameters a < 1 */
        e = 1.0 + exp_m1 * a;
        for(;;) {
            p = e * randomUniform();
            if (p >= 1.0) {
                x = -log((e - p) / a);
                if (randomExp() >= (1.0 - a) * log(x))
                    break;
            } else {
                x = exp(log(p) / a);
                if (randomExp() >= x)
                    break;
            }
        }
        return scale * x;
    }
	
    /* --- a >= 1 : GD algorithm --- */
	
    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a != aa) {
        aa = a;
        s2 = a - 0.5;
        s = sqrt(s2);
        d = sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
	 x = (s,1/2) -normal deviate. */
	
    /* immediate acceptance (i) */
    t = randomNorm();
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0)
        return scale * ret_val;
	
    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = randomUniform();
    if (d * u <= t * t * t)
        return scale * ret_val;
	
    /* Step 4: recalculations of q0, b, si, c if necessary */
	
    if (a != aaa) {
        aaa = a;
        r = 1.0 / a;
        q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
               + q2) * r + q1) * r;
		
        /* Approximation depending on size of parameter a */
        /* The constants in the expressions for b, si and c */
        /* were established by numerical experiments */
		
        if (a <= 3.686) {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
        } else if (a <= 13.022) {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
        } else {
            b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
        }
    }
    /* Step 5: no quotient test if x not positive */
	
    if (x > 0.0) {
        /* Step 6: calculation of v and quotient q */
        v = t / (s + s);
        if (fabs(v) <= 0.25)
            q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                      + a3) * v + a2) * v + a1) * v;
        else
            q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
		
		
        /* Step 7: quotient acceptance (q) */
        if (log(1.0 - u) <= q)
            return scale * ret_val;
    }
	
    for(;;) {
        /* Step 8: e = standard exponential deviate
         *      u =  0,1 -uniform deviate
         *      t = (b,si)-float exponential (laplace) sample */
        e = randomExp();
        u = randomUniform();
        u = u + u - 1.0;
        if (u < 0.0)
            t = b - si * e;
        else
            t = b + si * e;
        /* Step  9:  rejection if t < tau(1) = -0.71874483771719 */
        if (t >= -0.71874483771719) {
            /* Step 10:  calculation of v and quotient q */
            v = t / (s + s);
            if (fabs(v) <= 0.25)
                q = q0 + 0.5 * t * t *
				((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
				  + a2) * v + a1) * v;
            else
                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
            /* Step 11:  hat acceptance (h) */
            /* (if q not positive go to step 8) */
            if (q > 0.0) {
                w = expm1(q);
                /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                /* if t is rejected sample again at step 8 */
                if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
                    break;
            }
        }
    } /* for(;;) .. until  `t' is accepted */
    x = s + 0.5 * t;
    return scale * x * x;
}

/**************************************************************************************************/
//
//	essentially swiped from http://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
//
/**************************************************************************************************/

vector<float> RandomNumberGenerator::randomDirichlet(vector<float> alphas){

	int nAlphas = (int)alphas.size();
	vector<float> dirs(nAlphas, 0.0000);
	
	float sum = 0.0000;
	for(int i=0;i<nAlphas;i++){
		dirs[i] = randomGamma(alphas[i]);
		sum += dirs[i];
	}
	
	for(int i=0;i<nAlphas;i++){
		dirs[i] /= sum;
	}
	
	return dirs;
	
}

/**************************************************************************************************/

#endif
