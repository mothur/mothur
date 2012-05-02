#ifndef SPLINE
#define SPLINE


/*
 *  spline.h
 *  Mothur
 *
 *  Created by westcott on 12/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

//This class was translated into c++ using c and fortran code from www.koders.com as a reference.

#include "mothurout.h"

class Spline {
	
	public:
		Spline() { m = MothurOut::getInstance(); }
		~Spline() {}
		
		int sbart(double *penalt, double *dofoff,
			   double *xs, double *ys, double *ws, double *ssw,
			   int *n, double *knot, int *nk, double *coef,
			   double *sz, double *lev, double *crit, int *icrit,
			   double *spar, int *ispar, int *iter, double *lspar,
			   double *uspar, double *tol, double *eps, int *isetup,
			   int *ld4, int *ldnk, int *ier);
		
	private:
		MothurOut* m;
};

#endif


