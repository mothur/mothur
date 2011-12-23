
#include "spline.h"


extern"C" {
	void sgram_(double *sg0, double *sg1, double *sg2,
				double *sg3, double* knot, int* nk);
	void stxwx_(double *xs, double *ys, double *ws, int *n,
				double *knot, int *nk, double *xwy, double *hs0, double *hs1, double *hs2, double *hs3);
	void sslvrg_(double *penalt, double *dofoff, double *xs, double *ys, double *ws, double *ssw, int *n,
				 double *knot, int *nk,
				 double *coef, double *sz, double *lev, double *crit, int *icrit, double *lspar, double *xwy,
				 double *hs0, double *hs1, double *hs2, double *hs3,
				 double *sg0, double *sg1, double *sg2, double *sg3, double *abd,
				 double *p1ip, double *p2ip, int *ld4, int *ldnk, int *ier);
}
/***********************************************************/
// used as reference - http://www.koders.com/fortran/fid8AA63B49CF22F0138E9B3DBDC405696F4A62C1CF.aspx
//  http://www.koders.com/c/fidD995301A8A5549CE0361F4E7FFDFD3CDC4B4E4A3.aspx
/* A Cubic B-spline Smoothing routine.
 
 // sbart.f -- translated by f2c (version 20010821).
 // ------- and f2c-clean,v 1.9 2000/01/13
 //
 // According to the GAMFIT sources, this was derived from code by
 // Finbarr O'Sullivan.
 //
 
 
 The algorithm minimises:
 
 (1/n) * sum ws(i)^2 * (ys(i)-sz(i))^2 + lambda* int ( s"(x) )^2 dx
 
 lambda is a function of the spar which is assumed to be between 0 and 1
 
 INPUT
 -----
 penalt			A penalty > 1 to be used in the gcv criterion
 dofoff			either `df.offset' for GCV or `df' (to be matched).
 n				number of data points
 ys(n)			vector of length n containing the observations
 ws(n)			vector containing the weights given to each data point
 NB: the code alters the values here.
 xs(n)			vector containing the ordinates of the observations
 ssw			`centered weighted sum of y^2'
 nk				number of b-spline coefficients to be estimated
 nk <= n+2
 knot(nk+4)		vector of knot points defining the cubic b-spline basis.
 To obtain full cubic smoothing splines one might
 have (provided the xs-values are strictly increasing)
 spar			penalised likelihood smoothing parameter
 ispar			indicating if spar is supplied (ispar=1) or to be estimated
 lspar, uspar	lower and upper values for spar search;  0.,1. are good values
 tol, eps		used in Golden Search routine
 isetup			setup indicator [initially 0
 icrit			indicator saying which cross validation score is to be computed
 0: none ;  1: GCV ;  2: CV ;  3: 'df matching'
 ld4			the leading dimension of abd (ie ld4=4)
 ldnk			the leading dimension of p2ip (not referenced)
 
 OUTPUT
 ------
 coef(nk)		vector of spline coefficients
 sz(n)			vector of smoothed z-values
 lev(n)			vector of leverages
 crit			either ordinary or generalized CV score
 spar			if ispar != 1
 lspar			== lambda (a function of spar and the design)
 iter			number of iterations needed for spar search (if ispar != 1)
 ier			error indicator
 ier = 0 ___  everything fine
 ier = 1 ___  spar too small or too big
 problem in cholesky decomposition
 
 Working arrays/matrix
 xwy				X'Wy
 hs0,hs1,hs2,hs3	the diagonals of the X'WX matrix
 sg0,sg1,sg2,sg3	the diagonals of the Gram matrix SIGMA
 abd (ld4,nk)		[ X'WX + lambda*SIGMA ] in diagonal form
 p1ip(ld4,nk)		inner products between columns of L inverse
 p2ip(ldnk,nk)		all inner products between columns of L inverse
 where  L'L = [X'WX + lambda*SIGMA]  NOT REFERENCED
 */

int Spline::sbart
(double *penalt, double *dofoff,
 double *xs, double *ys, double *ws, double *ssw,
 int *n, double *knot, int *nk, double *coef,
 double *sz, double *lev, double *crit, int *icrit,
 double *spar, int *ispar, int *iter, double *lspar,
 double *uspar, double *tol, double *eps, int *isetup,
 int *ld4, int *ldnk, int *ier)
{
	try{
	/* A Cubic B-spline Smoothing routine.
	 
	 The algorithm minimises:
	 
	 (1/n) * sum ws(i)^2 * (ys(i)-sz(i))^2 + lambda* int ( s(x) )^2 dx
	 
	 lambda is a function of the spar which is assumed to be between 0 and 1
	 
	 INPUT
	 -----
	 penalt	A penalty > 1 to be used in the gcv criterion
	 dofoff	either `df.offset' for GCV or `df' (to be matched).
	 n		number of data points
	 ys(n)	vector of length n containing the observations
	 ws(n)	vector containing the weights given to each data point
	 NB: the code alters the values here.
	 xs(n)	vector containing the ordinates of the observations
	 ssw          `centered weighted sum of y^2
	 nk		number of b-spline coefficients to be estimated
	 nk <= n+2
	 knot(nk+4)	vector of knot points defining the cubic b-spline basis.
	 To obtain full cubic smoothing splines one might
	 have (provided the xs-values are strictly increasing)
	 spar		penalised likelihood smoothing parameter
	 ispar	indicating if spar is supplied (ispar=1) or to be estimated
	 lspar, uspar lower and upper values for spar search;  0.,1. are good values
	 tol, eps	used in Golden Search routine
	 isetup	setup indicator [initially 0
	 icrit	indicator saying which cross validation score is to be computed
	 0: none ;  1: GCV ;  2: CV ;  3: 'df matching'
	 ld4		the leading dimension of abd (ie ld4=4)
	 ldnk		the leading dimension of p2ip (not referenced)
	 
	 OUTPUT
	 ------
	 coef(nk)	vector of spline coefficients
	 sz(n)	vector of smoothed z-values
	 lev(n)	vector of leverages
	 crit		either ordinary or generalized CV score
	 spar         if ispar != 1
	 lspar         == lambda (a function of spar and the design)
	 iter		number of iterations needed for spar search (if ispar != 1)
	 ier		error indicator
	 ier = 0 ___  everything fine
	 ier = 1 ___  spar too small or too big
	 problem in cholesky decomposition
	 
	 Working arrays/matrix
	 xwy			XWy
	 hs0,hs1,hs2,hs3	the diagonals of the XWX matrix
	 sg0,sg1,sg2,sg3	the diagonals of the Gram matrix SIGMA
	 abd (ld4,nk)		[ XWX + lambda*SIGMA ] in diagonal form
	 p1ip(ld4,nk)		inner products between columns of L inverse
	 p2ip(ldnk,nk)	all inner products between columns of L inverse
	 where  LL = [XWX + lambda*SIGMA]  NOT REFERENCED
	 */
	
#define CRIT(FX) (*icrit == 3 ? FX - 3. : FX)
		/* cancellation in (3 + eps) - 3, but still...informative */
		
#define BIG_f (1e100)
		
		/* c_Gold is the squared inverse of the golden ratio */
		static const double c_Gold = 0.381966011250105151795413165634;
		/* == (3. - sqrt(5.)) / 2. */
		
		/* Local variables */
		static double ratio;/* must be static (not needed in R) */
		
		double a, b, d, e, p, q, r, u, v, w, x;
		double ax, fu, fv, fw, fx, bx, xm;
		double t1, t2, tol1, tol2;
		double* xwy = new double[*nk];
		double* hs0 = new double[*nk];
		double* hs1 = new double[*nk];
		double* hs2 = new double[*nk];
		double* hs3 = new double[*nk];
		double* sg0 = new double[*nk];
		double* sg1 = new double[*nk];
		double* sg2 = new double[*nk];
		double* sg3 = new double[*nk];
		double* abd = new double[*nk*(*ld4)]; 
		double* p1ip = new double[*nk*(*ld4)]; 
		double* p2ip = new double[*nk];	
		int i, maxit;

		/* unnecessary initializations to keep  -Wall happy */
		d = 0.; fu = 0.; u = 0.;
		ratio = 1.;
		
		/*  Compute SIGMA, X' W X, X' W z, trace ratio, s0, s1.
		 
		 SIGMA	-> sg0,sg1,sg2,sg3
		 X' W X	-> hs0,hs1,hs2,hs3
		 X' W Z	-> xwy
		 */
		
		/* trevor fixed this 4/19/88
		 * Note: sbart, i.e. stxwx() and sslvrg() {mostly, not always!}, use
		 *	 the square of the weights; the following rectifies that */
		for (i = 0; i < *n; ++i)
			if (ws[i] > 0.)
				ws[i] = sqrt(ws[i]);
		
		if (*isetup == 0) {
			/* SIGMA[i,j] := Int  B''(i,t) B''(j,t) dt  {B(k,.) = k-th B-spline} */
			sgram_(sg0, sg1, sg2, sg3, knot, nk);
			stxwx_(xs, ys, ws, n,
							knot, nk,
							xwy,
							hs0, hs1, hs2, hs3);
			/* Compute ratio :=  tr(X' W X) / tr(SIGMA) */
			t1 = t2 = 0.;
			for (i = 3 - 1; i < (*nk - 3); ++i) {
				t1 += hs0[i];
				t2 += sg0[i];
			}
			ratio = t1 / t2;
			*isetup = 1;
		}
		/*     Compute estimate */
		
		if (*ispar == 1) { /* Value of spar supplied */
			*lspar = ratio * pow(16.0, (*spar * 6.0 - 2.0));
			sslvrg_(penalt, dofoff, xs, ys, ws, ssw, n,
							 knot, nk,
							 coef, sz, lev, crit, icrit, lspar, xwy,
							 hs0, hs1, hs2, hs3,
							 sg0, sg1, sg2, sg3, abd,
							 p1ip, p2ip, ld4, ldnk, ier);
			/* got through check 2 */
			return 0;
		}
		
		/* ELSE ---- spar not supplied --> compute it ! ---------------------------
		 
		 Use Forsythe Malcom and Moler routine to MINIMIZE criterion
		 f denotes the value of the criterion
		 
		 an approximation	x  to the point where	f  attains a minimum  on
		 the interval  (ax,bx)  is determined.
		 */
		ax = *lspar;
		bx = *uspar;
		
		/* INPUT
		 
		 ax	 left endpoint of initial interval
		 bx	 right endpoint of initial interval
		 f	 function subprogram which evaluates  f(x)  for any  x
         in the interval  (ax,bx)
		 tol	 desired length of the interval of uncertainty of the final
         result ( >= 0 )
		 
		 OUTPUT
		 
		 fmin	 abcissa approximating the point where	f  attains a minimum
		 */
		
		/*
		 The method used is a combination of  golden  section  search  and
		 successive parabolic interpolation.	convergence is never much slower
		 than	 that  for  a  fibonacci search.  if  f	 has a continuous second
		 derivative which is positive at the minimum (which is not  at  ax  or
		 bx),	 then  convergence  is	superlinear, and usually of the order of
		 about  1.324....
		 the function  f  is never evaluated at two points closer together
		 than	 eps*abs(fmin) + (tol/3), where eps is	approximately the square
		 root	 of  the  relative  machine  precision.	  if   f   is a unimodal
		 function and the computed values of	 f   are  always  unimodal  when
		 separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
		 the abcissa of the global minimum of	 f  on the interval  ax,bx  with
		 an error less than  3*eps*abs(fmin) + tol.  if   f	is not unimodal,
		 then fmin may approximate a local, but perhaps non-global, minimum to
		 the same accuracy.
		 this function subprogram is a slightly modified	version	 of  the
		 algol  60 procedure	localmin  given in richard brent, algorithms for
		 minimization without derivatives, prentice - hall, inc. (1973).
		 
		 Double	 a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
		 Double	 fu,fv,fw,fx,x
		 */
		
		/*  eps is approximately the square root of the relative machine
		 precision.
		 
		 -	 eps = 1e0
		 - 10	 eps = eps/2e0
		 -	 tol1 = 1e0 + eps
		 -	 if (tol1 > 1e0) go to 10
		 -	 eps = sqrt(eps)
		 R Version <= 1.3.x had
		 eps = .000244     ( = sqrt(5.954 e-8) )
		 -- now eps is passed as argument
		 */
		
		/* initialization */
		
		maxit = *iter;
		*iter = 0;
		a = ax;
		b = bx;
		v = a + c_Gold * (b - a);
		w = v;
		x = v;
		e = 0.;
		*spar = x;
		*lspar = ratio * pow(16.0, (*spar * 6.0 - 2.0));
		sslvrg_(penalt, dofoff, xs, ys, ws, ssw, n,
						 knot, nk,
						 coef, sz, lev, crit, icrit, lspar, xwy,
				hs0, hs1, hs2, hs3,
				sg0, sg1, sg2, sg3, abd,
				p1ip, p2ip, ld4, ldnk, ier);
		fx = *crit;
		fv = fx;
		fw = fx;
		
		/* main loop
		 --------- */
		while(*ier == 0) { /* L20: */
			
			if (m->control_pressed) { return 0; }
			
			xm = (a + b) * .5;
			tol1 = *eps * fabs(x) + *tol / 3.;
			tol2 = tol1 * 2.;
			++(*iter);
			
			
			/* Check the (somewhat peculiar) stopping criterion: note that
			 the RHS is negative as long as the interval [a,b] is not small:*/
			if (fabs(x - xm) <= tol2 - (b - a) * .5 || *iter > maxit)
				goto L_End;
			
			
			/* is golden-section necessary */
			
			if (fabs(e) <= tol1 ||
				/*  if had Inf then go to golden-section */
				fx >= BIG_f || fv >= BIG_f || fw >= BIG_f) goto L_GoldenSect;
			
			/* Fit Parabola */
			
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = (q - r) * 2.;
			if (q > 0.)
				p = -p;
			q = fabs(q);
			r = e;
			e = d;
			
			/* is parabola acceptable?  Otherwise do golden-section */
			
			if (fabs(p) >= fabs(.5 * q * r) ||
				q == 0.)
			/* above line added by BDR;
			 * [the abs(.) >= abs() = 0 should have branched..]
			 * in FTN: COMMON above ensures q is NOT a register variable */
				
				goto L_GoldenSect;
			
			if (p <= q * (a - x) ||
				p >= q * (b - x))			goto L_GoldenSect;
			
			
			
			/* Parabolic Interpolation step */
			d = p / q;
			u = x + d;
			
			/* f must not be evaluated too close to ax or bx */
			if ((u - a < tol2 || b - u < tol2)) {
				d = abs(tol1) * sgn(xm - x);
			}	
			
			goto L50;
			/*------*/
			
		L_GoldenSect: /* a golden-section step */
			
			if (x >= xm)    e = a - x;
			else/* x < xm*/ e = b - x;
			d = c_Gold * e;
			
			
		L50:
			u = x + ((fabs(d) >= tol1) ? d : (abs(tol1)*sgn(d)));
			/*  tol1 check : f must not be evaluated too close to x */
			
			*spar = u;
			*lspar = ratio * pow(16.0, (*spar * 6.0 - 2.0));
			sslvrg_(penalt, dofoff, xs, ys, ws, ssw, n,
							 knot, nk,
							 coef, sz, lev, crit, icrit, lspar, xwy,
					hs0, hs1, hs2, hs3,
					sg0, sg1, sg2, sg3, abd,
					p1ip, p2ip, ld4, ldnk, ier);
			fu = *crit;
			
			if(isnan(fu)) {
				fu = 2. * BIG_f;
			}
			
			/*  update  a, b, v, w, and x */
			
			if (fu <= fx) {
				if (u >= x) a = x; else b = x;
				
				v = w; fv = fw;
				w = x; fw = fx;
				x = u; fx = fu;
			}
			else {
				if (u < x)  a = u; else b = u;
				
				if (fu <= fw || w == x) {		        /* L70: */
					v = w; fv = fw;
					w = u; fw = fu;
				} else if (fu <= fv || v == x || v == w) {	/* L80: */
					v = u; fv = fu;
				}
			}
		}/* end main loop -- goto L20; */
		
	L_End:
		
		*spar = x;
		*crit = fx;
				
		//free memory
		delete [] xwy;
		delete [] hs0;
		delete [] hs1;
		delete [] hs2;
		delete [] hs3;
		delete [] sg0;
		delete [] sg1;
		delete [] sg2;
		delete [] sg3;
		delete [] abd;
		delete [] p1ip;
		delete [] p2ip;
		
		return 0;
		/* sbart */

	}catch(exception& e) {
		m->errorOut(e, "Spline", "sbart");
		exit(1);
	}
}	

/***********************************************************/
