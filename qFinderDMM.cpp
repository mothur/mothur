//
//  qFinderDMM.cpp
//  pds_dmm
//
//  Created by Patrick Schloss on 11/8/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include "qFinderDMM.h"
#include "linearalgebra.h"

#define EPSILON numeric_limits<double>::epsilon()

/**************************************************************************************************/

qFinderDMM::qFinderDMM(vector<vector<int> > cm, int p): countMatrix(cm), numPartitions(p){
    try {
        m = MothurOut::getInstance();
        numSamples = (int)countMatrix.size();
        numOTUs = (int)countMatrix[0].size();
        
        
        kMeans();
        //    cout << "done kMeans" << endl;
        
        optimizeLambda();
        
        
        //    cout << "done optimizeLambda" << endl;
        
        double change = 1.0000;
        currNLL = 0.0000;
        
        int iter = 0;
        
        while(change > 1.0e-6 && iter < 100){
            
            //        cout << "Calc_Z: ";
            calculatePiK();
            
            optimizeLambda();
            
            //        printf("Iter:%d\t",iter);
            
            for(int i=0;i<numPartitions;i++){
                weights[i] = 0.0000;
                for(int j=0;j<numSamples;j++){
                    weights[i] += zMatrix[i][j];
                }
                //            printf("w_%d=%.3f\t",i,weights[i]);
                
            }
            
            double nLL = getNegativeLogLikelihood();
            
            change = abs(nLL - currNLL);
            
            currNLL = nLL;
            
            //        printf("NLL=%.5f\tDelta=%.4e\n",currNLL, change);
            
            iter++;
        }
        
        error.resize(numPartitions);
        
        logDeterminant = 0.0000;
        
        LinearAlgebra l;
        
        for(currentPartition=0;currentPartition<numPartitions;currentPartition++){
            
            error[currentPartition].assign(numOTUs, 0.0000);
            
            if(currentPartition > 0){
                logDeterminant += (2.0 * log(numSamples) - log(weights[currentPartition]));
            }
            vector<vector<double> > hessian = getHessian();
            vector<vector<double> > invHessian = l.getInverse(hessian);
            
            for(int i=0;i<numOTUs;i++){
                logDeterminant += log(abs(hessian[i][i]));
                error[currentPartition][i] = invHessian[i][i];
            }
            
        }
        
        int numParameters = numPartitions * numOTUs + numPartitions - 1;
        laplace = currNLL + 0.5 * logDeterminant - 0.5 * numParameters * log(2.0 * 3.14159);
        bic = currNLL + 0.5 * log(numSamples) * numParameters;
        aic = currNLL + numParameters;
    }
	catch(exception& e) {
		m->errorOut(e, "qFinderDMM", "qFinderDMM");
		exit(1);
	}
}

/**************************************************************************************************/

void qFinderDMM::printZMatrix(string fileName, vector<string> sampleName){
    try {
        ofstream printMatrix;
        m->openOutputFile(fileName, printMatrix); //(fileName.c_str());
        printMatrix.setf(ios::fixed, ios::floatfield);
        printMatrix.setf(ios::showpoint);
        
        for(int i=0;i<numPartitions;i++){   printMatrix << "\tPartition_" << i+1;   }   printMatrix << endl;
        
        for(int i=0;i<numSamples;i++){
            printMatrix << sampleName[i];
            for(int j=0;j<numPartitions;j++){
                printMatrix << setprecision(4) << '\t' << zMatrix[j][i];
            }
            printMatrix << endl;
        }
        printMatrix.close();
    }
	catch(exception& e) {
		m->errorOut(e, "qFinderDMM", "printZMatrix");
		exit(1);
	}
}

/**************************************************************************************************/

void qFinderDMM::printRelAbund(string fileName, vector<string> otuNames){
    try {
        ofstream printRA;
        m->openOutputFile(fileName, printRA); //(fileName.c_str());
        printRA.setf(ios::fixed, ios::floatfield);
        printRA.setf(ios::showpoint);
        
        vector<double> totals(numPartitions, 0.0000);
        for(int i=0;i<numPartitions;i++){
            for(int j=0;j<numOTUs;j++){
                totals[i] += exp(lambdaMatrix[i][j]);
            }
        }
        
        printRA << "Taxon";
        for(int i=0;i<numPartitions;i++){
            printRA << "\tPartition_" << i+1 << '_' << setprecision(4) << totals[i];
            printRA << "\tPartition_" << i+1 <<"_LCI" << "\tPartition_" << i+1 << "_UCI";
        }
        printRA << endl;
        
        for(int i=0;i<numOTUs;i++){
            
            if (m->control_pressed) { break; }
            
            printRA << otuNames[i];
            for(int j=0;j<numPartitions;j++){
                
                if(error[j][i] >= 0.0000){
                    double std = sqrt(error[j][i]);
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i]) / totals[j];
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i] - 2.0 * std) / totals[j];
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i] + 2.0 * std) / totals[j];
                }
                else{
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i]) / totals[j];
                    printRA << '\t' << "NA";
                    printRA << '\t' << "NA";
                }
            }
            printRA << endl;
        }
        
        printRA.close();
    }
	catch(exception& e) {
		m->errorOut(e, "qFinderDMM", "printRelAbund");
		exit(1);
	}
}

/**************************************************************************************************/

// these functions for bfgs2 solver were lifted from the gnu_gsl source code...

/* Find a minimum in x=[0,1] of the interpolating quadratic through
 * (0,f0) (1,f1) with derivative fp0 at x=0.  The interpolating
 * polynomial is q(x) = f0 + fp0 * z + (f1-f0-fp0) * z^2
 */

static double
interp_quad (double f0, double fp0, double f1, double zl, double zh)
{
    double fl = f0 + zl*(fp0 + zl*(f1 - f0 -fp0));
    double fh = f0 + zh*(fp0 + zh*(f1 - f0 -fp0));
    double c = 2 * (f1 - f0 - fp0);       /* curvature */
    
    double zmin = zl, fmin = fl;
    
    if (fh < fmin) { zmin = zh; fmin = fh; } 
    
    if (c > 0)  /* positive curvature required for a minimum */
    {
        double z = -fp0 / c;      /* location of minimum */
        if (z > zl && z < zh) {
            double f = f0 + z*(fp0 + z*(f1 - f0 -fp0));
            if (f < fmin) { zmin = z; fmin = f; };
        }
    }
    
    return zmin;
}

/**************************************************************************************************/

/* Find a minimum in x=[0,1] of the interpolating cubic through
 * (0,f0) (1,f1) with derivatives fp0 at x=0 and fp1 at x=1.
 *
 * The interpolating polynomial is:
 *
 * c(x) = f0 + fp0 * z + eta * z^2 + xi * z^3
 *
 * where eta=3*(f1-f0)-2*fp0-fp1, xi=fp0+fp1-2*(f1-f0). 
 */

double cubic (double c0, double c1, double c2, double c3, double z){
    return c0 + z * (c1 + z * (c2 + z * c3));
}

/**************************************************************************************************/

void check_extremum (double c0, double c1, double c2, double c3, double z,
                     double *zmin, double *fmin){
    /* could make an early return by testing curvature >0 for minimum */
    
    double y = cubic (c0, c1, c2, c3, z);
    
    if (y < *fmin)  
    {
        *zmin = z;  /* accepted new point*/
        *fmin = y;
    }
}

/**************************************************************************************************/

int gsl_poly_solve_quadratic (double a, double b, double c, 
                              double *x0, double *x1)
{
    
    double disc = b * b - 4 * a * c;
    
    if (a == 0) /* Handle linear case */
    {
        if (b == 0)
        {
            return 0;
        }
        else
        {
            *x0 = -c / b;
            return 1;
        };
    }
    
    if (disc > 0)
    {
        if (b == 0)
        {
            double r = fabs (0.5 * sqrt (disc) / a);
            *x0 = -r;
            *x1 =  r;
        }
        else
        {
            double sgnb = (b > 0 ? 1 : -1);
            double temp = -0.5 * (b + sgnb * sqrt (disc));
            double r1 = temp / a ;
            double r2 = c / temp ;
            
            if (r1 < r2) 
            {
                *x0 = r1 ;
                *x1 = r2 ;
            } 
            else 
            {
                *x0 = r2 ;
                *x1 = r1 ;
            }
        }
        return 2;
    }
    else if (disc == 0) 
    {
        *x0 = -0.5 * b / a ;
        *x1 = -0.5 * b / a ;
        return 2 ;
    }
    else
    {
        return 0;
    }
   
}

/**************************************************************************************************/

double interp_cubic (double f0, double fp0, double f1, double fp1, double zl, double zh){
    double eta = 3 * (f1 - f0) - 2 * fp0 - fp1;
    double xi = fp0 + fp1 - 2 * (f1 - f0);
    double c0 = f0, c1 = fp0, c2 = eta, c3 = xi;
    double zmin, fmin;
    double z0, z1;
    
    zmin = zl; fmin = cubic(c0, c1, c2, c3, zl);
    check_extremum (c0, c1, c2, c3, zh, &zmin, &fmin);
    
    {
        int n = gsl_poly_solve_quadratic (3 * c3, 2 * c2, c1, &z0, &z1);
        
        if (n == 2)  /* found 2 roots */
        {
            if (z0 > zl && z0 < zh) 
                check_extremum (c0, c1, c2, c3, z0, &zmin, &fmin);
            if (z1 > zl && z1 < zh) 
                check_extremum (c0, c1, c2, c3, z1, &zmin, &fmin);
        }
        else if (n == 1)  /* found 1 root */
        {
            if (z0 > zl && z0 < zh) 
                check_extremum (c0, c1, c2, c3, z0, &zmin, &fmin);
        }
    }
    
    return zmin;
}

/**************************************************************************************************/

double interpolate (double a, double fa, double fpa,
                    double b, double fb, double fpb, double xmin, double xmax){
    /* Map [a,b] to [0,1] */
    double z, alpha, zmin, zmax;
    
    zmin = (xmin - a) / (b - a);
    zmax = (xmax - a) / (b - a);
    
    if (zmin > zmax)
    {
        double tmp = zmin;
        zmin = zmax;
        zmax = tmp;
    };
    
    if(!isnan(fpb) ){
        z = interp_cubic (fa, fpa * (b - a), fb, fpb * (b - a), zmin, zmax);
    }
    else{
        z = interp_quad(fa, fpa * (b - a), fb, zmin, zmax);
    }

    
    alpha = a + z * (b - a);
    
    return alpha;
}

/**************************************************************************************************/

int qFinderDMM::lineMinimizeFletcher(vector<double>& x, vector<double>& p, double f0, double df0, double alpha1, double& alphaNew, double& fAlpha, vector<double>& xalpha, vector<double>& gradient ){
    try {
        
        double rho = 0.01;
        double sigma = 0.10;
        double tau1 = 9.00;
        double tau2 = 0.05;
        double tau3 = 0.50;
        
        double alpha = alpha1;
        double alpha_prev = 0.0000;
        
        xalpha.resize(numOTUs, 0.0000);
        
        double falpha_prev = f0;
        double dfalpha_prev = df0;
        
        double a = 0.0000;          double b = alpha;
        double fa = f0;             double fb = 0.0000;
        double dfa = df0;           double dfb = 0.0/0.0;
        
        int iter = 0;
        int maxIters = 100;
        while(iter++ < maxIters){
            if (m->control_pressed) { break; }
            
            for(int i=0;i<numOTUs;i++){
                xalpha[i] = x[i] + alpha * p[i];
            }
            
            fAlpha = negativeLogEvidenceLambdaPi(xalpha);
            
            if(fAlpha > f0 + alpha * rho * df0 || fAlpha >= falpha_prev){
                a = alpha_prev;         b = alpha;
                fa = falpha_prev;       fb = fAlpha;
                dfa = dfalpha_prev;     dfb = 0.0/0.0;
                break;
            }
            
            negativeLogDerivEvidenceLambdaPi(xalpha, gradient);
            double dfalpha = 0.0000;
            for(int i=0;i<numOTUs;i++){ dfalpha += gradient[i] * p[i]; }
            
            if(abs(dfalpha) <= -sigma * df0){
                alphaNew = alpha;
                return 1;
            }
            
            if(dfalpha >= 0){
                a = alpha;                  b = alpha_prev;
                fa = fAlpha;                fb = falpha_prev;
                dfa = dfalpha;              dfb = dfalpha_prev;
                break;
            }
            
            double delta = alpha - alpha_prev;
            
            double lower = alpha + delta;
            double upper = alpha + tau1 * delta;
            
            double alphaNext = interpolate(alpha_prev, falpha_prev, dfalpha_prev, alpha, fAlpha, dfalpha, lower, upper);
            
            alpha_prev = alpha;
            falpha_prev = fAlpha;
            dfalpha_prev = dfalpha;
            alpha = alphaNext;
        }
        
        iter = 0;
        while(iter++ < maxIters){
            if (m->control_pressed) { break; }
            
            double delta = b - a;
            
            double lower = a + tau2 * delta;
            double upper = b - tau3 * delta;
            
            alpha = interpolate(a, fa, dfa, b, fb, dfb, lower, upper);
            
            for(int i=0;i<numOTUs;i++){
                xalpha[i] = x[i] + alpha * p[i];
            }
            
            fAlpha = negativeLogEvidenceLambdaPi(xalpha);
            
            if((a - alpha) * dfa <= EPSILON){
                return 0;
            }
            
            if(fAlpha > f0 + rho * alpha * df0 || fAlpha >= fa){
                b = alpha;
                fb = fAlpha;
                dfb = 0.0/0.0;
            }
            else{
                double dfalpha = 0.0000;
                
                negativeLogDerivEvidenceLambdaPi(xalpha, gradient);
                dfalpha = 0.0000;
                for(int i=0;i<numOTUs;i++){ dfalpha += gradient[i] * p[i]; }
                
                if(abs(dfalpha) <= -sigma * df0){
                    alphaNew = alpha;
                    return 1;
                }
                
                if(((b-a >= 0 && dfalpha >= 0) || ((b-a) <= 0.000 && dfalpha <= 0))){
                    b = a;      fb = fa;        dfb = dfa;
                    a = alpha;  fa = fAlpha;    dfa = dfalpha;
                }
                else{
                    a = alpha;
                    fa = fAlpha;
                    dfa = dfalpha;
                }
            }
            
            
        }
        
        return 1;
    }
	catch(exception& e) {
		m->errorOut(e, "qFinderDMM", "lineMinimizeFletcher");
		exit(1);
	}
}

/**************************************************************************************************/

int qFinderDMM::bfgs2_Solver(vector<double>& x){
    try{
//        cout << "bfgs2_Solver" << endl;
        int bfgsIter = 0;
        double step = 1.0e-6;
        double delta_f = 0.0000;//f-f0;

        vector<double> gradient;
        double f = negativeLogEvidenceLambdaPi(x);
        
//        cout << "after negLE" << endl;
        
        negativeLogDerivEvidenceLambdaPi(x, gradient);

//        cout << "after negLDE" << endl;

        vector<double> x0 = x;
        vector<double> g0 = gradient;

        double g0norm = 0;
        for(int i=0;i<numOTUs;i++){
            g0norm += g0[i] * g0[i];
        }
        g0norm = sqrt(g0norm);

        vector<double> p = gradient;
        double pNorm = 0;
        for(int i=0;i<numOTUs;i++){
            p[i] *= -1 / g0norm;
            pNorm += p[i] * p[i];
        }
        pNorm = sqrt(pNorm);
        double df0 = -g0norm;

        int maxIter = 5000;
        
//        cout << "before while" << endl;
        
        while(g0norm > 0.001 && bfgsIter++ < maxIter){
            if (m->control_pressed) {  return 0; }

            double f0 = f;
            vector<double> dx(numOTUs, 0.0000);
            
            double alphaOld, alphaNew;

            if(pNorm == 0 || g0norm == 0 || df0 == 0){
                dx.assign(numOTUs, 0.0000);
                break;
            }
            if(delta_f < 0){
                double delta = max(-delta_f, 10 * EPSILON * abs(f0));
                alphaOld = min(1.0, 2.0 * delta / (-df0));
            }
            else{
                alphaOld = step;
            }
            
            int success = lineMinimizeFletcher(x0, p, f0, df0, alphaOld, alphaNew, f, x, gradient);
            
            if(!success){
                x = x0;
                break;   
            }
            
            delta_f = f - f0;
            
            vector<double> dx0(numOTUs);
            vector<double> dg0(numOTUs);
            
            for(int i=0;i<numOTUs;i++){
                dx0[i] = x[i] - x0[i];
                dg0[i] = gradient[i] - g0[i];
            }
            
            double dxg = 0;
            double dgg = 0;
            double dxdg = 0;
            double dgnorm = 0;
            
            for(int i=0;i<numOTUs;i++){
                dxg += dx0[i] * gradient[i];
                dgg += dg0[i] * gradient[i];
                dxdg += dx0[i] * dg0[i];
                dgnorm += dg0[i] * dg0[i];
            }
            dgnorm = sqrt(dgnorm);
            
            double A, B;
            
            if(dxdg != 0){
                B = dxg / dxdg;
                A = -(1.0 + dgnorm*dgnorm /dxdg) * B + dgg / dxdg;            
            }
            else{
                B = 0;
                A = 0;
            }
            
            for(int i=0;i<numOTUs;i++){     p[i] = gradient[i] - A * dx0[i] - B * dg0[i];   }
            
            x0 = x;
            g0 = gradient;
            

            double pg = 0;
            pNorm = 0.0000;
            g0norm = 0.0000;
            
            for(int i=0;i<numOTUs;i++){
                pg += p[i] * gradient[i];
                pNorm += p[i] * p[i];
                g0norm += g0[i] * g0[i];
            }
            pNorm = sqrt(pNorm);
            g0norm = sqrt(g0norm);
            
            double dir = (pg >= 0.0) ? -1.0 : +1.0;

            for(int i=0;i<numOTUs;i++){ p[i] *= dir / pNorm;    }
            
            pNorm = 0.0000;
            df0 = 0.0000;
            for(int i=0;i<numOTUs;i++){
                pNorm += p[i] * p[i];       
                df0 += p[i] * g0[i];
            }
            
            pNorm = sqrt(pNorm);

        }
//        cout << "bfgsIter:\t" << bfgsIter << endl;

        return bfgsIter;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "bfgs2_Solver");
        exit(1);
    }
}

/**************************************************************************************************/

//can we get these psi/psi1 calculations into their own math class?
//psi calcualtions swiped from gsl library...

static const double psi_cs[23] = {
    -.038057080835217922,
    .491415393029387130, 
    -.056815747821244730,
    .008357821225914313,
    -.001333232857994342,
    .000220313287069308,
    -.000037040238178456,
    .000006283793654854,
    -.000001071263908506,
    .000000183128394654,
    -.000000031353509361,
    .000000005372808776,
    -.000000000921168141,
    .000000000157981265,
    -.000000000027098646,
    .000000000004648722,
    -.000000000000797527,
    .000000000000136827,
    -.000000000000023475,
    .000000000000004027,
    -.000000000000000691,
    .000000000000000118,
    -.000000000000000020
};

static double apsi_cs[16] = {    
    -.0204749044678185,
    -.0101801271534859,
    .0000559718725387,
    -.0000012917176570,
    .0000000572858606,
    -.0000000038213539,
    .0000000003397434,
    -.0000000000374838,
    .0000000000048990,
    -.0000000000007344,
    .0000000000001233,
    -.0000000000000228,
    .0000000000000045,
    -.0000000000000009,
    .0000000000000002,
    -.0000000000000000 
};    

/**************************************************************************************************/

double qFinderDMM::cheb_eval(const double seriesData[], int order, double xx){
    try {
        double d = 0.0000;
        double dd = 0.0000;
        
        double x2 = xx * 2.0000;
        
        for(int j=order;j>=1;j--){
            if (m->control_pressed) {  return 0; }
            double temp = d;
            d = x2 * d - dd + seriesData[j];
            dd = temp;
        }
        
        d = xx * d - dd + 0.5 * seriesData[0];
        return d;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "cheb_eval");
        exit(1);
    }
}

/**************************************************************************************************/

double qFinderDMM::psi(double xx){
    try {
        double psiX = 0.0000;
        
        if(xx < 1.0000){
            
            double t1 = 1.0 / xx;
            psiX = cheb_eval(psi_cs, 22, 2.0*xx-1.0);
            psiX = -t1 + psiX;
            
        }
        else if(xx < 2.0000){
            
            const double v = xx - 1.0;
            psiX = cheb_eval(psi_cs, 22, 2.0*v-1.0);
            
        }
        else{
            const double t = 8.0/(xx*xx)-1.0;
            psiX = cheb_eval(apsi_cs, 15, t);
            psiX += log(xx) - 0.5/xx;
        }
        
        return psiX;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "psi");
        exit(1);
    }
}

/**************************************************************************************************/

/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
    1.00000000000000000000000000000,
    0.083333333333333333333333333333,
    -0.00138888888888888888888888888889,
    0.000033068783068783068783068783069,
    -8.2671957671957671957671957672e-07,
    2.0876756987868098979210090321e-08,
    -5.2841901386874931848476822022e-10,
    1.3382536530684678832826980975e-11,
    -3.3896802963225828668301953912e-13,
    8.5860620562778445641359054504e-15,
    -2.1748686985580618730415164239e-16,
    5.5090028283602295152026526089e-18,
    -1.3954464685812523340707686264e-19,
    3.5347070396294674716932299778e-21,
    -8.9535174270375468504026113181e-23
};

/**************************************************************************************************/

double qFinderDMM::psi1(double xx){
    try {
        
        /* Euler-Maclaurin summation formula
         * [Moshier, p. 400, with several typo corrections]
         */
        
        double s = 2.0000;
        const int jmax = 12;
        const int kmax = 10;
        int j, k;
        const double pmax  = pow(kmax + xx, -s);
        double scp = s;
        double pcp = pmax / (kmax + xx);
        double value = pmax*((kmax+xx)/(s-1.0) + 0.5);
        
        for(k=0; k<kmax; k++) {
            if (m->control_pressed) {  return 0; }
            value += pow(k + xx, -s);
        }
        
        for(j=0; j<=jmax; j++) {
            if (m->control_pressed) {  return 0; }
            double delta = hzeta_c[j+1] * scp * pcp;
            value += delta;
            
            if(fabs(delta/value) < 0.5*EPSILON) break;
            
            scp *= (s+2*j+1)*(s+2*j+2);
            pcp /= (kmax + xx)*(kmax + xx);
        }
        
        return value;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "psi1");
        exit(1);
    }
}

/**************************************************************************************************/

double qFinderDMM::negativeLogEvidenceLambdaPi(vector<double>& x){
    try{
        vector<double> sumAlphaX(numSamples, 0.0000);
        
        double logEAlpha = 0.0000;
        double sumLambda = 0.0000;
        double sumAlpha = 0.0000;
        double logE = 0.0000;
        double nu = 0.10000;
        double eta = 0.10000;
        
        double weight = 0.00000;
        for(int i=0;i<numSamples;i++){
            weight += zMatrix[currentPartition][i];
        }
        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) {  return 0; }
            double lambda = x[i];
            double alpha = exp(x[i]);
            logEAlpha += lgamma(alpha);
            sumLambda += lambda;
            sumAlpha += alpha;
            
            for(int j=0;j<numSamples;j++){
                double X = countMatrix[j][i];
                double alphaX = alpha + X;
                sumAlphaX[j] += alphaX;
                
                logE -= zMatrix[currentPartition][j] * lgamma(alphaX);
            }
        }
        
        logEAlpha -= lgamma(sumAlpha);

        for(int i=0;i<numSamples;i++){
            logE += zMatrix[currentPartition][i] * lgamma(sumAlphaX[i]);
        }

        return logE + weight * logEAlpha + nu * sumAlpha - eta * sumLambda;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "negativeLogEvidenceLambdaPi");
        exit(1);
    }
}

/**************************************************************************************************/

void qFinderDMM::negativeLogDerivEvidenceLambdaPi(vector<double>& x, vector<double>& df){
    try{
//        cout << "\tstart negativeLogDerivEvidenceLambdaPi" << endl;
        
        vector<double> storeVector(numSamples, 0.0000);
        vector<double> derivative(numOTUs, 0.0000);
        vector<double> alpha(numOTUs, 0.0000);
        
        double store = 0.0000;
        double nu = 0.1000;
        double eta = 0.1000;
        
        double weight = 0.0000;
        for(int i=0;i<numSamples;i++){
            weight += zMatrix[currentPartition][i];
        }

        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) {  return; }
//            cout << "start i loop" << endl;
//            
//            cout << i << '\t' << alpha[i] << '\t' << x[i] << '\t' << exp(x[i]) << '\t' << store << endl;
            
            alpha[i] = exp(x[i]);
            store += alpha[i];
            
//            cout << "before derivative" << endl;
            
            derivative[i] = weight * psi(alpha[i]);

//            cout << "after derivative" << endl;

//            cout << i << '\t' << alpha[i] << '\t' << psi(alpha[i]) << '\t' << derivative[i] << endl;

            for(int j=0;j<numSamples;j++){
                double X = countMatrix[j][i];
                double alphaX = X + alpha[i];
                
                derivative[i] -= zMatrix[currentPartition][j] * psi(alphaX);
                storeVector[j] += alphaX;
            }
//            cout << "end i loop" << endl;
        }

        double sumStore = 0.0000;
        for(int i=0;i<numSamples;i++){
            sumStore += zMatrix[currentPartition][i] * psi(storeVector[i]);
        }
        
        store = weight * psi(store);
        
        df.resize(numOTUs, 0.0000);
        
        for(int i=0;i<numOTUs;i++){
            df[i] = alpha[i] * (nu + derivative[i] - store + sumStore) - eta;
//            cout << i << '\t' << df[i] << endl;
        }
//        cout << df.size() << endl;
//        cout << "\tend negativeLogDerivEvidenceLambdaPi" << endl;
    }
    catch(exception& e){
         m->errorOut(e, "qFinderDMM", "negativeLogDerivEvidenceLambdaPi");
        exit(1);
    }
}

/**************************************************************************************************/

double qFinderDMM::getNegativeLogEvidence(vector<double>& lambda, int group){
    try {
        double sumAlpha = 0.0000;
        double sumAlphaX = 0.0000;
        double sumLnGamAlpha = 0.0000;
        double logEvidence = 0.0000;
        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) {  return 0; }
            double alpha = exp(lambda[i]);
            double X = countMatrix[group][i];
            double alphaX = alpha + X;
            
            sumLnGamAlpha += lgamma(alpha);
            sumAlpha += alpha;
            sumAlphaX += alphaX;
            
            logEvidence -= lgamma(alphaX);
        }
        
        sumLnGamAlpha -= lgamma(sumAlpha);
        logEvidence += lgamma(sumAlphaX);
        
        return logEvidence + sumLnGamAlpha;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "getNegativeLogEvidence");
        exit(1);
    }
}

/**************************************************************************************************/

void qFinderDMM::kMeans(){
    try {
        
        vector<vector<double> > relativeAbundance(numSamples);
        vector<vector<double> > alphaMatrix;
        
        alphaMatrix.resize(numPartitions);
        lambdaMatrix.resize(numPartitions);
        for(int i=0;i<numPartitions;i++){
            alphaMatrix[i].assign(numOTUs, 0);
            lambdaMatrix[i].assign(numOTUs, 0);
        }
        
        //get relative abundance
        for(int i=0;i<numSamples;i++){
            if (m->control_pressed) {  return; }
            int groupTotal = 0;
            
            relativeAbundance[i].assign(numOTUs, 0.0);
            
            for(int j=0;j<numOTUs;j++){
                groupTotal += countMatrix[i][j];
            }
            for(int j=0;j<numOTUs;j++){
                relativeAbundance[i][j] = countMatrix[i][j] / (double)groupTotal;
            }
        }
        
        //randomly assign samples into partitions
        zMatrix.resize(numPartitions);
        for(int i=0;i<numPartitions;i++){
            zMatrix[i].assign(numSamples, 0);
        }
        
        for(int i=0;i<numSamples;i++){
            zMatrix[rand()%numPartitions][i] = 1;
        }
        
        double maxChange = 1;
        int maxIters = 1000;
        int iteration = 0;
        
        weights.assign(numPartitions, 0);
        
        while(maxChange > 1e-6 && iteration < maxIters){
            
            if (m->control_pressed) {  return; }
            //calcualte average relative abundance
            maxChange = 0.0000;
            for(int i=0;i<numPartitions;i++){
                
                double normChange = 0.0;
                
                weights[i] = 0;
                
                for(int j=0;j<numSamples;j++){
                    weights[i] += (double)zMatrix[i][j];
                }
                
                vector<double> averageRelativeAbundance(numOTUs, 0);
                for(int j=0;j<numOTUs;j++){
                    for(int k=0;k<numSamples;k++){
                        averageRelativeAbundance[j] += zMatrix[i][k] * relativeAbundance[k][j];
                    }
                }
                
                for(int j=0;j<numOTUs;j++){
                    averageRelativeAbundance[j] /= weights[i];
                    double difference = averageRelativeAbundance[j] - alphaMatrix[i][j];
                    normChange += difference * difference;
                    alphaMatrix[i][j] = averageRelativeAbundance[j];
                }
                
                normChange = sqrt(normChange);
                
                if(normChange > maxChange){ maxChange = normChange; }
            }
            
            
            //calcualte distance between each sample in partition adn the average relative abundance
            for(int i=0;i<numSamples;i++){
                if (m->control_pressed) {  return; }
                
                double normalizationFactor = 0;
                vector<double> totalDistToPartition(numPartitions, 0);
                
                for(int j=0;j<numPartitions;j++){
                    for(int k=0;k<numOTUs;k++){
                        double difference = alphaMatrix[j][k] - relativeAbundance[i][k];
                        totalDistToPartition[j] += difference * difference;
                    }
                    totalDistToPartition[j] = sqrt(totalDistToPartition[j]);
                    normalizationFactor += exp(-50.0 * totalDistToPartition[j]);
                }
                
                
                for(int j=0;j<numPartitions;j++){
                    zMatrix[j][i] = exp(-50.0 * totalDistToPartition[j]) / normalizationFactor;
                }
                
            }
            
            iteration++;
            //        cout << "K means: " << iteration << '\t' << maxChange << endl;
            
        }
        
        //    cout << "Iter:-1";
        for(int i=0;i<numPartitions;i++){
            weights[i] = 0.0000;
            
            for(int j=0;j<numSamples;j++){
                weights[i] += zMatrix[i][j];
            }
            //        printf("\tw_%d=%.3f", i, weights[i]);
        }
        //    cout << endl;
        
        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) {  return; }
            for(int j=0;j<numPartitions;j++){
                if(alphaMatrix[j][i] > 0){
                    lambdaMatrix[j][i] = log(alphaMatrix[j][i]);
                }
                else{
                    lambdaMatrix[j][i] = -10.0;
                }
            }
        }
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "kMeans");
        exit(1);
    }
}

/**************************************************************************************************/

void qFinderDMM::optimizeLambda(){    
    try {
        for(currentPartition=0;currentPartition<numPartitions;currentPartition++){
            if (m->control_pressed) {  return; }
            bfgs2_Solver(lambdaMatrix[currentPartition]);
        }
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "optimizeLambda");
        exit(1);
    }
}

/**************************************************************************************************/

void qFinderDMM::calculatePiK(){
    try {
        vector<double> store(numPartitions);
        
        for(int i=0;i<numSamples;i++){
            if (m->control_pressed) {  return; }
            double sum = 0.0000;
            double minNegLogEvidence =numeric_limits<double>::max();
            
            for(int j=0;j<numPartitions;j++){
                double negLogEvidenceJ = getNegativeLogEvidence(lambdaMatrix[j], i);
                
                if(negLogEvidenceJ < minNegLogEvidence){
                    minNegLogEvidence = negLogEvidenceJ;
                }
                store[j] = negLogEvidenceJ;
            }
            
            for(int j=0;j<numPartitions;j++){
                if (m->control_pressed) {  return; }
                zMatrix[j][i] = weights[j] * exp(-(store[j] - minNegLogEvidence));
                sum += zMatrix[j][i];
            }
            
            for(int j=0;j<numPartitions;j++){
                zMatrix[j][i] /= sum;
            }
            
        }
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "calculatePiK");
        exit(1);
    }
    
}

/**************************************************************************************************/

double qFinderDMM::getNegativeLogLikelihood(){
    try {
        double eta = 0.10000;
        double nu = 0.10000;
        
        vector<double> pi(numPartitions, 0.0000);
        vector<double> logBAlpha(numPartitions, 0.0000);
        
        double doubleSum = 0.0000;
        
        for(int i=0;i<numPartitions;i++){
            if (m->control_pressed) {  return 0; }
            double sumAlphaK = 0.0000;
            
            pi[i] = weights[i] / (double)numSamples;
            
            for(int j=0;j<numOTUs;j++){
                double alpha = exp(lambdaMatrix[i][j]);
                sumAlphaK += alpha;
                
                logBAlpha[i] += lgamma(alpha);
            }
            logBAlpha[i] -= lgamma(sumAlphaK);
        }
        
        for(int i=0;i<numSamples;i++){
            if (m->control_pressed) {  return 0; }
            
            double probability = 0.0000;
            double factor = 0.0000;
            double sum = 0.0000;
            vector<double> logStore(numPartitions, 0.0000);
            double offset = -numeric_limits<double>::max();
            
            for(int j=0;j<numOTUs;j++){
                sum += countMatrix[i][j];
                factor += lgamma(countMatrix[i][j] + 1.0000);
            }
            factor -= lgamma(sum + 1.0);
            
            for(int k=0;k<numPartitions;k++){
                
                double sumAlphaKX = 0.0000;
                double logBAlphaX = 0.0000;
                
                for(int j=0;j<numOTUs;j++){
                    double alphaX = exp(lambdaMatrix[k][j]) + (double)countMatrix[i][j];
                    
                    sumAlphaKX += alphaX;
                    logBAlphaX += lgamma(alphaX);
                }
                
                logBAlphaX -= lgamma(sumAlphaKX);
                
                logStore[k] = logBAlphaX - logBAlpha[k] - factor;
                if(logStore[k] > offset){
                    offset = logStore[k];
                }
                
            }
            
            for(int k=0;k<numPartitions;k++){
                probability += pi[k] * exp(-offset + logStore[k]);
            }
            doubleSum += log(probability) + offset;
            
        }
        
        double L5 = - numOTUs * numPartitions * lgamma(eta);
        double L6 = eta * numPartitions * numOTUs * log(nu);
        
        double alphaSum, lambdaSum;
        alphaSum = lambdaSum = 0.0000;
        
        for(int i=0;i<numPartitions;i++){
            for(int j=0;j<numOTUs;j++){
                if (m->control_pressed) {  return 0; }
                alphaSum += exp(lambdaMatrix[i][j]);
                lambdaSum += lambdaMatrix[i][j];
            }
        }
        alphaSum *= -nu;
        lambdaSum *= eta;
        
        return (-doubleSum - L5 - L6 - alphaSum - lambdaSum);
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "getNegativeLogLikelihood");
        exit(1);
    }


}

/**************************************************************************************************/

vector<vector<double> > qFinderDMM::getHessian(){
    try {
        vector<double> alpha(numOTUs, 0.0000);
        double alphaSum = 0.0000;
        
        vector<double> pi = zMatrix[currentPartition];
        vector<double> psi_ajk(numOTUs, 0.0000);
        vector<double> psi_cjk(numOTUs, 0.0000);
        vector<double> psi1_ajk(numOTUs, 0.0000);
        vector<double> psi1_cjk(numOTUs, 0.0000);
        
        for(int j=0;j<numOTUs;j++){
            
            if (m->control_pressed) {  break; }
            
            alpha[j] = exp(lambdaMatrix[currentPartition][j]);
            alphaSum += alpha[j];
            
            for(int i=0;i<numSamples;i++){
                double X = (double) countMatrix[i][j];
                
                psi_ajk[j] += pi[i] * psi(alpha[j]);
                psi1_ajk[j] += pi[i] * psi1(alpha[j]);
                
                psi_cjk[j] += pi[i] * psi(alpha[j] + X);
                psi1_cjk[j] += pi[i] * psi1(alpha[j] + X);
            }
        }
        
        
        double psi_Ck = 0.0000;
        double psi1_Ck = 0.0000;
        
        double weight = 0.0000;
        
        for(int i=0;i<numSamples;i++){
            if (m->control_pressed) {  break; }
            weight += pi[i];
            double sum = 0.0000;
            for(int j=0;j<numOTUs;j++){     sum += alpha[j] + countMatrix[i][j];    }
            
            psi_Ck += pi[i] * psi(sum);
            psi1_Ck += pi[i] * psi1(sum);
        }
        
        double psi_Ak = weight * psi(alphaSum);
        double psi1_Ak = weight * psi1(alphaSum);
        
        vector<vector<double> > hessian(numOTUs);
        for(int i=0;i<numOTUs;i++){ hessian[i].assign(numOTUs, 0.0000); }
        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) {  break; }
            double term1 = -alpha[i] * (- psi_ajk[i] + psi_Ak + psi_cjk[i] - psi_Ck);
            double term2 = -alpha[i] * alpha[i] * (-psi1_ajk[i] + psi1_Ak + psi1_cjk[i] - psi1_Ck);
            double term3 = 0.1 * alpha[i];
            
            hessian[i][i] = term1 + term2 + term3;
            
            for(int j=0;j<i;j++){   
                hessian[i][j] = - alpha[i] * alpha[j] * (psi1_Ak - psi1_Ck);
                hessian[j][i] = hessian[i][j];
            }
        }
        
        return hessian;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "getHessian");
        exit(1);
    }
}

/**************************************************************************************************/
