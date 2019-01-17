//
//  qFinderDMM.cpp
//  pds_dmm
//
//  Created by Patrick Schloss on 11/8/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include "qFinderDMM.h"



/**************************************************************************************************/

qFinderDMM::qFinderDMM(vector<vector<int> > cm, int p) : CommunityTypeFinder() {
    try {
        
        numPartitions = p;
        countMatrix = cm;
        numSamples = (int)countMatrix.size();
        numOTUs = (int)countMatrix[0].size();
        
        findkMeans();
        optimizeLambda();
        
        double change = 1.0000;
        currNLL = 0.0000;
        
        int iter = 0;
        
        while(change > 1.0e-6 && iter < 100){
            
            calculatePiK();
            
            optimizeLambda();
            
            for(int i=0;i<numPartitions;i++){
                weights[i] = 0.0000;
                for(int j=0;j<numSamples;j++){
                    weights[i] += zMatrix[i][j];
                }
               
            }
            
            double nLL = getNegativeLogLikelihood();
            
            change = abs(nLL - currNLL);
            
            currNLL = nLL;
           
            iter++;
        }
        error.resize(numPartitions);
        
        logDeterminant = 0.0000;
        
        LinearAlgebra l;
        
        for(currentPartition=0;currentPartition<numPartitions;currentPartition++){
            
            error[currentPartition].assign(numOTUs, 0.0000);
            
            if (m->getDebug()) { m->mothurOut("current partition = " + toString(currentPartition) + "\n"); }
            
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
void qFinderDMM::printFitData(ofstream& out){
    try {
        out << setprecision (2) << numPartitions << '\t'  << getNLL() << '\t' << getLogDet() << '\t' <<  getBIC() << '\t' << getAIC() << '\t' << laplace << endl;
        return;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "printFitData");
        exit(1);
    }
}
/**************************************************************************************************/
void qFinderDMM::printFitData(ostream& out, double minLaplace){
    try {
        if(laplace < minLaplace){
            out << setprecision (2) << numPartitions << '\t'  << getNLL() << '\t' << getLogDet() << '\t' <<  getBIC() << '\t' << getAIC() << '\t' << laplace << "***" << endl;
        }else {
            out << setprecision (2) << numPartitions << '\t'  << getNLL() << '\t' << getLogDet() << '\t' <<  getBIC() << '\t' << getAIC() << '\t' << laplace << endl;
        }
        
        m->mothurOutJustToLog(toString(numPartitions) + '\t' + toString(getNLL()) + '\t' + toString(getLogDet()) + '\t');
        m->mothurOutJustToLog(toString(getBIC()) + '\t' + toString(getAIC()) + '\t' + toString(laplace));

        return;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "printFitData");
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
    Utils util;
    double disc = b * b - 4 * a * c;
    
    if (util.isEqual(a, 0))/* Handle linear case */
    {
        if (util.isEqual(b, 0))
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
        if (util.isEqual(b, 0))
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
    else if (util.isEqual(disc, 0))
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
            if (m->getControl_pressed()) { break; }
            
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
            if (m->getControl_pressed()) { break; }
            
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

        int bfgsIter = 0;
        double step = 1.0e-6;
        double delta_f = 0.0000;//f-f0;

        vector<double> gradient;
        double f = negativeLogEvidenceLambdaPi(x);
    
        negativeLogDerivEvidenceLambdaPi(x, gradient);

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
        
        
        while(g0norm > 0.001 && bfgsIter++ < maxIter){
            if (m->getControl_pressed()) {  return 0; }

            double f0 = f;
            vector<double> dx(numOTUs, 0.0000);
            
            double alphaOld, alphaNew;

            if(util.isEqual(pNorm, 0) || util.isEqual(g0norm, 0) || util.isEqual(df0, 0)){
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
            
            if(!util.isEqual(dxdg, 0)){
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


        return bfgsIter;
    }
    catch(exception& e){
        m->errorOut(e, "qFinderDMM", "bfgs2_Solver");
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
            if (m->getControl_pressed()) {  return 0; }
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
            if (m->getControl_pressed()) {  return; }
            
            alpha[i] = exp(x[i]);
            store += alpha[i];
            
            derivative[i] = weight * psi(alpha[i]);

            for(int j=0;j<numSamples;j++){
                double X = countMatrix[j][i];
                double alphaX = X + alpha[i];
                
                derivative[i] -= zMatrix[currentPartition][j] * psi(alphaX);
                storeVector[j] += alphaX;
            }
        }

        double sumStore = 0.0000;
        for(int i=0;i<numSamples;i++){
            sumStore += zMatrix[currentPartition][i] * psi(storeVector[i]);
        }
        
        store = weight * psi(store);
        
        df.resize(numOTUs, 0.0000);
        
        for(int i=0;i<numOTUs;i++){
            df[i] = alpha[i] * (nu + derivative[i] - store + sumStore) - eta;
        }
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
            if (m->getControl_pressed()) {  return 0; }
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

void qFinderDMM::optimizeLambda(){    
    try {
        for(currentPartition=0;currentPartition<numPartitions;currentPartition++){
            if (m->getControl_pressed()) {  return; }
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
            if (m->getControl_pressed()) {  return; }
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
                if (m->getControl_pressed()) {  return; }
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
            if (m->getControl_pressed()) {  return 0; }
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
            if (m->getControl_pressed()) {  return 0; }
            
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
                if (m->getControl_pressed()) {  return 0; }
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
