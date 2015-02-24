/*
 *  jacknife.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "jackknife.h"

/***********************************************************************/
void Jackknife::getAMatrix(void){
	try {
		vector<vector<double> > B = m->binomial(maxOrder);

		aMat.resize(maxOrder+1);

		for(int i=0;i<=maxOrder;i++){
			aMat[i].resize(maxOrder+1);
			for(int j=1;j<=maxOrder;j++){
			
				aMat[i][j] = 1 + B[i][j] * (int)(pow(-1.0,j+1));
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Jackknife", "getAMatrix");
		exit(1);
	}
}

/**************************************************************************************************/

double Jackknife::CN(double z){
	try {
		if(z>6.0)	{	return 0.0;		}
		if(z<-6.0)	{	return 0.0;		}
	
		const double b1=  0.31938153;
		const double b2= -0.356563782;
		const double b3=  1.781477937;
		const double b4= -1.821255978;
		const double b5=  1.330274429;
		const double p=   0.2316419;
		const double c2=  0.3989423;
	
		double a=abs(z);
		double t=1.0/(1.0+a*p);
		double b=c2*exp((-z)*(z/2.0));
		double n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
		n = 2*b*n;
		return n;
	}
	catch(exception& e) {
		m->errorOut(e, "Jackknife", "CN");
		exit(1);
	}
}

/***********************************************************************/

EstOutput Jackknife::getValues(SAbundVector* rank){
	try {
		//EstOutput jackData(3,0);
		data.resize(3,0);
	
		double jack, jacklci, jackhci;
	
		int maxRank = (double)rank->getMaxRank();
		int S = rank->getNumBins();

		double N[maxOrder+1];
		double variance[maxOrder+1];
		double p[maxOrder+1];
	
		int k = 0;

		for(int i=0;i<=maxOrder;i++){
			N[i]=0.0000;
			variance[i]=0.0000;
			for(int j=1;j<=maxRank;j++){
				if(j<=i){
					N[i] += aMat[i][j]*rank->get(j);
					variance[i] += aMat[i][j]*aMat[i][j]*rank->get(j);
				}
				else{
					N[i] += rank->get(j);
					variance[i] += rank->get(j);
				}
			}
			variance[i] = variance[i]-N[i];
			double var = 0.0000;
			if(i>0){
				for(int j=1;j<=maxRank;j++){
					if(j<=i){	var += rank->get(j)*pow((aMat[i][j]-aMat[i-1][j]),2.0);	}
					else	{	var += 0.0000;	}
				}
				var -= ((N[i]-N[i-1])*(N[i]-N[i-1]))/S;
				var = var * S / (S-1);
				double T = (N[i]-N[i-1])/sqrt(var);
				if(T<=0.00){	p[i-1] = 1.00000;		}
				else{			p[i-1] = CN(T);			}
			
				if(p[i-1]>=0.05){
					k = i-1;
					break;
				}
			}
			if(i == maxOrder){	k=1;	}
		}

		double ci = 0;
	
		if(k>1){
			double c = (0.05-p[k-1])/(p[k]-p[k-1]);
			ci = 0.0000;
			jack = c*N[k]+(1-c)*N[k-1];
			for(int j=1;j<=maxRank;j++){
				if(j<=k){	ci += rank->get(j)*pow((c*aMat[k][j]+(1-c)*aMat[k-1][j]),2.0);	}
				else	{	ci += rank->get(j);	}
			}
			ci = 1.96 * sqrt(ci - jack);
		}
		else if(k==1){
			jack = N[1];
			ci = 1.96*sqrt(variance[1]);
		}else{
			jack = 0.0;
			ci = 0.0;
		}
	
		jacklci = jack-ci;
		jackhci = jack+ci;
		
		data[0] = jack;
		data[1] = jacklci;
		data[2] = jackhci;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
	
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Jackknife", "getValues");
		exit(1);
	}
}

/***********************************************************************/
