/*
 *  ace.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "ace.h"

/***********************************************************************/

EstOutput Ace::getValues(SAbundVector* rank) {
	try {
		data.resize(3,0);
		double ace, acelci, acehci;
	
		double nrare = 0;
		double srare = 0;
		double sabund = 0;
	
		double Cace, term1, gamace;
		double numsum = 0;
	
		double maxRank = (double)rank->getMaxRank();
	
		for(int i=1;i<=maxRank;i++){
			if(i<=abund){
				srare += rank->get(i);
				nrare += i*rank->get(i);
				numsum += (i-1)*i*rank->get(i);
			}
			else if(i>abund)	{sabund += rank->get(i);}
		}
		double sobs = srare + sabund;
	
		if (nrare == 0){ Cace = 0.0000; }
		else { Cace = 1.0000 -(double)rank->get(1)/(double)nrare; }
	
		double denom = Cace * (double)(nrare * (nrare-1));
	
		if(denom <= 0.0){	term1=0.0000;	} else {	term1 = (double)(srare * numsum)/(double)denom - 1.0;	}
		if(term1 >= 0.0){	gamace = term1;	} else {	gamace = 0.0;											}
		
		if(gamace >= 0.64){
			gamace = gamace * (1 + (nrare * (1 - Cace) * numsum) / denom);
			if(gamace<0){			gamace = 0;			}
		}
	
		if(Cace == 0.0){
			ace = 0.00;}//ace
		else{
			ace = (double)sabund+((double)srare+(double)rank->get(1)*gamace)/Cace;//ace
		}
	
		/*		
		The following code was obtained from Anne Chao for calculating the SE for her ACE estimator			
		My modification was to reset the frequencies so that a singleton is found in rank[1] insted
		of in rank[0], etc.
		
		I have also added the forumlae to calculate the 95% confidence intervals.
		*/
	
		double j,D_s=0,nn=0,ww=0;
		int Max_Index=rank->getMaxRank()+1;
		double pp, temp1, temp2;
		vector<double> Part_N_Part_F(Max_Index+1,0.0);
    
		for (j=1; j<Max_Index; j++) if(j<=abund) D_s += rank->get(j);
		for (j=1; j<Max_Index; j++){
			if(j<=abund){
				nn += rank->get(j) * j;
				ww += rank->get(j) * j * ( j - 1);
			}
		}
		double C_hat = 1.-rank->get(1)/double(nn);
		double Gamma = ( D_s * ww) / ( C_hat * nn * ( nn - 1.)) - 1.; 
		temp1 = double(nn - rank->get(1));
		temp2 = double(nn - 1.); 
	
		if ( Gamma > 0.){
			Part_N_Part_F[1] =  ( D_s + nn) * ( 1. + rank->get(1) * ww / temp1 / temp2) / temp1 + nn * D_s * ww * ( temp1 - 1.) /
			( temp1 * temp1 * temp2 * temp2) - ( nn + rank->get(1)) / temp1;
			for ( j=2; j<=Max_Index; j++){
				if(j<=abund){
					Part_N_Part_F[j] = ( nn * temp1 - j * rank->get(1) * D_s) / temp1 / temp1 * ( 1. + rank->get(1) * ww / temp1 / temp2)
					+ j * rank->get(1) * D_s * nn * ( ( j - 1.) * temp1 * temp2 - ww * ( temp1 + temp2))
					/ temp1 / temp1 / temp1 / temp2 / temp2 + j * rank->get(1) * rank->get(1) / temp1 / temp1;
				}
			}
		}
		else{
			Part_N_Part_F[1] = ( nn + D_s ) / temp1;
			for ( j=2; j<=Max_Index; j++){
				if(j<=abund){
					Part_N_Part_F[j-1] = ( nn * temp1 - j * rank->get(1) * D_s ) / temp1 / temp1;
				}
			}
		}
		if(Max_Index>abund){
			for ( j=abund+1; j<=Max_Index; j++){				
				Part_N_Part_F[j-1] = 1.;	
			}
		} 
		for ( temp1=0., temp2=0., j=0; j<Max_Index; j++) {
			pp = Part_N_Part_F[j];
			temp1 += pp * rank->get(j);
			temp2 += pp * pp * rank->get(j);
		}
	
		double se = temp2 - temp1 * temp1 / ace;
	
		if(toString(se) == "nan"){
			acelci = ace;
			acehci = ace;
		}	
		else if(ace==0.000){
			acelci = ace;
			acehci = ace;
		}
		else if(ace==sobs){
			double ci = 1.96*pow(se,0.5);
			acelci = ace-ci;					//ace lci
			acehci = ace+ci;					//ace hci
		}else{
			double denom = pow(ace-sobs,2);
			double c = exp(1.96*pow((log(1+se/denom)),0.5));
			acelci = sobs+(ace-sobs)/c;			//ace lci 
			acehci = sobs+(ace-sobs)*c;			//ace hci
		}
		
		data[0] = ace;
		data[1] = acelci;
		data[2] = acehci;	
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Ace", "getValues");
		exit(1);
	}
}

/***********************************************************************/
