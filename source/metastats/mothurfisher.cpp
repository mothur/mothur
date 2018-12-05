/*
 *  mothurfisher.cpp
 *  Mothur
 *
 *  Created by westcott on 7/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

//translated to c++ using source code http://www.langsrud.com/stat/fisher.htm as a reference

#include "mothurfisher.h"
/***********************************************************/
double MothurFisher::fexact(double n11_, double n12_, double n21_, double n22_, string o)  {
	try {
		sleft = 0.0; sright = 0.0; sless = 0.0; slarg = 0.0;
        otuLabel = o;
		
		if(n11_<0) n11_ *= -1;
		if(n12_<0) n12_ *= -1;
		if(n21_<0) n21_ *= -1;
		if(n22_<0) n22_ *= -1; 
        
		double n1_ = n11_+n12_;
		double n_1 = n11_+n21_;
		double n   = n11_ +n12_ +n21_ +n22_;
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: fisher:fexact n11_, n1_, n_1, n " + toString(n11_) + " " + toString(n1_) + " " + toString(n_1) + " " + toString(n) + " \n"); }
		exact(n11_,n1_,n_1,n);
		double twotail = sleft+sright;
		
		if(twotail>1) twotail=1;
		double result = twotail;
		return result; 
		
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "fexact");
		exit(1);
	}	
}
/***********************************************************/
double MothurFisher::lngamm(double z) {
	// Reference: "Lanczos, C. 'A precision approximation 
	// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
	// Translation of  Alan Miller's FORTRAN-implementation
	// See http://lib.stat.cmu.edu/apstat/245
	try {
		double x = 0;
		x += 0.1659470187408462e-06/(z+7);
		x += 0.9934937113930748e-05/(z+6);
		x -= 0.1385710331296526    /(z+5);
		x += 12.50734324009056     /(z+4);
		x -= 176.6150291498386     /(z+3);
		x += 771.3234287757674     /(z+2);
		x -= 1259.139216722289     /(z+1);
		x += 676.5203681218835     /(z);
		x += 0.9999999999995183;
		
		return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
		
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "lngamm");
		exit(1);
	}	
}

/***********************************************************/
double MothurFisher::lnfact(double n){
	try {
		if(n <= 1) return(0);
		return(lngamm(n+1));
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "lnfact");
		exit(1);
	}	
}
/***********************************************************/
double MothurFisher::lnbico(double n, double k){
	try {
		return(lnfact(n)-lnfact(k)-lnfact(n-k));
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "lnbico");
		exit(1);
	}
}
/***********************************************************/
double MothurFisher::hyper_323(double n11, double n1_, double n_1, double n){
	try {
		return(exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1)));
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "hyper_323");
		exit(1);
	}
}
/***********************************************************/
double MothurFisher::myhyper(double n11){
	try {
		double hyper0Result = hyper0(n11,0,0,0);
		return hyper0Result;
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "myhyper");
		exit(1);
	}
}
/***********************************************************/
double MothurFisher::hyper0(double n11i, double n1_i, double n_1i, double ni) {
	try {
		if (!( !util.isEqual(n1_i, 0) && !util.isEqual(n_1i,0) && !util.isEqual(ni, 0) )) {
			if(!(((int)n11i % 10) == 0)){
				if(util.isEqual(n11i,sn11+1))
				{
					sprob *= ((sn1_-sn11)/(n11i))*((sn_1-sn11)/(n11i+sn-sn1_-sn_1));
					sn11 = n11i;
					return sprob;
				}
				if(util.isEqual(n11i,sn11-1))
				{
					sprob *= ((sn11)/(sn1_-n11i))*((sn11+sn-sn1_-sn_1)/(sn_1-n11i));
					sn11 = n11i;
					return sprob;
				}
			}
			sn11 = n11i;
		}else{
			sn11 = n11i;
			sn1_=n1_i;
			sn_1=n_1i;
			sn=ni;
		}
		
		sprob = hyper_323(sn11,sn1_,sn_1,sn);
		return sprob;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "hyper0");
		exit(1);
	}
}
/***********************************************************/
double MothurFisher::exact(double n11, double n1_, double n_1, double n){
	try {
		double p,i,j,prob;
		double max=n1_;
		if(n_1<max) max=n_1;
		double min = n1_+n_1-n;
		if(min<0) min=0;
		if(util.isEqual(min,max))
		{
			sless = 1;
			sright= 1;
			sleft = 1;
			slarg = 1;
			return 1;
		}
		prob=hyper0(n11,n1_,n_1,n);
		sleft=0;
		p=myhyper(min);
		for(i=min+1; p<0.99999999*prob; i++)
		{
			sleft += p;
			p=myhyper(i);
            if (i > max) {
                m->mothurOut("[WARNING]: i value too high. Take a closer look at the pvalue for " + otuLabel + ".\n");
                break;
            }
		}
		i--;
		if(p<1.00000001*prob) sleft += p;
		else i--;
		sright=0;
		p=myhyper(max);
		for(j=max-1; p<0.99999999*prob; j--)
		{
			sright += p;
			p=myhyper(j);
            if (j < 0) {
                m->mothurOut("[WARNING]: j value too low. Take a closer look at the pvalue for " + otuLabel + ".\n");
                break;
            }
		}
		j++;
		if(p<1.00000001*prob) sright += p;
		else j++;
		if(abs(i-n11)<abs(j-n11)) 
		{
			sless = sleft;
			slarg = 1 - sleft + prob;
		} 
		else 
		{
			sless = 1 - sright + prob;
			slarg = sright;
		}
		return prob;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurFisher", "exact");
		exit(1);
	}
}
/***********************************************************/



