/*
 *  linearalgebra.cpp
 *  mothur
 *
 *  Created by westcott on 1/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "linearalgebra.h"
#include "wilcox.h"

#define PI 3.1415926535897932384626433832795

// This class references functions used from "Numerical Recipes in C++" //

/*********************************************************************************************************************************/
inline double SQR(const double a)
{
    return a*a;
}
/*********************************************************************************************************************************/

inline double SIGN(const double a, const double b)
{
    return b>=0 ? (a>=0 ? a:-a) : (a>=0 ? -a:a);
}
/*********************************************************************************************************************************/
//NUmerical recipes pg. 245 - Returns the complementary error function erfc(x) with fractional error everywhere less than 1.2 × 10−7.
double LinearAlgebra::erfcc(double x){
    try {
        double t,z,ans;
        z=fabs(x);
        t=1.0/(1.0+0.5*z); 
        
        ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
            t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
            t*(-0.82215223+t*0.17087277))))))))); 
        
        //cout << "in erfcc " << t << '\t' << ans<< endl;
        return (x >= 0.0 ? ans : 2.0 - ans);
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "betai");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 232
double LinearAlgebra::betai(const double a, const double b, const double x) {
    try {
        double bt;
        double result = 0.0;
        
        if (x < 0.0 || x > 1.0) { m->mothurOut("[ERROR]: bad x in betai.\n"); m->control_pressed = true; return 0.0; }
        
        if (x == 0.0 || x == 1.0)  { bt = 0.0; }
        else { bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));  }
        
        if (x < (a+1.0) / (a+b+2.0)) { result = bt*betacf(a,b,x)/a; }
        else { result = 1.0-bt*betacf(b,a,1.0-x)/b; }
        
        return result;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "betai");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 219
double LinearAlgebra::gammln(const double xx) {
    try {
        int j;
        double x,y,tmp,ser;
        static const double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
            -1.231739572450155,0.120858003e-2,-0.536382e-5};
        
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<6;j++) {
            ser += cof[j]/++y;
        }
        return -tmp+log(2.5066282746310005*ser/x);
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "gammln");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 223
double LinearAlgebra::gammp(const double a, const double x) {
    try {
        double gamser,gammcf,gln;
        
        if (x < 0.0 || a <= 0.0) { m->mothurOut("[ERROR]: Invalid arguments in routine GAMMP\n"); m->control_pressed = true; return 0.0;}
        if (x < (a+1.0)) {
            gser(gamser,a,x,gln);
            return gamser;
        } else {
            gcf(gammcf,a,x,gln);
            return 1.0-gammcf;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "gammp");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 223
/*double LinearAlgebra::gammq(const double a, const double x) {
    try {
        double gamser,gammcf,gln;
        
        if (x < 0.0 || a <= 0.0) { m->mothurOut("[ERROR]: Invalid arguments in routine GAMMQ\n"); m->control_pressed = true; return 0.0; }
        if (x < (a+1.0)) {
            gser(gamser,a,x,gln);
            return 1.0-gamser;
        } else {
            gcf(gammcf,a,x,gln);
            return gammcf;
        }   
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "gammq");
		exit(1);
	}
}
*********************************************************************************************************************************/
//Numerical Recipes pg. 224
double LinearAlgebra::gcf(double& gammcf, const double a, const double x, double& gln){
    try {
        const int ITMAX=100;
        const double EPS=numeric_limits<double>::epsilon();
        const double FPMIN=numeric_limits<double>::min()/EPS;
        int i;
        double an,b,c,d,del,h;
        
        gln=gammln(a);
        b=x+1.0-a;
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=ITMAX;i++) {
            an = -i*(i-a);
            b += 2.0;
            d=an*d+b;
            if (fabs(d) < FPMIN) { d=FPMIN; }
            c=b+an/c;
            if (fabs(c) < FPMIN) { c=FPMIN; }
            d=1.0/d;
            del=d*c;
            h *= del;
            if (fabs(del-1.0) <= EPS) break;
        }
        if (i > ITMAX)  { m->mothurOut("[ERROR]: " + toString(a) + " too large, ITMAX=100 too small in gcf\n"); m->control_pressed = true; }
        gammcf=exp(-x+a*log(x)-gln)*h;
        
        return 0.0;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "gcf");
		exit(1);
	}

}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 223
double LinearAlgebra::gser(double& gamser, const double a, const double x, double& gln) {
    try {
        int n;
        double sum,del,ap;
        const double EPS = numeric_limits<double>::epsilon();
        
        gln=gammln(a);
        if (x <= 0.0) { 
            if (x < 0.0) {  m->mothurOut("[ERROR]: x less than 0 in routine GSER\n"); m->control_pressed = true;  }
            gamser=0.0; return 0.0;
        } else {
            ap=a;
            del=sum=1.0/a;
            for (n=0;n<100;n++) {
                ++ap;
                del *= x/ap;
                sum += del;
                if (fabs(del) < fabs(sum)*EPS) {
                    gamser=sum*exp(-x+a*log(x)-gln);
                    return 0.0;
                }
            }
            
            m->mothurOut("[ERROR]: a too large, ITMAX too small in routine GSER\n");
            return 0.0;
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "gser");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Numerical Recipes pg. 233
double LinearAlgebra::betacf(const double a, const double b, const double x) {
    try {
        const int MAXIT = 100;
        const double EPS = numeric_limits<double>::epsilon();
        const double FPMIN = numeric_limits<double>::min() / EPS;
        int m1, m2;
        double aa, c, d, del, h, qab, qam, qap;
        
        qab=a+b;
        qap=a+1.0;
        qam=a-1.0;
        c=1.0;
        d=1.0-qab*x/qap;
        if (fabs(d) < FPMIN) d=FPMIN;
        d=1.0/d;
        h=d;
        for (m1=1;m1<=MAXIT;m1++) {
            m2=2*m1;
            aa=m1*(b-m1)*x/((qam+m2)*(a+m2));
            d=1.0+aa*d;
            if (fabs(d) < FPMIN) d=FPMIN;
            c=1.0+aa/c;
            if (fabs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            h *= d*c;
            aa = -(a+m1)*(qab+m1)*x/((a+m2)*(qap+m2));
            d=1.0+aa*d;
            if (fabs(d) < FPMIN) d=FPMIN;
            c=1.0+aa/c;
            if (fabs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            del=d*c;
            h *= del;
            if (fabs(del-1.0) < EPS) break;
        }
        
        if (m1 > MAXIT) { m->mothurOut("[ERROR]: a or b too big or MAXIT too small in betacf."); m->mothurOutEndLine(); m->control_pressed = true; }
        return h;
        
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "betacf");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//[3][4] * [4][5] - columns in first must match rows in second, returns matrix[3][5]
vector<vector<double> > LinearAlgebra::matrix_mult(vector<vector<double> > first, vector<vector<double> > second){
	try {
		vector<vector<double> > product;
		
		int first_rows = first.size();
		int first_cols = first[0].size();
		int second_cols = second[0].size();
		
		product.resize(first_rows);
		for(int i=0;i<first_rows;i++){
			product[i].resize(second_cols);
		}
		
		for(int i=0;i<first_rows;i++){
			for(int j=0;j<second_cols;j++){
				
				if (m->control_pressed) { return product; }
					
				product[i][j] = 0.0;
				for(int k=0;k<first_cols;k++){
					product[i][j] += first[i][k] * second[k][j];
				}
			}
		}
		
		return product;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "matrix_mult");
		exit(1);
	}
	
}
/*********************************************************************************************************************************/

vector<vector<double> > LinearAlgebra::transpose(vector<vector<double> >matrix){
	try {
		vector<vector<double> > trans; trans.resize(matrix[0].size());
        for (int i = 0; i < trans.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) { trans[i].push_back(matrix[j][i]); }
        }
 				
		return trans;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "transpose");
		exit(1);
	}
	
}
/*********************************************************************************************************************************/

void LinearAlgebra::recenter(double offset, vector<vector<double> > D, vector<vector<double> >& G){
	try {
		int rank = D.size();
		
		vector<vector<double> > A(rank);
		vector<vector<double> > C(rank);
		for(int i=0;i<rank;i++){
			A[i].resize(rank);
			C[i].resize(rank);
		}
		
		double scale = -1.0000 / (double) rank;
		
		for(int i=0;i<rank;i++){
			A[i][i] = 0.0000;
			C[i][i] = 1.0000 + scale;
			for(int j=i+1;j<rank;j++){
				A[i][j] = A[j][i] = -0.5 * D[i][j] * D[i][j] + offset;
				C[i][j] = C[j][i] = scale;
			}
		}
		
		A = matrix_mult(C,A);
		G = matrix_mult(A,C);
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "recenter");
		exit(1);
	}
	
}
/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

int LinearAlgebra::tred2(vector<vector<double> >& a, vector<double>& d, vector<double>& e){
	try {
		double scale, hh, h, g, f;
		
		int n = a.size();
		
		d.resize(n);
		e.resize(n);
		
		for(int i=n-1;i>0;i--){
			int l=i-1;
			h = scale = 0.0000;
			if(l>0){
				for(int k=0;k<l+1;k++){
					scale += fabs(a[i][k]);
				}
				if(scale == 0.0){
					e[i] = a[i][l];
				}
				else{
					for(int k=0;k<l+1;k++){
						a[i][k] /= scale;
						h += a[i][k] * a[i][k];
					}
					f = a[i][l];
					g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
					e[i] = scale * g;
					h -= f * g;
					a[i][l] = f - g;
					f = 0.0;
					for(int j=0;j<l+1;j++){
						a[j][i] = a[i][j] / h;
						g = 0.0;
						for(int k=0;k<j+1;k++){
							g += a[j][k] * a[i][k];
						}
						for(int k=j+1;k<l+1;k++){
							g += a[k][j] * a[i][k];
						}
						e[j] = g / h;
						f += e[j] * a[i][j];
					}
					hh = f / (h + h);
					for(int j=0;j<l+1;j++){
						f = a[i][j];
						e[j] = g = e[j] - hh * f;
						for(int k=0;k<j+1;k++){
							a[j][k] -= (f * e[k] + g * a[i][k]);
						}
					}
				}
			}
			else{
				e[i] = a[i][l];
			}
			
			d[i] = h;
		}
		
		d[0] = 0.0000;
		e[0] = 0.0000;
		
		for(int i=0;i<n;i++){
			int l = i;
			if(d[i] != 0.0){
				for(int j=0;j<l;j++){
					g = 0.0000;
					for(int k=0;k<l;k++){
						g += a[i][k] * a[k][j];
					}
					for(int k=0;k<l;k++){
						a[k][j] -= g * a[k][i];
					}
				}
			}
			d[i] = a[i][i];
			a[i][i] = 1.0000;
			for(int j=0;j<l;j++){
				a[j][i] = a[i][j] = 0.0;
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "tred2");
		exit(1);
	}
	
}
/*********************************************************************************************************************************/

double LinearAlgebra::pythag(double a, double b)	{	return(pow(a*a+b*b,0.5));	}

/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

int LinearAlgebra::qtli(vector<double>& d, vector<double>& e, vector<vector<double> >& z) {
	try {
		int myM, i, iter;
		double s, r, p, g, f, dd, c, b;
		
		int n = d.size();
		for(int i=1;i<=n;i++){
			e[i-1] = e[i];
		}
		e[n-1] = 0.0000;
		
		for(int l=0;l<n;l++){
			iter = 0;
			do {
				for(myM=l;myM<n-1;myM++){
					dd = fabs(d[myM]) + fabs(d[myM+1]);
					if(fabs(e[myM])+dd == dd) break;
				}
				if(myM != l){
					if(iter++ == 3000) cerr << "Too many iterations in tqli\n";
					g = (d[l+1]-d[l]) / (2.0 * e[l]);
					r = pythag(g, 1.0);
					g = d[myM] - d[l] + e[l] / (g + SIGN(r,g));
					s = c = 1.0;
					p = 0.0000;
					for(i=myM-1;i>=l;i--){
						f = s * e[i];
						b = c * e[i];
						e[i+1] = (r=pythag(f,g));
						if(r==0.0){
							d[i+1] -= p;
							e[myM] = 0.0000;
							break;
						}
						s = f / r;
						c = g / r;
						g = d[i+1] - p;
						r = (d[i] - g) * s + 2.0 * c * b;
						d[i+1] = g + ( p = s * r);
						g = c * r - b;
						for(int k=0;k<n;k++){
							f = z[k][i+1];
							z[k][i+1] = s * z[k][i] + c * f;
							z[k][i] = c * z[k][i] - s * f;
						}
					}
					if(r == 0.00 && i >= l) continue;
					d[l] -= p;
					e[l] = g;
					e[myM] = 0.0;
				}
			} while (myM != l);
		}
		
		int k;
		for(int i=0;i<n;i++){
			p=d[k=i];
			for(int j=i;j<n;j++){
				if(d[j] >= p){
					p=d[k=j];
				}
			}
			if(k!=i){
				d[k]=d[i];
				d[i]=p;
				for(int j=0;j<n;j++){
					p=z[j][i];
					z[j][i] = z[j][k];
					z[j][k] = p;
				}
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "qtli");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//groups by dimension
vector< vector<double> > LinearAlgebra::calculateEuclidianDistance(vector< vector<double> >& axes, int dimensions){
	try {
		//make square matrix
		vector< vector<double> > dists; dists.resize(axes.size());
		for (int i = 0; i < dists.size(); i++) {  dists[i].resize(axes.size(), 0.0); }
		
		if (dimensions == 1) { //one dimension calc = abs(x-y)
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					dists[i][j] = abs(axes[i][0] - axes[j][0]);
					dists[j][i] = dists[i][j];
				}
			}
			
		}else if (dimensions > 1) { //two dimension calc = sqrt ((x1 - y1)^2 + (x2 - y2)^2)...
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					double sum = 0.0;
					for (int k = 0; k < dimensions; k++) {
						sum += ((axes[i][k] - axes[j][k]) * (axes[i][k] - axes[j][k]));
					}
					
					dists[i][j] = sqrt(sum);
					dists[j][i] = dists[i][j];
				}
			}
			
		}
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calculateEuclidianDistance");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//returns groups by dimensions from dimensions by groups
vector< vector<double> > LinearAlgebra::calculateEuclidianDistance(vector< vector<double> >& axes){
	try {
		//make square matrix
		vector< vector<double> > dists; dists.resize(axes[0].size());
		for (int i = 0; i < dists.size(); i++) {  dists[i].resize(axes[0].size(), 0.0); }
		
		if (axes.size() == 1) { //one dimension calc = abs(x-y)
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					dists[i][j] = abs(axes[0][i] - axes[0][j]);
					dists[j][i] = dists[i][j];
				}
			}
			
		}else if (axes.size() > 1) { //two dimension calc = sqrt ((x1 - y1)^2 + (x2 - y2)^2)...
			
			for (int i = 0; i < dists[0].size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					double sum = 0.0;
					for (int k = 0; k < axes.size(); k++) {
						sum += ((axes[k][i] - axes[k][j]) * (axes[k][i] - axes[k][j]));
					}
					
					dists[i][j] = sqrt(sum);
					dists[j][i] = dists[i][j];
				}
			}
			
		}
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calculateEuclidianDistance");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//assumes both matrices are square and the same size
double LinearAlgebra::calcPearson(vector< vector<double> >& euclidDists, vector< vector<double> >& userDists){
	try {
		
		//find average for - X
		int count = 0;
		float averageEuclid = 0.0; 
		for (int i = 0; i < euclidDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				averageEuclid += euclidDists[i][j];  
				count++;
			}
		}
		averageEuclid = averageEuclid / (float) count;   
			
		//find average for - Y
		count = 0;
		float averageUser = 0.0; 
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				averageUser += userDists[i][j]; 
				count++;
			}
		}
		averageUser = averageUser / (float) count;  

		double numerator = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < euclidDists.size(); i++) {
			
			for (int k = 0; k < i; k++) { //just lt dists
				
				float Yi = userDists[i][k];
				float Xi = euclidDists[i][k];
				
				numerator += ((Xi - averageEuclid) * (Yi - averageUser));
				denomTerm1 += ((Xi - averageEuclid) * (Xi - averageEuclid));
				denomTerm2 += ((Yi - averageUser) * (Yi - averageUser));
			}
		}
		
		double denom = (sqrt(denomTerm1) * sqrt(denomTerm2));
		double r = numerator / denom;
		
		//divide by zero error
		if (isnan(r) || isinf(r)) { r = 0.0; }
		
		return r;
		
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcPearson");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//assumes both matrices are square and the same size
double LinearAlgebra::calcSpearman(vector< vector<double> >& euclidDists, vector< vector<double> >& userDists){
	try {
		double r; 
		
		//format data
		map<float, int> tableX; 
		map<float, int>::iterator itTable;
		vector<spearmanRank> scores; 
		
		for (int i = 0; i < euclidDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				spearmanRank member(toString(scores.size()), euclidDists[i][j]);
				scores.push_back(member);  
				
				//count number of repeats
				itTable = tableX.find(euclidDists[i][j]);
				if (itTable == tableX.end()) { 
					tableX[euclidDists[i][j]] = 1;
				}else {
					tableX[euclidDists[i][j]]++;
				}
			}
		}
		
		//sort scores
		sort(scores.begin(), scores.end(), compareSpearman); 

		//calc LX
		double Lx = 0.0; 
		for (itTable = tableX.begin(); itTable != tableX.end(); itTable++) {
			double tx = (double) itTable->second;
			Lx += ((pow(tx, 3.0) - tx) / 12.0);
		}
		
		//find ranks of xi
		map<string, float> rankEuclid;
		vector<spearmanRank> ties;
		int rankTotal = 0;
		for (int j = 0; j < scores.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(scores[j]);
			
			if (j != (scores.size()-1)) { // you are not the last so you can look ahead
				if (scores[j].score != scores[j+1].score) { // you are done with ties, rank them and continue
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankEuclid[ties[k].name] = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					rankEuclid[ties[k].name] = thisrank;
				}
			}
		}
		
		
		//format data
		map<float, int> tableY; 
		scores.clear(); 
		
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				spearmanRank member(toString(scores.size()), userDists[i][j]);
				scores.push_back(member);  
				
				//count number of repeats
				itTable = tableY.find(userDists[i][j]);
				if (itTable == tableY.end()) { 
					tableY[userDists[i][j]] = 1;
				}else {
					tableY[userDists[i][j]]++;
				}
			}
		}
		
		//sort scores
		sort(scores.begin(), scores.end(), compareSpearman); 
		
		//calc LX
		double Ly = 0.0; 
		for (itTable = tableY.begin(); itTable != tableY.end(); itTable++) {
			double ty = (double) itTable->second;
			Ly += ((pow(ty, 3.0) - ty) / 12.0);
		}
		
		//find ranks of yi
		map<string, float> rankUser;
		ties.clear();
		rankTotal = 0;
		for (int j = 0; j < scores.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(scores[j]);
			
			if (j != (scores.size()-1)) { // you are not the last so you can look ahead
				if (scores[j].score != scores[j+1].score) { // you are done with ties, rank them and continue
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankUser[ties[k].name] = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					rankUser[ties[k].name] = thisrank;
				}
			}
		}
		
			
		double di = 0.0;	
		int count = 0;
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < i; j++) {
			
				float xi = rankEuclid[toString(count)];
				float yi = rankUser[toString(count)];
			
				di += ((xi - yi) * (xi - yi));
				
				count++;
			}
		}
		
		double n = (double) count;
		
		double SX2 = ((pow(n, 3.0) - n) / 12.0) - Lx;
		double SY2 = ((pow(n, 3.0) - n) / 12.0) - Ly;
		
		r = (SX2 + SY2 - di) / (2.0 * sqrt((SX2*SY2)));
		
		//divide by zero error
		if (isnan(r) || isinf(r)) { r = 0.0; }
	
		return r;
		
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcSpearman");
		exit(1);
	}
}

/*********************************************************************************************************************************/
//assumes both matrices are square and the same size
double LinearAlgebra::calcKendall(vector< vector<double> >& euclidDists, vector< vector<double> >& userDists){
	try {
		double r;
		
		//format data
		vector<spearmanRank> scores; 
		for (int i = 0; i < euclidDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				spearmanRank member(toString(scores.size()), euclidDists[i][j]);
				scores.push_back(member);
			}
		}
			
		//sort scores
		sort(scores.begin(), scores.end(), compareSpearman); 	
		
		//find ranks of xi
		map<string, float> rankEuclid;
		vector<spearmanRank> ties;
		int rankTotal = 0;
		for (int j = 0; j < scores.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(scores[j]);
			
			if (j != (scores.size()-1)) { // you are not the last so you can look ahead
				if (scores[j].score != scores[j+1].score) { // you are done with ties, rank them and continue
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankEuclid[ties[k].name] = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					rankEuclid[ties[k].name] = thisrank;
				}
			}
		}
		
		vector<spearmanRank> scoresUser; 
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < i; j++) {
				spearmanRank member(toString(scoresUser.size()), userDists[i][j]);
				scoresUser.push_back(member);  
			}
		}
		
		//sort scores
		sort(scoresUser.begin(), scoresUser.end(), compareSpearman); 	
		
		//find ranks of yi
		map<string, float> rankUser;
		ties.clear();
		rankTotal = 0;
		for (int j = 0; j < scoresUser.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(scoresUser[j]);
			
			if (j != (scoresUser.size()-1)) { // you are not the last so you can look ahead
				if (scoresUser[j].score != scoresUser[j+1].score) { // you are done with ties, rank them and continue
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankUser[ties[k].name] = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					rankUser[ties[k].name] = thisrank;
				}
			}
		}
		
		int numCoor = 0;
		int numDisCoor = 0;
		
		//order user ranks
		vector<spearmanRank> user; 
		for (int l = 0; l < scores.size(); l++) {   
			spearmanRank member(scores[l].name, rankUser[scores[l].name]);
			user.push_back(member);
		}
				
		int count = 0;
		for (int l = 0; l < scores.size(); l++) {
					
			int numWithHigherRank = 0;
			int numWithLowerRank = 0;
			float thisrank = user[l].score;
					
			for (int u = l+1; u < scores.size(); u++) {
				if (user[u].score > thisrank) { numWithHigherRank++; }
				else if (user[u].score < thisrank) { numWithLowerRank++; }
				count++;
			}
					
			numCoor += numWithHigherRank;
			numDisCoor += numWithLowerRank;
		}
				
		r = (numCoor - numDisCoor) / (float) count;
		
		//divide by zero error
		if (isnan(r) || isinf(r)) { r = 0.0; }
		
		return r;
		
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcKendall");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcKendall(vector<double>& x, vector<double>& y, double& sig){
	try {
		if (x.size() != y.size()) { m->mothurOut("[ERROR]: vector size mismatch."); m->mothurOutEndLine(); return 0.0; }
		
		//format data
		vector<spearmanRank> xscores; 
		for (int i = 0; i < x.size(); i++) {
			spearmanRank member(toString(i), x[i]);
			xscores.push_back(member);  
		}
		
		//sort xscores
		sort(xscores.begin(), xscores.end(), compareSpearman);
		
		//convert scores to ranks of x
		vector<spearmanRank*> ties;
		int rankTotal = 0;
		for (int j = 0; j < xscores.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(&(xscores[j]));
				
			if (j != xscores.size()-1) { // you are not the last so you can look ahead
				if (xscores[j].score != xscores[j+1].score) { // you are done with ties, rank them and continue
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						(*ties[k]).score = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					(*ties[k]).score = thisrank;
				}
			}
		}
		
			
		//format data
		vector<spearmanRank> yscores;
		for (int j = 0; j < y.size(); j++) {
			spearmanRank member(toString(j), y[j]);
			yscores.push_back(member);
		}
		
		//sort yscores
		sort(yscores.begin(), yscores.end(), compareSpearman);
		
		//convert to ranks
		map<string, float> rank;
		vector<spearmanRank> yties;
		rankTotal = 0;
		for (int j = 0; j < yscores.size(); j++) {
			rankTotal += (j+1);
			yties.push_back(yscores[j]);
				
			if (j != yscores.size()-1) { // you are not the last so you can look ahead
				if (yscores[j].score != yscores[j+1].score) { // you are done with ties, rank them and continue
					for (int k = 0; k < yties.size(); k++) {
						float thisrank = rankTotal / (float) yties.size();
						rank[yties[k].name] = thisrank;
					}
					yties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				for (int k = 0; k < yties.size(); k++) {
					float thisrank = rankTotal / (float) yties.size();
					rank[yties[k].name] = thisrank;
				}
			}
		}
			
			
		int numCoor = 0;
		int numDisCoor = 0;
		
		//associate x and y
		vector<spearmanRank> otus; 
		for (int l = 0; l < xscores.size(); l++) {   
			spearmanRank member(xscores[l].name, rank[xscores[l].name]);
			otus.push_back(member);
		}
				
		int count = 0;
		for (int l = 0; l < xscores.size(); l++) {
					
			int numWithHigherRank = 0;
			int numWithLowerRank = 0;
			float thisrank = otus[l].score;
					
			for (int u = l+1; u < xscores.size(); u++) {
				if (otus[u].score > thisrank) { numWithHigherRank++; }
				else if (otus[u].score < thisrank) { numWithLowerRank++; }
				count++;
			}
					
			numCoor += numWithHigherRank;
			numDisCoor += numWithLowerRank;
		}
				
		double p = (numCoor - numDisCoor) / (float) count;
		
		sig = calcKendallSig(x.size(), p);
		
		return p;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcKendall");
		exit(1);
	}
}
double LinearAlgebra::ran0(int& idum)
{
    const int IA=16807,IM=2147483647,IQ=127773;
    const int IR=2836,MASK=123459876;
    const double AM=1.0/double(IM);
    int k;
    double ans;
    
    idum ^= MASK;
    k=idum/IQ;
    idum=IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    ans=AM*idum;
    idum ^= MASK;
    return ans;
}

double LinearAlgebra::ran1(int &idum)
{
	const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	const int NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static int iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;
    
	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double LinearAlgebra::ran2(int &idum)
{
	const int IM1=2147483563,IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	static int idum2=123456789,iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;
    
	if (idum <= 0) {
		idum=(idum==0 ? 1 : -idum);
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double LinearAlgebra::ran3(int &idum)
{
	static int inext,inextp;
	static int iff=0;
	const int MBIG=1000000000,MSEED=161803398,MZ=0;
	const double FAC=(1.0/MBIG);
	static vector<int> ma(56);
	int i,ii,k,mj,mk;
    
	if (idum < 0 || iff == 0) {
		iff=1;
		mj=labs(MSEED-labs(idum));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < int(MZ)) mk += MBIG;
			mj=ma[ii];
		}
		for (k=0;k<4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < int(MZ)) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < int(MZ)) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

double LinearAlgebra::ran4(int &idum)
{
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static const unsigned long jflone = 0x00004080;
	static const unsigned long jflmsk = 0xffff007f;
#else
	static const unsigned long jflone = 0x3f800000;
	static const unsigned long jflmsk = 0x007fffff;
#endif
	unsigned long irword,itemp,lword;
	static int idums = 0;
    
	if (idum < 0) {
		idums = -idum;
		idum=1;
	}
	irword=idum;
	lword=idums;
	psdes(lword,irword);
	itemp=jflone | (jflmsk & irword);
	++idum;
	return (*(float *)&itemp)-1.0;
}

void LinearAlgebra::psdes(unsigned long &lword, unsigned long &irword)
{
	const int NITER=4;
	static const unsigned long c1[NITER]={
		0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static const unsigned long c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
	unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
    
	for (i=0;i<NITER;i++) {
		ia=(iswap=irword) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		irword=lword ^ (((ia = (ib >> 16) |
                          ((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		lword=iswap;
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcKendallSig(double n, double r){
    try {
        
        double sig = 0.0;
        double svar=(4.0*n+10.0)/(9.0*n*(n-1.0)); 
        double z= r/sqrt(svar); 
        sig=erfcc(fabs(z)/1.4142136);

		if (isnan(sig) || isinf(sig)) { sig = 0.0; }
        
        return sig;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcKendallSig");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcKruskalWallis(vector<spearmanRank>& values, double& pValue){
	try {
        double H;
        set<string> treatments;
        
        //rank values
        sort(values.begin(), values.end(), compareSpearman);
        vector<spearmanRank*> ties;
        int rankTotal = 0;
        vector<int> TIES;
        for (int j = 0; j < values.size(); j++) {
            treatments.insert(values[j].name);
            rankTotal += (j+1);
            ties.push_back(&(values[j]));
            
            if (j != values.size()-1) { // you are not the last so you can look ahead
                if (values[j].score != values[j+1].score) { // you are done with ties, rank them and continue
                    if (ties.size() > 1) { TIES.push_back(ties.size()); }
                    for (int k = 0; k < ties.size(); k++) {
                        double thisrank = rankTotal / (double) ties.size();
                        (*ties[k]).score = thisrank;
                    }
                    ties.clear();
                    rankTotal = 0;
                }
            }else { // you are the last one
                if (ties.size() > 1) { TIES.push_back(ties.size()); }
                for (int k = 0; k < ties.size(); k++) {
                    double thisrank = rankTotal / (double) ties.size();
                    (*ties[k]).score = thisrank;
                }
            }
        }
        
        
        // H = 12/(N*(N+1)) * (sum Ti^2/n) - 3(N+1)
        map<string, double> sums;
        map<string, double> counts;
        for (set<string>::iterator it = treatments.begin(); it != treatments.end(); it++) { sums[*it] = 0.0; counts[*it] = 0; }
        
        for (int j = 0; j < values.size(); j++) {
            sums[values[j].name] += values[j].score;
            counts[values[j].name]+= 1.0;
        }
        
        double middleTerm = 0.0;
        for (set<string>::iterator it = treatments.begin(); it != treatments.end(); it++) {
            middleTerm += ((sums[*it]*sums[*it])/counts[*it]);
        }
        
        double firstTerm = 12 / (double) (values.size()*(values.size()+1));
        double lastTerm = 3 * (values.size()+1);
        
        H = firstTerm * middleTerm - lastTerm;
       
        //adjust for ties
        if (TIES.size() != 0) {
            double sum = 0.0;
            for (int j = 0; j < TIES.size(); j++) { sum += ((TIES[j]*TIES[j]*TIES[j])-TIES[j]); }
            double result = 1.0 - (sum / (double) ((values.size()*values.size()*values.size())-values.size()));
            H /= result;
        }
        
        if (isnan(H) || isinf(H)) { H = 0; }
        
        //Numerical Recipes pg221
        pValue = 1.0 - (gammp(((treatments.size()-1)/(double)2.0), H/2.0));
        
        return H;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcKruskalWallis");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::normalvariate(double mean, double standardDeviation) {
    try {
        double u1 = ((double)(rand()) + 1.0 )/( (double)(RAND_MAX) + 1.0);
        double u2 = ((double)(rand()) + 1.0 )/( (double)(RAND_MAX) + 1.0);
        //double r = sqrt( -2.0*log(u1) );
        //double theta = 2.0*PI*u2;
        //cout << cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1)) << endl;
        return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "normalvariate");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//thanks http://www.johndcook.com/cpp_phi.html
double LinearAlgebra::pnorm(double x){
    try {
        // constants
        double a1 =  0.254829592;
        double a2 = -0.284496736;
        double a3 =  1.421413741;
        double a4 = -1.453152027;
        double a5 =  1.061405429;
        double p  =  0.3275911;
        
        // Save the sign of x
        int sign = 1;
        if (x < 0)
            sign = -1;
        x = fabs(x)/sqrt(2.0);
        
        // A&S formula 7.1.26
        double t = 1.0/(1.0 + p*x);
        double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
        
        return 0.5*(1.0 + sign*y);
        
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "pnorm");
		exit(1);
	}
}

/*********************************************************************************************************************************/
double LinearAlgebra::calcWilcoxon(vector<double>& x, vector<double>& y, double& sig){
	try {		
		double W = 0.0;
        sig = 0.0;
        
        vector<spearmanRank> ranks;
        for (int i = 0; i < x.size(); i++) {
            if (m->control_pressed) { return W; }
            spearmanRank member("x", x[i]);
            ranks.push_back(member);
        }
        
        for (int i = 0; i < y.size(); i++) {
            if (m->control_pressed) { return W; }
            spearmanRank member("y", y[i]);
            ranks.push_back(member);
        }
        
        //sort values
		sort(ranks.begin(), ranks.end(), compareSpearman);
		
		//convert scores to ranks of x
		vector<spearmanRank*> ties;
		int rankTotal = 0;
        vector<int> TIES;
		for (int j = 0; j < ranks.size(); j++) {
            if (m->control_pressed) { return W; }
			rankTotal += (j+1);
			ties.push_back(&(ranks[j]));
            
			if (j != ranks.size()-1) { // you are not the last so you can look ahead
				if (ranks[j].score != ranks[j+1].score) { // you are done with ties, rank them and continue
                    if (ties.size() > 1) { TIES.push_back(ties.size()); }
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						(*ties[k]).score = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
                if (ties.size() > 1) { TIES.push_back(ties.size()); }
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					(*ties[k]).score = thisrank;
				}
			}
		}
        
        //from R wilcox.test function
        //STATISTIC <- sum(r[seq_along(x)]) - n.x * (n.x + 1)/2
        double sumRanks = 0.0;
        for (int i = 0; i < ranks.size(); i++) {
            if (m->control_pressed) { return W; }
            if (ranks[i].name == "x") { sumRanks += ranks[i].score; }
        }
        
        W = sumRanks - x.size() * ((double)(x.size() + 1)) / 2.0;
        
        //exact <- (n.x < 50) && (n.y < 50)
        bool findExact = false;
        if ((x.size() < 50) && (y.size() < 50)) { findExact = true; }
        
        
        if (findExact && (TIES.size() == 0)) { //find exact and no ties
            //PVAL <- switch(alternative, two.sided = {
            //p <- if (STATISTIC > (n.x * n.y/2))
            PWilcox wilcox;
            double pval = 0.0;
            if (W > ((double)x.size()*y.size()/2.0)) {
                //pwilcox(STATISTIC-1, n.x, n.y, lower.tail = FALSE)
                pval = wilcox.pwilcox(W-1, x.size(), y.size(), false);
            }else {
                //pwilcox(STATISTIC,n.x, n.y)
                pval = wilcox.pwilcox(W, x.size(), y.size(), true);
            }
            sig = 2.0 * pval;
            if (1.0 < sig) { sig = 1.0; }
        }else {
            //z <- STATISTIC - n.x * n.y/2
            double z = W - (double)(x.size() * y.size()/2.0);
            //NTIES <- table(r)
            double sum = 0.0;
            for (int j = 0; j < TIES.size(); j++) { sum += ((TIES[j]*TIES[j]*TIES[j])-TIES[j]); }
           
            //SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) -
                                            //sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y -
                                                                            //1))))
            double sigma = 0.0;
            double firstTerm = (double)(x.size() * y.size()/12.0);
            double secondTerm = (double)(x.size() + y.size() + 1) - sum / (double)((x.size() + y.size()) * (x.size() + y.size() - 1));
            sigma = sqrt(firstTerm * secondTerm);
            
            //CORRECTION <- switch(alternative, two.sided = sign(z) * 0.5, greater = 0.5, less = -0.5)
            double CORRECTION = 0.0;
            if (z < 0) { CORRECTION = -1.0; }
            else if (z > 0) { CORRECTION = 1.0; }
            CORRECTION *= 0.5;
            
            z = (z - CORRECTION)/sigma;
            
            //PVAL <- switch(alternative,  two.sided = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
            sig = pnorm(z);
            if ((1.0-sig) < sig) { sig = 1.0 - sig; }
            sig *= 2;
        }
        
        return W;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcWilcoxon");
		exit(1);
	}
}

/*********************************************************************************************************************************/
double LinearAlgebra::choose(double n, double k){
	try {
        n = floor(n + 0.5);
        k = floor(k + 0.5);
        
        double lchoose = gammln(n + 1.0) - gammln(k + 1.0) - gammln(n - k + 1.0);
        
        return (floor(exp(lchoose) + 0.5));
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "choose");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcSpearman(vector<double>& x, vector<double>& y, double& sig){
	try {
		if (x.size() != y.size()) { m->mothurOut("[ERROR]: vector size mismatch."); m->mothurOutEndLine(); return 0.0; }
		
		//format data
        double sf = 0.0; //f^3 - f where f is the number of ties in x;
        double sg = 0.0; //f^3 - f where f is the number of ties in y;
		map<float, int> tableX; 
		map<float, int>::iterator itTable;
		vector<spearmanRank> xscores; 
		
		for (int i = 0; i < x.size(); i++) {
			spearmanRank member(toString(i), x[i]);
			xscores.push_back(member);  
				
			//count number of repeats
			itTable = tableX.find(x[i]);
			if (itTable == tableX.end()) { 
				tableX[x[i]] = 1;
			}else {
				tableX[x[i]]++;
			}
		}
		
		
		//calc LX
		double Lx = 0.0;
		for (itTable = tableX.begin(); itTable != tableX.end(); itTable++) {
			double tx = (double) itTable->second;
			Lx += ((pow(tx, 3.0) - tx) / 12.0);
		}
		
		
		//sort x
		sort(xscores.begin(), xscores.end(), compareSpearman);
		
		//convert scores to ranks of x
		//convert to ranks
		map<string, float> rankx;
		vector<spearmanRank> xties;
		int rankTotal = 0;
		for (int j = 0; j < xscores.size(); j++) {
			rankTotal += (j+1);
			xties.push_back(xscores[j]);
			
			if (j != xscores.size()-1) { // you are not the last so you can look ahead
				if (xscores[j].score != xscores[j+1].score) { // you are done with ties, rank them and continue
					for (int k = 0; k < xties.size(); k++) {
						float thisrank = rankTotal / (float) xties.size();
						rankx[xties[k].name] = thisrank;
					}
                    int t = xties.size();
                    sf += (t*t*t-t);
					xties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				for (int k = 0; k < xties.size(); k++) {
					float thisrank = rankTotal / (float) xties.size();
					rankx[xties[k].name] = thisrank;
				}
			}
		}		
			
		//format x
		vector<spearmanRank> yscores;
		map<float, int> tableY;
		for (int j = 0; j < y.size(); j++) {
			spearmanRank member(toString(j), y[j]);
			yscores.push_back(member);
				
			itTable = tableY.find(member.score);
			if (itTable == tableY.end()) { 
				tableY[member.score] = 1;
			}else {
				tableY[member.score]++;
			}
				
		}
			
		//calc Ly
		double Ly = 0.0;
		for (itTable = tableY.begin(); itTable != tableY.end(); itTable++) {
			double ty = (double) itTable->second;
			Ly += ((pow(ty, 3.0) - ty) / 12.0);
		}
			
		sort(yscores.begin(), yscores.end(), compareSpearman);
			
		//convert to ranks
		map<string, float> rank;
		vector<spearmanRank> yties;
		rankTotal = 0;
		for (int j = 0; j < yscores.size(); j++) {
			rankTotal += (j+1);
			yties.push_back(yscores[j]);
			
			if (j != yscores.size()-1) { // you are not the last so you can look ahead
				if (yscores[j].score != yscores[j+1].score) { // you are done with ties, rank them and continue
					for (int k = 0; k < yties.size(); k++) {
						float thisrank = rankTotal / (float) yties.size();
						rank[yties[k].name] = thisrank;
					}
                    int t = yties.size();
                    sg += (t*t*t-t);
					yties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				for (int k = 0; k < yties.size(); k++) {
					float thisrank = rankTotal / (float) yties.size();
					rank[yties[k].name] = thisrank;
				}
			}
		}
		
		double di = 0.0;
		for (int k = 0; k < x.size(); k++) {
					
			float xi = rankx[toString(k)];
			float yi = rank[toString(k)];
					
			di += ((xi - yi) * (xi - yi));
		}
				
		double p = 0.0;
				
		double n = (double) x.size();
		double SX2 = ((pow(n, 3.0) - n) / 12.0) - Lx;
		double SY2 = ((pow(n, 3.0) - n) / 12.0) - Ly;
				
		p = (SX2 + SY2 - di) / (2.0 * sqrt((SX2*SY2)));
		
		//Numerical Recipes 646
        sig = calcSpearmanSig(n, sf, sg, di);
		
		return p;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcSpearman");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcSpearmanSig(double n, double sf, double sg, double d){
    try {
        
        double sig = 0.0;
        double probrs = 0.0;
        double en=n;
        double en3n=en*en*en-en;
        double aved=en3n/6.0-(sf+sg)/12.0;
        double fac=(1.0-sf/en3n)*(1.0-sg/en3n);
        double vard=((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac;
        double zd=(d-aved)/sqrt(vard);
        double probd=erfcc(fabs(zd)/1.4142136);
        double rs=(1.0-(6.0/en3n)*(d+(sf+sg)/12.0))/sqrt(fac);
        fac=(rs+1.0)*(1.0-rs);
        if (fac > 0.0) {
            double t=rs*sqrt((en-2.0)/fac);
            double df=en-2.0;
            probrs=betai(0.5*df,0.5,df/(df+t*t));
        }else {
            probrs = 0.0;
        }
        
        //smaller of probd and probrs is sig
        sig = probrs;
        if (probd < probrs) { sig = probd; }
        
		if (isnan(sig) || isinf(sig)) { sig = 0.0; }
		
        return sig;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcSpearmanSig");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcPearson(vector<double>& x, vector<double>& y, double& sig){
	try {
		if (x.size() != y.size()) { m->mothurOut("[ERROR]: vector size mismatch."); m->mothurOutEndLine(); return 0.0; }
		
		//find average X
		float averageX = 0.0; 
		for (int i = 0; i < x.size(); i++) { averageX += x[i];  }
		averageX = averageX / (float) x.size(); 
		
		//find average Y
		float sumY = 0.0;
		for (int j = 0; j < y.size(); j++) { sumY += y[j]; }
		float Ybar = sumY / (float) y.size();
			
		double r = 0.0;
		double numerator = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
				
		for (int j = 0; j < x.size(); j++) {
			float Yi = y[j];
			float Xi = x[j];
					
			numerator += ((Xi - averageX) * (Yi - Ybar));
			denomTerm1 += ((Xi - averageX) * (Xi - averageX));
			denomTerm2 += ((Yi - Ybar) * (Yi - Ybar));
		}
				
		double denom = (sqrt(denomTerm1) * sqrt(denomTerm2));
				
		r = numerator / denom;
		
		//Numerical Recipes pg.644
        sig = calcPearsonSig(x.size(), r);
		
		return r;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcPearson");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double LinearAlgebra::calcPearsonSig(double n, double r){
    try {
        
        double sig = 0.0;
        const double TINY = 1.0e-20;
        double z = 0.5*log((1.0+r+TINY)/(1.0-r+TINY)); //Fisher's z transformation
    
        //code below was giving an error in betacf with sop files
        //int df = n-2;
        //double t = r*sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)));
        //sig = betai(0.5+df, 0.5, df/(df+t*t));
        
        //Numerical Recipes says code below gives approximately the same result
        sig = erfcc(fabs(z*sqrt(n-1.0))/1.4142136);
		if (isnan(sig) || isinf(sig)) { sig = 0.0; }
		
        return sig;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calcPearsonSig");
		exit(1);
	}
}
/*********************************************************************************************************************************/

vector<vector<double> > LinearAlgebra::getObservedEuclideanDistance(vector<vector<double> >& relAbundData){
    try {

        int numSamples = relAbundData.size();
        int numOTUs = relAbundData[0].size();
        
        vector<vector<double> > dMatrix(numSamples);
        for(int i=0;i<numSamples;i++){
            dMatrix[i].resize(numSamples);
        }
        
        for(int i=0;i<numSamples;i++){
            for(int j=0;j<numSamples;j++){
                
                if (m->control_pressed) { return dMatrix; }
                
                double d = 0;
                for(int k=0;k<numOTUs;k++){
                    d += pow((relAbundData[i][k] - relAbundData[j][k]), 2.0000);
                }
                dMatrix[i][j] = pow(d, 0.50000);
                dMatrix[j][i] = dMatrix[i][j];
                
            }
        }
        return dMatrix;
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "getObservedEuclideanDistance");
		exit(1);
	}
}

/*********************************************************************************************************************************/
vector<double> LinearAlgebra::solveEquations(vector<vector<double> > A, vector<double> b){
    try {
        int length = (int)b.size();
        vector<double> x(length, 0);
        vector<int> index(length);
        for(int i=0;i<length;i++){  index[i] = i;   }
        double d;
        
        ludcmp(A, index, d);  if (m->control_pressed) { return b; }
        lubksb(A, index, b);
        
        return b;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "solveEquations");
		exit(1);
	}
}
/*********************************************************************************************************************************/
vector<float> LinearAlgebra::solveEquations(vector<vector<float> > A, vector<float> b){
    try {
        int length = (int)b.size();
        vector<double> x(length, 0);
        vector<int> index(length);
        for(int i=0;i<length;i++){  index[i] = i;   }
        float d;
        
        ludcmp(A, index, d);  if (m->control_pressed) { return b; }
        lubksb(A, index, b);
        
        return b;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "solveEquations");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void LinearAlgebra::ludcmp(vector<vector<double> >& A, vector<int>& index, double& d){
    try {
        double tiny = 1e-20;
        
        int n = (int)A.size();
        vector<double> vv(n, 0.0);
        double temp;
        int imax;
        
        d = 1.0;
        
        for(int i=0;i<n;i++){
            double big = 0.0;
            for(int j=0;j<n;j++){   if((temp=fabs(A[i][j])) > big ) big=temp;  }
            if(big==0.0){   m->mothurOut("Singular matrix in routine ludcmp\n");    }
            vv[i] = 1.0/big;
        }
        
        for(int j=0;j<n;j++){
            if (m->control_pressed) { break; }
            for(int i=0;i<j;i++){
                double sum = A[i][j];
                for(int k=0;k<i;k++){   sum -= A[i][k] * A[k][j];   }
                A[i][j] = sum;
            }
            
            double big = 0.0;
            for(int i=j;i<n;i++){
                double sum = A[i][j];
                for(int k=0;k<j;k++){   sum -= A[i][k] * A[k][j];   }
                A[i][j] = sum;
                double dum;
                if((dum = vv[i] * fabs(sum)) >= big){
                    big = dum;
                    imax = i;
                }
            }
            if(j != imax){
                for(int k=0;k<n;k++){
                    double dum = A[imax][k];
                    A[imax][k] = A[j][k];
                    A[j][k] = dum;
                }
                d = -d;
                vv[imax] = vv[j];
            }
            index[j] = imax;
            
            if(A[j][j] == 0.0){ A[j][j] = tiny; }
            
            if(j != n-1){
                double dum = 1.0/A[j][j];
                for(int i=j+1;i<n;i++){ A[i][j] *= dum; }
            }
        }
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "ludcmp");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void LinearAlgebra::lubksb(vector<vector<double> >& A, vector<int>& index, vector<double>& b){
    try {
        //if(m->debug){   m->mothurOut("lubksb\n");    }
        double total;
        int n = (int)A.size();
        int ii = 0;
        
        for(int i=0;i<n;i++){
            //if(m->debug){   m->mothurOut("i loop " + toString(i) + "\n");    }
            if (m->control_pressed) { break; }
            int ip = index[i];
            total = b[ip];
            b[ip] = b[i];
            
            if (ii != 0) {
                for(int j=ii-1;j<i;j++){
                    //if(m->debug){   m->mothurOut("j loop " + toString(j) + "\n");    }
                    total -= A[i][j] * b[j];
                }
            }
            else if(total != 0){  ii = i+1;   }
            b[i] = total;
        }
        for(int i=n-1;i>=0;i--){
            //if(m->debug){   m->mothurOut("i loop " + toString(i) + "\n");    }
            total = b[i];
            for(int j=i+1;j<n;j++){ total -= A[i][j] * b[j];   }
            b[i] = total / A[i][i];
            //if (A[i][i] == 0) { cout << "ohno!!" << endl; }
        }
        //if(m->debug){   m->mothurOut("end lubksb\n");    }
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "lubksb");
		exit(1);
	}
}
/*********************************************************************************************************************************/

void LinearAlgebra::ludcmp(vector<vector<float> >& A, vector<int>& index, float& d){
    try {
        //if(m->debug){   m->mothurOut("ludcmp\n");    }
        double tiny = 1e-20;
        
        int n = (int)A.size();
        vector<float> vv(n, 0.0);
        double temp;
        int imax;
        
        d = 1.0;
        
        for(int i=0;i<n;i++){
            //if(m->debug){   m->mothurOut("i loop " + toString(i) + "\n");    }
            float big = 0.0;
            for(int j=0;j<n;j++){   if((temp=fabs(A[i][j])) > big ) big=temp;  }
            if(big==0.0){   m->mothurOut("Singular matrix in routine ludcmp\n");    }
            vv[i] = 1.0/big;
        }
        
        for(int j=0;j<n;j++){
            if (m->control_pressed) { break; }
            //if(m->debug){   m->mothurOut("j loop " + toString(j) + "\n");    }
            for(int i=0;i<j;i++){
                //if(m->debug){   m->mothurOut("i loop " + toString(i) + "\n");    }
                float sum = A[i][j];
                for(int k=0;k<i;k++){   sum -= A[i][k] * A[k][j];   }
                A[i][j] = sum;
            }
            
            float big = 0.0;
            for(int i=j;i<n;i++){
                //if(m->debug){   m->mothurOut("j loop " + toString(j) + "\n");    }
                float sum = A[i][j];
                for(int k=0;k<j;k++){   sum -= A[i][k] * A[k][j];   }
                A[i][j] = sum;
                float dum;
                if((dum = vv[i] * fabs(sum)) >= big){
                    big = dum;
                    imax = i;
                }
            }
            if(j != imax){
                for(int k=0;k<n;k++){
                    float dum = A[imax][k];
                    A[imax][k] = A[j][k];
                    A[j][k] = dum;
                }
                d = -d;
                vv[imax] = vv[j];
            }
            index[j] = imax;
            
            if(A[j][j] == 0.0){ A[j][j] = tiny; }
            
            if(j != n-1){
                float dum = 1.0/A[j][j];
                for(int i=j+1;i<n;i++){ A[i][j] *= dum; }
            }
        }
        //if(m->debug){   m->mothurOut("end ludcmp\n");    }
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "ludcmp");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void LinearAlgebra::lubksb(vector<vector<float> >& A, vector<int>& index, vector<float>& b){
    try {
        float total;
        int n = (int)A.size();
        int ii = 0;
        
        for(int i=0;i<n;i++){
            if (m->control_pressed) { break; }
            int ip = index[i];
            total = b[ip];
            b[ip] = b[i];
            
            if (ii != 0) {
                for(int j=ii-1;j<i;j++){
                    total -= A[i][j] * b[j];
                }
            }
            else if(total != 0){  ii = i+1;   }
            b[i] = total;
        }
        for(int i=n-1;i>=0;i--){
            total = b[i];
            for(int j=i+1;j<n;j++){ total -= A[i][j] * b[j];  }
            b[i] = total / A[i][i];
        }
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "lubksb");
		exit(1);
	}
}

/*********************************************************************************************************************************/

vector<vector<double> > LinearAlgebra::getInverse(vector<vector<double> > matrix){
    try {
        int n = (int)matrix.size();
        
        vector<vector<double> > inverse(n);
        for(int i=0;i<n;i++){   inverse[i].assign(n, 0.0000);   }
        
        vector<double> column(n, 0.0000);
        vector<int> index(n, 0);
        double dummy;
        
        ludcmp(matrix, index, dummy);
        
        for(int j=0;j<n;j++){
            if (m->control_pressed) { break; }
            
            column.assign(n, 0);
            
            column[j] = 1.0000;
            
            lubksb(matrix, index, column);
            
            for(int i=0;i<n;i++){   inverse[i][j] = column[i];  }
        }
        
        return inverse;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "getInverse");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//modelled R lda function - MASS:::lda.default
vector< vector<double> > LinearAlgebra::lda(vector< vector<double> >& a, vector<string> groups, vector< vector<double> >& means, bool& ignore) {
    try {
        
        set<string> uniqueGroups;
        for (int i = 0; i < groups.size(); i++) { uniqueGroups.insert(groups[i]); }
        int numGroups = uniqueGroups.size();
        
        map<string, int> quickIndex; //className to index. hoping to save counts, proportions and means in vectors to save time. This map will allow us to know index 0 in counts refers to group1.
        int count = 0;
        for (set<string>::iterator it = uniqueGroups.begin(); it != uniqueGroups.end(); it++) { quickIndex[*it] = count; count++; }
        
        int numSampled = groups.size(); //number of sampled groups
        int numOtus = a.size(); //number of flagged bins
        
        //counts <- as.vector(table(g)) //number of samples from each class in random sampling
        vector<int> counts; counts.resize(numGroups, 0);
        for (int i = 0; i < groups.size(); i++) {
            counts[quickIndex[groups[i]]]++;
        }
        
        vector<double> proportions; proportions.resize(numGroups, 0.0);
        for (int i = 0; i < numGroups; i++) {  proportions[i] = counts[i] / (double) numSampled; }
        
        means.clear(); //means[0] -> means[0][0] average for [group0][OTU0].
        means.resize(numGroups); for (int i = 0; i < means.size(); i++) { means[i].resize(numOtus, 0.0); }
        for (int j = 0; j < numSampled; j++) { //total for each class for each OTU
            for (int i = 0; i < numOtus; i++) { means[quickIndex[groups[j]]][i] += a[i][j]; }
        }
        //average for each class for each OTU
        for (int j = 0; j < numGroups; j++) { for (int i = 0; i < numOtus; i++) { means[j][i] /= counts[j]; }  }
        
        //randCov <- x - group.means[g, ]
        vector< vector<double> > randCov; //randCov[0][0] -> (random sample value0 for OTU0 - average for samples group in OTU0). example OTU0, random sample 0.01 from class early. average of class early for OTU0 is 0.005. randCov[0][0] = (0.01-0.005)
        for (int i = 0; i < numOtus; i++) { //for each flagged OTU
            vector<double> tempRand;
            for (int j = 0; j < numSampled; j++) { tempRand.push_back(a[i][j] - means[quickIndex[groups[j]]][i]);  }
            randCov.push_back(tempRand);
        }
        
        //find variance and std for each OTU
        //f1 <- sqrt(diag(var(x - group.means[g, ])))
        vector<double> stdF1;
        vector<double> ave;
        for (int i = 0; i < numOtus; i++) {
            stdF1.push_back(0.0);
            ave.push_back(m->getAverage(randCov[i]));
        }
        
        for (int i = 0; i < numOtus; i++) {
            for (int j = 0; j < numSampled; j++) { stdF1[i] += ((randCov[i][j] - ave[i]) * (randCov[i][j] - ave[i]));  }
        }
        
        //fac <- 1/(n - ng)
        double fac = 1 / (double) (numSampled-numGroups);
        
        for (int i = 0; i < stdF1.size(); i++) {
            stdF1[i] /= (double) (numSampled-1);
            stdF1[i] = sqrt(stdF1[i]);
        }
        
        vector< vector<double> > scaling; //[numOTUS][numOTUS]
        for (int i = 0; i < numOtus; i++) {
            vector<double> temp;
            for (int j = 0; j < numOtus; j++) {
                if (i == j) { temp.push_back(1.0/stdF1[i]); }
                else { temp.push_back(0.0); }
                
            }
            scaling.push_back(temp);
        }
        /*
         cout << "scaling = " << endl;
         for (int i = 0; i < scaling.size(); i++) {
         for (int j = 0; j < scaling[i].size(); j++) { cout << scaling[i][j] << '\t'; }
         cout << endl;
         }*/
        
        //X <- sqrt(fac) * ((x - group.means[g, ]) %*% scaling)
        vector< vector<double> > X = randCov; //[numOTUS][numSampled]
        //((x - group.means[g, ]) %*% scaling)
        //matrix multiplication of randCov and scaling
        LinearAlgebra linear;
        X = linear.matrix_mult(scaling, randCov); //[numOTUS][numOTUS] * [numOTUS][numSampled] = [numOTUS][numSampled]
        fac = sqrt(fac);
        
        for (int i = 0; i < X.size(); i++) {
            for (int j = 0; j < X[i].size(); j++) { X[i][j] *= fac;  }
        }
        
        vector<double> d;
        vector< vector<double> > v;
        vector< vector<double> > Xcopy; //X = [numOTUS][numSampled]
        bool transpose = false; //svd requires rows < columns, so if they are not then I need to transpose and look for the results in v.
        if (X.size() < X[0].size()) { Xcopy = linear.transpose(X); transpose=true; }
        else                        { Xcopy = X;                    }
        linear.svd(Xcopy, d, v); //Xcopy gets the results we want for v below, because R's version is [numSampled][numOTUS]
        
        /*cout << "Xcopy = " << endl;
        for (int i = 0; i < Xcopy.size(); i++) {
            for (int j = 0; j < Xcopy[i].size(); j++) { cout << Xcopy[i][j] << '\t'; }
            cout << endl;
        }
        cout << "v = " << endl;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) { cout << v[i][j] << '\t'; }
            cout << endl;
        }
         */
        
        int rank = 0;
        set<int> goodColumns;
        //cout << "d = " << endl;
        for (int i = 0; i < d.size(); i++) {  if (d[i] > 0.0000000001) { rank++; goodColumns.insert(i); } } //cout << d[i] << endl;
        
        if (rank == 0) {
            ignore=true; //m->mothurOut("[ERROR]: rank = 0: variables are numerically const\n"); m->control_pressed = true;
            return scaling; }
        
        //scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank], , rank)
        //X.s$v[, 1L:rank] = columns in Xcopy that correspond to "good" d values
        //diag(1/X.s$d[1L:rank], , rank) = matrix size rank * rank where the diagonal is 1/"good" dvalues
        /*example:
         d
         [1] 3.721545e+00 3.034607e+00 2.296649e+00 7.986927e-16 6.922408e-16
         [6] 5.471102e-16
         
         $v
         [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
         [1,]  0.31122175  0.10944725  0.20183340 -0.30136820  0.60786235 -0.13537095
         [2,] -0.29563726 -0.20568893  0.11233366 -0.05073289  0.48234270  0.21965978
         ...
         
         [1] "X.s$v[, 1L:rank]"
         [,1]        [,2]        [,3]
         [1,]  0.31122175  0.10944725  0.20183340
         [2,] -0.29563726 -0.20568893  0.11233366
         ...
         [1] "1/X.s$d[1L:rank]"
         [1] 0.2687056 0.3295320 0.4354170
         
         [1] "diag(1/X.s$d[1L:rank], , rank)"
         [,1]     [,2]     [,3]
         [1,] 0.2687056 0.000000 0.000000
         [2,] 0.0000000 0.329532 0.000000
         [3,] 0.0000000 0.000000 0.435417
         */
        if (transpose) {
            Xcopy = linear.transpose(v);
            /*
            cout << "Xcopy = " << endl;
            for (int i = 0; i < Xcopy.size(); i++) {
                for (int j = 0; j < Xcopy[i].size(); j++) { cout << Xcopy[i][j] << '\t'; }
                cout << endl;
            }*/
        }
        v.clear(); //store "good" columns - X.s$v[, 1L:rank]
        v.resize(Xcopy.size()); //[numOTUS]["good" columns]
        for (set<int>::iterator it = goodColumns.begin(); it != goodColumns.end(); it++) {
            for (int i = 0; i < Xcopy.size(); i++) {
                v[i].push_back(Xcopy[i][*it]);
            }
        }
        
        vector< vector<double> > diagRanks; diagRanks.resize(rank);
        for (int i = 0; i < rank; i++) { diagRanks[i].resize(rank, 0.0); }
        count = 0;
        for (set<int>::iterator it = goodColumns.begin(); it != goodColumns.end(); it++) {  diagRanks[count][count] = 1.0 / d[*it]; count++; }
        
        scaling = linear.matrix_mult(linear.matrix_mult(scaling, v), diagRanks); //([numOTUS][numOTUS]*[numOTUS]["good" columns]) = [numOTUS]["good" columns] then ([numOTUS]["good" columns] * ["good" columns]["good" columns] = scaling = [numOTUS]["good" columns]
        
        /*cout << "scaling = " << endl;
        for (int i = 0; i < scaling.size(); i++) {
            for (int j = 0; j < scaling[i].size(); j++) { cout << scaling[i][j] << '\t'; }
            cout << endl;
        }*/
        
        //Note: linear.matrix_mult [1][numGroups] * [numGroups][numOTUs] - columns in first must match rows in second, returns matrix[1][numOTUs]
        vector< vector<double> > prior; prior.push_back(proportions);
        vector< vector<double> >  xbar = linear.matrix_mult(prior, means);
        vector<double> xBar = xbar[0]; //length numOTUs
        
        /*cout << "xbar" << endl;
        for (int j = 0; j < numOtus; j++) {  cout << xBar[j] <<'\t'; }cout <<  endl;*/
        //fac <- 1/(ng - 1)
        fac = 1 / (double) (numGroups-1);
        //scale(group.means, center = xbar, scale = FALSE) %*% scaling
        vector< vector<double> > scaledMeans = means; //[numGroups][numOTUs]
        for (int i = 0; i < numGroups; i++) {
            for (int j = 0; j < numOtus; j++) {  scaledMeans[i][j] -= xBar[j]; }
        }
        scaledMeans = linear.matrix_mult(scaledMeans, scaling); //[numGroups][numOTUS]*[numOTUS]["good"columns] = [numGroups]["good"columns]
        
        
        //sqrt((n * prior) * fac)
        vector<double> temp = proportions; //[numGroups]
        for (int i = 0; i < temp.size(); i++) { temp[i] *= numSampled * fac; temp[i] = sqrt(temp[i]);  }
        
        //X <- sqrt((n * prior) * fac) * (scale(group.means, center = xbar, scale = FALSE) %*% scaling)
        //X <- temp * scaledMeans
        X.clear(); X = scaledMeans; //[numGroups]["good"columns]
        for (int i = 0; i < X.size(); i++) {
            for (int j = 0; j < X[i].size(); j++) {  X[i][j] *= temp[j];  }
        }
        /*
        cout << "X = " << endl;
        for (int i = 0; i < X.size(); i++) {
            for (int j = 0; j < X[i].size(); j++) { cout << X[i][j] << '\t'; }
            cout << endl;
        }
        */
        
        d.clear(); v.clear();
        //we want to transpose so results are in Xcopy, but if that makes rows > columns then we don't since svd requires rows < cols.
        transpose=false;
        if (X.size() > X[0].size()) {   Xcopy = X;  transpose=true;     }
        else                        {   Xcopy = linear.transpose(X);    }
        linear.svd(Xcopy, d, v); //Xcopy gets the results we want for v below
        /*cout << "Xcopy = " << endl;
        for (int i = 0; i < Xcopy.size(); i++) {
            for (int j = 0; j < Xcopy[i].size(); j++) { cout << Xcopy[i][j] << '\t'; }
            cout << endl;
        }
        
        cout << "v = " << endl;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) { cout << v[i][j] << '\t'; }
            cout << endl;
        }
        
        cout << "d = " << endl;
        for (int i = 0; i < d.size(); i++) { cout << d[i] << endl; }*/
        
        //rank <- sum(X.s$d > tol * X.s$d[1L])
        //X.s$d[1L] = larger value in d vector
        double largeD = m->max(d);
        rank = 0; goodColumns.clear();
        for (int i = 0; i < d.size(); i++) { if (d[i] > (0.0000000001*largeD)) { rank++; goodColumns.insert(i); } }
        
        if (rank == 0) {
            ignore=true;//m->mothurOut("[ERROR]: rank = 0: class means are numerically identical.\n"); m->control_pressed = true;
            return scaling; }
        
        if (transpose) { Xcopy = linear.transpose(v);  }
        //scaling <- scaling %*% X.s$v[, 1L:rank] - scaling * "good" columns
        v.clear(); //store "good" columns - X.s$v[, 1L:rank]
        v.resize(Xcopy.size()); //Xcopy = ["good"columns][numGroups]
        for (set<int>::iterator it = goodColumns.begin(); it != goodColumns.end(); it++) {
            for (int i = 0; i < Xcopy.size(); i++) {
                v[i].push_back(Xcopy[i][*it]);
            }
        }
        
        scaling = linear.matrix_mult(scaling, v); //[numOTUS]["good" columns] * ["good"columns][new "good" columns]
        
        /*cout << "scaling = " << endl;
        for (int i = 0; i < scaling.size(); i++) {
            for (int j = 0; j < scaling[i].size(); j++) { cout << scaling[i][j] << '\t'; }
            cout << endl;
        }*/
        ignore=false;
        return scaling;
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "lda");
		exit(1);
	}
}
/*********************************************************************************************************************************/
//Singular value decomposition (SVD) - adapted from http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp
/*
 * svdcomp - SVD decomposition routine.
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
 */

int LinearAlgebra::svd(vector< vector<double> >& a, vector<double>& w, vector< vector<double> >& v) {
    try {
        int flag, i, its, j, jj, k, l, nm;
        double c, f, h, s, x, y, z;
        double anorm = 0.0, g = 0.0, scale = 0.0;

        int numRows = a.size(); if (numRows == 0) { return 0; }
        int numCols = a[0].size();
        w.resize(numCols, 0.0);
        v.resize(numCols); for (int i = 0; i < numCols; i++) { v[i].resize(numRows, 0.0); }
    
        vector<double> rv1; rv1.resize(numCols, 0.0);
        if (numRows < numCols){  m->mothurOut("[ERROR]: numRows < numCols\n"); m->control_pressed = true; return 0; }

        /* Householder reduction to bidiagonal form */
        for (i = 0; i < numCols; i++)
        {
            /* left-hand reduction */
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < numRows)
            {
                for (k = i; k < numRows; k++)
                    scale += fabs((double)a[k][i]);
                if (scale)
                {
                    for (k = i; k < numRows; k++)
                    {
                        a[k][i] = (double)((double)a[k][i]/scale);
                        s += ((double)a[k][i] * (double)a[k][i]);
                    }
                    f = (double)a[i][i];
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    a[i][i] = (double)(f - g);
                    if (i != numCols - 1)
                    {
                        for (j = l; j < numCols; j++)
                        {
                            for (s = 0.0, k = i; k < numRows; k++)
                                s += ((double)a[k][i] * (double)a[k][j]);
                            f = s / h;
                            for (k = i; k < numRows; k++)
                                a[k][j] += (double)(f * (double)a[k][i]);
                        }
                    }
                    for (k = i; k < numRows; k++)
                        a[k][i] = (double)((double)a[k][i]*scale);
                }
            }
            w[i] = (double)(scale * g);
            
            /* right-hand reduction */
            g = s = scale = 0.0;
            if (i < numRows && i != numCols - 1)
            {
                for (k = l; k < numCols; k++)
                    scale += fabs((double)a[i][k]);
                if (scale)
                {
                    for (k = l; k < numCols; k++)
                    {
                        a[i][k] = (double)((double)a[i][k]/scale);
                        s += ((double)a[i][k] * (double)a[i][k]);
                    }
                    f = (double)a[i][l];
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    a[i][l] = (double)(f - g);
                    for (k = l; k < numCols; k++)
                        rv1[k] = (double)a[i][k] / h;
                    if (i != numRows - 1)
                    {
                        for (j = l; j < numRows; j++)
                        {
                            for (s = 0.0, k = l; k < numCols; k++)
                                s += ((double)a[j][k] * (double)a[i][k]);
                            for (k = l; k < numCols; k++)
                                a[j][k] += (double)(s * rv1[k]);
                        }
                    }
                    for (k = l; k < numCols; k++)
                        a[i][k] = (double)((double)a[i][k]*scale);
                }
            }
            anorm = max(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
        }
        
        /* accumulate the right-hand transformation */
        for (i = numCols - 1; i >= 0; i--)
        {
            if (i < numCols - 1)
            {
                if (g)
                {
                    for (j = l; j < numCols; j++)
                        v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                    for (j = l; j < numCols; j++)
                    {
                        for (s = 0.0, k = l; k < numCols; k++)
                            s += ((double)a[i][k] * (double)v[k][j]);
                        for (k = l; k < numCols; k++)
                            v[k][j] += (double)(s * (double)v[k][i]);
                    }
                }
                for (j = l; j < numCols; j++)
                    v[i][j] = v[j][i] = 0.0;
            }
            v[i][i] = 1.0;
            g = rv1[i];
            l = i;
        }
        
        /* accumulate the left-hand transformation */
        for (i = numCols - 1; i >= 0; i--)
        {
            l = i + 1;
            g = (double)w[i];
            if (i < numCols - 1)
                for (j = l; j < numCols; j++)
                    a[i][j] = 0.0;
            if (g)
            {
                g = 1.0 / g;
                if (i != numCols - 1)
                {
                    for (j = l; j < numCols; j++)
                    {
                        for (s = 0.0, k = l; k < numRows; k++)
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = (s / (double)a[i][i]) * g;
                        for (k = i; k < numRows; k++)
                            a[k][j] += (double)(f * (double)a[k][i]);
                    }
                }
                for (j = i; j < numRows; j++)
                    a[j][i] = (double)((double)a[j][i]*g);
            }
            else
            {
                for (j = i; j < numRows; j++)
                    a[j][i] = 0.0;
            }
            ++a[i][i];
        }
        
        /* diagonalize the bidiagonal form */
        for (k = numCols - 1; k >= 0; k--)
        {                             /* loop over singular values */
            for (its = 0; its < 30; its++)
            {                         /* loop over allowed iterations */
                flag = 1;
                for (l = k; l >= 0; l--)
                {                     /* test for splitting */
                    nm = l - 1;
                    if (fabs(rv1[l]) + anorm == anorm)
                    {
                        flag = 0;
                        break;
                    }
                    if (fabs((double)w[nm]) + anorm == anorm)
                        break;
                }
                if (flag)
                {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i <= k; i++)
                    {
                        f = s * rv1[i];
                        if (fabs(f) + anorm != anorm)
                        {
                            g = (double)w[i];
                            h = pythag(f, g);
                            w[i] = (double)h;
                            h = 1.0 / h;
                            c = g * h;
                            s = (- f * h);
                            for (j = 0; j < numRows; j++)
                            {
                                y = (double)a[j][nm];
                                z = (double)a[j][i];
                                a[j][nm] = (double)(y * c + z * s);
                                a[j][i] = (double)(z * c - y * s);
                            }
                        }
                    }
                }
                z = (double)w[k];
                if (l == k)
                {                  /* convergence */
                    if (z < 0.0)
                    {              /* make singular value nonnegative */
                        w[k] = (double)(-z);
                        for (j = 0; j < numCols; j++)
                            v[j][k] = (-v[j][k]);
                    }
                    break;
                }
                if (its >= 30) {
                    m->mothurOut("No convergence after 30,000! iterations \n"); m->control_pressed = true;
                    return(0);
                }
                
                /* shift from bottom 2 x 2 minor */
                x = (double)w[l];
                nm = k - 1;
                y = (double)w[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = pythag(f, 1.0);
                f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                
                /* next QR transformation */
                c = s = 1.0;
                for (j = l; j <= nm; j++)
                {
                    i = j + 1;
                    g = rv1[i];
                    y = (double)w[i];
                    h = s * g;
                    g = c * g;
                    z = pythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y = y * c;
                    for (jj = 0; jj < numCols; jj++)
                    {
                        x = (double)v[jj][j];
                        z = (double)v[jj][i];
                        v[jj][j] = (float)(x * c + z * s);
                        v[jj][i] = (float)(z * c - x * s);
                    }
                    z = pythag(f, h);
                    w[j] = (float)z;
                    if (z) 
                    {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = (c * g) + (s * y);
                    x = (c * y) - (s * g);
                    for (jj = 0; jj < numRows; jj++)
                    {
                        y = (double)a[jj][j];
                        z = (double)a[jj][i];
                        a[jj][j] = (double)(y * c + z * s);
                        a[jj][i] = (double)(z * c - y * s);
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k] = (double)x;
            }
        }
        
        return(0);
        
    }
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "svd");
		exit(1);
	}
}


