/*
 *  linearalgebra.cpp
 *  mothur
 *
 *  Created by westcott on 1/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "linearalgebra.h"

/*********************************************************************************************************************************/

inline double SIGN(const double a, const double b)
{
    return b>=0 ? (a>=0 ? a:-a) : (a>=0 ? -a:a);
}
/*********************************************************************************************************************************/

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
double LinearAlgebra::calcPearson(vector< vector<double> >& euclidDists, vector< vector<double> >& userDists){
	try {
		
	/*	euclidDists.clear();
		userDists.clear();
		
		euclidDists.resize(1);
		userDists.resize(1);
		
		userDists[0].push_back(0.3070833);
		userDists[0].push_back(0.3244475);
		userDists[0].push_back(0.6055993);
		userDists[0].push_back(0.3372481);
		userDists[0].push_back(0.9151715);
		userDists[0].push_back(0.6182255);
		userDists[0].push_back(0.7748142);
		userDists[0].push_back(0.08554735);
		userDists[0].push_back(0.6343481);
		userDists[0].push_back(0.4049274);
		
		euclidDists[0].push_back(0.3342815);
		euclidDists[0].push_back(0.3173829);
		euclidDists[0].push_back(0.6852404);
		euclidDists[0].push_back(0.7819186);
		euclidDists[0].push_back(0.5705242);
		euclidDists[0].push_back(0.8007263);
		euclidDists[0].push_back(0.8561724);
		euclidDists[0].push_back(0.4901089);
		euclidDists[0].push_back(0.7027247);
		euclidDists[0].push_back(0.7669696);*/
		
		
		//find average for - X
		int count = 0;
		vector<float> averageEuclid; averageEuclid.resize(euclidDists.size(), 0.0);
		for (int i = 0; i < euclidDists.size(); i++) {
			for (int j = 0; j < euclidDists[i].size(); j++) {
				averageEuclid[i] += euclidDists[i][j];  
				count++;
			}
		}
		for (int i = 0; i < averageEuclid.size(); i++) {  averageEuclid[i] = averageEuclid[i] / (float) count;   }
			
		//find average for - Y
		count = 0;
		vector<float> averageUser; averageUser.resize(userDists.size(), 0.0);
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < userDists[i].size(); j++) {
				averageUser[i] += userDists[i][j]; 
				count++;
			}
		}
		for (int i = 0; i < averageUser.size(); i++) {  averageUser[i] = averageUser[i] / (float) count;  }

		double numerator = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < euclidDists.size(); i++) {
			
			for (int k = 0; k < euclidDists[i].size(); k++) {
				
				float Yi = userDists[i][k];
				float Xi = euclidDists[i][k];
				
				numerator += ((Xi - averageEuclid[i]) * (Yi - averageUser[i]));
				denomTerm1 += ((Xi - averageEuclid[i]) * (Xi - averageEuclid[i]));
				denomTerm2 += ((Yi - averageUser[i]) * (Yi - averageUser[i]));
			}
		}
		
		double denom = (sqrt(denomTerm1) * sqrt(denomTerm2));
		double r = numerator / denom;

		return r;
		
	}
	catch(exception& e) {
		m->errorOut(e, "LinearAlgebra", "calculateEuclidianDistance");
		exit(1);
	}
}
/*********************************************************************************************************************************/


