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
					
			for (int u = l; u < scores.size(); u++) {
				if (user[u].score > thisrank) { numWithHigherRank++; }
				else if (user[u].score < thisrank) { numWithLowerRank++; }
				count++;
			}
					
			numCoor += numWithHigherRank;
			numDisCoor += numWithLowerRank;
		}
				
		//comparing to yourself
		count -= userDists.size();
				
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


