#ifndef MOTHUR_METASTATS
#define MOTHUR_METASTATS

/*
 *  mothurmetastats.h
 *  Mothur
 *
 *  Created by westcott on 7/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "mothurout.h"

class MothurMetastats {
	
	public:
		MothurMetastats(double, int); //threshold, numPermutations
		~MothurMetastats();
	
		int runMetastats(string, vector< vector<double> >&, int); //outputFileName, data, secondGroupingStart
	
	private:
		MothurOut* m;
		int row, column, numPermutations;
		double threshold;
	
		int start(vector<double>&, int, vector<double>&, vector< vector<double> >&); //Find the initial values for the matrix
		int meanvar(vector<double>&, int, vector<double>&);
		int testp(vector<double>&, vector<double>&, vector<double>&, int, vector<double>&, vector<double>&);
		int permute_matrix(vector<double>&, vector<double>&, int, vector<double>&, vector<double>&, vector<double>&);
		int permute_array(vector<int>&);
		int calc_twosample_ts(vector<double>&, int, vector<double>&, vector<double>&, vector<double>&);



	
};
	
	//void testp(double *permuted_ttests,int *B,double *permuted,double 
		//	   *Imatrix,int *nc,int *nr,int *g,double *Tinitial,double *ps);
	//void permute_matrix(double *Imatrix,int *nc,int *nr,double 
					//	*permuted,int *g,double *trial_ts,double *Tinitial,double 
					//	*counter);
	//void permute_array(int *array, int n);
	//void calc_twosample_ts(double *Pmatrix,int *g,int *nc,int *nr,double 
	//					   *Ts,double *Tinitial,double *counter1);
	//void meanvar(double *pmatrix,int *g,int *nr,int *nc,double *storage);
	//void start(double *Imatrix,int *g,int *nr,int *nc,double *testing,
	//		   double storage[][9]);
	
	//int metastat_main (char*, int, int, double, int, double**, int);
	

#endif

