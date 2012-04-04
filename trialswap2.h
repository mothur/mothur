#ifndef TRIALSWAP2
#define TRIALSWAP2

/*
 *  trialswap2.h
 *  Mothur
 *
 *  Created by Kathryn Iverson on June 27, 2011.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"


class TrialSwap2 {
    
public:
	TrialSwap2(){  m = MothurOut::getInstance(); };
    ~TrialSwap2(){};
    
    double calc_pvalue_lessthan (vector<double>, double);
    double calc_pvalue_greaterthan (vector<double>, double);
    int swap_checkerboards (vector<vector<int> > &);
    int calc_combo (vector<vector<int> > &);
    double calc_vratio (vector<int>, vector<int>);
    int calc_checker (vector<vector<int> > &,vector<int>);
    double calc_c_score (vector<vector<int> > &,vector<int>);
    
    int sim1 (vector<vector<int> > &);
    void sim2(vector<vector<int> >&);
    int sim2plus(vector<int>, vector<vector<int> > &);
    void sim3(vector<vector<int> > &);
    int sim4(vector<int>, vector<int>, vector<vector<int> > &);
    int sim5(vector<int>, vector<int>, vector<vector<int> > &);
    int sim6(vector<int>, vector<vector<int> > &);
    int sim7(vector<int>, vector<vector<int> > &);
    int sim8(vector<int>, vector<int>, vector<vector<int> > &);
    int transpose_matrix (vector<vector<int> > &, vector<vector<int> > &);
    int update_row_col_totals(vector<vector<int> > &, vector<int>&, vector<int>&);

    
private:
    MothurOut* m;
    
    double t_test (double, int, double, vector<double>);
    int print_matrix(vector<vector<int> > &, int, int);
    
    

};

#endif


