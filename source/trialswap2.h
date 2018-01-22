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
#include "utils.hpp"


class TrialSwap2 {
    
public:
    TrialSwap2(){ m = MothurOut::getInstance(); };
    ~TrialSwap2(){};
    
    double calc_pvalue_lessthan (vector<double>, double);
    double calc_pvalue_greaterthan (vector<double>, double);
    int swap_checkerboards (vector<vector<int> > &, int, int);
    int calc_combo (int, int, vector<vector<int> > &);
    double calc_vratio (int, int, vector<int>, vector<int>);
    int calc_checker (vector<vector<int> > &, vector<int>, int, int);
    double calc_c_score (vector<vector<int> > &, vector<int>, int, int);
    double get_zscore (double, double, double);
    double getSD (int, vector<double>, double);
    
    
private:
    MothurOut* m;
    Utils util;
    
    double t_test (double, int, double, vector<double>);
    int print_matrix(vector<vector<int> > &, int, int);
    
    
        
};
#endif


