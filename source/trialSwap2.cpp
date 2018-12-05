#include "trialswap2.h"


//The sum_of_squares, havel_hakimi and calc_c_score algorithms have been adapted from I. Miklos and J. Podani. 2004. Randomization of presence-absence matrices: comments and new algorithms. Ecology 85:86-92.


double TrialSwap2::calc_c_score (vector<vector<int> > &co_matrix, vector<int> rowtotal, int ncols, int nrows)
{
    try {
        double cscore = 0.0;
        double maxD;
        double D;
        double normcscore = 0.0;
        int nonzeros = 0;
        //int ncols = co_matrix[0].size(); int nrows = rowtotal.size();
        vector<vector<double> > s; s.resize(nrows);
        for (int i = 0; i < nrows; i++) { s[i].resize(nrows,0.0); }//only fill half the matrix
        
        
        for(int i=0;i<nrows-1;i++)
        {
            
            for(int j=i+1;j<nrows;j++)
            {
                if (m->getControl_pressed()) { return 0; }
                for(int k=0;k<ncols;k++)
                {
                    if((co_matrix[i][k]==1)&&(co_matrix[j][k]==1)) //if both are 1s ie co-occurrence
                        s[i][j]++; //s counts co-occurrences
                }
                
                //rowtotal[i] = A, rowtotal[j] = B, ncols = P, s[i][j] = J
                cscore += (rowtotal[i]-s[i][j])*(rowtotal[j]-s[i][j]);///(nrows*(nrows-1)/2);
                D = (rowtotal[i]-s[i][j])*(rowtotal[j]-s[i][j]);
                
                if(ncols < (rowtotal[i] + rowtotal[j]))
                {
                    maxD = (ncols-rowtotal[i])*(ncols-rowtotal[j]);
                }
                else
                {
                    maxD = rowtotal[i] * rowtotal[j];
                }
                
                if(!util.isEqual(maxD, 0))
                {
                    normcscore += D/maxD;
                    nonzeros++;
                }
            }
        }
        
        
        cscore = normcscore/(double)nonzeros;

        return cscore;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_c_score");
        exit(1);
    }
}
/**************************************************************************************************/
int TrialSwap2::calc_checker (vector<vector<int> > &co_matrix, vector<int> rowtotal, int ncols, int nrows)
{
    try {
        int cunits=0;
        //int s[nrows][ncols];
        //int ncols = co_matrix[0].size(); int nrows = rowtotal.size();
        vector<vector<int> > s; s.resize(nrows);
        for (int i = 0; i < nrows; i++) { s[i].resize(nrows,0); }//only fill half the matrix
        
        for(int i=0;i<nrows-1;i++)
        {
            for(int j=i+1;j<nrows;j++)
            {
                if (m->getControl_pressed()) { return 0; }
                //s[i][j]=0;
                for(int k=0;k<ncols;k++)
                {
                    
                    //iterates through the row and counts co-occurrences. The total number of co-occurrences for each row pair is kept in matrix s at location s[i][j].
                    if((co_matrix[i][k]==1)&&(co_matrix[j][k]==1)) //if both are 1s ie co-occurrence
                        s[i][j]++; //s counts co-occurrences
                    
                }
               
                if (s[i][j] == 0) {  cunits+=1; }
            }
        }
        
        return cunits;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_checker");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::calc_vratio (int nrows, int ncols, vector<int> rowtotal, vector<int> columntotal)
{
    try {
        //int nrows = rowtotal.size();
        //int ncols = columntotal.size();
        int sumCol = accumulate(columntotal.begin(), columntotal.end(), 0 );
        // int sumRow = accumulate(rowtotal.begin(), rowtotal.end(), 0 );
        
        double colAvg = (double) sumCol / (double) ncols;
        // double rowAvg = (double) sumRow / (double) nrows;
        
        double p = 0.0;
        
        // double totalRowVar = 0.0;
        double rowVar = 0.0;
        double colVar = 0.0;
        
        for(int i=0;i<nrows;i++)
        {
            if (m->getControl_pressed()) { return 0; }
            p = (double) rowtotal[i]/(double) ncols;
            rowVar += p * (1.0-p);
        }
        
        for(int i=0;i<ncols;i++)
        {
            if (m->getControl_pressed()) { return 0; }
            colVar += pow(((double) columntotal[i]-colAvg),2);
        }
        
        colVar = (1.0/(double)ncols) * colVar;
        
        return colVar/rowVar;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_vratio");
        exit(1);
    }
    
}
/**************************************************************************************************/
int TrialSwap2::calc_combo (int nrows, int ncols, vector<vector<int> > &nullmatrix)
{
    try {
        //need to transpose so we can compare rows (row-major order)
        //int tmpnrows = nrows;
        vector<vector<int> > tmpmatrix;
        
        vector<int> tmprow;
        if(!tmpmatrix.empty())
            tmpmatrix.clear();
        for (int i=0;i<ncols;i++)
        {
            for (int j=0;j<nrows;j++)
            {
                tmprow.push_back(nullmatrix[j][i]);
            }
            
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }
        
        int unique = 0;
        int match = 0;
        for(int j=0;j<ncols;j++)
        {
            match = 0;
            for(int i=j+1;i<=ncols;i++)
            {
                //comparing matrix rows
                if( (tmpmatrix[j] == tmpmatrix[i]))
                {
                    match++;
                    break;
                }
            }
            
            //on the last iteration of a previously matched row it will add itself because it doesn't match any following rows, so that combination is counted
            if (match == 0)
                unique++;
        }
        return unique;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_combo");
        exit(1);
    }
}
/**************************************************************************************************/
int TrialSwap2::swap_checkerboards (vector<vector<int> > &co_matrix, int ncols, int nrows)
{
    try {
        Utils util;
        //do 100 runs to make sure enough swaps are happening. This does NOT mean that there will be 1000 swaps, but that is the theoretical max.
        for(int a=0;a<1000;a++){
            int i, j, k, l;
            i = util.getRandomIndex(nrows-1);
            while((j = util.getRandomIndex(nrows-1) ) == i ) {;if (m->getControl_pressed()) { return 0; }}
            k = util.getRandomIndex(ncols-1);
            while((l = util.getRandomIndex(ncols-1)) == k ) {;if (m->getControl_pressed()) { return 0; }}

            if((co_matrix[i][k]*co_matrix[j][l]==1 && co_matrix[i][l]+co_matrix[j][k]==0)||(co_matrix[i][k]+co_matrix[j][l]==0 && co_matrix[i][l]*co_matrix[j][k]==1)) //checking for checkerboard value and swap
            {
                co_matrix[i][k]=1-co_matrix[i][k];
                co_matrix[i][l]=1-co_matrix[i][l];
                co_matrix[j][k]=1-co_matrix[j][k];
                co_matrix[j][l]=1-co_matrix[j][l];

            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "swap_checkerboards");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::calc_pvalue_greaterthan (vector<double> scorevec, double initialscore)
{
    try {
        int runs = scorevec.size();
        double p = 0.0;
        for( int i=0;i<runs;i++)
        {
            if (m->getControl_pressed()) { return 0; }
            if(scorevec[i]>=initialscore)
                p++;
        }
        return p/(double)runs;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_pvalue_greaterthan");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::calc_pvalue_lessthan (vector<double> scorevec, double initialscore)
{
    try {
        int runs = scorevec.size();
        double p = 0.0;
        for( int i=0;i<runs;i++)
        {
            if (m->getControl_pressed()) { return 0; }
            if(scorevec[i]<=initialscore)
                p++;
        }
        return p/(double)runs;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "calc_pvalue_lessthan");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::t_test (double initialscore, int runs, double nullMean, vector<double> scorevec)
{
    try {
        double t;
        double sampleSD;
        double sum = 0;
        
        for(int i=0;i<runs;i++)
        {
            if (m->getControl_pressed()) { return 0; }
            sum += pow((scorevec[i] - nullMean),2);
           
        }
        
        m->mothurOut("nullMean: " + toString(nullMean)); m->mothurOutEndLine();
        
        m->mothurOut("sum: " + toString(sum)); m->mothurOutEndLine();
        
        sampleSD = sqrt( (1/runs) * sum );
        
        m->mothurOut("samplSD: " + toString(sampleSD)); m->mothurOutEndLine();
        
        t = (nullMean - initialscore) / (sampleSD / sqrt(runs));
        
        return t;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "t_test");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::getSD (int runs, vector<double> scorevec, double nullMean)
{
    try{
        double sum = 0;
        for(int i=0;i<runs;i++)
            {
                if (m->getControl_pressed()) { return 0; }
                sum += pow((scorevec[i] - nullMean),2);
            }
        return sqrt( (1/double(runs)) * sum );
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "getSD");
        exit(1);
    }
}
/**************************************************************************************************/
double TrialSwap2::get_zscore (double sd, double nullMean, double initscore)
{
    try {
        return (initscore - nullMean) / sd;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "get_zscore");
        exit(1);
    }
}
/**************************************************************************************************/
int TrialSwap2::print_matrix(vector<vector<int> > &matrix, int nrows, int ncols)
{
    try {
        m->mothurOut("matrix:"); m->mothurOutEndLine();
        
        for (int i = 0; i < nrows; i++)
        {
            if (m->getControl_pressed()) { return 0; }
            for (int j = 0; j < ncols; j++)
            {
                m->mothurOut(toString(matrix[i][j]));
            }
            m->mothurOutEndLine();
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "print_matrix");
        exit(1);
    }
}
/**************************************************************************************************/





