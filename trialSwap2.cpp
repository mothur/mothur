#include "trialswap2.h"


//The sum_of_squares, havel_hakimi and calc_c_score algorithms have been adapted from I. Miklos and J. Podani. 2004. Randomization of presence-absence matrices: comments and new algorithms. Ecology 85:86-92.


/**************************************************************************************************/
int TrialSwap2::intrand(int n){
    try {
        double z;
        
        z = (double)random() * (double)n / (double)RAND_MAX;
        if(z>=n)
            z=n-1;
        if(z<0)
            z=0;
        return((int)floor(z));
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "intrand");
		exit(1);
	}
}
/**************************************************************************************************/
/* completely random matrix, all column and row totals are variable, matrix size is the same
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim1(vector<vector<int> > &co_matrix){ 
    try {
        vector<int> randRow;
        vector<vector<int> > tmpmatrix;
        int nrows = co_matrix.size();
        int ncols = co_matrix[0].size();
        
        //clear co_matrix
        //     for(i=0;i<nrows;i++)
        //     {
        //         co_matrix.clear();
        //     }
        
        //cout << "building matrix" << endl;
        for(int i=0;i<nrows;i++){
            if (m->control_pressed) { break; }
            
            for(int j=0;j<ncols;j++){
                double randNum = rand() / double(RAND_MAX);
                //cout << randNum << endl;
                
                if(randNum > 0.5) {
                    randRow.push_back(1);
                }else{
                    randRow.push_back(0);
                }
            }
            tmpmatrix.push_back(randRow);
            randRow.clear();
            //cout << endl;
        }
        co_matrix = tmpmatrix;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim1");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 *row sums fixed, columns equiprobable 
 */
void TrialSwap2::sim2(vector<vector<int> > &co_matrix)
{ 
    try {
        
        for(int i=0;i<co_matrix.size();i++)
        {
            if (m->control_pressed) { break; }
            random_shuffle( co_matrix[i].begin(), co_matrix[i].end() ); 
        }
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim2");
		exit(1);
	}
}
/**************************************************************************************************/
int TrialSwap2::sim2plus(vector<int> rowtotal, vector<vector<int> > &co_matrix)
{
    try {
        int nrows = co_matrix.size();
        int ncols = co_matrix[0].size();
        double cellprob = 1.0/ncols;
        vector<double> cellprobvec;
        vector<int> tmprow;
        vector<vector<int> > tmpmatrix;
        //double randNum;
        
        double start = 0.0;
        
        for(int i=0; i<ncols; i++)
        {
            if (m->control_pressed) { return 0; }
            cellprobvec.push_back(start + cellprob);
            start = cellprobvec[i];
        }
        
        for(int i=0; i<nrows; i++)
        {
            tmprow.assign(ncols, 0);
            
            while( accumulate( tmprow.begin(), tmprow.end(), 0 ) < rowtotal[i])
            {
                if (m->control_pressed) { return 0; }
                double randNum = rand() / double(RAND_MAX);
                //cout << randNum << endl;
                if(randNum <= cellprobvec[0])
                {
                    tmprow[0] = 1;
                    continue;
                }
                for(int j=1;j<ncols;j++)
                {
                    //cout << range[j] << endl;
                    if(randNum <= cellprobvec[j] && randNum > cellprobvec[j-1] && tmprow[j] != 1)
                    {
                        tmprow[j] = 1;
                    }
                }
            }
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }
        co_matrix = tmpmatrix;
        tmpmatrix.clear();
        cellprobvec.clear();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim2plus");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 * same as sim2 but using initmatrix which is the initial co-occurrence matrix before transposition
 * may have to be changed depending on what matrix 'seed' is used. One way to use is to transpose
 * every null matrix before using an index and use the random matrix as a seed for the next null.
 */
/**************************************************************************************************/
void TrialSwap2::sim3(vector<vector<int> > &initmatrix)
{
    try {
        for(int i=0;i<initmatrix.size();i++)
        {
            if (m->control_pressed) { break; }
            random_shuffle( initmatrix[i].begin(), initmatrix[i].end() ); 
        }
        
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim3");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 *
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim4(vector<int> columntotal, vector<int> rowtotal, vector<vector<int> > &co_matrix)
{   
    try {
        vector<double> colProb;
        vector<int> tmprow;//(ncols, 7);
        vector<vector<int> > tmpmatrix;
        vector<double> range;
        vector<double> randNums;
        int ncols = columntotal.size();
        int nrows = rowtotal.size();
        tmprow.clear();
        
        double colSum = accumulate( columntotal.begin(), columntotal.end(), 0 );
        //cout << "col sum: " << colSum << endl;
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            colProb.push_back(columntotal[i]/colSum);
        }
        
        double start = 0.0;
        
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            range.push_back(start + colProb[i]);
            start = range[i];
        }
        
        for(int i=0;i<nrows;i++)
        {
            tmprow.assign(ncols, 0);
            if (m->control_pressed) { return 0; }
            
            while ( accumulate( tmprow.begin(), tmprow.end(), 0 ) < rowtotal[i])
            {
                if (m->control_pressed) { return 0; }
                
                double randNum = rand() / double(RAND_MAX);
                if(randNum <= range[0])
                {
                    tmprow[0] = 1;
                    continue;
                }
                for(int j=1;j<ncols;j++)
                {
                    if(randNum <= range[j] && randNum > range[j-1] && tmprow[j] != 1)
                    {
                        tmprow[j] = 1;
                    }
                    
                }
            }
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }
        
        co_matrix = tmpmatrix;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim4");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 * inverse of sim4, MUST BE TRANSPOSED BEFORE CO-OCCURRENCE ANALYSIS
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim5(vector<int> initcolumntotal,vector<int> initrowtotal, vector<vector<int> > &initmatrix)
{
    try {
        vector<double> colProb;
        vector<int> tmprow;//(ncols, 7);
        vector<vector<int> > tmpmatrix;
        vector<double> range;
        vector<double> randNums;
        int ncols = initcolumntotal.size();
        int nrows = initrowtotal.size();
        
        tmprow.clear();
        
        double colSum = accumulate( initcolumntotal.begin(), initcolumntotal.end(), 0 );
        //cout << "col sum: " << colSum << endl;
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            colProb.push_back(initcolumntotal[i]/colSum);
        }
        
        double start = 0.0;
        
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            range.push_back(start + colProb[i]);
            start = range[i];
        }
        
        for(int i=0;i<nrows;i++)
        {
            tmprow.assign(ncols, 0);
            if (m->control_pressed) { return 0; }
            
            while ( accumulate( tmprow.begin(), tmprow.end(), 0 ) < initrowtotal[i])
            {
                if (m->control_pressed) { return 0; }
                
                double randNum = rand() / double(RAND_MAX);
                if(randNum <= range[0])
                {
                    tmprow[0] = 1;
                    continue;
                }
                for(int j=1;j<ncols;j++)
                {
                    if(randNum <= range[j] && randNum > range[j-1] && tmprow[j] != 1)
                    {
                        tmprow[j] = 1;
                    }
                    
                }
            }
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }
        
        initmatrix = tmpmatrix;
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim5");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 *
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim6(vector<int> columntotal, vector<vector<int> > &co_matrix)
{
    try {
        vector<vector<int> > tmpmatrix;
        vector<double> colProb;
        vector<int> tmprow;
        vector<double> range;
        int ncols = columntotal.size();
        int nrows = co_matrix.size();
        
        int colSum = accumulate( columntotal.begin(), columntotal.end(), 0 );
        
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            colProb.push_back(columntotal[i]/double (colSum));
        }
        
        double start = 0.0;
        
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
            range.push_back(start + colProb[i]);
            start = range[i];
        }
        
        for(int i=0;i<nrows;i++)
        {
            if (m->control_pressed) { return 0; }
            tmprow.assign(ncols, 0);
            int tmprowtotal;
            tmprowtotal = (rand() / double (RAND_MAX)) * 10;
            while ( tmprowtotal > ncols) {
                if (m->control_pressed) { return 0; }
                tmprowtotal = (rand() / double (RAND_MAX)) * 10;
            }
            //cout << tmprowtotal << endl;
            //cout << accumulate( tmprow.begin(), tmprow.end(), 0 ) << endl;
            
            while ( accumulate( tmprow.begin(), tmprow.end(), 0 ) < tmprowtotal)
            {
                if (m->control_pressed) { return 0; }
                double randNum = rand() / double(RAND_MAX);
                //cout << randNum << endl;
                if(randNum <= range[0])
                {
                    tmprow[0] = 1;
                    continue;
                }
                for(int j=1;j<ncols;j++)
                {
                    //cout << range[j] << endl;
                    if(randNum <= range[j] && randNum > range[j-1] && tmprow[j] != 1)
                    {
                        tmprow[j] = 1;
                    }
                    
                }
                
                
            }
            
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }
        
        co_matrix = tmpmatrix;
        tmpmatrix.clear();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim6");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 * MUST BE TRANSPOSED BEFORE CO-OCCURRENCE ANALYSIS
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim7(vector<int> initrowtotal, vector<vector<int> > &co_matrix)
{
    try {
        vector<vector<double> > probmatrix;
        vector<vector<int> > tmpmatrix;
        vector<double> colProb;
        vector<double> probrow;
        vector<int> tmprow;
        vector<double> range;
        double nc;
        int ncols = co_matrix[0].size(); int nrows = co_matrix.size(); 
        
        tmpmatrix.assign(nrows, vector<int>(ncols, 0.));
        
        int rowsum = accumulate( initrowtotal.begin(), initrowtotal.end(), 0 );
        
        nc = rowsum * ncols;
        //cout << nc << endl;
        
        //assign null matrix based on probabilities
        
        double start = 0.0; // don't reset start -- probs should be from 0-1 thoughout the entire matrix 
        
        for(int i=0;i<nrows;i++)
        {
            if (m->control_pressed) { return 0; }
            //cout << initrowtotal[i]/double(nc) << endl;
            double cellprob = initrowtotal[i]/double(nc);
            //cout << cellprob << endl;
            for(int j=0;j<ncols;j++)
            {
                
                probrow.push_back(start + cellprob);
                //cout << probrow[j] << endl;
                //cout << start << endl;
                start = start + cellprob;
            }
            probmatrix.push_back(probrow);
            probrow.clear();
        }
        
        
        //while(tmprowsum < rowsum)
        //for(int k=0;k<rowsum;k++)
        int k = 0;
        while(k < rowsum)
        {
            if (m->control_pressed) { return 0; }
        done:
            //cout << k << endl;
            //tmprowsum = accumulate( tmprowtotal.begin(), tmprowtotal.end(), 0 );
            double randNum = rand() / double(RAND_MAX);
            //cout << randNum << "+" << endl;
            //special case for the first entry
            if(randNum <= probmatrix[0][0] && tmpmatrix[0][0] != 1)
            {
                tmpmatrix[0][0] = 1;
                k++;
                //cout << k << endl;
                continue;
            }
            
            
            for(int i=0;i<nrows;i++)
            {
                if (m->control_pressed) { return 0; }
                for(int j=0;j<ncols;j++)
                {
                    //cout << probmatrix[i][j] << endl;
                    if(randNum <= probmatrix[i][j] && randNum > probmatrix[i][j-1] && tmpmatrix[i][j] != 1)
                    {
                        tmpmatrix[i][j] = 1;
                        k++;
                        //cout << k << endl;
                        goto done;
                    }
                    //else
                    //k = k-1;
                }
                
            }
            
        }
        
        co_matrix = tmpmatrix;
        return 0;
    //build probibility matrix
    /* for(int i=0;i<nrows;i++)
     {
     for(int j=0;j<ncols;j++)
     {
     probrow.push_back(rowtotal[i]/nc);
     }
     probmatrix.pushback(probrow);
     probrow.clear;
     }
     */
    
    /* int colSum = accumulate( initcolumntotal.begin(), initcolumntotal.end(), 0 );
        
        for(int i=0;i<ncols;i++)
        {
            colProb.push_back(initcolumntotal[i]/double (colSum));
        }
        
        double start = 0.0;
        
        for(int i=0;i<ncols;i++)
        {
            range.push_back(start + colProb[i]);
            start = range[i];
        }
        
        for(int i=0;i<nrows;i++)
        {
            tmprow.assign(ncols, 0);
            int tmprowtotal;
            tmprowtotal = (rand() / double (RAND_MAX)) * 10;
            while ( tmprowtotal > ncols)
                tmprowtotal = (rand() / double (RAND_MAX)) * 10;
            //cout << tmprowtotal << endl;
            //cout << accumulate( tmprow.begin(), tmprow.end(), 0 ) << endl;
            
            while ( accumulate( tmprow.begin(), tmprow.end(), 0 ) < tmprowtotal)
            {
                double randNum = rand() / double(RAND_MAX);
                //cout << randNum << endl;
                if(randNum <= range[0])
                {
                    tmprow[0] = 1;
                    continue;
                }
                for(int j=1;j<ncols;j++)
                {
                    //cout << range[j] << endl;
                    if(randNum <= range[j] && randNum > range[j-1] && tmprow[j] != 1)
                    {
                        tmprow[j] = 1;
                    }
                }
            }
            
            tmpmatrix.push_back(tmprow);
            tmprow.clear();
        }

        initmatrix = tmpmatrix;
     */
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim7");
		exit(1);
	}
}
/**************************************************************************************************/
/*
 *
 *
 *
 */
/**************************************************************************************************/
int TrialSwap2::sim8(vector<int> columntotal, vector<int> rowtotal, vector<vector<int> > &co_matrix)
{   
    try {
        double prob; 
        double start = 0.0;
        int ncols = columntotal.size(); int nrows = rowtotal.size(); 
        double probarray[nrows * ncols];
        double randnum;
        int grandtotal; 
        int total = 0;
        
        //double colSum = accumulate( columntotal.begin(), columntotal.end(), 0 );
        double rowSum = accumulate( rowtotal.begin(), rowtotal.end(), 0 );
        
        if (m->control_pressed) { return 0; }
        
        //cout << "rowsum: " << rowSum << endl;
        
        grandtotal = rowSum;
        
        //create probability matrix with each site being between 0 and 1
        for (int i=0;i<nrows;i++) {
            if (m->control_pressed) { return 0; }
            for (int j=0;j<ncols;j++) {
                prob = (rowtotal[i] * columntotal[j])/(rowSum*rowSum);
                if (prob == 0.0)
                    probarray[ncols * i + j] = -1;
                else
                    probarray[ncols * i + j] = start + prob;
                //probmatrixrow.pushback(start + prob);
                start += prob;
            }
        }
        //cout << "prbarray" << endl;
        //for(int i=0;i<(nrows*ncols);i++)
        //cout << probarray[i] << " ";
        //cout << endl;
        
        //generate random muber between 0 and 1 and interate through probarray until found
        while (total < grandtotal)  {
            if (m->control_pressed) { return 0; }
            randnum = rand() / double(RAND_MAX);
            //cout << "rand num: " << randnum << endl;
            if((randnum <= probarray[0]) && (probarray[0] != 2) ) {
                probarray[0] = 2;
                total++;
                continue;
            }
            for(int i=1;i<(nrows*ncols);i++) {
                if (m->control_pressed) { return 0; }
                if((randnum <= probarray[i]) && (randnum > probarray[i-1]) && (probarray[i] != 2) ) {
                    probarray[i] = 2;
                    total++;
                    break;
                }
                else
                    continue;
            }
        }
        //cout << "prbarray" << endl;
        //for(int i=0;i<(nrows*ncols);i++)
        //cout << probarray[i] << " ";
        //cout << endl;
        for(int i=0;i<nrows;i++) {
            if (m->control_pressed) { return 0; }
            for(int j=0;j<ncols;j++) {
                if(probarray[ncols * i + j] == 2)
                    co_matrix[i][j] = 1;
                else
                    co_matrix[i][j] = 0;
            }
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sim8");
		exit(1);
	}
}
/**************************************************************************************************/
int TrialSwap2::havel_hakimi(vector<int> rowtotal,vector<int> columntotal,vector<vector<int> > &co_matrix) 
{
    try {
        int nrows = co_matrix.size();
        int ncols = co_matrix[0].size();
        int i,j,k;
        vector<int> r1; r1.resize(nrows,0);
        vector<int> c;  c.resize(ncols,0);
        vector<int> c1; c1.resize(ncols,0);
       
        
        for(i=0;i<nrows;i++) {
            for(j=0;j<ncols;j++) {
                co_matrix[i][j]=0;
            }
        }
        for(i=0;i<nrows;i++) {
            r1[i]=1; 
        }
           
        for(i=0;i<ncols;i++)
        {
            c[i]=columntotal[i];
            c1[i]=i;
        }
        
        for(k=0;k<nrows;k++)
        {
            if (m->control_pressed) { return 0; }
            i=intrand(nrows);
            while(r1[i]==0) {
                if (m->control_pressed) { return 0; }
                i=intrand(nrows);
            }
            r1[i]=0;
            sho(c,c1,ncols);
            for(j=0;j<rowtotal[i];j++)
            {
                if (m->control_pressed) { return 0; }
                co_matrix[i][c1[j]]=1;
                c[j]--;
                if(c[j]<0)
                    m->mothurOut("Uhh! " + toString(c1[j]) + "\n");
            }
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "havel_hakimi");
		exit(1);
	}
}
/**************************************************************************************************/
int TrialSwap2::sho(vector<int> c, vector<int> c1, int k)
{
    try {
        int i,j,temp;
        
        for(j=k-1;j>0;j--)
        {
            if (m->control_pressed) { return 0; }
            for(i=0;i<j;i++)
            {
                if(c[i]<c[i+1])
                {
                    temp=c[i];
                    c[i]=c[i+1];
                    c[i+1]=temp;
                    temp=c1[i];
                    c1[i]=c1[i+1];
                    c1[i+1]=temp;
                }
            }
        }
        for(j=1;j<1000;j++)
        {
            if (m->control_pressed) { return 0; }
            i=intrand(k-1);
            if(c[i]==c[i+1])
            {
                temp=c[i];
                c[i]=c[i+1];
                c[i+1]=temp;
                temp=c1[i];
                c1[i]=c1[i+1];
                c1[i+1]=temp;
            }
        }
        return(0);
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "sho");
		exit(1);
	}
}
/**************************************************************************************************/
double TrialSwap2::calc_c_score (vector<vector<int> > &co_matrix,vector<int>  rowtotal)
{
    try {
        double cscore = 0.0;
        double maxD;
        double D;
        double normcscore = 0.0;
        int nonzeros = 0;
        int ncols = co_matrix[0].size(); int nrows = rowtotal.size(); 
        vector<vector<double> > s(nrows, vector<double>(nrows,0.0)); //only fill half the matrix
        
        for(int i=0;i<nrows-1;i++)
        {
            
            for(int j=i+1;j<nrows;j++)
            {
                if (m->control_pressed) { return 0; }
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
                
                if(maxD != 0)
                {
                    normcscore += D/maxD;
                    nonzeros++;    
                }            
            }
        }
        
        cscore = cscore/(double)(nrows*(nrows-1)/2);
        //cout << "normalized c score: " << normcscore/nonzeros << endl;
        
        return cscore;
    }
	catch(exception& e) {
		m->errorOut(e, "TrialSwap2", "calc_c_score");
		exit(1);
	}
}
/**************************************************************************************************/
int TrialSwap2::calc_checker (vector<vector<int> > &co_matrix, vector<int>  rowtotal)
{
    try {
        int cunits=0;
        //int s[nrows][ncols];
        int ncols = co_matrix[0].size(); int nrows = rowtotal.size(); 
        vector<vector<int> > s(nrows, vector<int>(nrows,0)); //only fill half the matrix
        
        for(int i=0;i<nrows-1;i++)
        {
            for(int j=i+1;j<nrows;j++)
            {
                if (m->control_pressed) { return 0; }
                //s[i][j]=0;
                for(int k=0;k<ncols;k++)
                {
                    //cout << s[i][j] << endl;
                    //iterates through the row and counts co-occurrences. The total number of co-occurrences for each row pair is kept in matrix s at location s[i][j].
                    if((co_matrix[i][k]==1)&&(co_matrix[j][k]==1)) //if both are 1s ie co-occurrence
                        s[i][j]++; //s counts co-occurrences
                    
                }
                //cout << "rowtotal: " << rowtotal[i] << endl;
                //cout << "co-occurrences: " << s[i][j] << endl;
                //cunits+=(rowtotal[i]-s[i][j])*(rowtotal[j]-s[i][j]);
                if (s[i][j] == 0)
                {
                    cunits+=1;
                }
                //cunits+=s[i][j];
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
double TrialSwap2::calc_vratio (vector<int> rowtotal, vector<int> columntotal)
{
    try {
        int nrows = rowtotal.size();
        int ncols = columntotal.size();
        int sumCol = accumulate(columntotal.begin(), columntotal.end(), 0 );
       // int sumRow = accumulate(rowtotal.begin(), rowtotal.end(), 0 );
        
        double colAvg = (double) sumCol / (double) ncols;
 //       double rowAvg = (double) sumRow / (double) nrows;
        
        double p = 0.0;
        
 //       double totalRowVar = 0.0;
        double rowVar = 0.0;
        double colVar = 0.0;
        
        for(int i=0;i<nrows;i++)
        {
            if (m->control_pressed) { return 0; }
            p = (double) rowtotal[i]/(double) ncols;
            rowVar += p * (1.0-p);
        } 
        
        for(int i=0;i<ncols;i++)
        {
            if (m->control_pressed) { return 0; }
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
int TrialSwap2::calc_combo (vector<vector<int> > &initmatrix)
{
    try {
        int initrows = initmatrix.size();
        int unique = 0;
        int match = 0;
        int matches = 0;
        for(int i=0;i<initrows;i++)
        {
            match = 0;
            for(int j=i+1;j<=initrows;j++)
            {
                if (m->control_pressed) { return 0; }
                if( (initmatrix[i] == initmatrix[j])) 
                {
                    match++;
                    matches++;
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
int TrialSwap2::swap_checkerboards (vector<vector<int> > &co_matrix)
{
    try {
        int ncols = co_matrix[0].size(); int nrows = co_matrix.size(); 
        int i, j, k, l;
        i=intrand(nrows);
        while((j = intrand(nrows) ) == i ) {;if (m->control_pressed) { return 0; }}
        k=intrand(ncols);
        while((l = intrand(ncols) ) == k ) {;if (m->control_pressed) { return 0; }}
        //cout << co_matrix[i][k] << " " << co_matrix[j][l] << endl;
        //cout << co_matrix[i][l] << " " << co_matrix[j][k] << endl;
        //cout << co_matrix[i][l] << " " << co_matrix[j][k] << endl;
        //cout << co_matrix[i][l] << " " << co_matrix[j][k] << endl;
        if((co_matrix[i][k]*co_matrix[j][l]==1 && co_matrix[i][l]+co_matrix[j][k]==0)||(co_matrix[i][k]+co_matrix[j][l]==0 && co_matrix[i][l]*co_matrix[j][k]==1)) //checking for checkerboard value and swap
        {
            co_matrix[i][k]=1-co_matrix[i][k];
            co_matrix[i][l]=1-co_matrix[i][l];
            co_matrix[j][k]=1-co_matrix[j][k];
            co_matrix[j][l]=1-co_matrix[j][l];
            //cout << "swapped!" << endl;
        }
        //cout << "i: " << i << " j: " << j << " k: " << " l: " << l << endl;
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
            if (m->control_pressed) { return 0; }
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
            if (m->control_pressed) { return 0; }
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
            if (m->control_pressed) { return 0; }
            sum += pow((scorevec[i] - nullMean),2);
            //cout << "scorevec[" << i << "]" << scorevec[i] << endl;
        }
        
        m->mothurOut("nullMean: " + toString(nullMean)); m->mothurOutEndLine();
        
        m->mothurOut("sum: " + toString(sum));  m->mothurOutEndLine();
        
        sampleSD = sqrt( (1/runs) * sum );
        
        m->mothurOut("samplSD: " + toString(sampleSD));  m->mothurOutEndLine();
        
        t = (nullMean - initialscore) / (sampleSD / sqrt(runs));
        
        return t;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "t_test");
        exit(1);
    }
}
/**************************************************************************************************/
int TrialSwap2::print_matrix(vector<vector<int> > &matrix, int nrows, int ncols)
{
    try {
         m->mothurOut("matrix:");  m->mothurOutEndLine();
        
        for (int i = 0; i < nrows; i++)
        {
            if (m->control_pressed) { return 0; }
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
int TrialSwap2::transpose_matrix (vector<vector<int> > &initmatrix, vector<vector<int> > &co_matrix)//, int nrows, int nocols)
{    
    try {
        int ncols = initmatrix.size(); int nrows = initmatrix[0].size(); 
        int tmpnrows = nrows;
        //vector<vector<int> > tmpvec;
        vector<int> tmprow;
        if(!co_matrix.empty())
            co_matrix.clear();
        for (int i=0;i<nrows;i++)
        {       
            if (m->control_pressed) { return 0; }
            for (int j=0;j<ncols;j++)
            {
                tmprow.push_back(initmatrix[j][i]);
            }
            /*if (accumulate( tmprow.begin(), tmprow.end(), 0 ) == 0)
             {
             tmpnrows--;
             }
             else */
            co_matrix.push_back(tmprow);
            tmprow.clear();
        }
        nrows = tmpnrows;
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "transpose_matrix");
        exit(1);
    }
}
/**************************************************************************************************/
int TrialSwap2::update_row_col_totals(vector<vector<int> > &co_matrix, vector<int> &rowtotal, vector<int> &columntotal)
{
    try {
        //rowtotal.clear();
        //columntotal.clear();
        //generate (rowtotal.begin(), rowtotal.end(), 0);
        //generate (columntotal.begin(), columntotal.end(), 0);
        int nrows = co_matrix.size();
        int ncols = co_matrix[0].size();
        vector<int> tmpcolumntotal(ncols, 0);
        vector<int> tmprowtotal(nrows, 0);
        
        int rowcount = 0;
        
        for (int i = 0; i < nrows; i++)
        {
            if (m->control_pressed) { return 0; }
            for (int j = 0; j < ncols; j++)
            {
                if (co_matrix[i][j] == 1)
                {
                    rowcount++;
                    tmpcolumntotal[j]++;
                }           
            }    
            tmprowtotal[i] = rowcount;
            rowcount = 0;
        }
        columntotal = tmpcolumntotal;
        rowtotal = tmprowtotal;
        /*cout << "rowtotal: ";
        for(int i = 0; i<nrows; i++) { cout << rowtotal[i]; }
        cout << "  ";
        cout << " coltotal: ";
        for(int i = 0; i<ncols; i++) { cout << columntotal[i]; }
        cout << endl;*/
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TrialSwap2", "update_row_col_totals");
        exit(1);
    }
}
/**************************************************************************************************/
/*int main(int argc, char *argv[])
{
    srand (time(0));
    char* input_filename = argv[1];
    std::ifstream infile (input_filename);
   
    //skip the first line of headers
    getline(infile, line);
    //get the first line of data
    getline(infile, line);
    
    nrows = 0;
    ncols = 0;
    
    //int numspaces = 0;
    char nextChar;
    
    for (int i=0; i<int(line.length()); i++)
    {
      nextChar = line.at(i); // gets a character
      if (isspace(line[i]))
          ncols++;
    }
    
    ncols = ncols-3;
    
    cout << "number of OTUs: ";
    cout << ncols << endl;
    
    infile.close();
    
    std::ifstream infile2 (input_filename);
    
    //skip first line of headers
    getline(infile2, line);
    
    while (!infile2.eof())
    { 
        getline(infile2, line);
        if (!line.empty())
            nrows++;
    }
    
    cout << "number of sites: ";
    cout << nrows << endl;
    
    infile2.close();
    
    std::ifstream infile3 (input_filename);
    
    //skip first line
    getline(infile3, line);
    
    //variables that depend on info from initial matrix
    vector<vector<int> > co_matrix;//[nrows][ncols];
    vector<vector<int> > initmatrix;
    vector<int> tmprow;
    vector<double> stats;
    int tmpnrows = nrows;
    
    for (int row1=0; row1<nrows; row1++) // first line was skipped when counting, so we can start from 0
    {
        //ignore first 3 cols in each row, data starts on the 3th col
        for (int i = 0; i < 3; i++)
            infile3 >> tmp;

        for (int col=0; col<ncols; col++) 
        {
            infile3 >> tmp;
            //cout << tmp << endl;
            if (atoi(tmp.c_str()) > 0)
                tmprow.push_back(1);
            else
                tmprow.push_back(0);        
        }
        if (accumulate( tmprow.begin(), tmprow.end(), 0 ) == 0)
        {
            tmpnrows--;
        }
        else
            initmatrix.push_back(tmprow);
        //add the row to the matrix
        //initmatrix.push_back(tmprow);
        tmprow.clear();
        //cout << tmprow << endl;
    }   
    
    infile3.close();
    nrows = tmpnrows;
    
    //print init matrix    
    /* cout << "original matrix:" << endl;

    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            cout << initmatrix[i][j];            
        }    
        cout << endl;
    } */
    
        //for (i=0;i<ncols;i++)
        //cout << "col "<< i<< ": " << columntotal[i] << endl;
    
    //co_matrix is now initmatrix and newmatrix is now co_matrix
    
    //remove cols where sum is 0
    
    //transpose matrix
   /* int newmatrows = ncols;
    int newmatcols = nrows;
    int initcols = ncols; //for the combo metric
    int initrows = nrows; //for the combo metric
    //swap for transposed matrix
    nrows = newmatrows;//ncols;
    ncols = newmatcols;//nrows;
    
    vector<int> columntotal(ncols, 0);
    vector<int> initcolumntotal(ncols, 0);
    vector<int> initrowtotal(nrows, 0);
    vector<int> rowtotal(nrows, 0);
    
    transpose_matrix(initmatrix,co_matrix);
    //remove degenerate rows and cols

    //cout << "transposed matrix:" << endl;
    int rowcount = 0;
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            if (co_matrix[i][j] == 1)
            {
                rowcount++;
                columntotal[j]++;
            }
            //cout << co_matrix[i][j];            
        }    
        //cout << " row total: " << rowcount << endl;
        //cout << endl;
        rowtotal[i] = rowcount;
        rowcount = 0;
    }
    
    initcolumntotal = rowtotal;
    initrowtotal = columntotal;
    
    cout << endl;    
    
    runs = atol(argv[2]);    
    int metric = atol(argv[3]);
    int nullModel = atol(argv[4]);
    double initscore;
    update_row_col_totals(co_matrix, rowtotal, columntotal, ncols, nrows);
    //do initial metric: checker, c score, v ratio or combo
    switch(metric) 
    {
        case 1:
            //c score
            initscore = calc_c_score(co_matrix, rowtotal);
            cout << "initial c score: " << initscore << endl;
            //print_matrix(co_matrix, nrows, ncols);
            break;
            
        case 2:
            //checker
            initscore = calc_checker(co_matrix, rowtotal);
            cout << "initial checker score: " << initscore << endl;
            break;
            
        case 3:
            //v ratio
            initscore = calc_vratio(nrows, ncols, rowtotal, columntotal);
            cout << "initial v ratio: " << initscore << endl;
            break;
            
        case 4:
            //combo
            initscore = calc_combo(initrows, initcols, initmatrix);
            cout << "initial combo score: " << initscore << endl;
            //set co_matrix equal to initmatrix because combo requires row comparisons
            co_matrix = initmatrix;
            break;
            
        case 5:
            //test!
            
            //print_matrix(co_matrix, nrows, ncols);
            //sim1(nrows, ncols, co_matrix);
            //sim2(nrows, ncols, co_matrix);
            //sim3(initrows, initcols, initmatrix);
            //sim4(columntotal, rowtotal, co_matrix);
            //sim5(initcolumntotal, initmatrix);
            //sim6(columntotal, co_matrix);
            //sim7(initcolumntotal, initmatrix);          
            sim8(columntotal, rowtotal, co_matrix);
            //print_matrix(initmatrix, initrows, initcols);
            //print_matrix(co_matrix, nrows, ncols);
            
            break;
            
        default:
            cout << "no metric selected!" << endl;
            return 1;
            
    }
      
    //matrix initialization
    //havel_hakimi(nrows, ncols, rowtotal, columntotal, co_matrix);
    //sum_of_square(nrows, ncols, rowtotal, columntotal, co_matrix);
    //co-matrix is now a random matrix with the same row and column totals as the initial matrix
    
    //null matrix burn in
    cout << "initializing null matrix...";
    for(int l=0;l<10000;l++)
    {
       //swap_checkerboards (co_matrix); 
       //if(l%10 == 0)        
        switch(nullModel)
        {
            case 1:
                //
                sim1(nrows, ncols, co_matrix);
                break;

            case 2:
                //sim2
                sim2(nrows, ncols, co_matrix);
                //sim2plus(nrows, ncols, initrowtotal, co_matrix);
                break;

            case 3:
                //sim3
                sim3(initrows, initcols, initmatrix);
                //transpose_matrix(initmatrix,co_matrix);
                co_matrix = initmatrix;
                break;

            case 4:
                //sim4
                sim4(columntotal, rowtotal, co_matrix);
                break;

            case 5:
                //sim5
                sim5(initcolumntotal, initrowtotal, initmatrix);
                transpose_matrix(initmatrix,co_matrix);
                //co_matrix = initmatrix;
                break;

            case 6:
                sim6(columntotal, co_matrix);
                break;

            case 7:
                //sim7(ncols, nrows, initrowtotal, co_matrix);          
                //transpose_matrix(initmatrix,co_matrix);
                //co_matrix = initmatrix;
                break;

            case 8:
                sim8(columntotal, rowtotal, co_matrix);
                break;

            case 9:
                //swap_checkerboards
                swap_checkerboards (co_matrix);
                break;

            default:
                cout << "no null model selected!" << endl;
                return 1;
        }
    }
    cout << "done!" << endl;
      
    //generate null matrices and calculate the metrics
    
    cout << "run: " << endl;
    for(int trial=0;trial<runs;trial++) //runs
    {
        printf("\b\b\b\b\b\b\b%7d",trial+1);
        fflush(stdout);
        
        switch(nullModel)
        {
            case 1: 
                //
                sim1(nrows, ncols, co_matrix);
                break;

            case 2:
                //sim2
                sim2(nrows, ncols, co_matrix);
                //for(int i=0;i<nrows;i++)
                    //cout << rowtotal[i] << " ";
                //sim2plus(nrows, ncols, initrowtotal, co_matrix);
                break;

            case 3:
                //sim3
                for(int i=0;i<nrows;i++)
                    cout  << " " << rowtotal[i];
                sim3(initrows, initcols, initmatrix);
                transpose_matrix(initmatrix,co_matrix);
                break;

            case 4:
                //sim4
                sim4(columntotal, rowtotal, co_matrix);
                break;

            case 5:
                //sim5
                sim5(initcolumntotal, initrowtotal, initmatrix);
                transpose_matrix(initmatrix,co_matrix);
                break;

            case 6:
                sim6(columntotal, co_matrix);
                break;

            case 7:
                sim7(ncols, nrows, initrowtotal, co_matrix);
                //print_matrix(co_matrix, nrows, ncols);
                //transpose_matrix(initmatrix,co_matrix);
                break;

            case 8:
                //sim8(initcolumntotal, initrowtotal, co_matrix);
                //initrow and initcol are flipped because of transposition. this is super ugly!
                sim8(initrowtotal, initcolumntotal, co_matrix);
                break;

            case 9:
                //swap_checkerboards
                swap_checkerboards (co_matrix);
                break;

            default:
                cout << "no null model selected!" << endl;
                return 1;
        }

        //cout << endl;
        //print_matrix(co_matrix, nrows, ncols);
        update_row_col_totals(co_matrix, rowtotal, columntotal, ncols, nrows);
        //cout << columntotal[1]<<endl;
        double tmp;
        switch(metric) 
        {
            case 1:
                //c score
                //swap_checkerboards(co_matrix);
                //cout <<  "%" << tmp << " nrows " << nrows << " ncols " << ncols << " rowtotals ";
                //for(int i = 0; i<nrows; i++) { cout << rowtotal[i]; }
                //cout << endl;
                tmp = calc_c_score(co_matrix, rowtotal);
                //cout << "%" << tmp << " ";
                stats.push_back(tmp);
                //cout << "%" << tmp << " ";
                //print_matrix(co_matrix, nrows, ncols);
                break;
                
            case 2:
                //checker
                //swap_checkerboards(co_matrix);
                //cout <<  "%" << tmp << " nrows " << nrows << " ncols " << ncols << " rowtotals ";
                //for(int i = 0; i<nrows; i++) { cout << rowtotal[i]; }
                //cout << endl;
                tmp = calc_checker(co_matrix, rowtotal);
                stats.push_back(tmp);
                //cout << "%" << tmp << endl;
                break;
                
            case 3:
                //v ratio
                stats.push_back(calc_vratio(nrows, ncols, rowtotal, columntotal));
                break;
                
            case 4:
                //combo
                stats.push_back(calc_combo(nrows, ncols, co_matrix));
                break;
                
            case 5:
                //test!
                break;
                
            default:
                cout << "no metric selected!" << endl;
                return 1;
                
        } 
        
    //print_matrix(co_matrix, nrows, ncols);
    //print_matrix(initmatrix, initrows, initcols);

    }
    
    cout << endl;
    
    double total = 0.0;
    for (int i=0; i<stats.size();i++)
    {
        total+=stats[i];
    }
        
    double nullMean = double (total/stats.size()); 

    cout << '\n' << "average metric score: " << nullMean << endl;

    //cout << "t-test: " << t_test(initscore, runs, nullMean, stats) << endl;
    
    if (metric == 1 || metric == 2)
        cout << "pvalue: " << calc_pvalue_greaterthan (stats, initscore) << endl;
    else
        cout << "pvalue: " << calc_pvalue_lessthan (stats, initscore) << endl;
         
    //print_matrix(co_matrix, nrows, ncols);
    
    return 0;
    
}*/
/**************************************************************************************************/


