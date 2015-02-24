/*
 *  mothurmetastats.cpp
 *  Mothur
 *
 *  Created by westcott on 7/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothurmetastats.h"
#include "mothurfisher.h"
#include "spline.h"

/***********************************************************/
MothurMetastats::MothurMetastats(double t, int n) {
	try {
		m = MothurOut::getInstance(); 
		threshold = t;
		numPermutations = n;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "MothurMetastats");
		exit(1);
	}	
}
/***********************************************************/
MothurMetastats::~MothurMetastats() {}
/***********************************************************/
 //main metastats function
int MothurMetastats::runMetastats(string outputFileName, vector< vector<double> >& data, int secGroupingStart) {
    try {
        row = data.size();		 //numBins
		column = data[0].size(); //numGroups in subset
        secondGroupingStart = secGroupingStart; //g number of samples in group 1
         
        vector< vector<double> > Pmatrix; Pmatrix.resize(row);
        for (int i = 0; i < row; i++) { Pmatrix[i].resize(column, 0.0);  } // the relative proportion matrix
        vector< vector<double> > C1; C1.resize(row);
        for (int i = 0; i < row; i++) { C1[i].resize(3, 0.0);  } // statistic profiles for class1 and class 2
        vector< vector<double> > C2; C2.resize(row);            // mean[1], variance[2], standard error[3] 
        for (int i = 0; i < row; i++) { C2[i].resize(3, 0.0);  } 
        vector<double> T_statistics; T_statistics.resize(row, 1); // a place to store the true t-statistics 
        vector<double> pvalues; pvalues.resize(row, 1); // place to store pvalues
        vector<double> qvalues; qvalues.resize(row, 1); // stores qvalues
       
        //*************************************
        //      convert to proportions
        //      generate Pmatrix
        //*************************************
        vector<double> totals; totals.resize(column, 0); // sum of columns
        //total[i] = total abundance for group[i]
		for (int i = 0; i < column; i++) {
			for (int j = 0; j < row; j++) {
				totals[i] += data[j][i];
			}
        }
        
        for (int i = 0; i < column; i++) {
			for (int j = 0; j < row; j++) {
				Pmatrix[j][i] = data[j][i]/totals[i];
			}
        }
        
        //#********************************************************************************
        //# ************************** STATISTICAL TESTING ********************************
        //#********************************************************************************
        
        if (column == 2){  //# then we have a two sample comparison
            //#************************************************************
            //#  generate p values fisher's exact test
            //#************************************************************
            double total1, total2; total1 = 0; total2 = 0;
			//total for first grouping
            for (int i = 0; i < secondGroupingStart; i++) { total1 += totals[i]; }
            
            //total for second grouping
            for (int i = secondGroupingStart; i < column; i++) { total2 += totals[i]; }
            
            vector<double> fish;	fish.resize(row, 0.0);
			vector<double> fish2;	fish2.resize(row, 0.0);
            
			for(int i = 0; i < row; i++){
				
				for(int j = 0; j < secondGroupingStart; j++)		{ fish[i] += data[i][j];	}
				for(int j = secondGroupingStart; j < column; j++)	{ fish2[i] += data[i][j];	}
				
				double f11, f12, f21, f22;
				f11 = fish[i];
				f12 = fish2[i];
				f21 = total1 - fish[i];
				f22 = total2 - fish2[i];
				
				MothurFisher fisher;
				double pre = fisher.fexact(f11, f12, f21, f22);
				if (pre > 0.999999999)	{ pre = 1.0; }
                
				if (m->control_pressed) { return 1; }
				
				pvalues[i] = pre;
			}
            
            //#*************************************
            //#  calculate q values from p values
            //#*************************************
            qvalues = calc_qvalues(pvalues);
            
        }else { //we have multiple subjects per population
            
            //#*************************************
            //#  generate statistics mean, var, stderr    
            //#*************************************
            for(int i = 0; i < row; i++){ // for each taxa
                //# find the mean of each group
                double g1Total = 0.0; double g2Total = 0.0;
                for (int j = 0; j < secondGroupingStart; j++)       {     g1Total += Pmatrix[i][j]; }
                C1[i][0] = g1Total/(double)(secondGroupingStart);
                for (int j = secondGroupingStart; j < column; j++)  {     g2Total += Pmatrix[i][j]; }
                C2[i][0] = g2Total/(double)(column-secondGroupingStart);
                
                 //# find the variance of each group
                double g1Var = 0.0; double g2Var = 0.0;
                for (int j = 0; j < secondGroupingStart; j++)       {     g1Var += pow((Pmatrix[i][j]-C1[i][0]), 2);  }
                C1[i][1] = g1Var/(double)(secondGroupingStart-1);
                for (int j = secondGroupingStart; j < column; j++)  {     g2Var += pow((Pmatrix[i][j]-C2[i][0]), 2);  }
                C2[i][1] = g2Var/(double)(column-secondGroupingStart-1);
                
                //# find the std error of each group -std err^2 (will change to std err at end)
                C1[i][2] = C1[i][1]/(double)(secondGroupingStart);    
                C2[i][2] = C2[i][1]/(double)(column-secondGroupingStart);
            }
            
            //#*************************************
            //#  two sample t-statistics
            //#*************************************
            for(int i = 0; i < row; i++){                  // # for each taxa
                double xbar_diff = C1[i][0] - C2[i][0]; 
                double denom = sqrt(C1[i][2] + C2[i][2]);
                T_statistics[i] = xbar_diff/denom;  // calculate two sample t-statistic
            }
            
           /*for (int i = 0; i < row; i++) {
                for (int j = 0; j < 3; j++) {
                    cout << "C1[" << i+1 << "," << j+1 << "]=" << C1[i][j] << ";" << endl;
                    cout << "C2[" << i+1 << "," << j+1 << "]=" << C2[i][j] << ";" << endl;
                }
                cout << "T_statistics[" << i+1 << "]=" << T_statistics[i] << ";" << endl;
            }
            
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < column; j++) {
                    cout << "Fmatrix[" << i+1 << "," << j+1 << "]=" << data[i][j] << ";" << endl;
                }
            }*/

            //#*************************************
            //# generate initial permuted p-values
            //#*************************************
            pvalues = permuted_pvalues(Pmatrix, T_statistics, data);
            
            //#*************************************
            //#  generate p values for sparse data 
            //#  using fisher's exact test
            //#*************************************
            double total1, total2; total1 = 0; total2 = 0;
			//total for first grouping
            for (int i = 0; i < secondGroupingStart; i++) { total1 += totals[i];  }
            
            //total for second grouping
            for (int i = secondGroupingStart; i < column; i++) { total2 += totals[i];  }
            
            vector<double> fish;	fish.resize(row, 0.0);
			vector<double> fish2;	fish2.resize(row, 0.0);
            
			for(int i = 0; i < row; i++){
				
				for(int j = 0; j < secondGroupingStart; j++)		{ fish[i] += data[i][j];	}
				for(int j = secondGroupingStart; j < column; j++)	{ fish2[i] += data[i][j];	}
                
                if ((fish[i] < secondGroupingStart) && (fish2[i] < (column-secondGroupingStart))) {
    
                    double f11, f12, f21, f22;
                    f11 = fish[i];
                    f12 = fish2[i];
                    f21 = total1 - fish[i];
                    f22 = total2 - fish2[i];
				
                    MothurFisher fisher;
                    double pre = fisher.fexact(f11, f12, f21, f22);
                    if (pre > 0.999999999)	{ pre = 1.0; }
                
                    if (m->control_pressed) { return 1; }
				
                    pvalues[i] = pre;
                }
			}

            //#*************************************
            //#  calculate q values from p values
            //#*************************************
            qvalues = calc_qvalues(pvalues);
            
            //#*************************************
            //#  convert stderr^2 to std error
            //#*************************************
            for(int i = 0; i < row; i++){
                C1[i][2] = sqrt(C1[i][2]);
                C2[i][2] = sqrt(C2[i][2]);
            }
        }
        
        // And now we write the files to a text file.
		struct tm *local;
		time_t t; t = time(NULL);
		local = localtime(&t);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
		out << "Local time and date of test: " << asctime(local) << endl;
		out << "# rows = " << row << ", # col = " << column << ", g = " << secondGroupingStart << endl << endl;
		out << numPermutations << " permutations" << endl << endl;	
		
		//output column headings - not really sure... documentation labels 9 columns, there are 10 in the output file
		//storage 0 = meanGroup1 - line 529, 1 = varGroup1 - line 532, 2 = err rate1 - line 534, 3 = mean of counts group1?? - line 291, 4 = meanGroup2 - line 536, 5 = varGroup2 - line 539, 6 = err rate2 - line 541, 7 = mean of counts group2?? - line 292, 8 = pvalues - line 293
		out << "OTU\tmean(group1)\tvariance(group1)\tstderr(group1)\tmean(group2)\tvariance(group2)\tstderr(group2)\tp-value\tq-value\n";
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { out.close(); return 0; }
			
            //if there are binlabels use them otherwise count.
			if (i < m->currentSharedBinLabels.size()) { out << m->currentSharedBinLabels[i] << '\t'; }
            else { out << (i+1) << '\t'; }
            
            out << C1[i][0] << '\t' << C1[i][1] << '\t' << C1[i][2] << '\t' << C2[i][0] << '\t' << C2[i][1] << '\t' << C2[i][2] << '\t' << pvalues[i] << '\t' << qvalues[i] << endl;
		}  
		
		out << endl << endl;
		out.close();
		


        return 0;
        
    }catch(exception& e) {
        m->errorOut(e, "MothurMetastats", "runMetastats");
        exit(1);
    }	
}
/***********************************************************/
vector<double> MothurMetastats::permuted_pvalues(vector< vector<double> >& Imatrix, vector<double>& tstats, vector< vector<double> >& Fmatrix) {
	try {
        //# matrix stores tstats for each taxa(row) for each permuted trial(column)
        vector<double> ps;  ps.resize(row, 0.0); //# to store the pvalues
        vector< vector<double> > permuted_ttests; permuted_ttests.resize(numPermutations);            
        for (int i = 0; i < numPermutations; i++) { permuted_ttests[i].resize(row, 0.0);  } 
 
        //# calculate null version of tstats using B permutations.
        for (int i = 0; i < numPermutations; i++) {   
            permuted_ttests[i] = permute_and_calc_ts(Imatrix);
        }
        
        //# calculate each pvalue using the null ts
        if ((secondGroupingStart) < 8 || (column-secondGroupingStart) < 8){
            vector< vector<double> > cleanedpermuted_ttests; cleanedpermuted_ttests.resize(numPermutations);  //# the array pooling just the frequently observed ts
            //# then pool the t's together!
            //# count how many high freq taxa there are
            int hfc = 1;
            for (int i = 0; i < row; i++) {                 // # for each taxa
                double group1Total = 0.0; double group2Total = 0.0;
                for(int j = 0; j < secondGroupingStart; j++)		{ group1Total += Fmatrix[i][j];	}
				for(int j = secondGroupingStart; j < column; j++)	{ group2Total += Fmatrix[i][j];	}
                
                if (group1Total >= secondGroupingStart || group2Total >= (column-secondGroupingStart)){ 
                    hfc++;
                    for (int j = 0; j < numPermutations; j++) {   cleanedpermuted_ttests[j].push_back(permuted_ttests[j][i]); }
                }
            }
            
            //#now for each taxa
            for (int i = 0; i < row; i++) { 
                //number of cleanedpermuted_ttests greater than tstat[i]
                int numGreater = 0;
                for (int j = 0; j < numPermutations; j++) {
                    for (int k = 0; k < hfc; k++) {
                        if (cleanedpermuted_ttests[j][k] > abs(tstats[i])) { numGreater++; }
                    }
                }
                
                ps[i] = (1/(double)(numPermutations*hfc))*numGreater;
            }
        }else{
            for (int i = 0; i < row; i++) { 
                //number of permuted_ttests[i] greater than tstat[i] //(sum(permuted_ttests[i,] > abs(tstats[i]))+1)
                int numGreater = 1;
                for (int j = 0; j < numPermutations; j++) { if (permuted_ttests[j][i] > abs(tstats[i])) { numGreater++; }   }
                ps[i] = (1/(double)(numPermutations+1))*numGreater;
            }
        }
        
        return ps;
        
    }catch(exception& e) {
        m->errorOut(e, "MothurMetastats", "permuted_pvalues");
        exit(1);
    }	
}
/***********************************************************/
vector<double> MothurMetastats::permute_and_calc_ts(vector< vector<double> >& Imatrix) {
	try {
        vector< vector<double> > permutedMatrix = Imatrix;
        
        //randomize columns, ie group abundances.
        map<int, int> randomMap;
        vector<int> randoms;
        for (int i = 0; i < column; i++) { randoms.push_back(i); }
        random_shuffle(randoms.begin(), randoms.end());
        for (int i = 0; i < randoms.size(); i++) {   randomMap[i] = randoms[i];   }
        
        //calc ts
        vector< vector<double> > C1; C1.resize(row);
        for (int i = 0; i < row; i++) { C1[i].resize(3, 0.0);  } // statistic profiles for class1 and class 2
        vector< vector<double> > C2; C2.resize(row);            // mean[1], variance[2], standard error[3] 
        for (int i = 0; i < row; i++) { C2[i].resize(3, 0.0);  } 
        vector<double> Ts; Ts.resize(row, 0.0); // a place to store the true t-statistics 

        //#*************************************
        //#  generate statistics mean, var, stderr    
        //#*************************************
        for(int i = 0; i < row; i++){ // for each taxa
            //# find the mean of each group
            double g1Total = 0.0; double g2Total = 0.0;
            for (int j = 0; j < secondGroupingStart; j++)       {     g1Total += permutedMatrix[i][randomMap[j]]; }
            C1[i][0] = g1Total/(double)(secondGroupingStart);
            for (int j = secondGroupingStart; j < column; j++)  {     g2Total += permutedMatrix[i][randomMap[j]]; }
            C2[i][0] = g2Total/(double)(column-secondGroupingStart);
            
            //# find the variance of each group
            double g1Var = 0.0; double g2Var = 0.0;
            for (int j = 0; j < secondGroupingStart; j++)       {     g1Var += pow((permutedMatrix[i][randomMap[j]]-C1[i][0]), 2);  }
            C1[i][1] = g1Var/(double)(secondGroupingStart-1);
            for (int j = secondGroupingStart; j < column; j++)  {     g2Var += pow((permutedMatrix[i][randomMap[j]]-C2[i][0]), 2);  }
            C2[i][1] = g2Var/(double)(column-secondGroupingStart-1);
            
            //# find the std error of each group -std err^2 (will change to std err at end)
            C1[i][2] = C1[i][1]/(double)(secondGroupingStart);    
            C2[i][2] = C2[i][1]/(double)(column-secondGroupingStart);
        }
        
        

        //#*************************************
        //#  two sample t-statistics
        //#*************************************
        for(int i = 0; i < row; i++){                  // # for each taxa
            double xbar_diff = C1[i][0] - C2[i][0]; 
            double denom = sqrt(C1[i][2] + C2[i][2]);
            Ts[i] = abs(xbar_diff/denom);  // calculate two sample t-statistic
        }

        return Ts;

        
    }catch(exception& e) {
        m->errorOut(e, "MothurMetastats", "permuted_ttests");
        exit(1);
    }	
}
/***********************************************************/
vector<double> MothurMetastats::calc_qvalues(vector<double>& pValues) {
	try {
        int numRows = pValues.size();
		vector<double> qvalues(numRows, 0.0);

		//fill lambdas with 0.00, 0.01, 0.02... 0.95
		vector<double> lambdas(96, 0);
		for (int i = 1; i < lambdas.size(); i++) { lambdas[i] = lambdas[i-1] + 0.01; }
		
		vector<double> pi0_hat(lambdas.size(), 0);
		
		//calculate pi0_hat
		for (int l = 0; l < lambdas.size(); l++){ // for each lambda value
			int count = 0;
			for (int i = 0; i < numRows; i++){ // for each p-value in order
				if (pValues[i] > lambdas[l]){ count++; }
			}
			pi0_hat[l] = count/(double)(numRows*(1.0-lambdas[l]));
            //cout << lambdas[l] << '\t' << count << '\t' << numRows*(1.0-lambdas[l]) << '\t' << pi0_hat[l] << endl;
		}
		
		double pi0 = smoothSpline(lambdas, pi0_hat, 3);
		
		//order p-values
		vector<double> ordered_qs = qvalues;
		vector<int> ordered_ps(pValues.size(), 0);
		for (int i = 1; i < ordered_ps.size(); i++) {  ordered_ps[i] = ordered_ps[i-1] + 1; }
		vector<double> tempPvalues = pValues;
		OrderPValues(0, numRows-1, tempPvalues, ordered_ps);
		
		ordered_qs[numRows-1] = min((pValues[ordered_ps[numRows-1]]*pi0), 1.0);
		for (int i = (numRows-2); i >= 0; i--){
			double p = pValues[ordered_ps[i]];
			double temp = p*numRows*pi0/(double)(i+1);

			ordered_qs[i] = min(temp, ordered_qs[i+1]);
		}
		
		//re-distribute calculated qvalues to appropriate rows
		for (int i = 0; i < numRows; i++){
			qvalues[ordered_ps[i]] = ordered_qs[i];
		}
		
		return qvalues;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "calc_qvalues");
		exit(1);
	}
}
/***********************************************************/
int MothurMetastats::OrderPValues(int low, int high, vector<double>& p, vector<int>& order) {
	try {
		
		if (low < high) {
			int i = low+1;
			int j = high;
			int pivot = (low+high) / 2;
			
			swapElements(low, pivot, p, order);  //puts pivot in final spot
			
			/* compare value */
			double key = p[low];
			
			/* partition */
			while(i <= j) {
				/* find member above ... */
				while((i <= high) && (p[i] <= key))	{  i++;  }  
				
				/* find element below ... */
				while((j >= low) && (p[j] > key))	{  j--;  } 
				
				if(i < j) {
					swapElements(i, j, p, order);
				}
			} 
			
			swapElements(low, j, p, order);
			
			/* recurse */
			OrderPValues(low, j-1, p, order);
			OrderPValues(j+1, high, p, order); 
		}		
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "OrderPValues");
		exit(1);
	}
}
/***********************************************************/
int MothurMetastats::swapElements(int i, int j, vector<double>& p, vector<int>& order) {
	try {
		
		double z = p[i];
		p[i] = p[j];
		p[j] = z;
		
		int temp = order[i];
		order[i] = order[j];
		order[j] = temp;
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "swapElements");
		exit(1);
	}
}
/***********************************************************/
double MothurMetastats::smoothSpline(vector<double>& x, vector<double>& y, int df) {
	try {
				
		double result = 0.0;
		int n = x.size();
		vector<double> w(n, 1);
		double* xb = new double[n];
		double* yb = new double[n];
		double* wb = new double[n];
		double yssw = 0.0;
		for (int i = 0; i < n; i++) {
			wb[i] = w[i];
			yb[i] = w[i]*y[i];
			yssw += (w[i] * y[i] * y[i]) - wb[i] * (yb[i] * yb[i]);
			xb[i] = (x[i] - x[0]) / (x[n-1] - x[0]);
		}
		
		vector<double> knot = sknot1(xb, n);
		int nk = knot.size() - 4;

		double low = -1.5; double high = 1.5; double tol = 1e-04; double eps = 2e-08; int maxit = 500;
		int ispar = 0; int icrit = 3; double dofoff = 3.0;
		double penalty = 1.0; 
		int ld4 = 4; int isetup = 0; int ldnk = 1; int ier; double spar = 1.0; double crit;
		
		double* knotb = new double[knot.size()];
		double* coef1 = new double[nk];
		double* lev1 = new double[n];
		double* sz1 = new double[n];
		for (int i = 0; i < knot.size(); i++) {	knotb[i] = knot[i];	}
		
		Spline spline;
		spline.sbart(&penalty, &dofoff, &xb[0], &yb[0], &wb[0], &yssw, &n, &knotb[0], &nk, &coef1[0], &sz1[0], &lev1[0], &crit,
				&icrit, &spar, &ispar, &maxit, &low, &high, &tol, &eps, &isetup, &ld4, &ldnk, &ier);
		
		result = coef1[nk-1];
		
		//free memory
		delete [] xb;
		delete [] yb;
		delete [] wb;
		delete [] knotb;
		delete [] coef1;
		delete [] lev1;
		delete [] sz1;
							
		return result;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "smoothSpline");
		exit(1);
	}
}
/***********************************************************/
vector<double> MothurMetastats::sknot1(double* x, int n) {
	try {
		vector<double> knots;
		int nk = nkn(n);
		
		//R equivalent - rep(x[1L], 3L)
		knots.push_back(x[0]); knots.push_back(x[0]); knots.push_back(x[0]);
		
		//generate a sequence of nk equally spaced values from 1 to n. R equivalent = seq.int(1, n, length.out = nk)
		vector<int> indexes = getSequence(0, n-1, nk);
		for (int i = 0; i < indexes.size(); i++) { knots.push_back(x[indexes[i]]);  } 
		
		//R equivalent - rep(x[n], 3L)
		knots.push_back(x[n-1]); knots.push_back(x[n-1]); knots.push_back(x[n-1]);
				
		return knots;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "sknot1");
		exit(1);
	}
}
/***********************************************************/
vector<int> MothurMetastats::getSequence(int start, int end, int length) {
	try {
		vector<int> sequence;
		double increment = (end-start) / (double) (length-1);
		
		sequence.push_back(start);
		for (int i = 1; i < length-1; i++) {
			sequence.push_back(int(i*increment));
		}
		sequence.push_back(end);
		
		return sequence;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "getSequence");
		exit(1);
	}
}	
/***********************************************************/
//not right, havent fully figured out the variable types yet...
int MothurMetastats::nkn(int n) {
	try {
		
		if (n < 50) { return n; }
		else {
			double a1 = log2(50);
			double a2 = log2(100);
			double a3 = log2(140);
			double a4 = log2(200);
			
			if (n < 200) {
				return (int)pow(2.0, (a1 + (a2 - a1) * (n - 50)/(float)150));
			}else if (n < 800) {
				return (int)pow(2.0, (a2 + (a3 - a2) * (n - 200)/(float)600));
			}else if (n < 3200) {
				return (int)pow(2.0, (a3 + (a4 - a3) * (n - 800)/(float)2400));
			}else {
				return (int)pow((double)(200 + (n - 3200)), 0.2);
			}
		}
	
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "nkn");
		exit(1);
	}
}
/***********************************************************/





