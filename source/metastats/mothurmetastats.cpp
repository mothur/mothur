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
				double pre = fisher.fexact(f11, f12, f21, f22, m->currentSharedBinLabels[i]);
				if (pre > 0.999999999)	{ pre = 1.0; }
                
				if (m->control_pressed) { return 1; }
				
				pvalues[i] = pre;
			}
            
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
            
            if (m->debug) {
                for (int i = 0; i < row; i++) {
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
                }
            }
            //#*************************************
            //# generate initial permuted p-values
            //#*************************************
            pvalues = permuted_pvalues(Pmatrix, T_statistics, data);
            
            if (m->debug) {  for (int i = 0; i < row; i++) { m->mothurOut("[DEBUG]: " + m->currentSharedBinLabels[i] + " pvalue = " + toString(pvalues[i]) + "\n"); } }
            
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
                    if (m->debug) {   m->mothurOut("[DEBUG]: about to run fisher for Otu " + m->currentSharedBinLabels[i] + " F11, F12, F21, F22 = " + toString(f11) + " " + toString(f12) + " " + toString(f21) + " " + toString(f22) + " " + "\n"); }
                    double pre = fisher.fexact(f11, f12, f21, f22, m->currentSharedBinLabels[i]);
                    if (m->debug) {   m->mothurOut("[DEBUG]: about to completed fisher for Otu " + m->currentSharedBinLabels[i] + " pre = " + toString(pre) + "\n"); }
                    
                    if (pre > 0.999999999)	{ pre = 1.0; }
                
                    if (m->control_pressed) { return 1; }
				
                    pvalues[i] = pre;
                }
			}
            
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
		out << "OTU\tmean(group1)\tvariance(group1)\tstderr(group1)\tmean(group2)\tvariance(group2)\tstderr(group2)\tp-value\n";
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { out.close(); return 0; }
			
            //if there are binlabels use them otherwise count.
			if (i < m->currentSharedBinLabels.size()) { out << m->currentSharedBinLabels[i] << '\t'; }
            else { out << (i+1) << '\t'; }
            
            out << C1[i][0] << '\t' << C1[i][1] << '\t' << C1[i][2] << '\t' << C2[i][0] << '\t' << C2[i][1] << '\t' << C2[i][2] << '\t' << pvalues[i] <<  endl;
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
        m->mothurRandomShuffle(randoms);
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





