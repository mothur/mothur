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
MothurMetastats::~MothurMetastats() {
	try {
		
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "~MothurMetastats");
		exit(1);
	}	
}
/***********************************************************/
//main metastats function
int MothurMetastats::runMetastats(string outputFileName, vector< vector<double> >& data, int secondGroupingStart) {
	try {
		int bflag = 0;
		row = data.size();		 //numBins
		column = data[0].size(); //numGroups in subset
		int size = row*column;
		
		//consistent with original, but this should never be true
		if ((secondGroupingStart >= column) || (secondGroupingStart <= 0)) { m->mothurOut("[ERROR]: Check your g value."); m->mothurOutEndLine(); return 0; }
		
		//Initialize the matrices
		vector<double> pmatrix; pmatrix.resize(size, 0.0);
		vector<double> permuted; permuted.resize(size, 0.0);
		vector< vector<double> > storage; storage.resize(row);
		for (int i = 0; i < storage.size(); i++) { storage[i].resize(9, 0.0); }
		
		//Produces the sum of each column
		vector<double> total; total.resize(column, 0.0);
		vector<double> ratio; ratio.resize(column, 0.0);
		double total1 = 0.0; double total2 = 0.0;
		
		//total[i] = total abundance for group[i]
		for (int i = 0; i < column; i++) {
			for (int j = 0; j < row; j++) {
				total[i] += data[j][i];
			}
		}
		
		//total for first grouping
		for (int i = 0; i < (secondGroupingStart-1); i++) { total1 += total[i]; }
		
		//total for second grouping
		for (int i = (secondGroupingStart-1); i < column; i++) { total2 += total[i]; }
		
		//Creates the ratios by first finding the minimum of totals
		double min = total[0];
		for (int i = 0; i < total.size(); i++) {
			 if (total[i] < min) { min = total[i]; }
		}
		
		//sanity check
		if (min <= 0.0) { m->mothurOut("[ERROR]: the sum of one of the columns <= 0."); m->mothurOutEndLine(); return 0; }
		
		//Ratio time...
		for(int i = 0; i < ratio.size(); i++){  ratio[i] = total[i] / min; }
		
		//Change matrix into an array as received by R for compatibility - kept to be consistent with original
		int count = 0;
		for(int i = 0; i < column; i++){
			for(int j = 0; j < row; j++){
				pmatrix[count]=data[j][i];
				count++;
			}
		}
		
		if(row == 1){
			for (int i =0; i < column; i++){ pmatrix[i] /= ratio[i]; }
		}else {
			count = 0; int j=-1;
			
			for (int i=0; i < size; i++) {
				if (count % row == 0) { j++; }
				pmatrix[i] /= ratio[j];
				count++; 
			}   
		}
		
		vector<double> permuted_ttests; permuted_ttests.resize(row, 0.0);
		vector<double> pvalues;			pvalues.resize(row, 0.0);
		vector<double> tinitial;		tinitial.resize(row, 0.0);
		
		if (m->control_pressed) { return 1; }
		
		//Find the initial values for the matrix.
		start(pmatrix, secondGroupingStart, tinitial, storage);
		
		if (m->control_pressed) { return 1; }
		
		// Start the calculations.
		if ( (column == 2) || ((secondGroupingStart-1) < 8) || ((column-secondGroupingStart+1) < 8) ){ 
			
			vector<double> fish;	fish.resize(row, 0.0);
			vector<double> fish2;	fish2.resize(row, 0.0);
			
			for(int i = 0; i < row; i++){
				
				for(int j = 0; j < (secondGroupingStart-1); j++)		{ fish[i] += data[i][j];	}
				for(int j = (secondGroupingStart-1); j < column; j++)	{ fish2[i] += data[i][j];	}
				
				//vector<double> tempData; tempData.resize(4, 0.0);
				double f11, f12, f21, f22;
				f11 = fish[i];
				f12 = fish2[i];
				f21 = total1 - fish[i];
				f22 = total2 - fish2[i];
				
				double pre = 0.0;
				
				MothurFisher fisher;
				pre = fisher.fexact(f11, f12, f21, f22);
				
				if (m->control_pressed) { return 1; }
				
				if (pre > 0.999999999)	{ pre = 1.0; }
				storage[i][8] = pre;
				pvalues[i] = pre;
			}
			
		}else {
	
			testp(permuted_ttests, permuted, pmatrix, secondGroupingStart, tinitial, pvalues);
			
			if (m->control_pressed) { return 1; }
			
			// Checks to make sure the matrix isn't sparse.
			vector<double> sparse;		sparse.resize(row, 0.0);
			vector<double> sparse2;		sparse2.resize(row, 0.0);
			
			int c = 0;
			
			for(int i = 0; i < row; i++){
				
				for(int j = 0; j < (secondGroupingStart-1); j++){	sparse[i] += data[i][j]; }
				if(sparse[i] < (double)(secondGroupingStart-1)){	c++; }
				
				// ?<= for col
				for(int j = (secondGroupingStart-1); j < column; j++){  sparse2[i] += data[i][j]; }
				if( (sparse2[i] < (double)(column-secondGroupingStart+1))) { c++; }
				
				if (c == 2) {
					c=0;
					double f11,f12,f21,f22;
					
					f11=sparse[i];  sparse[i]=0;
					f12=sparse2[i];  sparse2[i]=0;
					f21 = total1 - f11;
					f22 = total2 - f12;
					
					double pre = 0.0;
					
					MothurFisher fisher;
					pre = fisher.fexact(f11, f12, f21, f22);
					
					if (m->control_pressed) { return 1; }
					
					if (pre > 0.999999999){
						pre = 1.0;
					}
					
					storage[i][8] = pre;
					pvalues[i] = pre;
				}				
			}
			
			bflag = 1;
		}

		// Calculates the mean of counts (not normalized)
		vector< vector<double> > temp; temp.resize(row);
		for (int i = 0; i < temp.size(); i++) { temp[i].resize(2, 0.0); }
		
		for (int j = 0; j < row; j++){
			if (m->control_pressed) { return 1; }
			
			for (int i = 1; i <= (secondGroupingStart-1); i++){ temp[j][0] += data[j][i-1]; }
			temp[j][0] /= (double)(secondGroupingStart-1);
			
			for(int i = secondGroupingStart; i <= column; i++){ temp[j][1] += data[j][i-1]; }
			temp[j][1] /= (double)(column-secondGroupingStart+1);
		}
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { return 1; }
			
			storage[i][3]=temp[i][0];
			storage[i][7]=temp[i][1];
			storage[i][8]=pvalues[i];
		}
		
		// BACKUP checks
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		for (int i = 0; i < row; i++){
			
			if (m->control_pressed) { return 1; }
			
			if(pvalues[i] < threshold){
				m->mothurOut("Feature " + toString((i+1)) + " is significant, p = "); 
				cout << pvalues[i];
				m->mothurOutJustToLog(toString(pvalues[i])); m->mothurOutEndLine();
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
		if (bflag == 1){ out << numPermutations << " permutations" << endl << endl;	}
		
		//output column headings - not really sure... documentation labels 9 columns, there are 10 in the output file
		//storage 0 = meanGroup1 - line 529, 1 = varGroup1 - line 532, 2 = err rate1 - line 534, 3 = mean of counts group1?? - line 291, 4 = meanGroup2 - line 536, 5 = varGroup2 - line 539, 6 = err rate2 - line 541, 7 = mean of counts group2?? - line 292, 8 = pvalues - line 293
		out << "OTU\tmean(group1)\tvariance(group1)\tstderr(group1)\tmean_of_counts(group1)\tmean(group2)\tvariance(group2)\tstderr(group2)\tmean_of_counts(group1)\tp-value\n";
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { out.close(); return 0; }
			
			out << (i+1);
			for(int j = 0; j < 9; j++){ out << '\t' << storage[i][j]; }
			out << endl;
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
//Find the initial values for the matrix
int MothurMetastats::start(vector<double>& Imatrix, int secondGroupingStart, vector<double>& initial, vector< vector<double> >& storage) {
	try {
		
		int a = row; a*=4;
		
		double xbardiff = 0.0; double denom = 0.0;
		vector<double> store;	store.resize(a, 0.0);
		vector<double> tool;	tool.resize(a, 0.0);
		vector< vector<double> > C1; C1.resize(row);
		for (int i = 0; i < C1.size(); i++) { C1[i].resize(3, 0.0); }
		vector< vector<double> > C2; C2.resize(row);
		for (int i = 0; i < C2.size(); i++) { C2[i].resize(3, 0.0); }
		
		meanvar(Imatrix, secondGroupingStart, store);
		
		if (m->control_pressed) { return 0; }
		
		//copy store into tool
		tool = store;
		
		for (int i = 0; i < row; i++){
			C1[i][0]=tool[i]; //mean group 1
			storage[i][0]=C1[i][0];
			C1[i][1]=tool[i+row+row]; // var group 1
			storage[i][1]=C1[i][1];
			C1[i][2]=C1[i][1]/(secondGroupingStart-1);
			storage[i][2]=sqrt(C1[i][2]);
			
			C2[i][0]=tool[i+row]; // mean group 2
			storage[i][4]=C2[i][0];    
			C2[i][1]=tool[i+row+row+row]; // var group 2 
			storage[i][5]=C2[i][1];        
			C2[i][2]=C2[i][1]/(column-secondGroupingStart+1);
			storage[i][6]=sqrt(C2[i][2]);   
		}
		
		if (m->control_pressed) { return 0; }
		
		for (int i = 0; i < row; i++){
			xbardiff = C1[i][0]-C2[i][0];
			denom = sqrt(C1[i][2]+C2[i][2]);
			initial[i]=fabs(xbardiff/denom);
		}	

		return 0; 
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "start");
		exit(1);
	}	
}
/***********************************************************/
int MothurMetastats::meanvar(vector<double>& pmatrix, int secondGroupingStart, vector<double>& store) {
	try {
		vector<double> temp;	temp.resize(row, 0.0);
		vector<double> temp2;	temp2.resize(row, 0.0);
		vector<double> var;		var.resize(row, 0.0);
		vector<double> var2;	var2.resize(row, 0.0);
		
		double a = secondGroupingStart-1;
		double b = column - a;
		int m = a * row;
		int n = row * column;
		
		for (int i = 0; i < m; i++)		{ temp[i%row] += pmatrix[i];	}
		for (int i = 0; i < n; i++)		{ temp2[i%row]+= pmatrix[i];	}
		for (int i = 0; i < row; i++)	{ temp2[i] -= temp[i];		}
		for (int i = 0; i <= row-1;i++)	{
			store[i] = temp[i]/a;
			store[i+row]=temp2[i]/b;
		}
		
		//That completes the mean calculations.
		
		for (int i = 0; i < m; i++)		{ var[i%row] += pow((pmatrix[i]-store[i%row]),2);		}
		for (int i = m; i < n; i++)		{ var2[i%row]+= pow((pmatrix[i]-store[(i%row)+row]),2); }
		for (int i = 0; i <= row-1; i++){
			store[i+2*row]=var[i]/(a-1);
			store[i+3*row]=var2[i]/(b-1);
		}
		
		// That completes var calculations.
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "meanvar");
		exit(1);
	}	
}
/***********************************************************/
int MothurMetastats::testp(vector<double>& permuted_ttests, vector<double>& permuted, vector<double>& Imatrix, int secondGroupingStart, vector<double>& Tinitial, vector<double>& ps) {
	try {
		
		vector<double> Tvalues;		Tvalues.resize(row, 0.0);
		vector<double> counter;		counter.resize(row, 0.0);
		int a, b, n;
		
		a = numPermutations;
		b = row;
		n = a*b;
		
		for (int j = 1; j <= row; j++)	{	
			if (m->control_pressed) { return 0; }
			permute_matrix(Imatrix, permuted, secondGroupingStart, Tvalues, Tinitial, counter);	
		}
		
		for(int j = 0; j < row; j++)	{	
			if (m->control_pressed) { return 0; }
			ps[j] = ((counter[j]+1)/(double)(a+1)); 
		}
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "testp");
		exit(1);
	}	
}	
/***********************************************************/
int MothurMetastats::permute_matrix(vector<double>& Imatrix, vector<double>& permuted, int secondGroupingStart, vector<double>& trial_ts, vector<double>& Tinitial, vector<double>& counter1){
	try {
	
		vector<int> y; y.resize(column, 0);
		for (int i = 1; i <= column; i++){ y[i-1] = i; }
		
		permute_array(y); 
		
		int f = 0; int c = 0; int k = 0;
		for (int i = 0; i < column; i++){
			
			if (m->control_pressed) { return 0; }
			
			f = y[i]; //column number
			c = 1;
			c *= (f-1);
			c *= row;
			if (f == 1){ c = 0; } // starting value position in the Imatrix
			
			for(int j = 1; j <= row; j++){
				permuted[k] = Imatrix[c];
				c++; k++;
			}
		}
		
		calc_twosample_ts(permuted, secondGroupingStart, trial_ts, Tinitial, counter1);
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "permute_matrix");
		exit(1);
	}	
}
/***********************************************************/
int MothurMetastats::permute_array(vector<int>& array) {
	try {
		static int seeded = 0;
		
		if (! seeded) {
			seeded = 1;
			srand(time(NULL));
		}
		
		for (int i = 0; i < array.size(); i++) {
			if (m->control_pressed) { return 0; }
			
			int selection = rand() % (array.size() - i);
			int tmp = array[i + selection];
			array[i + selection] = array[i];
			array[i] = tmp;
		}
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "permute_array");
		exit(1);
	}	
}
/***********************************************************/
int MothurMetastats::calc_twosample_ts(vector<double>& Pmatrix, int secondGroupingStart, vector<double>& Ts, vector<double>& Tinitial, vector<double>& counter) {
	try {
		int a = row * 4;
		
		vector< vector<double> > C1; C1.resize(row);
		for (int i = 0; i < C1.size(); i++) { C1[i].resize(3, 0.0); }
		vector< vector<double> > C2; C2.resize(row);
		for (int i = 0; i < C2.size(); i++) { C2[i].resize(3, 0.0); }
		vector<double> storage; storage.resize(a, 0.0);
		vector<double> tool;	tool.resize(a, 0.0);
		double xbardiff = 0.0; double denom = 0.0;
		
		meanvar(Pmatrix, secondGroupingStart, storage);
		
		for(int i = 0;i <= (a-1); i++) {	
			if (m->control_pressed) { return 0; }
			tool[i] = storage[i];	
		}
		
		for (int i = 0; i < row; i++){
			if (m->control_pressed) { return 0; }
			C1[i][0]=tool[i];
			C1[i][1]=tool[i+row+row];
			C1[i][2]=C1[i][1]/(secondGroupingStart-1);
			
			C2[i][0]=tool[i+row];
			C2[i][1]=tool[i+row+row+row]; // var group 2 
			C2[i][2]=C2[i][1]/(column-secondGroupingStart+1);
		}
		
		for (int i = 0; i < row; i++){
			if (m->control_pressed) { return 0; }
			xbardiff = C1[i][0]-C2[i][0];
			denom = sqrt(C1[i][2]+C2[i][2]);
			Ts[i]=fabs(xbardiff/denom);
			if (fabs(Ts[i])>(fabs(Tinitial[i])+.0000000000001)){ //13th place
				counter[i]++;
			}
		}
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "MothurMetastats", "calc_twosample_ts");
		exit(1);
	}
}
/***********************************************************/


