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
		for (int i = 0; i < secondGroupingStart; i++) { total1 += total[i]; }
		
		//total for second grouping
		for (int i = secondGroupingStart; i < column; i++) { total2 += total[i]; }
		
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
		if ( (column == 2) || (secondGroupingStart < 8) || ((column-secondGroupingStart) < 8) ){ 
			
			vector<double> fish;	fish.resize(row, 0.0);
			vector<double> fish2;	fish2.resize(row, 0.0);
			
			for(int i = 0; i < row; i++){
				
				for(int j = 0; j < secondGroupingStart; j++)		{ fish[i] += data[i][j];	}
				for(int j = secondGroupingStart; j < column; j++)	{ fish2[i] += data[i][j];	}
				
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
				
				for(int j = 0; j < secondGroupingStart; j++)	{	sparse[i] += data[i][j];	}
				if(sparse[i] < (double)secondGroupingStart)		{	c++;						}
				
				// ?<= for col
				for(int j = secondGroupingStart; j < column; j++)		{  sparse2[i] += data[i][j]; }
				if( (sparse2[i] < (double)(column-secondGroupingStart)))	{ c++;						 }
				
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
			
			for (int i = 0; i < secondGroupingStart; i++){ temp[j][0] += data[j][i]; }
			temp[j][0] /= (double)secondGroupingStart;
			
			for(int i = secondGroupingStart; i < column; i++){ temp[j][1] += data[j][i]; }
			temp[j][1] /= (double)(column-secondGroupingStart);
		}
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { return 1; }
			
			storage[i][3]=temp[i][0];
			storage[i][7]=temp[i][1];
			storage[i][8]=pvalues[i];
		}
		
		vector<double> qvalues = calc_qvalues(pvalues);
		
		// BACKUP checks
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		for (int i = 0; i < row; i++){
			
			if (m->control_pressed) { return 1; }
			
			if(qvalues[i] < threshold){
				m->mothurOut("Feature " + toString((i+1)) + " is significant, q = "); 
				cout << qvalues[i];
				m->mothurOutJustToLog(toString(qvalues[i])); m->mothurOutEndLine();
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
		out << "OTU\tmean(group1)\tvariance(group1)\tstderr(group1)\tmean_of_counts(group1)\tmean(group2)\tvariance(group2)\tstderr(group2)\tmean_of_counts(group1)\tp-value\tq-value\n";
		
		for(int i = 0; i < row; i++){
			if (m->control_pressed) { out.close(); return 0; }
			
			out << (i+1);
			for(int j = 0; j < 9; j++){ out << '\t' << storage[i][j]; }
			out << '\t' << qvalues[i];
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
			C1[i][2]=C1[i][1]/(secondGroupingStart);
			storage[i][2]=sqrt(C1[i][2]);
			
			C2[i][0]=tool[i+row]; // mean group 2
			storage[i][4]=C2[i][0];    
			C2[i][1]=tool[i+row+row+row]; // var group 2 
			storage[i][5]=C2[i][1];        
			C2[i][2]=C2[i][1]/(column-secondGroupingStart);
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
		
		double a = secondGroupingStart;
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
			C1[i][2]=C1[i][1]/(secondGroupingStart);
			
			C2[i][0]=tool[i+row];
			C2[i][1]=tool[i+row+row+row]; // var group 2 
			C2[i][2]=C2[i][1]/(column-secondGroupingStart);
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
vector<double> MothurMetastats::calc_qvalues(vector<double>& pValues) {
	try {
		
       /* cout << "x <- c(" << pValues[0];
        for (int l = 1; l < pValues.size(); l++){
            cout << ", " << pValues[l];
        }
        cout << ")\n";*/
        
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
			pi0_hat[l] = count/(double)(numRows*(1-lambdas[l]));
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





