
/*
 *  pcacommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "pcacommand.h"

//**********************************************************************************************************************

PCACommand::PCACommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip","lt"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; abort = true; }	
			else {	filename = phylipfile;  }
			
			//columnfile = validParameter.validFile(parameters, "column", true);
			//if (columnfile == "not open") { abort = true; }	
			//else if (columnfile == "not found") { columnfile = ""; }
			//else {  format = "column";	}
			
			//namefile = validParameter.validFile(parameters, "name", true);
			//if (namefile == "not open") { abort = true; }	
			//else if (namefile == "not found") { namefile = ""; }
			
			
			//error checking on files	
			if (phylipfile == "")	{ mothurOut("You must provide a distance file before running the pca command."); mothurOutEndLine(); abort = true; }		
			//if ((phylipfile == "") && (columnfile == ""))	{ mothurOut("You must provide a distance file before running the pca command."); mothurOutEndLine(); abort = true; }
			//else if ((phylipfile != "") && (columnfile != "")) { mothurOut("You may not use both the column and the phylip parameters."); mothurOutEndLine(); abort = true; }
			
			//if (columnfile != "") {
			//	if (namefile == "") {  mothurOut("You need to provide a namefile if you are going to use the column format."); mothurOutEndLine(); abort = true; }
			//}
			
			string temp = validParameter.validFile(parameters, "lt", false);				if (temp == "not found") { temp = "false"; }
			bool lt = isTrue(temp);
			
			if (lt)		{  matrix = 2;	}
			else		{  matrix = 1;  }


		}

	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "PCACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void PCACommand::help(){
	try {
	
		mothurOut("The pca command..."); mothurOutEndLine();
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
PCACommand::~PCACommand(){}
//**********************************************************************************************************************
int PCACommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);
		
		vector<string> names;
		vector<vector<double> > D;
	
		fbase = filename;
		if(fbase.find_last_of(".")!=string::npos){
			fbase.erase(fbase.find_last_of(".")+1); 
		}
		else{
			fbase += ".";
		}
		read(filename, matrix, names, D);
   	
		double offset = 0.0000;
		vector<double> d;
		vector<double> e;
		vector<vector<double> > G = D;
		vector<vector<double> > copy_G;
		int rank = D.size();
		
		cout << "\nProcessing...\n";
		
		for(int count=0;count<2;count++){
			recenter(offset, D, G);   
			tred2(G, d, e);
			qtli(d, e, G);
			offset = d[d.size()-1];
			if(offset > 0.0) break;
		} 
		
		
		output(fbase, names, G, d);
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "execute");
		exit(1);
	}
}
/*********************************************************************************************************************************/

inline double SIGN(const double a, const double b)
{
    return b>=0 ? (a>=0 ? a:-a) : (a>=0 ? -a:a);
}

/*********************************************************************************************************************************/

void PCACommand::get_comment(istream& f, char begin, char end){
	try {
		char d=f.get();
		while(d != end){	d = f.get();	}
		d = f.peek();
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "get_comment");
		exit(1);
	}
}	

/*********************************************************************************************************************************/

void PCACommand::read_mega(istream& f, int square_m, vector<string>& name_list, vector<vector<double> >& d){
	try {
		get_comment(f, '#', '\n');
		
		char test = f.peek();
		
		while(test == '!'){								//get header comments
			get_comment(f, '!', ';');
			while(isspace(test=f.get()))		{;}
			f.putback(test);
			test = f.peek();
		}
		while(test != '\n'){							//get sequence names
			get_comment(f, '[', ']');
			char d = f.get();
			d = f.get();
			if(d == '#'){
				string name;
				f >> name;
				name_list.push_back(name);
				while(isspace(test=f.get()))		{;}
				f.putback(test);
			}
			else{
				break;
			}
		}
		int rank = name_list.size();
		d.resize(rank);
		for(int i=0;i<rank;i++){		d[i].resize(rank);		}
		
		d[0][0] = 0.0000;
		get_comment(f, '[', ']');
		for(int i=1;i<rank;i++){
			get_comment(f, '[', ']');
			d[i][i]=0.0000;
			for(int j=0;j<i;j++){
				f >> d[i][j];
				if (d[i][j] == -0.0000)
					d[i][j] = 0.0000;
				d[j][i]=d[i][j];
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "read_mega");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void PCACommand::read_phylip(istream& f, int square_m, vector<string>& name_list, vector<vector<double> >& d){
	try {
		//     int count1=0;
		//     int count2=0;
		
		int rank;
		f >> rank;
		
		name_list.resize(rank);
		d.resize(rank);
		if(square_m == 1){
			for(int i=0;i<rank;i++)
				d[i].resize(rank);
			for(int i=0;i<rank;i++) {
				f >> name_list[i];
				//			cout << i << "\t" << name_list[i] << endl;
				for(int j=0;j<rank;j++) {
					f >> d[i][j];
					if (d[i][j] == -0.0000)
						d[i][j] = 0.0000;
				}
			}
		}
		else if(square_m == 2){
			for(int i=0;i<rank;i++){
				d[i].resize(rank);
			}
			d[0][0] = 0.0000;
			f >> name_list[0];
			for(int i=1;i<rank;i++){
				f >> name_list[i];
				d[i][i]=0.0000;
				for(int j=0;j<i;j++){
					f >> d[i][j];
					if (d[i][j] == -0.0000)
						d[i][j] = 0.0000;
					d[j][i]=d[i][j];
				}
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "read_phylip");
		exit(1);
	}

}

/*********************************************************************************************************************************/

void PCACommand::read(string fname, int m, vector<string>& names, vector<vector<double> >& D){
	try {
		ifstream f(fname.c_str());
		if(!f) {
			cerr << "Error: Could not open " << fname << endl;
			exit(1);
		}
		char test = f.peek();
		
		if(test == '#'){
			read_mega(f, m, names, D);
		}
		else{
			read_phylip(f, m, names, D);
		}
		
		int rank = D.size();
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "read");
		exit(1);
	}
}

/*********************************************************************************************************************************/
double PCACommand::pythag(double a, double b){
    return(pow(a*a+b*b,0.5));
}
/*********************************************************************************************************************************/

void PCACommand::matrix_mult(vector<vector<double> > first, vector<vector<double> > second, vector<vector<double> >& product){
	try {
		int first_rows = first.size();
		int first_cols = first[0].size();
		int second_cols = second[0].size();
		
		product.resize(first_rows);
		for(int i=0;i<first_rows;i++){
			product[i].resize(second_cols);
		}
		
		for(int i=0;i<first_rows;i++){
			for(int j=0;j<second_cols;j++){
				product[i][j] = 0.0;
				for(int k=0;k<first_cols;k++){
					product[i][j] += first[i][k] * second[k][j];
				}
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "matrix_mult");
		exit(1);
	}

}

/*********************************************************************************************************************************/

void PCACommand::recenter(double offset, vector<vector<double> > D, vector<vector<double> >& G){
	try {
		int rank = D.size();
		
		vector<vector<double> > A(rank);
		vector<vector<double> > C(rank);
		for(int i=0;i<rank;i++){
			A[i].resize(rank);
			C[i].resize(rank);
		}
		
		double scale = -1.0000 / (double) rank;
		
		for(int i=0;i<rank;i++){
			A[i][i] = 0.0000;
			C[i][i] = 1.0000 + scale;
			for(int j=i+1;j<rank;j++){
				A[i][j] = A[j][i] = -0.5 * D[i][j] * D[i][j] + offset;
				C[i][j] = C[j][i] = scale;
			}
		}
		
		matrix_mult(C,A,A);
		matrix_mult(A,C,G);
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "recenter");
		exit(1);
	}

}

/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

void PCACommand::tred2(vector<vector<double> >& a, vector<double>& d, vector<double>& e){
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
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "tred2");
		exit(1);
	}

}

/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

void PCACommand::qtli(vector<double>& d, vector<double>& e, vector<vector<double> >& z) {
	try {
		int m, i, iter;
		double s, r, p, g, f, dd, c, b;
		
		int n = d.size();
		for(int i=1;i<=n;i++){
			e[i-1] = e[i];
		}
		e[n-1] = 0.0000;
		
		for(int l=0;l<n;l++){
			iter = 0;
			do {
				for(m=l;m<n-1;m++){
					dd = fabs(d[m]) + fabs(d[m+1]);
					if(fabs(e[m])+dd == dd) break;
				}
				if(m != l){
					if(iter++ == 30) cerr << "Too many iterations in tqli\n";
					g = (d[l+1]-d[l]) / (2.0 * e[l]);
					r = pythag(g, 1.0);
					g = d[m] - d[l] + e[l] / (g + SIGN(r,g));
					s = c = 1.0;
					p = 0.0000;
					for(i=m-1;i>=l;i--){
						f = s * e[i];
						b = c * e[i];
						e[i+1] = (r=pythag(f,g));
						if(r==0.0){
							d[i+1] -= p;
							e[m] = 0.0000;
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
					e[m] = 0.0;
				}
			} while (m != l);
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
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "qtli");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void PCACommand::output(string fnameRoot, vector<string> name_list, vector<vector<double> > G, vector<double> d) {
	try {
		int rank = name_list.size();
		double dsum = 0.0000;
		for(int i=0;i<rank;i++){
			dsum += d[i];
			for(int j=0;j<rank;j++){
				if(d[j] >= 0)	{	G[i][j] *= pow(d[j],0.5);	}
				else			{	G[i][j] = 0.00000;			}
			}
		}
		
		ofstream pcaData((fnameRoot+"pca").c_str(), ios::trunc);
		pcaData.setf(ios::fixed, ios::floatfield);
		pcaData.setf(ios::showpoint);	
		
		ofstream pcaLoadings((fnameRoot+"pca.loadings").c_str(), ios::trunc);
		pcaLoadings.setf(ios::fixed, ios::floatfield);
		pcaLoadings.setf(ios::showpoint);	
		
		pcaLoadings << "axis\tloading\n";
		for(int i=0;i<rank;i++){
			pcaLoadings << i+1 << '\t' << d[i] * 100.0 / dsum << endl;
		}
		
		pcaData << "SeqName";
		for(int i=0;i<rank;i++){
			pcaData << '\t' << "axis" << i+1;
		}
		pcaData << endl;
		
		for(int i=0;i<rank;i++){
			pcaData << name_list[i] << '\t';
			for(int j=0;j<rank;j++){
				pcaData << G[i][j] << '\t';
			}
			pcaData << endl;
		}
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "output");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void PCACommand::print_matrix(vector<vector<double> > A) {
	try {
		int rank = A.size();
		for(int i=0;i<rank;i++){
			for(int j=0;j<rank;j++){
				cout << A[i][j] << "  ";
			}
			cout << endl;
		}
	}
	catch(exception& e) {
		errorOut(e, "PCACommand", "print_matrix");
		exit(1);
	}
}
/*********************************************************************************************************************************/

