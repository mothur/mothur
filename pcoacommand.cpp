
/*
 *  pcacommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "pcoacommand.h"

//**********************************************************************************************************************
vector<string> PCOACommand::getValidParameters(){	
	try {
		string Array[] =  {"phylip", "metric","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
PCOACommand::PCOACommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["pcoa"] = tempOutNames;
		outputTypes["loadings"] = tempOutNames;
		outputTypes["corr"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "PCOACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PCOACommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PCOACommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

PCOACommand::PCOACommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip","metric","outputdir", "inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["pcoa"] = tempOutNames;
			outputTypes["loadings"] = tempOutNames;
			outputTypes["corr"] = tempOutNames;
			
			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; abort = true; }	
			else {	filename = phylipfile;  }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(phylipfile); //if user entered a file with a path then preserve it	
			}
			
			//error checking on files	
			if (phylipfile == "")	{ m->mothurOut("You must provide a distance file before running the pcoa command."); m->mothurOutEndLine(); abort = true; }		
		
			string temp = validParameter.validFile(parameters, "metric", false);	if (temp == "not found"){	temp = "T";				}
			metric = m->isTrue(temp); 
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "PCOACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void PCOACommand::help(){
	try {
	
		m->mothurOut("The pcoa command parameters are phylip and metric"); m->mothurOutEndLine();
		m->mothurOut("The phylip parameter allows you to enter your distance file."); m->mothurOutEndLine();
		m->mothurOut("The metric parameter allows indicate you if would like the pearson correlation coefficient calculated. Default=True"); m->mothurOutEndLine();
		m->mothurOut("Example pcoa(phylip=yourDistanceFile).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourDistanceFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
PCOACommand::~PCOACommand(){}
//**********************************************************************************************************************
int PCOACommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);
		
		vector<string> names;
		vector<vector<double> > D;
	
		fbase = outputDir + m->getRootName(m->getSimpleName(filename));
		
		read(filename, names, D);
		
		if (m->control_pressed) { return 0; }
   	
		double offset = 0.0000;
		vector<double> d;
		vector<double> e;
		vector<vector<double> > G = D;
		vector<vector<double> > copy_G;
		//int rank = D.size();
		
		m->mothurOut("\nProcessing...\n\n");
		
		for(int count=0;count<2;count++){
			recenter(offset, D, G);		if (m->control_pressed) { return 0; }
			tred2(G, d, e);				if (m->control_pressed) { return 0; }
			qtli(d, e, G);				if (m->control_pressed) { return 0; }
			offset = d[d.size()-1];
			if(offset > 0.0) break;
		} 
		
		if (m->control_pressed) { return 0; }
		
		output(fbase, names, G, d);
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		
		if (metric) {   
			
			for (int i = 1; i < 4; i++) {
							
				vector< vector<double> > EuclidDists = calculateEuclidianDistance(G, i); //G is the pcoa file
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
				
				double corr = calcPearson(EuclidDists, D); //G is the pcoa file, D is the users distance matrix
				
				m->mothurOut("Pearson's coefficient using " + toString(i) + " axis: " + toString(corr)); m->mothurOutEndLine();
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			}
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "execute");
		exit(1);
	}
}
/*********************************************************************************************************************************/
vector< vector<double> > PCOACommand::calculateEuclidianDistance(vector< vector<double> >& axes, int dimensions){
	try {
		//make square matrix
		vector< vector<double> > dists; dists.resize(axes.size());
		for (int i = 0; i < dists.size(); i++) {  dists[i].resize(axes.size(), 0.0); }
			
		if (dimensions == 1) { //one dimension calc = abs(x-y)
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					dists[i][j] = abs(axes[i][0] - axes[j][0]);
					dists[j][i] = dists[i][j];
				}
			}
			
		}else if (dimensions == 2) { //two dimension calc = sqrt ((x1 - y1)^2 + (x2 - y2)^2)
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					double firstDim = ((axes[i][0] - axes[j][0]) * (axes[i][0] - axes[j][0]));
					double secondDim = ((axes[i][1] - axes[j][1]) * (axes[i][1] - axes[j][1]));

					dists[i][j] = sqrt((firstDim + secondDim));
					dists[j][i] = dists[i][j];
				}
			}
			
		}else if (dimensions == 3) { //two dimension calc = sqrt ((x1 - y1)^2 + (x2 - y2)^2 + (x3 - y3)^2)
			
			for (int i = 0; i < dists.size(); i++) {
				
				if (m->control_pressed) { return dists; }
				
				for (int j = 0; j < i; j++) {
					double firstDim = ((axes[i][0] - axes[j][0]) * (axes[i][0] - axes[j][0]));
					double secondDim = ((axes[i][1] - axes[j][1]) * (axes[i][1] - axes[j][1]));
					double thirdDim = ((axes[i][2] - axes[j][2]) * (axes[i][2] - axes[j][2]));
					
					dists[i][j] = sqrt((firstDim + secondDim + thirdDim));
					dists[j][i] = dists[i][j];
				}
			}
			
		}else { m->mothurOut("[ERROR]: too many dimensions, aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "calculateEuclidianDistance");
		exit(1);
	}
}
/*********************************************************************************************************************************/
double PCOACommand::calcPearson(vector< vector<double> >& euclidDists, vector< vector<double> >& userDists){
	try {
		
		//find average for - X
		vector<float> averageEuclid; averageEuclid.resize(euclidDists.size(), 0.0);
		for (int i = 0; i < euclidDists.size(); i++) {
			for (int j = 0; j < euclidDists[i].size(); j++) {
				averageEuclid[i] += euclidDists[i][j];  
			}
		}
		for (int i = 0; i < averageEuclid.size(); i++) {  averageEuclid[i] = averageEuclid[i] / (float) euclidDists.size();   }
		
		//find average for - Y
		vector<float> averageUser; averageUser.resize(userDists.size(), 0.0);
		for (int i = 0; i < userDists.size(); i++) {
			for (int j = 0; j < userDists[i].size(); j++) {
				averageUser[i] += userDists[i][j];  
			}
		}
		for (int i = 0; i < averageUser.size(); i++) {  averageUser[i] = averageUser[i] / (float) userDists.size();  }
		
		double numerator = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < euclidDists.size(); i++) {
			
			for (int k = 0; k < i; k++) {
				
				float Yi = userDists[i][k];
				float Xi = euclidDists[i][k];
					
				numerator += ((Xi - averageEuclid[k]) * (Yi - averageUser[k]));
				denomTerm1 += ((Xi - averageEuclid[k]) * (Xi - averageEuclid[k]));
				denomTerm2 += ((Yi - averageUser[k]) * (Yi - averageUser[k]));
			}
		}
		
		double denom = (sqrt(denomTerm1) * sqrt(denomTerm2));
		double r = numerator / denom;
		
		return r;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "calculateEuclidianDistance");
		exit(1);
	}
}
/*********************************************************************************************************************************/

inline double SIGN(const double a, const double b)
{
    return b>=0 ? (a>=0 ? a:-a) : (a>=0 ? -a:a);
}

/*********************************************************************************************************************************/

void PCOACommand::get_comment(istream& f, char begin, char end){
	try {
		char d=f.get();
		while(d != end){	d = f.get();	}
		d = f.peek();
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "get_comment");
		exit(1);
	}
}	

/*********************************************************************************************************************************/

int PCOACommand::read_phylip(istream& f, int square_m, vector<string>& name_list, vector<vector<double> >& d){
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
					if (m->control_pressed) { return 0; }
					
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
					if (m->control_pressed) { return 0; }
					f >> d[i][j];
					if (d[i][j] == -0.0000)
						d[i][j] = 0.0000;
					d[j][i]=d[i][j];
				}
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "read_phylip");
		exit(1);
	}

}

/*********************************************************************************************************************************/

void PCOACommand::read(string fname, vector<string>& names, vector<vector<double> >& D){
	try {
		ifstream f;
		m->openInputFile(fname, f);
			
		//check whether matrix is square
		char d;
		int q = 1;
		int numSeqs;
		string name;
		
		f >> numSeqs >> name; 
		
		while((d=f.get()) != EOF){
			
			//is d a number meaning its square
			if(isalnum(d)){ 
				q = 1; 
				break; 
			}
			
			//is d a line return meaning its lower triangle
			if(d == '\n'){
				q = 2;
				break;
			}
		}
		f.close();
		
		//reopen to get back to beginning
		m->openInputFile(fname, f);			
		read_phylip(f, q, names, D);
	}
		catch(exception& e) {
		m->errorOut(e, "PCOACommand", "read");
		exit(1);
	}
}

/*********************************************************************************************************************************/

double PCOACommand::pythag(double a, double b)	{	return(pow(a*a+b*b,0.5));	}

/*********************************************************************************************************************************/

void PCOACommand::matrix_mult(vector<vector<double> > first, vector<vector<double> > second, vector<vector<double> >& product){
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
		m->errorOut(e, "PCOACommand", "matrix_mult");
		exit(1);
	}

}

/*********************************************************************************************************************************/

void PCOACommand::recenter(double offset, vector<vector<double> > D, vector<vector<double> >& G){
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
		m->errorOut(e, "PCOACommand", "recenter");
		exit(1);
	}

}

/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

void PCOACommand::tred2(vector<vector<double> >& a, vector<double>& d, vector<double>& e){
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
		m->errorOut(e, "PCOACommand", "tred2");
		exit(1);
	}

}

/*********************************************************************************************************************************/

//  This function is taken from Numerical Recipes in C++ by Press et al., 2nd edition, pg. 479

void PCOACommand::qtli(vector<double>& d, vector<double>& e, vector<vector<double> >& z) {
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
		m->errorOut(e, "PCOACommand", "qtli");
		exit(1);
	}
}

/*********************************************************************************************************************************/

void PCOACommand::output(string fnameRoot, vector<string> name_list, vector<vector<double> >& G, vector<double> d) {
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
		
		ofstream pcaData((fnameRoot+"pcoa").c_str(), ios::trunc);
		pcaData.setf(ios::fixed, ios::floatfield);
		pcaData.setf(ios::showpoint);	
		outputNames.push_back(fnameRoot+"pcoa");
		outputTypes["pcoa"].push_back(fnameRoot+"pcoa");
		
		ofstream pcaLoadings((fnameRoot+"pcoa.loadings").c_str(), ios::trunc);
		pcaLoadings.setf(ios::fixed, ios::floatfield);
		pcaLoadings.setf(ios::showpoint);
		outputNames.push_back(fnameRoot+"pcoa.loadings");
		outputTypes["loadings"].push_back(fnameRoot+"pcoa.loadings");	
		
		pcaLoadings << "axis\tloading\n";
		for(int i=0;i<rank;i++){
			pcaLoadings << i+1 << '\t' << d[i] * 100.0 / dsum << endl;
		}
		
		pcaData << "group";
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
		m->errorOut(e, "PCOACommand", "output");
		exit(1);
	}
}

/*********************************************************************************************************************************/

