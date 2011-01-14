/*
 *  nmdscommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/11/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "nmdscommand.h"
#include "readphylipvector.h"

//**********************************************************************************************************************
vector<string> NMDSCommand::getValidParameters(){	
	try {
		string Array[] =  {"phylip","axes","dimension","maxiters","step","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
NMDSCommand::NMDSCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["nmds"] = tempOutNames;
		outputTypes["stress"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "NMDSCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> NMDSCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> NMDSCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

NMDSCommand::NMDSCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip","axes","dimension","maxiters","step","outputdir", "inputdir"};
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
				
				it = parameters.find("axes");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["axes"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["nmds"] = tempOutNames;
			outputTypes["stress"] = tempOutNames;
			
			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; m->mothurOut("You must provide a distance file before running the nmds command."); m->mothurOutEndLine(); abort = true; }	
			
			axesfile = validParameter.validFile(parameters, "axes", true);
			if (axesfile == "not open") { axesfile = ""; abort = true; }
			else if (axesfile == "not found") { axesfile = "";  }				
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(phylipfile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "dimension", false);	if (temp == "not found") {	temp = "2";	}
			convert(temp, dimension);
			
			temp = validParameter.validFile(parameters, "maxiters", false);	if (temp == "not found") {	temp = "1000";	}
			convert(temp, maxIters);
			
			temp = validParameter.validFile(parameters, "step", false);	if (temp == "not found") {	temp = "0.2";	}
			convert(temp, step);
			
			temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found") {	temp = "2";	}
			convert(temp, cutoff); 
			cutoff /= 100.0;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "NMDSCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void NMDSCommand::help(){
	try {
		
		m->mothurOut("The nmds command parameters are phylip, axes, dimension, maxiters, cutoff and step."); m->mothurOutEndLine();
		m->mothurOut("The phylip parameter allows you to enter your distance file."); m->mothurOutEndLine();
		m->mothurOut("The axes parameter allows you to enter a file containing a starting configuration."); m->mothurOutEndLine();
		m->mothurOut("The dimension parameter allows you to select how many dimensions to use. Default=2"); m->mothurOutEndLine();
		m->mothurOut("The maxiters parameter allows you to select the maximum number of iters to try. Default=1000"); m->mothurOutEndLine();
		m->mothurOut("The cutoff parameter allows you to select set an acceptable percentage of magnitude. Default=2, meaning when magnitude of g reaches 2% of it's starting value the process will stop."); m->mothurOutEndLine();
		m->mothurOut("The step parameter allows you to set a starting step. Default=0.2"); m->mothurOutEndLine();
		m->mothurOut("Example nmds(phylip=yourDistanceFile).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourDistanceFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
NMDSCommand::~NMDSCommand(){}
//**********************************************************************************************************************
int NMDSCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);
		
		vector<string> names;
		vector<seqDist> matrix; //seqDist = int, int, float - index of seq1 in names, index of seq2 in names, their distance
		
		//read in phylip file
		ReadPhylipVector readFile(phylipfile);
		names = readFile.read(matrix);
		if (m->control_pressed) { return 0; }
	
		//randomly generate the starting configuration - step 2
		vector< vector<double> > axes;
		if (axesfile == "") {	axes = generateStartingConfiguration(names.size());		}
		else				{	axes = readAxes(names);									}
		if (m->control_pressed) { return 0; }
		
		//sort matrix from smallest distance to largest - step 5
		sort(matrix.begin(), matrix.end(), compareSequenceDistance);
		
		bool stable = false;
		int count = 0;
		vector<double> previousStresses;
		vector< vector<double> > previousGradient = axes;
		double initialMagnitude;
		m->mothurOutEndLine(); m->mothurOut("Iter\tStress\tMagnitude"); m->mothurOutEndLine();
		while ((count != maxIters) && (!stable)) {
			count++;
			
			//normalize axes - step 3
			normalizeConfiguration(axes, names.size());
			if (m->control_pressed) { return 0; }
			
			//calculate Euclidean distances - step 4
			vector< vector<double> > euclid = linearCalc.calculateEuclidianDistance(axes);
			if (m->control_pressed) { return 0; }
			
			//order euclid elements in same order as matrix - step 6
			//if there are ties in the matrix we want to arrange the euclid distances in the best way so we do not to add unnecessary stress
			vector<seqDist> eDists;
			vector<seqDist> ties;
			for (int i = 0; i < matrix.size(); i++) {
				
				seqDist temp(matrix[i].seq1, matrix[i].seq2, euclid[matrix[i].seq1][matrix[i].seq2]);
				ties.push_back(temp);
				
				if (i != matrix.size()-1) { // you are not the last so you can look ahead
					if (matrix[i].dist != matrix[i+1].dist) { // you are done with ties, sort and save them, then continue
						sort(ties.begin(), ties.end(), compareSequenceDistance);
						for (int k = 0; k < ties.size(); k++) {	eDists.push_back(ties[k]);	}
						ties.clear();
					}
				}else { // you are the last one
					sort(ties.begin(), ties.end(), compareSequenceDistance);
					for (int k = 0; k < ties.size(); k++) {	eDists.push_back(ties[k]);	}
				}
			}
			
			for (int i = 0; i < euclid.size(); i++) {  euclid[i].clear(); } euclid.clear();
			if (m->control_pressed) { return 0; }
			
			//find D - from step 7
			vector<seqDist> D = satisfyMonotonicity(eDists);
			if (m->control_pressed) { return 0; }
			
			//calculate the raw stress and normalize it - steps 8 and 9
			double rawStress;
			double stress = calculateStress(eDists, D, rawStress);
			previousStresses.push_back(stress);
			if (stress == 0) { m->mothurOut("Stress reached zero after " + toString(count) + " iters, stopping."); m->mothurOutEndLine(); break; }
			if (m->control_pressed) { return 0; }
			
			//calculate stress gradient - step 10
			vector< vector<double> > stressGradient = calculateStressGradientVector(eDists, D, rawStress, stress, axes);
			if (m->control_pressed) { return 0; }
			
			//calculate magnitude
			double magnitude = calculateMagnitude(stressGradient);	
			if (count == 1) { initialMagnitude = magnitude; }
			if (m->control_pressed) { return 0; }
			
			//save gradient before adjusting config.
			previousGradient = stressGradient;
			
			if ((count % 100) == 0) { m->mothurOut(toString(count) + "\t" + toString(previousStresses[previousStresses.size()-1]) + "\t" + toString(magnitude)); m->mothurOutEndLine(); }

			//are we done - we are done if percentage of magnitude compared to initial magnitude is less than cutoff
			double percentage = magnitude / initialMagnitude;
			if (percentage < cutoff) { stable = true; }
			else {
			
				//calculate new step size
				step = calculateStep(previousGradient, stressGradient, previousStresses);
				cout << "count = " << count << '\t' << step << endl;
				if (m->control_pressed) { return 0; }
			
				//find new config.
				axes = calculateNewConfiguration(magnitude, axes, stressGradient);
				if (m->control_pressed) { return 0; }
			}
		}
		
		if (m->control_pressed) { return 0; }
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(phylipfile)) + "nmds";
		string stressFileName = outputDir + m->getRootName(m->getSimpleName(phylipfile)) + "stress.nmds";
		outputNames.push_back(outputFileName); outputTypes["nmds"].push_back(outputFileName);
		outputNames.push_back(stressFileName); outputTypes["stress"].push_back(stressFileName);
		
		output(outputFileName, stressFileName, previousGradient, previousStresses, names);
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
//generate random config
vector< vector<double> > NMDSCommand::generateStartingConfiguration(int numNames) {
	try {
		vector< vector<double> > axes;  axes.resize(dimension);
		for (int i = 0; i < axes.size(); i++) {  axes[i].resize(numNames); }
		
		//generate random number between -1 and 1, precision 6
		for (int i = 0; i < axes.size(); i++) {
			for (int j = 0; j < axes[i].size(); j++) {
				
				if (m->control_pressed) { return axes; }
				
				//generate random int between 0 and 99999
				int myrand = (int)((float)(rand()) / ((RAND_MAX / 99998) + 1));
				
				//generate random sign
				int mysign = (int)((float)(rand()) / ((RAND_MAX / 99998) + 1));
				
				//if mysign is even then sign = positive, else sign = negative
				if ((mysign % 2) == 0) { mysign = 1.0; }
				else { mysign = -1.0; }
				
				axes[i][j] = mysign * myrand / (float) 100000;
			}
		}

		return axes;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "generateStartingConfiguration");
		exit(1);
	}
}
//**********************************************************************************************************************
//normalize configuration
int NMDSCommand::normalizeConfiguration(vector< vector<double> >& axes, int numNames) {
	try {
		vector<double> averageAxes; averageAxes.resize(dimension, 0.0);
		
		//find average
		for (int i = 0; i < axes.size(); i++) {
			for (int j = 0; j < axes[i].size(); j++) {	averageAxes[i] += axes[i][j];	}
			
			averageAxes[i] /= (float) numNames;
		}
		
		//normalize axes
		double sumDenom = 0.0;
		for (int i = 0; i < axes.size(); i++) {
			for (int j = 0; j < axes[i].size(); j++) {
				sumDenom += ((axes[i][j] - averageAxes[i]) * (axes[i][j] - averageAxes[i]));
			}
		}
		
		double denom = sqrt((sumDenom / (float) (axes.size() * numNames)));
		
		for (int i = 0; i < axes.size(); i++) {
			for (int j = 0; j < axes[i].size(); j++) {
				axes[i][j] = (axes[i][j] - averageAxes[i]) / denom;
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "normalizeConfiguration");
		exit(1);
	}
}
//**********************************************************************************************************************
//adjust eDists so that it creates monotonically increasing series of succesive values that increase or stay the same, but never decrease
vector<seqDist> NMDSCommand::satisfyMonotonicity(vector<seqDist> eDists) {
	try {
		
		vector<seqDist> D = eDists; 
		
		for (int i = 0; i < (D.size()-1); i++) {
			
			if (m->control_pressed) { return D; }
			
			//is the distance in i+1 smaller than i, if yes then adjust
			if (D[i+1].dist < D[i].dist) {  D[i+1].dist = D[i].dist;  }
		}
		
		return D;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "satisfyMonotonicity");
		exit(1);
	}
}
//**********************************************************************************************************************
//find raw stress, and normalize using
double NMDSCommand::calculateStress(vector<seqDist>& eDists, vector<seqDist>& D, double& rawStress) {
	try {
		double normStress = 0.0;
		double denom = 0.0;
		rawStress = 0.0;
		
		//find raw stress
		for (int i = 0; i < D.size(); i++) {
			
			if (m->control_pressed) { return normStress; }
			
			rawStress += ((eDists[i].dist - D[i].dist) * (eDists[i].dist - D[i].dist));
			denom += (eDists[i].dist * eDists[i].dist);
		}
		
		//normalize stress
		if (rawStress != 0.0) {
			normStress = 100 * sqrt((rawStress / denom));
		}
		
		return normStress;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "calculateStress");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<double> > NMDSCommand::calculateStressGradientVector(vector<seqDist>& eDists, vector<seqDist>& D, double rawStress, double stress, vector< vector<double> >& axes) {
	try {
		vector< vector<double> > gradient; gradient.resize(dimension);
		for (int i = 0; i < gradient.size(); i++) { gradient[i].resize(axes[0].size(), 0.0); }
	
		double sumDij = 0.0;
		for (int i = 0; i < eDists.size(); i++) {  sumDij += (eDists[i].dist * eDists[i].dist); }
		
		for (int i = 0; i < eDists.size(); i++) {
			
			for (int j = 0; j < dimension; j++) {
			
				if (m->control_pressed) { return gradient; }
				
				double firstTerm1 = (stress / rawStress) * (eDists[i].dist - D[i].dist);
				double firstTerm2 = eDists[i].dist * (stress / sumDij);
				double firstTerm = firstTerm1 - firstTerm2;
				
				double secondTerm = (axes[j][eDists[i].seq1] - axes[j][eDists[i].seq2]) / eDists[i].dist; 
				
				double results = (firstTerm * secondTerm);
				
				gradient[j][eDists[i].seq1] += results;
				gradient[j][eDists[i].seq2] -= results;
			}
		}
		
		return gradient;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "calculateStressGradientVector");
		exit(1);
	}
}
//**********************************************************************************************************************
double NMDSCommand::calculateMagnitude(vector< vector<double> >& gradient) {
	try {
		double magnitude = 0.0;
		
		double sum = 0.0;
		for (int i = 0; i < gradient.size(); i++) {
			for (int j = 0; j < gradient[i].size(); j++) {
				sum += (gradient[i][j] * gradient[i][j]);
			}
		}
		
		magnitude = sqrt(((1.0/(float)gradient[0].size()) * sum));
		
		return magnitude;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "calculateMagnitude");
		exit(1);
	}
}
//**********************************************************************************************************************
//described in Kruskal paper page 121 + 122
double NMDSCommand::calculateStep(vector< vector<double> >& prevGrad, vector< vector<double> >& grad, vector<double>& prevStress) {
	try {
		double newStep = step;
		
		//calc the cos theta
		double sumNum = 0.0;
		double sumDenom1 = 0.0;
		double sumDenom2 = 0.0;
		for (int i = 0; i < prevGrad.size(); i++) {
			for (int j = 0; j < prevGrad[i].size(); j++) {
				sumDenom1 += (grad[i][j] * grad[i][j]);
				sumDenom2 += (prevGrad[i][j] * prevGrad[i][j]);
				sumNum += (grad[i][j] * prevGrad[i][j]);
			}
		}
		
		double cosTheta = sumNum / (sqrt(sumDenom1) * sqrt(sumDenom2));
		cosTheta *= cosTheta;
	
		//calc angle factor
		double angle = pow(4.0, cosTheta);
	
		//calc 5 step ratio
		double currentStress = prevStress[prevStress.size()-1];
		double lastStress = prevStress[0];
		if (prevStress.size() > 1) {  lastStress = prevStress[prevStress.size()-2];		}
		double fivePrevStress = prevStress[0];
		if (prevStress.size() > 5) {  fivePrevStress = prevStress[prevStress.size()-6]; }
			
		double fiveStepRatio = min(1.0, (currentStress / fivePrevStress));
		
		//calc relaxation factor
		double relaxation = 1.3 / (1.0 + pow(fiveStepRatio, 5.0));
		
		//calc good luck factor
		double goodLuck = min(1.0, (currentStress / lastStress));
		
		//calc newStep
		cout << "\ncos = " << cosTheta << " step = " << step << " angle = " << angle << " relaxation = " << relaxation << " goodluck = " << goodLuck << endl;
		newStep = step * angle * relaxation * goodLuck;
		
		return newStep;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "calculateStep");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<double> > NMDSCommand::calculateNewConfiguration(double magnitude, vector< vector<double> >& axes, vector< vector<double> >& gradient) {
	try {
		
		vector< vector<double> > newAxes = axes;
		
		for (int i = 0; i < newAxes.size(); i++) {
			
			if (m->control_pressed) { return newAxes; }
			
			for (int j = 0; j < newAxes[i].size(); j++) {
				newAxes[i][j] = axes[i][j] + ((step / magnitude) * gradient[i][j]);
			}
		}
		
		return newAxes;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "calculateNewConfiguration");
		exit(1);
	}
}
//**********************************************************************************************************************
int NMDSCommand::output(string outputFileName, string stressFileName, vector< vector<double> >& config, vector<double>& stresses, vector<string>& names) {
	try {
		
		ofstream out, out2;
		m->openOutputFile(outputFileName, out);
		m->openOutputFile(stressFileName, out2);
		
		//output headers
		out << "group\t";
		for (int i = 0; i < dimension; i++) { out << "axis" << (i+1) << '\t'; }
		out << endl;
		
		out2 << "Iter\tStress" << endl;
		
		//output nmds file
		for (int i = 0; i < config[0].size(); i++) {
			
			if (m->control_pressed) { out.close(); out2.close(); return 0; }
			
			out << names[i] << '\t';
			
			for (int j = 0; j < config.size(); j++) {
				out << config[j][i] << '\t';
			}
			
			out << endl;
		}
		out.close();
		
		//output stress file
		for (int j = 0; j < stresses.size(); j++) {
			if (m->control_pressed) { out2.close(); return 0; }
			
			out2 << (j+1) << '\t' << stresses[j] << endl;
		}
		out2.close();
		
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "output");
		exit(1);
	}
}
/*****************************************************************/
vector< vector<double> > NMDSCommand::readAxes(vector<string> names){
	try {
		vector< vector<double> > axes;  
		
		ifstream in;
		m->openInputFile(axesfile, in);
		
		string headerLine = m->getline(in); m->gobble(in);
		
		//count the number of axis you are reading
		bool done = false;
		int count = 0;
		while (!done) {
			int pos = headerLine.find("axis");
			if (pos != string::npos) {
				count++;
				headerLine = headerLine.substr(pos+4);
			}else { done = true; }
		}
		
		if (dimension > count) { m->mothurOut("You requested " + toString(dimension) + " axes, but your file only includes " + toString(count) + ". Using " + toString(count) + "."); m->mothurOutEndLine(); dimension = count; }
		
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); return axes; }
			
			string group = "";
			in >> group; m->gobble(in);
			
			bool ignore = false;
			if (!m->inUsersGroups(group, names)) { ignore = true; m->mothurOut(group + " is in your axes file and not in your distance file, ignoring."); m->mothurOutEndLine(); }
			
			vector<double> thisGroupsAxes;
			for (int i = 0; i < count; i++) {
				float temp = 0.0;
				in >> temp; 
				
				//only save the axis we want
				if (i < dimension) {  thisGroupsAxes.push_back(temp); }
			}
			
			if (!ignore) {	axes.push_back(thisGroupsAxes); }
			
			m->gobble(in);
		}
		in.close();
		
		//sanity check
		if (names.size() != axes.size()) { m->mothurOut("[ERROR]: your axes file does not match your distance file, aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		return axes;
	}
	catch(exception& e) {
		m->errorOut(e, "NMDSCommand", "readAxes");	
		exit(1);
	}
}
//**********************************************************************************************************************


