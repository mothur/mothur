/*
 *  matrixoutputcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "matrixoutputcommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "sharedjackknife.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"

//**********************************************************************************************************************
vector<string> MatrixOutputCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pcalc("calc", "Multiple", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-whittaker-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-hamming-structchi2-gower-memchi2-memchord-memeuclidean-mempearson", "jclass-thetayc", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter poutput("output", "Multiple", "lt-square", "lt", "", "", "",false,false); parameters.push_back(poutput);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MatrixOutputCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The dist.shared command parameters are shared, groups, calc, output, processors and label.  shared is a required, unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n";
		helpString += "The group names are separated by dashes. The label parameter allows you to select what distance levels you would like distance matrices created for, and is also separated by dashes.\n";
		helpString += "The dist.shared command should be in the following format: dist.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels).\n";
		helpString += "The output parameter allows you to specify format of your distance matrix. Options are lt, and square. The default is lt.\n";
		helpString += "Example dist.shared(groups=A-B-C, calc=jabund-sorabund).\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The default value for calc is jclass and thetayc.\n";
		helpString += validCalculator.printCalc("matrix");
		helpString += "The dist.shared command outputs a .dist file for each calculator you specify at each distance you choose.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
MatrixOutputCommand::MatrixOutputCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["phylip"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "MatrixOutputCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MatrixOutputCommand::MatrixOutputCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters  = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["phylip"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not found") { 			
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else { m->setSharedFile(sharedfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(sharedfile); //if user entered a file with a path then preserve it	
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "lt"; }
			if ((output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are lt and square. I will use lt."); m->mothurOutEndLine(); output = "lt"; }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}

			if (abort == false) {
			
				ValidCalculators validCalculator;
				
				int i;
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("matrix", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedsobs") { 
							matrixCalculators.push_back(new SharedSobsCS());
						}else if (Estimators[i] == "sharedchao") { 
							matrixCalculators.push_back(new SharedChao1());
						}else if (Estimators[i] == "sharedace") { 
							matrixCalculators.push_back(new SharedAce());
						}else if (Estimators[i] == "jabund") { 	
							matrixCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							matrixCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							matrixCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							matrixCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							matrixCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							matrixCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							matrixCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							matrixCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "kstest") { 
							matrixCalculators.push_back(new KSTest());
						}else if (Estimators[i] == "sharednseqs") { 
							matrixCalculators.push_back(new SharedNSeqs());
						}else if (Estimators[i] == "ochiai") { 
							matrixCalculators.push_back(new Ochiai());
						}else if (Estimators[i] == "anderberg") { 
							matrixCalculators.push_back(new Anderberg());
						}else if (Estimators[i] == "kulczynski") { 
							matrixCalculators.push_back(new Kulczynski());
						}else if (Estimators[i] == "kulczynskicody") { 
							matrixCalculators.push_back(new KulczynskiCody());
						}else if (Estimators[i] == "lennon") { 
							matrixCalculators.push_back(new Lennon());
						}else if (Estimators[i] == "morisitahorn") { 
							matrixCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							matrixCalculators.push_back(new BrayCurtis());
						}else if (Estimators[i] == "whittaker") { 
							matrixCalculators.push_back(new Whittaker());
						}else if (Estimators[i] == "odum") { 
							matrixCalculators.push_back(new Odum());
						}else if (Estimators[i] == "canberra") { 
							matrixCalculators.push_back(new Canberra());
						}else if (Estimators[i] == "structeuclidean") { 
							matrixCalculators.push_back(new StructEuclidean());
						}else if (Estimators[i] == "structchord") { 
							matrixCalculators.push_back(new StructChord());
						}else if (Estimators[i] == "hellinger") { 
							matrixCalculators.push_back(new Hellinger());
						}else if (Estimators[i] == "manhattan") { 
							matrixCalculators.push_back(new Manhattan());
						}else if (Estimators[i] == "structpearson") { 
							matrixCalculators.push_back(new StructPearson());
						}else if (Estimators[i] == "soergel") { 
							matrixCalculators.push_back(new Soergel());
						}else if (Estimators[i] == "spearman") { 
							matrixCalculators.push_back(new Spearman());
						}else if (Estimators[i] == "structkulczynski") { 
							matrixCalculators.push_back(new StructKulczynski());
						}else if (Estimators[i] == "speciesprofile") { 
							matrixCalculators.push_back(new SpeciesProfile());
						}else if (Estimators[i] == "hamming") { 
							matrixCalculators.push_back(new Hamming());
						}else if (Estimators[i] == "structchi2") { 
							matrixCalculators.push_back(new StructChi2());
						}else if (Estimators[i] == "gower") { 
							matrixCalculators.push_back(new Gower());
						}else if (Estimators[i] == "memchi2") { 
							matrixCalculators.push_back(new MemChi2());
						}else if (Estimators[i] == "memchord") { 
							matrixCalculators.push_back(new MemChord());
						}else if (Estimators[i] == "memeuclidean") { 
							matrixCalculators.push_back(new MemEuclidean());
						}else if (Estimators[i] == "mempearson") { 
							matrixCalculators.push_back(new MemPearson());
						}
					}
				}
				
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "MatrixOutputCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

MatrixOutputCommand::~MatrixOutputCommand(){}

//**********************************************************************************************************************

int MatrixOutputCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
			
		//if the users entered no valid calculators don't execute command
		if (matrixCalculators.size() == 0) { m->mothurOut("No valid calculators."); m->mothurOutEndLine();  return 0; }
			
		input = new InputData(sharedfile, "sharedfile");
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
					
		if (lookup.size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command."); m->mothurOutEndLine(); delete input; for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } return 0;}
		
		numGroups = lookup.size();
		lines.resize(processors);
		for (int i = 0; i < processors; i++) {
			lines[i].start = int (sqrt(float(i)/float(processors)) * numGroups);
			lines[i].end = int (sqrt(float(i+1)/float(processors)) * numGroups);
		}	
		
		if (m->control_pressed) { delete input; for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } m->clearGroups(); return 0;  }
				
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
		
			if (m->control_pressed) { outputTypes.clear(); delete input; for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } m->clearGroups(); return 0;  }
		
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}

			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		if (m->control_pressed) { outputTypes.clear(); delete input; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } m->clearGroups(); return 0;  }

		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it);  
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		if (m->control_pressed) { outputTypes.clear(); delete input;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } m->clearGroups(); return 0;  }

		//run last label if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {  delete lookup[i]; }  } 
			lookup = input->getSharedRAbundVectors(lastLabel);

			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			process(lookup);
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
		if (m->control_pressed) { outputTypes.clear();  delete input;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } m->clearGroups(); return 0;  }
		
		//reset groups parameter
		m->clearGroups();  
		
		//set phylip file as new current phylipfile
		string current = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setPhylipFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();


		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void MatrixOutputCommand::printSims(ostream& out, vector< vector<float> >& simMatrix) {
	try {
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//output num seqs
		out << simMatrix.size() << endl;
		
		if (output == "lt") {
			for (int m = 0; m < simMatrix.size(); m++)	{
				out << lookup[m]->getGroup() << '\t';
				for (int n = 0; n < m; n++)	{
					out << simMatrix[m][n] << '\t'; 
				}
				out << endl;
			}
		}else{
			for (int m = 0; m < simMatrix.size(); m++)	{
				out << lookup[m]->getGroup() << '\t';
				for (int n = 0; n < simMatrix[m].size(); n++)	{
					out << simMatrix[m][n] << '\t'; 
				}
				out << endl;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "printSims");
		exit(1);
	}
}
/***********************************************************/
int MatrixOutputCommand::process(vector<SharedRAbundVector*> thisLookup){
	try {
		EstOutput data;
		vector<SharedRAbundVector*> subset;
		vector< vector<seqDist> > calcDists; calcDists.resize(matrixCalculators.size()); //one for each calc, this will be used to make .dist files
		
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			driver(thisLookup, 0, numGroups, calcDists);
		}else{
			int process = 1;
			vector<int> processIDS;
			
			//loop through and create all the processes you want
			while (process != processors) {
				int pid = fork();
				
				if (pid > 0) {
					processIDS.push_back(pid); 
					process++;
				}else if (pid == 0){
					driver(thisLookup, lines[process].start, lines[process].end, calcDists);   
					
					string tempdistFileName = m->getRootName(m->getSimpleName(sharedfile)) + toString(getpid()) + ".dist";
					ofstream outtemp;
					m->openOutputFile(tempdistFileName, outtemp);
						
					for (int i = 0; i < calcDists.size(); i++) {
						outtemp << calcDists[i].size() << endl;
							
						for (int j = 0; j < calcDists[i].size(); j++) {
							outtemp << calcDists[i][j].seq1 << '\t' << calcDists[i][j].seq2 << '\t' << calcDists[i][j].dist << endl;
						}
					}
					outtemp.close();
									
					exit(0);
				}else { 
					m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
					for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
					exit(0);
				}
			}
			
			//parent do your part
			driver(thisLookup, lines[0].start, lines[0].end, calcDists);   
						
			//force parent to wait until all the processes are done
			for (int i = 0; i < processIDS.size(); i++) {
				int temp = processIDS[i];
				wait(&temp);
			}
			
			for (int i = 0; i < processIDS.size(); i++) {
				string tempdistFileName = m->getRootName(m->getSimpleName(sharedfile)) + toString(processIDS[i]) +  ".dist";
				ifstream intemp;
				m->openInputFile(tempdistFileName, intemp);
					
				for (int k = 0; k < calcDists.size(); k++) {
					int size = 0;
					intemp >> size; m->gobble(intemp);
						
					for (int j = 0; j < size; j++) {
						int seq1 = 0;
						int seq2 = 0;
						float dist = 1.0;
							
						intemp >> seq1 >> seq2 >> dist;   m->gobble(intemp);
							
						seqDist tempDist(seq1, seq2, dist);
						calcDists[k].push_back(tempDist);
					}
				}
				intemp.close();
				m->mothurRemove(tempdistFileName);
			}
			
		}
#else
		driver(thisLookup, 0, numGroups, calcDists);
#endif
		
		for (int i = 0; i < calcDists.size(); i++) {
			if (m->control_pressed) { break; }
				
			//initialize matrix
			vector< vector<float> > matrix; //square matrix to represent the distance
			matrix.resize(thisLookup.size());
			for (int k = 0; k < thisLookup.size(); k++) {  matrix[k].resize(thisLookup.size(), 0.0); }
				
			for (int j = 0; j < calcDists[i].size(); j++) {
				int row = calcDists[i][j].seq1;
				int column = calcDists[i][j].seq2;
				float dist = calcDists[i][j].dist;
					
				matrix[row][column] = dist;
				matrix[column][row] = dist;
			}
			
			string distFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + matrixCalculators[i]->getName() + "." + thisLookup[0]->getLabel()  + "." + output + ".dist";
			outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
			ofstream outDist;
			m->openOutputFile(distFileName, outDist);
			outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
			
			printSims(outDist, matrix);
			
			outDist.close();
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
int MatrixOutputCommand::driver(vector<SharedRAbundVector*> thisLookup, int start, int end, vector< vector<seqDist> >& calcDists) { 
	try {
		
		vector<SharedRAbundVector*> subset;
		for (int k = start; k < end; k++) { // pass cdd each set of groups to compare
			
			for (int l = 0; l < k; l++) {
				
				if (k != l) { //we dont need to similiarity of a groups to itself
					subset.clear(); //clear out old pair of sharedrabunds
					//add new pair of sharedrabunds
					subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
					
					for(int i=0;i<matrixCalculators.size();i++) {
						
						//if this calc needs all groups to calculate the pair load all groups
						if (matrixCalculators[i]->getNeedsAll()) { 
							//load subset with rest of lookup for those calcs that need everyone to calc for a pair
							for (int w = 0; w < thisLookup.size(); w++) {
								if ((w != k) && (w != l)) { subset.push_back(thisLookup[w]); }
							}
						}
						
						vector<double> tempdata = matrixCalculators[i]->getValues(subset); //saves the calculator outputs
						
						if (m->control_pressed) { return 1; }
						
						seqDist temp(l, k, tempdata[0]);
						calcDists[i].push_back(temp);
					}
				}
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MatrixOutputCommand", "driver");
		exit(1);
	}
}
/***********************************************************/



