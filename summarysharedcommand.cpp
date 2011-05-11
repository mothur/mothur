/*
 *  summarysharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarysharedcommand.h"
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
vector<string> SummarySharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pdistance("distance", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pdistance);
		CommandParameter pcalc("calc", "Multiple", "sharedchao-sharedsobs-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-whittaker-sharednseqs-ochiai-anderberg-skulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-structchi2-hamming-gower-memchi2-memchord-memeuclidean-mempearson", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pall("all", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pall);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummarySharedCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The summary.shared command parameters are shared, label, calc, distance, processors and all.  shared is required if there is no current sharedfile.\n";
		helpString += "The summary.shared command should be in the following format: \n";
		helpString += "summary.shared(label=yourLabel, calc=yourEstimators, groups=yourGroups).\n";
		helpString += "Example summary.shared(label=unique-.01-.03, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan).\n";
		helpString +=  validCalculator.printCalc("sharedsummary");
		helpString += "The default value for calc is sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The distance parameter allows you to indicate you would like a distance file created for each calculator for each label, default=f.\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The all parameter is used to specify if you want the estimate of all your groups together.  This estimate can only be made for sharedsobs and sharedchao calculators. The default is false.\n";
		helpString += "If you use sharedchao and run into memory issues, set all to false. \n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		helpString += "Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SummarySharedCommand::SummarySharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "SummarySharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

SummarySharedCommand::SummarySharedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			map<string, string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
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
			
			//get shared file
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
					
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			else { 
				 if (calc == "default")  {  calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
			}
			
			string temp = validParameter.validFile(parameters, "all", false);				if (temp == "not found") { temp = "false"; }
			all = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "distance", false);					if (temp == "not found") { temp = "false"; }
			createPhylip = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors); 
			
			if (abort == false) {
			
				ValidCalculators validCalculator;
				int i;
				
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("sharedsummary", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedsobs") { 
							sumCalculators.push_back(new SharedSobsCS());
						}else if (Estimators[i] == "sharedchao") { 
							sumCalculators.push_back(new SharedChao1());
						}else if (Estimators[i] == "sharedace") { 
							sumCalculators.push_back(new SharedAce());
						}else if (Estimators[i] == "jabund") { 	
							sumCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							sumCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							sumCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							sumCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							sumCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							sumCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							sumCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							sumCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "kstest") { 
							sumCalculators.push_back(new KSTest());
						}else if (Estimators[i] == "sharednseqs") { 
							sumCalculators.push_back(new SharedNSeqs());
						}else if (Estimators[i] == "ochiai") { 
							sumCalculators.push_back(new Ochiai());
						}else if (Estimators[i] == "anderberg") { 
							sumCalculators.push_back(new Anderberg());
						}else if (Estimators[i] == "kulczynski") { 
							sumCalculators.push_back(new Kulczynski());
						}else if (Estimators[i] == "kulczynskicody") { 
							sumCalculators.push_back(new KulczynskiCody());
						}else if (Estimators[i] == "lennon") { 
							sumCalculators.push_back(new Lennon());
						}else if (Estimators[i] == "morisitahorn") { 
							sumCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							sumCalculators.push_back(new BrayCurtis());
						}else if (Estimators[i] == "whittaker") { 
							sumCalculators.push_back(new Whittaker());
						}else if (Estimators[i] == "odum") { 
							sumCalculators.push_back(new Odum());
						}else if (Estimators[i] == "canberra") { 
							sumCalculators.push_back(new Canberra());
						}else if (Estimators[i] == "structeuclidean") { 
							sumCalculators.push_back(new StructEuclidean());
						}else if (Estimators[i] == "structchord") { 
							sumCalculators.push_back(new StructChord());
						}else if (Estimators[i] == "hellinger") { 
							sumCalculators.push_back(new Hellinger());
						}else if (Estimators[i] == "manhattan") { 
							sumCalculators.push_back(new Manhattan());
						}else if (Estimators[i] == "structpearson") { 
							sumCalculators.push_back(new StructPearson());
						}else if (Estimators[i] == "soergel") { 
							sumCalculators.push_back(new Soergel());
						}else if (Estimators[i] == "spearman") { 
							sumCalculators.push_back(new Spearman());
						}else if (Estimators[i] == "structkulczynski") { 
							sumCalculators.push_back(new StructKulczynski());
						}else if (Estimators[i] == "speciesprofile") { 
							sumCalculators.push_back(new SpeciesProfile());
						}else if (Estimators[i] == "hamming") { 
							sumCalculators.push_back(new Hamming());
						}else if (Estimators[i] == "structchi2") { 
							sumCalculators.push_back(new StructChi2());
						}else if (Estimators[i] == "gower") { 
							sumCalculators.push_back(new Gower());
						}else if (Estimators[i] == "memchi2") { 
							sumCalculators.push_back(new MemChi2());
						}else if (Estimators[i] == "memchord") { 
							sumCalculators.push_back(new MemChord());
						}else if (Estimators[i] == "memeuclidean") { 
							sumCalculators.push_back(new MemEuclidean());
						}else if (Estimators[i] == "mempearson") { 
							sumCalculators.push_back(new MemPearson());
						}
					}
				}
				
				mult = false;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "SummarySharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SummarySharedCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ofstream outputFileHandle, outAll;
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "shared.summary";
		
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }
		//check if any calcs can do multiples
		else{
			if (all){ 
				for (int i = 0; i < sumCalculators.size(); i++) {
					if (sumCalculators[i]->getMultiple() == true) { mult = true; }
				}
			}
		}
			
		input = new InputData(sharedfile, "sharedfile");
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
	
		/******************************************************/
		//output headings for files
		/******************************************************/
		//output estimator names as column headers
		m->openOutputFile(outputFileName, outputFileHandle);
		outputFileHandle << "label" <<'\t' << "comparison" << '\t'; 
		for(int i=0;i<sumCalculators.size();i++){
			outputFileHandle << '\t' << sumCalculators[i]->getName();
			if (sumCalculators[i]->getCols() == 3) {   outputFileHandle << "\t" << sumCalculators[i]->getName() << "_lci\t" << sumCalculators[i]->getName() << "_hci";  }
		}
		outputFileHandle << endl;
		outputFileHandle.close();
		
		//create file and put column headers for multiple groups file
		string outAllFileName = ((m->getRootName(sharedfile)) + "sharedmultiple.summary");
		if (mult == true) {
			m->openOutputFile(outAllFileName, outAll);
			outputNames.push_back(outAllFileName);
			
			outAll << "label" <<'\t' << "comparison" << '\t'; 
			for(int i=0;i<sumCalculators.size();i++){
				if (sumCalculators[i]->getMultiple() == true) { 
					outAll << '\t' << sumCalculators[i]->getName();
				}
			}
			outAll << endl;
			outAll.close();
		}
		
		if (lookup.size() < 2) { 
			m->mothurOut("I cannot run the command without at least 2 valid groups."); 
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			
			//close files and clean up
			remove(outputFileName.c_str());
			if (mult == true) { remove(outAllFileName.c_str());  }
			return 0;
		//if you only have 2 groups you don't need a .sharedmultiple file
		}else if ((lookup.size() == 2) && (mult == true)) { 
			mult = false;
			remove(outAllFileName.c_str());
			outputNames.pop_back();
		}
		
		if (m->control_pressed) {
			if (mult) {  remove(outAllFileName.c_str());  }
			remove(outputFileName.c_str()); 
			delete input;
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
			m->Groups.clear(); 
			return 0;
		}
		/******************************************************/
		
		
		/******************************************************/
		//comparison breakup to be used by different processes later
		numGroups = m->Groups.size();
		lines.resize(processors);
		for (int i = 0; i < processors; i++) {
			lines[i].start = int (sqrt(float(i)/float(processors)) * numGroups);
			lines[i].end = int (sqrt(float(i+1)/float(processors)) * numGroups);
		}		
		/******************************************************/
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
			
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {
				if (mult) {  remove(outAllFileName.c_str());  }
				remove(outputFileName.c_str()); 
				delete input; 
				for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
				for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
				m->Groups.clear(); 
				return 0;
			}

		
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup, outputFileName, outAllFileName);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
					
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					process(lookup, outputFileName, outAllFileName);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					//restore real lastlabel to save below
					lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();			
				
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		if (m->control_pressed) {
			if (mult) { remove(outAllFileName.c_str());  }
			remove(outputFileName.c_str()); 
			delete input; 
			for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
			m->Groups.clear(); 
			return 0;
		}

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
		
		//run last label if you need to
		if (needToRun == true)  {
				for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	} } 
				lookup = input->getSharedRAbundVectors(lastLabel);

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup, outputFileName, outAllFileName);
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
				
		//reset groups parameter
		m->Groups.clear();  
		
		for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		delete input;  
		
		if (m->control_pressed) {
			remove(outAllFileName.c_str());  
			remove(outputFileName.c_str()); 
			return 0;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		if (mult) { m->mothurOut(outAllFileName); m->mothurOutEndLine();	outputTypes["summary"].push_back(outAllFileName); }
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	} outputTypes["summary"].push_back(outputFileName);
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
int SummarySharedCommand::process(vector<SharedRAbundVector*> thisLookup, string sumFileName, string sumAllFileName) {
	try {
			vector< vector<seqDist> > calcDists;  //vector containing vectors that contains the summary results for each group compare
			calcDists.resize(sumCalculators.size()); //one for each calc, this will be used to make .dist files
				
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					driver(thisLookup, 0, numGroups, sumFileName+".temp", sumAllFileName+".temp", calcDists);
					m->appendFiles((sumFileName + ".temp"), sumFileName);
					remove((sumFileName + ".temp").c_str());
					if (mult) {
						m->appendFiles((sumAllFileName + ".temp"), sumAllFileName);
						remove((sumAllFileName + ".temp").c_str());
					}
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
							driver(thisLookup, lines[process].start, lines[process].end, sumFileName + toString(getpid()) + ".temp", sumAllFileName + toString(getpid()) + ".temp", calcDists);   
							
							//only do this if you want a distance file
							if (createPhylip) {
								string tempdistFileName = m->getRootName(m->getSimpleName(sumFileName)) + toString(getpid()) + ".dist";
								ofstream outtemp;
								m->openOutputFile(tempdistFileName, outtemp);
								
								for (int i = 0; i < calcDists.size(); i++) {
									outtemp << calcDists[i].size() << endl;
									
									for (int j = 0; j < calcDists[i].size(); j++) {
										outtemp << calcDists[i][j].seq1 << '\t' << calcDists[i][j].seq2 << '\t' << calcDists[i][j].dist << endl;
									}
								}
								outtemp.close();
							}
							
							exit(0);
						}else { 
							m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
							for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
							exit(0);
						}
					}
					
					//parent do your part
					driver(thisLookup, lines[0].start, lines[0].end, sumFileName + toString(getpid()) + ".temp", sumAllFileName + toString(getpid()) + ".temp", calcDists);   
					m->appendFiles((sumFileName + toString(getpid()) + ".temp"), sumFileName);
					remove((sumFileName + toString(getpid()) + ".temp").c_str());
					if (mult) { m->appendFiles((sumAllFileName + toString(getpid()) + ".temp"), sumAllFileName); }
						
					//force parent to wait until all the processes are done
					for (int i = 0; i < processIDS.size(); i++) {
						int temp = processIDS[i];
						wait(&temp);
					}
					
					for (int i = 0; i < processIDS.size(); i++) {
						m->appendFiles((sumFileName + toString(processIDS[i]) + ".temp"), sumFileName);
						remove((sumFileName + toString(processIDS[i]) + ".temp").c_str());
						if (mult) {	remove((sumAllFileName + toString(processIDS[i]) + ".temp").c_str());	}
						
						if (createPhylip) {
							string tempdistFileName = m->getRootName(m->getSimpleName(sumFileName)) + toString(processIDS[i]) +  ".dist";
							ifstream intemp;
							m->openInputFile(tempdistFileName, intemp);
							
							for (int i = 0; i < calcDists.size(); i++) {
								int size = 0;
								intemp >> size; m->gobble(intemp);
									
								for (int j = 0; j < size; j++) {
									int seq1 = 0;
									int seq2 = 0;
									float dist = 1.0;
									
									intemp >> seq1 >> seq2 >> dist;   m->gobble(intemp);
									
									seqDist tempDist(seq1, seq2, dist);
									calcDists[i].push_back(tempDist);
								}
							}
							intemp.close();
							remove(tempdistFileName.c_str());
						}
					}

				}
			#else
				driver(thisLookup, 0, numGroups, (sumFileName + ".temp"), (sumAllFileName + ".temp"), calcDists);
				m->appendFiles((sumFileName + ".temp"), sumFileName);
				remove((sumFileName + ".temp").c_str());
				if (mult) {
					m->appendFiles((sumAllFileName + ".temp"), sumAllFileName);
					remove((sumAllFileName + ".temp").c_str());
				}
			#endif
			
			if (createPhylip) {
				for (int i = 0; i < calcDists.size(); i++) {
					if (m->control_pressed) { break; }
				
					string distFileName = outputDir + m->getRootName(m->getSimpleName(sumFileName)) + sumCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".dist";
					outputNames.push_back(distFileName);
					ofstream outDist;
					m->openOutputFile(distFileName, outDist);
					outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
					
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
					
					//output to file
					outDist << thisLookup.size() << endl;
					for (int r=0; r<thisLookup.size(); r++) { 
						//output name
						string name = thisLookup[r]->getGroup();
						if (name.length() < 10) { //pad with spaces to make compatible
							while (name.length() < 10) {  name += " ";  }
						}
						outDist << name << '\t';
					
						//output distances
						for (int l = 0; l < r; l++) {	outDist  << matrix[r][l] << '\t';  }
						outDist << endl;
					}
					
					outDist.close();
				}
			}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
int SummarySharedCommand::driver(vector<SharedRAbundVector*> thisLookup, int start, int end, string sumFile, string sumAllFile, vector< vector<seqDist> >& calcDists) { 
	try {
		
		//loop through calculators and add to file all for all calcs that can do mutiple groups
		if (mult == true) {
			ofstream outAll;
			m->openOutputFile(sumAllFile, outAll);
			
			//output label
			outAll << thisLookup[0]->getLabel() << '\t';
			
			//output groups names
			string outNames = "";
			for (int j = 0; j < thisLookup.size(); j++) {
				outNames += thisLookup[j]->getGroup() +  "-";
			}
			outNames = outNames.substr(0, outNames.length()-1); //rip off extra '-';
			outAll << outNames << '\t';
			
			for(int i=0;i<sumCalculators.size();i++){
				if (sumCalculators[i]->getMultiple() == true) { 
					sumCalculators[i]->getValues(thisLookup);
					
					if (m->control_pressed) { outAll.close(); return 1; }
					
					outAll << '\t';
					sumCalculators[i]->print(outAll);
				}
			}
			outAll << endl;
			outAll.close();
		}
		
		ofstream outputFileHandle;
		m->openOutputFile(sumFile, outputFileHandle);
		
		vector<SharedRAbundVector*> subset;
		for (int k = start; k < end; k++) { // pass cdd each set of groups to compare

			for (int l = 0; l < k; l++) {
				
				outputFileHandle << thisLookup[0]->getLabel() << '\t';
				
				subset.clear(); //clear out old pair of sharedrabunds
				//add new pair of sharedrabunds
				subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
				
				//sort groups to be alphanumeric
				if (thisLookup[k]->getGroup() > thisLookup[l]->getGroup()) {
					outputFileHandle << (thisLookup[l]->getGroup() +'\t' + thisLookup[k]->getGroup()) << '\t'; //print out groups
				}else{
					outputFileHandle << (thisLookup[k]->getGroup() +'\t' + thisLookup[l]->getGroup()) << '\t'; //print out groups
				}
				
				for(int i=0;i<sumCalculators.size();i++) {
					
					//if this calc needs all groups to calculate the pair load all groups
					if (sumCalculators[i]->getNeedsAll()) { 
						//load subset with rest of lookup for those calcs that need everyone to calc for a pair
						for (int w = 0; w < thisLookup.size(); w++) {
							if ((w != k) && (w != l)) { subset.push_back(thisLookup[w]); }
						}
					}
					
					vector<double> tempdata = sumCalculators[i]->getValues(subset); //saves the calculator outputs
					
					if (m->control_pressed) { outputFileHandle.close(); return 1; }
					
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
					
					seqDist temp(l, k, (1.0 - tempdata[0]));
					calcDists[i].push_back(temp);
				}
				outputFileHandle << endl;
			}
		}
		
		outputFileHandle.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/


