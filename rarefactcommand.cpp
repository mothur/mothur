/*
 *  rarefactcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactcommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "simpsoneven.h"
#include "heip.h"
#include "smithwilson.h"
#include "invsimpson.h"
#include "npshannon.h"
#include "shannoneven.h"
#include "shannon.h"
#include "jackknife.h"
#include "coverage.h"

//**********************************************************************************************************************
vector<string> RareFactCommand::getValidParameters(){	
	try {
		string Array[] =  {"iters","freq","label","calc","abund","processors","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactCommand::getRequiredFiles(){	
	try {
		string AlignArray[] =  {"shared","list","rabund","sabund","or"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
RareFactCommand::RareFactCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["rarefaction"] = tempOutNames;
		outputTypes["r_chao"] = tempOutNames;
		outputTypes["r_ace"] = tempOutNames;
		outputTypes["r_jack"] = tempOutNames;
		outputTypes["r_shannon"] = tempOutNames;
		outputTypes["r_shannoneven"] = tempOutNames;
		outputTypes["r_heip"] = tempOutNames;
		outputTypes["r_smithwilson"] = tempOutNames;
		outputTypes["r_npshannon"] = tempOutNames;
		outputTypes["r_simpson"] = tempOutNames;
		outputTypes["r_simpsoneven"] = tempOutNames;
		outputTypes["r_invsimpson"] = tempOutNames;
		outputTypes["r_bootstrap"] = tempOutNames;
		outputTypes["r_coverage"] = tempOutNames;
		outputTypes["r_nseqs"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "RareFactCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RareFactCommand::RareFactCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		Estimators.clear();
				
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); delete validCalculator; abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"iters","freq","label","calc","abund","processors","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["rarefaction"] = tempOutNames;
			outputTypes["r_chao"] = tempOutNames;
			outputTypes["r_ace"] = tempOutNames;
			outputTypes["r_jack"] = tempOutNames;
			outputTypes["r_shannon"] = tempOutNames;
			outputTypes["r_shannoneven"] = tempOutNames;
			outputTypes["r_heip"] = tempOutNames;
			outputTypes["r_smithwilson"] = tempOutNames;
			outputTypes["r_npshannon"] = tempOutNames;
			outputTypes["r_simpson"] = tempOutNames;
			outputTypes["r_simpsoneven"] = tempOutNames;
			outputTypes["r_invsimpson"] = tempOutNames;
			outputTypes["r_bootstrap"] = tempOutNames;
			outputTypes["r_coverage"] = tempOutNames;
			outputTypes["r_nseqs"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}

			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "") && (globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { m->mothurOut("You must read a list, sabund, rabund or shared file before you can use the rarefact.single command."); m->mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if(label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs";  }
			else { 
				 if (calc == "default")  {  calc = "sobs";  }
			}
			m->splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "abund", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, nIters); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "RareFactCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void RareFactCommand::help(){
	try {
		m->mothurOut("The rarefaction.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION.\n");
		m->mothurOut("The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		m->mothurOut("The rarefaction.single command parameters are label, iters, freq, calc, processors and abund.  No parameters are required. \n");
		m->mothurOut("The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n");
		m->mothurOut("The processors parameter allows you to specify the number of processors to use. The default is 1.\n");
		m->mothurOut("The rarefaction.single command should be in the following format: \n");
		m->mothurOut("rarefaction.single(label=yourLabel, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n");
		m->mothurOut("Example rarefaction.single(label=unique-.01-.03, iters=10000, freq=10, calc=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson).\n");
		m->mothurOut("The default values for iters is 1000, freq is 100, and calc is rarefaction which calculates the rarefaction curve for the observed richness.\n");
		validCalculator->printCalc("rarefaction", cout);
		m->mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

RareFactCommand::~RareFactCommand(){}

//**********************************************************************************************************************

int RareFactCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		string hadShared = "";
		if ((globaldata->getFormat() != "sharedfile")) { inputFileNames.push_back(globaldata->inputFileName);  }
		else { hadShared = globaldata->getSharedFile(); inputFileNames = parseSharedFile(globaldata->getSharedFile());  globaldata->setFormat("rabund");  }
				
		if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } return 0; }
		
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p]));
			globaldata->inputFileName = inputFileNames[p];
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } return 0; }
			
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
			int i;
			validCalculator = new ValidCalculators();
			
			
			for (i=0; i<Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("rarefaction", Estimators[i]) == true) { 
					if (Estimators[i] == "sobs") { 
						rDisplays.push_back(new RareDisplay(new Sobs(), new ThreeColumnFile(fileNameRoot+"rarefaction")));
						outputNames.push_back(fileNameRoot+"rarefaction"); outputTypes["rarefaction"].push_back(fileNameRoot+"rarefaction");
					}else if (Estimators[i] == "chao") { 
						rDisplays.push_back(new RareDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"r_chao")));
						outputNames.push_back(fileNameRoot+"r_chao"); outputTypes["r_chao"].push_back(fileNameRoot+"r_chao");
					}else if (Estimators[i] == "ace") { 
						if(abund < 5)
							abund = 10;
						rDisplays.push_back(new RareDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"r_ace")));
						outputNames.push_back(fileNameRoot+"r_ace"); outputTypes["r_ace"].push_back(fileNameRoot+"r_ace");
					}else if (Estimators[i] == "jack") { 
						rDisplays.push_back(new RareDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"r_jack")));
						outputNames.push_back(fileNameRoot+"r_jack"); outputTypes["r_jack"].push_back(fileNameRoot+"r_jack");
					}else if (Estimators[i] == "shannon") { 
						rDisplays.push_back(new RareDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"r_shannon")));
						outputNames.push_back(fileNameRoot+"r_shannon"); outputTypes["r_shannon"].push_back(fileNameRoot+"r_shannon");
					}else if (Estimators[i] == "shannoneven") { 
						rDisplays.push_back(new RareDisplay(new ShannonEven(), new ThreeColumnFile(fileNameRoot+"r_shannoneven")));
						outputNames.push_back(fileNameRoot+"r_shannoneven"); outputTypes["r_shannoneven"].push_back(fileNameRoot+"r_shannoneven");
					}else if (Estimators[i] == "heip") { 
						rDisplays.push_back(new RareDisplay(new Heip(), new ThreeColumnFile(fileNameRoot+"r_heip")));
						outputNames.push_back(fileNameRoot+"r_heip"); outputTypes["r_heip"].push_back(fileNameRoot+"r_heip");
					}else if (Estimators[i] == "smithwilson") { 
						rDisplays.push_back(new RareDisplay(new SmithWilson(), new ThreeColumnFile(fileNameRoot+"r_smithwilson")));
						outputNames.push_back(fileNameRoot+"r_smithwilson"); outputTypes["r_smithwilson"].push_back(fileNameRoot+"r_smithwilson");
					}else if (Estimators[i] == "npshannon") { 
						rDisplays.push_back(new RareDisplay(new NPShannon(), new ThreeColumnFile(fileNameRoot+"r_npshannon")));
						outputNames.push_back(fileNameRoot+"r_npshannon"); outputTypes["r_npshannon"].push_back(fileNameRoot+"r_npshannon");
					}else if (Estimators[i] == "simpson") { 
						rDisplays.push_back(new RareDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"r_simpson")));
						outputNames.push_back(fileNameRoot+"r_simpson"); outputTypes["r_simpson"].push_back(fileNameRoot+"r_simpson");
					}else if (Estimators[i] == "simpsoneven") { 
						rDisplays.push_back(new RareDisplay(new SimpsonEven(), new ThreeColumnFile(fileNameRoot+"r_simpsoneven")));
						outputNames.push_back(fileNameRoot+"r_simpsoneven"); outputTypes["r_simpsoneven"].push_back(fileNameRoot+"r_simpsoneven");
					}else if (Estimators[i] == "invsimpson") { 
						rDisplays.push_back(new RareDisplay(new InvSimpson(), new ThreeColumnFile(fileNameRoot+"r_invsimpson")));
						outputNames.push_back(fileNameRoot+"r_invsimpson"); outputTypes["r_invsimpson"].push_back(fileNameRoot+"r_invsimpson");
					}else if (Estimators[i] == "bootstrap") { 
						rDisplays.push_back(new RareDisplay(new Bootstrap(), new ThreeColumnFile(fileNameRoot+"r_bootstrap")));
						outputNames.push_back(fileNameRoot+"r_bootstrap"); outputTypes["r_bootstrap"].push_back(fileNameRoot+"r_bootstrap");
					}else if (Estimators[i] == "coverage") { 
						rDisplays.push_back(new RareDisplay(new Coverage(), new ThreeColumnFile(fileNameRoot+"r_coverage")));
						outputNames.push_back(fileNameRoot+"r_coverage"); outputTypes["r_coverage"].push_back(fileNameRoot+"r_coverage");
					}else if (Estimators[i] == "nseqs") { 
						rDisplays.push_back(new RareDisplay(new NSeqs(), new ThreeColumnFile(fileNameRoot+"r_nseqs")));
						outputNames.push_back(fileNameRoot+"r_nseqs"); outputTypes["r_nseqs"].push_back(fileNameRoot+"r_nseqs");
					}
				}
			}
			
			
			//if the users entered no valid calculators don't execute command
			if (rDisplays.size() == 0) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	} delete validCalculator; return 0; }
			
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			order = globaldata->gorder;
			string lastLabel = order->getLabel();
			input = globaldata->ginput;
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	} delete validCalculator; delete read; delete input; globaldata->ginput = NULL; delete order; globaldata->gorder = NULL; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	} delete validCalculator; delete read; delete input; globaldata->ginput = NULL; delete order; globaldata->gorder = NULL; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }

				
				if(allLines == 1 || labels.count(order->getLabel()) == 1){
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					rCurve = new Rarefact(order, rDisplays, processors);
					rCurve->getCurve(freq, nIters);
					delete rCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
				}
				
				if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = order->getLabel();
					
					delete order;
					order = (input->getOrderVector(lastLabel));
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					rCurve = new Rarefact(order, rDisplays, processors);
					rCurve->getCurve(freq, nIters);
					delete rCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					//restore real lastlabel to save below
					order->setLabel(saveLabel);
				}
				
				lastLabel = order->getLabel();		
				
				delete order;
				order = (input->getOrderVector());
			}
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	} delete validCalculator; delete read; delete input; globaldata->ginput = NULL; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  return 0; }

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
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	} delete validCalculator; delete read; delete input; globaldata->ginput = NULL;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }

			//run last label if you need to
			if (needToRun == true)  {
				if (order != NULL) {	delete order;	}
				order = (input->getOrderVector(lastLabel));
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				rCurve = new Rarefact(order, rDisplays, processors);
				rCurve->getCurve(freq, nIters);
				delete rCurve;
				
				delete order;
			}
			
			
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
			rDisplays.clear();
			globaldata->gorder = NULL;
			delete input;  globaldata->ginput = NULL;
			delete read;
			delete validCalculator;
			
		}
		
		if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  }
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactCommand::parseSharedFile(string filename) {
	try {
		vector<string> filenames;
		
		map<string, ofstream*> filehandles;
		map<string, ofstream*>::iterator it3;
		
				
		//read first line
		read = new ReadOTUFile(filename);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		
		string sharedFileRoot = m->getRootName(filename);
		
		//clears file before we start to write to it below
		for (int i=0; i<lookup.size(); i++) {
			remove((sharedFileRoot + lookup[i]->getGroup() + ".rabund").c_str());
			filenames.push_back((sharedFileRoot + lookup[i]->getGroup() + ".rabund"));
		}
		
		ofstream* temp;
		for (int i=0; i<lookup.size(); i++) {
			temp = new ofstream;
			filehandles[lookup[i]->getGroup()] = temp;
			groups.push_back(lookup[i]->getGroup());
		}

		while(lookup[0] != NULL) {
		
			for (int i = 0; i < lookup.size(); i++) {
				RAbundVector rav = lookup[i]->getRAbundVector();
				m->openOutputFileAppend(sharedFileRoot + lookup[i]->getGroup() + ".rabund", *(filehandles[lookup[i]->getGroup()]));
				rav.print(*(filehandles[lookup[i]->getGroup()]));
				(*(filehandles[lookup[i]->getGroup()])).close();
			}
		
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		//free memory
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}
		delete read;
		delete input;
		globaldata->ginput = NULL;

		return filenames;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "parseSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************



