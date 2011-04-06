/*
 *  collectsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectsharedcommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
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
#include "sharednseqs.h"
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
vector<string> CollectSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "",false,false); parameters.push_back(pfreq);
		CommandParameter pcalc("calc", "Multiple", "sharedchao-sharedsobs-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-whittaker-sharednseqs-ochiai-anderberg-skulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-structchi2-hamming-gower-memchi2-memchord-memeuclidean-mempearson", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pall("all", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pall);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CollectSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The collect.shared command parameters are shared, label, freq, calc and groups.  shared is required if there is no current sharedfile. \n";
		helpString += "The collect.shared command should be in the following format: \n";
		helpString += "collect.shared(label=yourLabel, freq=yourFreq, calc=yourEstimators, groups=yourGroups).\n";
		helpString += "Example collect.shared(label=unique-.01-.03, freq=10, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan).\n";
		helpString += "The default values for freq is 100 and calc are sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan.\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += validCalculator.printCalc("shared");
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The all parameter is used to specify if you want the estimate of all your groups together.  This estimate can only be made for sharedsobs and sharedchao calculators. The default is false.\n";
		helpString += "If you use sharedchao and run into memory issues, set all to false. \n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		helpString += "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile).\n\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
CollectSharedCommand::CollectSharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sharedchao"] = tempOutNames;
		outputTypes["sharedsobs"] = tempOutNames;
		outputTypes["sharedace"] = tempOutNames;
		outputTypes["jabund"] = tempOutNames;
		outputTypes["sorabund"] = tempOutNames;
		outputTypes["jclass"] = tempOutNames;
		outputTypes["sorclass"] = tempOutNames;
		outputTypes["jest"] = tempOutNames;
		outputTypes["sorest"] = tempOutNames;
		outputTypes["thetayc"] = tempOutNames;
		outputTypes["thetan"] = tempOutNames;
		outputTypes["kstest"] = tempOutNames;
		outputTypes["whittaker"] = tempOutNames;
		outputTypes["sharednseqs"] = tempOutNames;
		outputTypes["ochiai"] = tempOutNames;
		outputTypes["anderberg"] = tempOutNames;
		outputTypes["skulczynski"] = tempOutNames;
		outputTypes["kulczynskicody"] = tempOutNames;
		outputTypes["lennon"] = tempOutNames;
		outputTypes["morisitahorn"] = tempOutNames;
		outputTypes["braycurtis"] = tempOutNames;
		outputTypes["odum"] = tempOutNames;
		outputTypes["canberra"] = tempOutNames;
		outputTypes["structeuclidean"] = tempOutNames;
		outputTypes["structchord"] = tempOutNames;
		outputTypes["hellinger"] = tempOutNames;
		outputTypes["manhattan"] = tempOutNames;
		outputTypes["structpearson"] = tempOutNames;
		outputTypes["soergel"] = tempOutNames;
		outputTypes["spearman"] = tempOutNames;
		outputTypes["structkulczynski"] = tempOutNames;
		outputTypes["structchi2"] = tempOutNames;
		outputTypes["speciesprofile"] = tempOutNames;
		outputTypes["hamming"] = tempOutNames;
		outputTypes["gower"] = tempOutNames;
		outputTypes["memchi2"] = tempOutNames;
		outputTypes["memchord"] = tempOutNames;
		outputTypes["memeuclidean"] = tempOutNames;
		outputTypes["mempearson"] = tempOutNames;
		
	}
	catch(exception& e) {
		m->errorOut(e, "CollectSharedCommand", "CollectSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CollectSharedCommand::CollectSharedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters=parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
	
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["sharedchao"] = tempOutNames;
			outputTypes["sharedsobs"] = tempOutNames;
			outputTypes["sharedace"] = tempOutNames;
			outputTypes["jabund"] = tempOutNames;
			outputTypes["sorabund"] = tempOutNames;
			outputTypes["jclass"] = tempOutNames;
			outputTypes["sorclass"] = tempOutNames;
			outputTypes["jest"] = tempOutNames;
			outputTypes["sorest"] = tempOutNames;
			outputTypes["thetayc"] = tempOutNames;
			outputTypes["thetan"] = tempOutNames;
			outputTypes["kstest"] = tempOutNames;
			outputTypes["whittaker"] = tempOutNames;
			outputTypes["sharednseqs"] = tempOutNames;
			outputTypes["ochiai"] = tempOutNames;
			outputTypes["anderberg"] = tempOutNames;
			outputTypes["skulczynski"] = tempOutNames;
			outputTypes["kulczynskicody"] = tempOutNames;
			outputTypes["lennon"] = tempOutNames;
			outputTypes["morisitahorn"] = tempOutNames;
			outputTypes["braycurtis"] = tempOutNames;
			outputTypes["odum"] = tempOutNames;
			outputTypes["canberra"] = tempOutNames;
			outputTypes["structeuclidean"] = tempOutNames;
			outputTypes["structchord"] = tempOutNames;
			outputTypes["hellinger"] = tempOutNames;
			outputTypes["manhattan"] = tempOutNames;
			outputTypes["structpearson"] = tempOutNames;
			outputTypes["soergel"] = tempOutNames;
			outputTypes["spearman"] = tempOutNames;
			outputTypes["structkulczynski"] = tempOutNames;
			outputTypes["speciesprofile"] = tempOutNames;
			outputTypes["structchi2"] = tempOutNames;
			outputTypes["hamming"] = tempOutNames;
			outputTypes["gower"] = tempOutNames;
			outputTypes["memchi2"] = tempOutNames;
			outputTypes["memchord"] = tempOutNames;
			outputTypes["memeuclidean"] = tempOutNames;
			outputTypes["mempearson"] = tempOutNames;
			
			
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
			// ...at some point should added some additional type checking..
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
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
			}
			m->Groups = Groups;
			
			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "all", false);				if (temp == "not found") { temp = "false"; }
			all = m->isTrue(temp);
						
			if (abort == false) {
				
				string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(sharedfile));
				
				ValidCalculators validCalculator;
				
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("shared", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedchao") { 
							cDisplays.push_back(new CollectDisplay(new SharedChao1(), new SharedOneColumnFile(fileNameRoot+"shared.chao")));
							outputNames.push_back(fileNameRoot+"shared.chao"); outputTypes["sharedchao"].push_back(fileNameRoot+"shared.chao");
						}else if (Estimators[i] == "sharedsobs") { 
							cDisplays.push_back(new CollectDisplay(new SharedSobsCS(), new SharedOneColumnFile(fileNameRoot+"shared.sobs")));
							outputNames.push_back(fileNameRoot+"shared.sobs"); outputTypes["sharedsobs"].push_back(fileNameRoot+"shared.sobs");
						}else if (Estimators[i] == "sharedace") { 
							cDisplays.push_back(new CollectDisplay(new SharedAce(), new SharedOneColumnFile(fileNameRoot+"shared.ace")));
							outputNames.push_back(fileNameRoot+"shared.ace"); outputTypes["sharedace"].push_back(fileNameRoot+"shared.ace");
						}else if (Estimators[i] == "jabund") { 	
							cDisplays.push_back(new CollectDisplay(new JAbund(), new SharedOneColumnFile(fileNameRoot+"jabund")));
							outputNames.push_back(fileNameRoot+"jabund"); outputTypes["jabund"].push_back(fileNameRoot+"jabund");
						}else if (Estimators[i] == "sorabund") { 
							cDisplays.push_back(new CollectDisplay(new SorAbund(), new SharedOneColumnFile(fileNameRoot+"sorabund")));
							outputNames.push_back(fileNameRoot+"sorabund"); outputTypes["sorabund"].push_back(fileNameRoot+"sorabund");
						}else if (Estimators[i] == "jclass") { 
							cDisplays.push_back(new CollectDisplay(new Jclass(), new SharedOneColumnFile(fileNameRoot+"jclass")));
							outputNames.push_back(fileNameRoot+"jclass"); outputTypes["jclass"].push_back(fileNameRoot+"jclass");
						}else if (Estimators[i] == "sorclass") { 
							cDisplays.push_back(new CollectDisplay(new SorClass(), new SharedOneColumnFile(fileNameRoot+"sorclass")));
							outputNames.push_back(fileNameRoot+"sorclass"); outputTypes["sorclass"].push_back(fileNameRoot+"sorclass");
						}else if (Estimators[i] == "jest") { 
							cDisplays.push_back(new CollectDisplay(new Jest(), new SharedOneColumnFile(fileNameRoot+"jest")));
							outputNames.push_back(fileNameRoot+"jest"); outputTypes["jest"].push_back(fileNameRoot+"jest");
						}else if (Estimators[i] == "sorest") { 
							cDisplays.push_back(new CollectDisplay(new SorEst(), new SharedOneColumnFile(fileNameRoot+"sorest")));
							outputNames.push_back(fileNameRoot+"sorest"); outputTypes["sorest"].push_back(fileNameRoot+"sorest");
						}else if (Estimators[i] == "thetayc") { 
							cDisplays.push_back(new CollectDisplay(new ThetaYC(), new SharedOneColumnFile(fileNameRoot+"thetayc")));
							outputNames.push_back(fileNameRoot+"thetayc"); outputTypes["thetayc"].push_back(fileNameRoot+"thetayc");
						}else if (Estimators[i] == "thetan") { 
							cDisplays.push_back(new CollectDisplay(new ThetaN(), new SharedOneColumnFile(fileNameRoot+"thetan")));
							outputNames.push_back(fileNameRoot+"thetan"); outputTypes["thetan"].push_back(fileNameRoot+"thetan");
						}else if (Estimators[i] == "kstest") { 
							cDisplays.push_back(new CollectDisplay(new KSTest(), new SharedOneColumnFile(fileNameRoot+"kstest")));
							outputNames.push_back(fileNameRoot+"kstest"); outputTypes["kstest"].push_back(fileNameRoot+"kstest");
						}else if (Estimators[i] == "whittaker") { 
							cDisplays.push_back(new CollectDisplay(new Whittaker(), new SharedOneColumnFile(fileNameRoot+"whittaker")));
							outputNames.push_back(fileNameRoot+"whittaker"); outputTypes["whittaker"].push_back(fileNameRoot+"whittaker");
						}else if (Estimators[i] == "sharednseqs") { 
							cDisplays.push_back(new CollectDisplay(new SharedNSeqs(), new SharedOneColumnFile(fileNameRoot+"shared.nseqs")));
							outputNames.push_back(fileNameRoot+"shared.nseqs"); outputTypes["shared.nseqs"].push_back(fileNameRoot+"shared.nseqs");
						}else if (Estimators[i] == "ochiai") { 
							cDisplays.push_back(new CollectDisplay(new Ochiai(), new SharedOneColumnFile(fileNameRoot+"ochiai")));
							outputNames.push_back(fileNameRoot+"ochiai"); outputTypes["ochiai"].push_back(fileNameRoot+"ochiai");
						}else if (Estimators[i] == "anderberg") { 
							cDisplays.push_back(new CollectDisplay(new Anderberg(), new SharedOneColumnFile(fileNameRoot+"anderberg")));
							outputNames.push_back(fileNameRoot+"anderberg"); outputTypes["anderberg"].push_back(fileNameRoot+"anderberg");
						}else if (Estimators[i] == "skulczynski") { 
							cDisplays.push_back(new CollectDisplay(new Kulczynski(), new SharedOneColumnFile(fileNameRoot+"kulczynski")));
							outputNames.push_back(fileNameRoot+"kulczynski"); outputTypes["kulczynski"].push_back(fileNameRoot+"kulczynski");
						}else if (Estimators[i] == "kulczynskicody") { 
							cDisplays.push_back(new CollectDisplay(new KulczynskiCody(), new SharedOneColumnFile(fileNameRoot+"kulczynskicody")));
							outputNames.push_back(fileNameRoot+"kulczynskicody"); outputTypes["kulczynskicody"].push_back(fileNameRoot+"kulczynskicody");
						}else if (Estimators[i] == "lennon") { 
							cDisplays.push_back(new CollectDisplay(new Lennon(), new SharedOneColumnFile(fileNameRoot+"lennon")));
							outputNames.push_back(fileNameRoot+"lennon"); outputTypes["lennon"].push_back(fileNameRoot+"lennon");
						}else if (Estimators[i] == "morisitahorn") { 
							cDisplays.push_back(new CollectDisplay(new MorHorn(), new SharedOneColumnFile(fileNameRoot+"morisitahorn")));
							outputNames.push_back(fileNameRoot+"morisitahorn"); outputTypes["morisitahorn"].push_back(fileNameRoot+"morisitahorn");
						}else if (Estimators[i] == "braycurtis") { 
							cDisplays.push_back(new CollectDisplay(new BrayCurtis(), new SharedOneColumnFile(fileNameRoot+"braycurtis")));
							outputNames.push_back(fileNameRoot+"braycurtis"); outputTypes["braycurtis"].push_back(fileNameRoot+"braycurtis");
						}else if (Estimators[i] == "odum") { 
							cDisplays.push_back(new CollectDisplay(new Odum(), new SharedOneColumnFile(fileNameRoot+"odum")));
							outputNames.push_back(fileNameRoot+"odum"); outputTypes["odum"].push_back(fileNameRoot+"odum");
						}else if (Estimators[i] == "canberra") { 
							cDisplays.push_back(new CollectDisplay(new Canberra(), new SharedOneColumnFile(fileNameRoot+"canberra")));
							outputNames.push_back(fileNameRoot+"canberra"); outputTypes["canberra"].push_back(fileNameRoot+"canberra");
						}else if (Estimators[i] == "structeuclidean") { 
							cDisplays.push_back(new CollectDisplay(new StructEuclidean(), new SharedOneColumnFile(fileNameRoot+"structeuclidean")));
							outputNames.push_back(fileNameRoot+"structeuclidean"); outputTypes["structeuclidean"].push_back(fileNameRoot+"structeuclidean");
						}else if (Estimators[i] == "structchord") { 
							cDisplays.push_back(new CollectDisplay(new StructChord(), new SharedOneColumnFile(fileNameRoot+"structchord")));
							outputNames.push_back(fileNameRoot+"structchord"); outputTypes["structchord"].push_back(fileNameRoot+"structchord");
						}else if (Estimators[i] == "hellinger") { 
							cDisplays.push_back(new CollectDisplay(new Hellinger(), new SharedOneColumnFile(fileNameRoot+"hellinger")));
							outputNames.push_back(fileNameRoot+"hellinger"); outputTypes["hellinger"].push_back(fileNameRoot+"hellinger");
						}else if (Estimators[i] == "manhattan") { 
							cDisplays.push_back(new CollectDisplay(new Manhattan(), new SharedOneColumnFile(fileNameRoot+"manhattan")));
							outputNames.push_back(fileNameRoot+"manhattan"); outputTypes["manhattan"].push_back(fileNameRoot+"manhattan");
						}else if (Estimators[i] == "structpearson") { 
							cDisplays.push_back(new CollectDisplay(new StructPearson(), new SharedOneColumnFile(fileNameRoot+"structpearson")));
							outputNames.push_back(fileNameRoot+"structpearson"); outputTypes["structpearson"].push_back(fileNameRoot+"structpearson");
						}else if (Estimators[i] == "soergel") { 
							cDisplays.push_back(new CollectDisplay(new Soergel(), new SharedOneColumnFile(fileNameRoot+"soergel")));
							outputNames.push_back(fileNameRoot+"soergel"); outputTypes["soergel"].push_back(fileNameRoot+"soergel");
						}else if (Estimators[i] == "spearman") { 
							cDisplays.push_back(new CollectDisplay(new Spearman(), new SharedOneColumnFile(fileNameRoot+"spearman")));
							outputNames.push_back(fileNameRoot+"spearman"); outputTypes["spearman"].push_back(fileNameRoot+"spearman");
						}else if (Estimators[i] == "structkulczynski") { 
							cDisplays.push_back(new CollectDisplay(new StructKulczynski(), new SharedOneColumnFile(fileNameRoot+"structkulczynski")));
							outputNames.push_back(fileNameRoot+"structkulczynski"); outputTypes["structkulczynski"].push_back(fileNameRoot+"structkulczynski");
						}else if (Estimators[i] == "speciesprofile") { 
							cDisplays.push_back(new CollectDisplay(new SpeciesProfile(), new SharedOneColumnFile(fileNameRoot+"speciesprofile")));
							outputNames.push_back(fileNameRoot+"speciesprofile"); outputTypes["speciesprofile"].push_back(fileNameRoot+"speciesprofile");
						}else if (Estimators[i] == "hamming") { 
							cDisplays.push_back(new CollectDisplay(new Hamming(), new SharedOneColumnFile(fileNameRoot+"hamming")));
							outputNames.push_back(fileNameRoot+"hamming"); outputTypes["hamming"].push_back(fileNameRoot+"hamming");
						}else if (Estimators[i] == "structchi2") { 
							cDisplays.push_back(new CollectDisplay(new StructChi2(), new SharedOneColumnFile(fileNameRoot+"structchi2")));
							outputNames.push_back(fileNameRoot+"structchi2"); outputTypes["structchi2"].push_back(fileNameRoot+"structchi2");
						}else if (Estimators[i] == "gower") { 
							cDisplays.push_back(new CollectDisplay(new Gower(), new SharedOneColumnFile(fileNameRoot+"gower")));
							outputNames.push_back(fileNameRoot+"gower"); outputTypes["gower"].push_back(fileNameRoot+"gower");
						}else if (Estimators[i] == "memchi2") { 
							cDisplays.push_back(new CollectDisplay(new MemChi2(), new SharedOneColumnFile(fileNameRoot+"memchi2")));
							outputNames.push_back(fileNameRoot+"memchi2"); outputTypes["memchi2"].push_back(fileNameRoot+"memchi2");
						}else if (Estimators[i] == "memchord") { 
							cDisplays.push_back(new CollectDisplay(new MemChord(), new SharedOneColumnFile(fileNameRoot+"memchord")));
							outputNames.push_back(fileNameRoot+"memchord"); outputTypes["memchord"].push_back(fileNameRoot+"memchord");
						}else if (Estimators[i] == "memeuclidean") { 
							cDisplays.push_back(new CollectDisplay(new MemEuclidean(), new SharedOneColumnFile(fileNameRoot+"memeuclidean")));
							outputNames.push_back(fileNameRoot+"memeuclidean"); outputTypes["memeuclidean"].push_back(fileNameRoot+"memeuclidean");
						}else if (Estimators[i] == "mempearson") { 
							cDisplays.push_back(new CollectDisplay(new MemPearson(), new SharedOneColumnFile(fileNameRoot+"mempearson")));
							outputNames.push_back(fileNameRoot+"mempearson"); outputTypes["mempearson"].push_back(fileNameRoot+"mempearson");
						}
						
					}
				}	
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "CollectSharedCommand", "CollectSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CollectSharedCommand::~CollectSharedCommand(){}
//**********************************************************************************************************************

int CollectSharedCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }
		for(int i=0;i<cDisplays.size();i++){	cDisplays[i]->setAll(all);	}	
	
		input = new InputData(sharedfile, "sharedfile");
		order = input->getSharedOrderVector();
		string lastLabel = order->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
			
		//set users groups
		SharedUtil* util = new SharedUtil();
		util->setGroups(m->Groups, m->namesOfGroups, "collect");
		delete util;

		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) { 
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}  outputTypes.clear();
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					delete order; delete input;
					m->Groups.clear();
					return 0;
			}

			if(allLines == 1 || labels.count(order->getLabel()) == 1){
			
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				cCurve->getSharedCurve(freq);
				delete cCurve;
			
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			//you have a label the user want that is smaller than this label and the last label has not already been processed
			if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = order->getLabel();
				
				delete order;
				order = input->getSharedOrderVector(lastLabel);
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				cCurve->getSharedCurve(freq);
				delete cCurve;
				
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
				
				//restore real lastlabel to save below
				order->setLabel(saveLabel);
			}
			
			
			lastLabel = order->getLabel();			
			
			//get next line to process
			delete order;
			order = input->getSharedOrderVector();
		}
		
		if (m->control_pressed) { 
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}   outputTypes.clear();
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					m->Groups.clear();
					delete input;
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
			if (order != NULL) {  delete order;  }
			order = input->getSharedOrderVector(lastLabel);
			
			m->mothurOut(order->getLabel()); m->mothurOutEndLine();
			cCurve = new Collect(order, cDisplays);
			cCurve->getSharedCurve(freq);
			delete cCurve;
			
			if (m->control_pressed) { 
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}  outputTypes.clear();
				for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
				delete order; 
				delete input;
				m->Groups.clear();
				return 0;
			}

			delete order;
		}
		
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}	
		
		//reset groups parameter
		m->Groups.clear(); 
		delete input;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectSharedCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
