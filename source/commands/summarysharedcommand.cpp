/*
 *  summarysharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarysharedcommand.h"
#include "subsample.h"

//**********************************************************************************************************************
vector<string> SummarySharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","phylip",false,false); parameters.push_back(psubsample);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
		CommandParameter pdistance("distance", "Boolean", "", "F", "", "", "","phylip",false,false); parameters.push_back(pdistance);
		CommandParameter pcalc("calc", "Multiple", "sharedchao-sharedsobs-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-whittaker-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-structchi2-hamming-gower-memchi2-memchord-memeuclidean-mempearson-jsd-rjsd", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan", "", "", "","",true,false,true); parameters.push_back(pcalc);
        CommandParameter poutput("output", "Multiple", "lt-square", "lt", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter pall("all", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pall);
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;  allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
        outputTypes["phylip"] = tempOutNames;
		
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
		helpString += "The summary.shared command parameters are " + getCommandParameters() + ".\nThe shared is required if there is no current sharedfile.\n";
		helpString += "The summary.shared command should be in the following format: \n";
		helpString += "summary.shared(shared=yourSharedFile).\n";
		helpString += "Example summary.shared(shared=final.opti_mcc.shared)\n";
		helpString +=  validCalculator.printCalc("sharedsummary");
        helpString += "The iters parameter allows you to choose the number of times you would like to run the subsample.\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group.\n";
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
        helpString += "The output parameter allows you to specify format of your distance matrix. Options are lt, and square. The default is lt.\n";
		helpString += "The default value for calc is sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The distance parameter allows you to indicate you would like a distance file created for each calculator for each label, default=f.\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The all parameter is used to specify if you want the estimate of all your groups together.  This estimate can only be made for sharedsobs and sharedchao calculators. The default is false.\n";
		helpString += "If you use sharedchao and run into memory issues, set all to false. \n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummarySharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],summary-[filename],[tag],summary"; } 
        else if (type == "phylip") {  pattern = "[filename],[calc],[distance],[outputtag],[tag2],dist"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SummarySharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SummarySharedCommand::SummarySharedCommand(string option) : Command()  {
	try {

		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required.\n"); abort = true; }
			}else { current->setSharedFile(sharedfile); }
			
			
			 
					if (outputdir == ""){    outputdir = util.hasPath(sharedfile);		}
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
					
			calc = validParameter.valid(parameters, "calc");			
			if (calc == "not found") { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			else { 
				 if (calc == "default")  {  calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			}
			util.splitAtDash(calc, Estimators);
			if (util.inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else {  util.splitAtDash(groups, Groups);
                    if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } } }
			
			string temp = validParameter.valid(parameters, "all");				if (temp == "not found") { temp = "false"; }
			all = util.isTrue(temp);
			
            temp = validParameter.valid(parameters, "iters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
            
            output = validParameter.valid(parameters, "output");
            if(output == "not found"){	output = "lt"; }
            else { createPhylip = true; }
			if ((output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are lt and square. I will use lt.\n"); output = "lt"; }
            
            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
            
            if (!subsample) { iters = 1; }
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "distance");					if (temp == "not found") { temp = "false"; }
			createPhylip = util.isTrue(temp);
            if (subsample) { createPhylip = true; }
            
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
            mult = false;
            numCalcs = 0;
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
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        ValidCalculators validCalculator;
        vector<Calculator*> sumCalculators;
        for (int i=0; i<Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("sharedsummary", Estimators[i]) ) {
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
                }else if (Estimators[i] == "jsd") {
                    sumCalculators.push_back(new JSD());
                }else if (Estimators[i] == "rjsd") {
                    sumCalculators.push_back(new RJSD());
                }
            }
        }
		
		ofstream outputFileHandle, outAll;
        map<string, string> variables; 
		variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
		string outputFileName = getOutputFileName("summary",variables);
		
		//if the users entered no valid calculators don't execute command
		if (Estimators.size() == 0) { return 0; }
		//check if any calcs can do multiples
		else{ if (all){  for (int i = 0; i < sumCalculators.size(); i++) { if (sumCalculators[i]->getMultiple() ) { mult = true; } } } }
			
		InputData input(sharedfile, "sharedfile", Groups);
		set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();

		/******************************************************/
		//output headings for files
		/******************************************************/
		//output estimator names as column headers
        if (!subsample) {
            util.openOutputFile(outputFileName, outputFileHandle);
            outputFileHandle << "label" <<'\t' << "comparison" << '\t';
            for(int i=0;i<sumCalculators.size();i++){
                outputFileHandle << '\t' << sumCalculators[i]->getName();
                if (sumCalculators[i]->getCols() == 3) {   outputFileHandle << "\t" << sumCalculators[i]->getName() << "_lci\t" << sumCalculators[i]->getName() << "_hci";  }
            }
            outputFileHandle << endl;
            outputFileHandle.close();
        }
		//create file and put column headers for multiple groups file
        variables["[tag]"]= "multiple";
		string outAllFileName = getOutputFileName("summary",variables);
		if (mult ) {
            if (!subsample) {
                util.openOutputFile(outAllFileName, outAll);
                outputNames.push_back(outAllFileName);
                
                outAll << "label" <<'\t' << "comparison" << '\t';
                for(int i=0;i<sumCalculators.size();i++){ if (sumCalculators[i]->getMultiple() ) {  outAll << '\t' << sumCalculators[i]->getName(); } }
                outAll << endl;
                outAll.close();
            }
		}
		
		if (lookup->size() < 2) {
			m->mothurOut("I cannot run the command without at least 2 valid groups."); 
			delete lookup;
			
            if (!subsample) {
                //close files and clean up
                util.mothurRemove(outputFileName);
                if (mult ) { util.mothurRemove(outAllFileName);  }
            }
			return 0;
		//if you only have 2 groups you don't need a .sharedmultiple file
		}else if ((lookup->size() == 2) && (mult )) {
			mult = false;
			util.mothurRemove(outAllFileName);
			outputNames.pop_back();
		}
		
		if (m->getControl_pressed()) { if (mult) {  util.mothurRemove(outAllFileName);  } util.mothurRemove(outputFileName); delete lookup; for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }   return 0; }
		/******************************************************/
        if (subsample) { 
            if (subsampleSize == -1) { //user has not set size, set size = smallest samples size
                subsampleSize = lookup->getNumSeqsSmallestGroup();
            }else {
                lookup->removeGroups(subsampleSize);
                Groups = lookup->getNamesGroups();
            }
            
            if (lookup->size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command.\n"); m->setControl_pressed(true);  return 0; }
        }

		
		/******************************************************/
		//comparison breakup to be used by different processes later
		numGroups = lookup->size();
        numCalcs = sumCalculators.size();
        for(int i=0;i<sumCalculators.size();i++){  sumCalculatorsNames.push_back(sumCalculators[i]->getName()); }
        
		lines.resize(processors);
		for (int i = 0; i < processors; i++) {
			lines[i].start = int (sqrt(float(i)/float(processors)) * numGroups);
			lines[i].end = int (sqrt(float(i+1)/float(processors)) * numGroups);
		}		
		/******************************************************/
		
        vector<string> currentLabels = lookup->getOTUNames();
        
        while (lookup != nullptr) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            process(lookup, outputFileName, outAllFileName, currentLabels); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }

		for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		
		if (m->getControl_pressed()) { util.mothurRemove(outAllFileName);   util.mothurRemove(outputFileName);  return 0; }
		
		m->mothurOut("\nOutput File Names:\n");
        if (!subsample) {
            m->mothurOut(outputFileName+"\n");
            if (mult) { m->mothurOut(outAllFileName+"\n"); outputTypes["summary"].push_back(outAllFileName); }
            outputTypes["summary"].push_back(outputFileName);
        }else { for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]+"\n"); 	} } m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
int SummarySharedCommand::printSims(ostream& out, vector< vector<double> >& simMatrix, vector<string> theseGroups) {
	try {
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//output num seqs
		out << simMatrix.size() << endl;
		
		if (output == "lt") {
			for (int b = 0; b < simMatrix.size(); b++)	{
				out << theseGroups[b];
				for (int n = 0; n < b; n++)	{
                    if (m->getControl_pressed()) { return 0; }
					out << '\t' << simMatrix[b][n];
				}
				out << endl;
			}
		}else{
			for (int b = 0; b < simMatrix.size(); m++)	{
				out << theseGroups[b];
				for (int n = 0; n < simMatrix[b].size(); n++)	{
                    if (m->getControl_pressed()) { return 0; }
					out << '\t' << simMatrix[b][n];
				}
				out << endl;
			}
		}
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummarySharedCommand", "printSims");
		exit(1);
	}
}
/**************************************************************************************************/
void driverSummaryShared(summarySharedData* params) {
    try {
        
        vector<Calculator*> sumCalculators;
        ValidCalculators validCalculator;
        for (int i=0; i<params->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("sharedsummary", params->Estimators[i]) ) {
                if (params->Estimators[i] == "sharedsobs") {
                    sumCalculators.push_back(new SharedSobsCS());
                }else if (params->Estimators[i] == "sharedchao") {
                    sumCalculators.push_back(new SharedChao1());
                }else if (params->Estimators[i] == "sharedace") {
                    sumCalculators.push_back(new SharedAce());
                }else if (params->Estimators[i] == "jabund") {
                    sumCalculators.push_back(new JAbund());
                }else if (params->Estimators[i] == "sorabund") {
                    sumCalculators.push_back(new SorAbund());
                }else if (params->Estimators[i] == "jclass") {
                    sumCalculators.push_back(new Jclass());
                }else if (params->Estimators[i] == "sorclass") {
                    sumCalculators.push_back(new SorClass());
                }else if (params->Estimators[i] == "jest") {
                    sumCalculators.push_back(new Jest());
                }else if (params->Estimators[i] == "sorest") {
                    sumCalculators.push_back(new SorEst());
                }else if (params->Estimators[i] == "thetayc") {
                    sumCalculators.push_back(new ThetaYC());
                }else if (params->Estimators[i] == "thetan") {
                    sumCalculators.push_back(new ThetaN());
                }else if (params->Estimators[i] == "kstest") {
                    sumCalculators.push_back(new KSTest());
                }else if (params->Estimators[i] == "sharednseqs") {
                    sumCalculators.push_back(new SharedNSeqs());
                }else if (params->Estimators[i] == "ochiai") {
                    sumCalculators.push_back(new Ochiai());
                }else if (params->Estimators[i] == "anderberg") {
                    sumCalculators.push_back(new Anderberg());
                }else if (params->Estimators[i] == "kulczynski") {
                    sumCalculators.push_back(new Kulczynski());
                }else if (params->Estimators[i] == "kulczynskicody") {
                    sumCalculators.push_back(new KulczynskiCody());
                }else if (params->Estimators[i] == "lennon") {
                    sumCalculators.push_back(new Lennon());
                }else if (params->Estimators[i] == "morisitahorn") {
                    sumCalculators.push_back(new MorHorn());
                }else if (params->Estimators[i] == "braycurtis") {
                    sumCalculators.push_back(new BrayCurtis());
                }else if (params->Estimators[i] == "whittaker") {
                    sumCalculators.push_back(new Whittaker());
                }else if (params->Estimators[i] == "odum") {
                    sumCalculators.push_back(new Odum());
                }else if (params->Estimators[i] == "canberra") {
                    sumCalculators.push_back(new Canberra());
                }else if (params->Estimators[i] == "structeuclidean") {
                    sumCalculators.push_back(new StructEuclidean());
                }else if (params->Estimators[i] == "structchord") {
                    sumCalculators.push_back(new StructChord());
                }else if (params->Estimators[i] == "hellinger") {
                    sumCalculators.push_back(new Hellinger());
                }else if (params->Estimators[i] == "manhattan") {
                    sumCalculators.push_back(new Manhattan());
                }else if (params->Estimators[i] == "structpearson") {
                    sumCalculators.push_back(new StructPearson());
                }else if (params->Estimators[i] == "soergel") {
                    sumCalculators.push_back(new Soergel());
                }else if (params->Estimators[i] == "spearman") {
                    sumCalculators.push_back(new Spearman());
                }else if (params->Estimators[i] == "structkulczynski") {
                    sumCalculators.push_back(new StructKulczynski());
                }else if (params->Estimators[i] == "speciesprofile") {
                    sumCalculators.push_back(new SpeciesProfile());
                }else if (params->Estimators[i] == "hamming") {
                    sumCalculators.push_back(new Hamming());
                }else if (params->Estimators[i] == "structchi2") {
                    sumCalculators.push_back(new StructChi2());
                }else if (params->Estimators[i] == "gower") {
                    sumCalculators.push_back(new Gower());
                }else if (params->Estimators[i] == "memchi2") {
                    sumCalculators.push_back(new MemChi2());
                }else if (params->Estimators[i] == "memchord") {
                    sumCalculators.push_back(new MemChord());
                }else if (params->Estimators[i] == "memeuclidean") {
                    sumCalculators.push_back(new MemEuclidean());
                }else if (params->Estimators[i] == "mempearson") {
                    sumCalculators.push_back(new MemPearson());
                }else if (params->Estimators[i] == "jsd") {
                    sumCalculators.push_back(new JSD());
                }else if (params->Estimators[i] == "rjsd") {
                    sumCalculators.push_back(new RJSD());
                }
            }
        }
        
        params->calcDists.resize(sumCalculators.size());
        
        //loop through calculators and add to file all for all calcs that can do mutiple groups
        if (params->mult && params->main) {
            ofstream outAll;
            
            if (!params->subsample) { //print names
                params->util.openOutputFile(params->sumAllFile, outAll);
                outAll << params->thisLookup[0]->getLabel() << '\t'; //output label
            
                //output groups names
                string outNames = "";
                for (int j = 0; j < params->thisLookup.size(); j++) { outNames += params->thisLookup[j]->getGroup() +  "-"; }
                
                outNames = outNames.substr(0, outNames.length()-1); //rip off extra '-';
                outAll << outNames << '\t';
            }
            
            for(int i=0;i<sumCalculators.size();i++){
                if (sumCalculators[i]->getMultiple() ) {
                    sumCalculators[i]->getValues(params->thisLookup);
                    
                    if (params->m->getControl_pressed()) { break; }
                    
                    if (!params->subsample) {
                        outAll << '\t';
                        sumCalculators[i]->print(outAll);
                    }
                }
            }
            if (!params->subsample) { outAll << endl; outAll.close(); }
        }
        
        ofstream outputFileHandle;
        if (!params->subsample) {  params->util.openOutputFile(params->sumFile, outputFileHandle); }
        
        vector<SharedRAbundVector*> subset;
        for (int k = params->start; k < params->end; k++) { // pass cdd each set of groups to compare
            
            if (params->m->getControl_pressed()) { break; }
            
            for (int l = 0; l < k; l++) {
                
                subset.clear(); //clear out old pair of sharedrabunds
                //add new pair of sharedrabunds
                subset.push_back(params->thisLookup[k]); subset.push_back(params->thisLookup[l]);
                
                if (!params->subsample) {
                    outputFileHandle << params->thisLookup[0]->getLabel() << '\t';
                    //sort groups to be alphanumeric
                    if (params->thisLookup[k]->getGroup() > params->thisLookup[l]->getGroup()) {
                        outputFileHandle << (params->thisLookup[l]->getGroup() +'\t' + params->thisLookup[k]->getGroup()) << '\t'; //print out groups
                    }else{
                        outputFileHandle << (params->thisLookup[k]->getGroup() +'\t' + params->thisLookup[l]->getGroup()) << '\t'; //print out groups
                    }
                }
                
                for(int i=0;i<sumCalculators.size();i++) {
                    
                    //if this calc needs all groups to calculate the pair load all groups
                    if (sumCalculators[i]->getNeedsAll()) {
                        //load subset with rest of lookup for those calcs that need everyone to calc for a pair
                        for (int w = 0; w < params->thisLookup.size(); w++) {
                            if ((w != k) && (w != l)) { subset.push_back(params->thisLookup[w]); }
                        }
                    }
                    
                    vector<double> tempdata = sumCalculators[i]->getValues(subset); //saves the calculator outputs
                    
                    if (params->m->getControl_pressed()) { break; }
                    
                    if (!params->subsample) {
                        outputFileHandle << '\t';
                        sumCalculators[i]->print(outputFileHandle);
                    }
                    seqDist temp(l, k, tempdata[0]);
                    params->calcDists[i].push_back(temp);
                }
                if (!params->subsample) { outputFileHandle << endl; }
            }
        }
        
        if (!params->subsample) {  outputFileHandle.close(); }
        for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "SummarySharedCommand", "driverSummaryShared");
        exit(1);
    }
}

/***********************************************************/
int SummarySharedCommand::runCalcs(SharedRAbundVectors*& thisItersLookup, string sumFileName, string sumAllFile, vector< vector<seqDist>  >& calcDists) {
    try{
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<summarySharedData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            // Allocate memory for thread data.
            string extension = toString(i+1) + ".temp";
            summarySharedData* dataBundle = new summarySharedData(sumFileName+extension, sumAllFile+extension, m, lines[i+1].start, lines[i+1].end, Estimators, thisItersLookup, false, mult, subsample);
            
            data.push_back(dataBundle);
            workerThreads.push_back(new std::thread(driverSummaryShared, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        //SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisItersLookup);
        string extension = toString(0) + ".temp";
        summarySharedData* dataBundle = new summarySharedData(sumFileName+extension, sumAllFile+extension, m, lines[0].start, lines[0].end, Estimators, thisItersLookup, true, mult, subsample);
        
        driverSummaryShared(dataBundle);
        for (int k = 0; k < calcDists.size(); k++) {
            int size = dataBundle->calcDists[k].size();
            for (int j = 0; j < size; j++) {    calcDists[k].push_back(dataBundle->calcDists[k][j]);    }
        }
        
        if (!subsample) {
            util.appendFiles((sumFileName + extension), sumFileName);
            util.mothurRemove((sumFileName + extension));
            if (mult) { util.appendFiles((sumAllFile + extension), sumAllFile);  util.mothurRemove((sumAllFile + extension)); }
        }
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            if (!subsample) {
                string extension = toString(i+1) + ".temp";
                util.appendFiles((sumFileName + extension), sumFileName);
                util.mothurRemove((sumFileName + extension));
            }
            
            if (createPhylip) {
                for (int k = 0; k < calcDists.size(); k++) {
                    int size = data[i]->calcDists[k].size();
                    for (int j = 0; j < size; j++) {    calcDists[k].push_back(data[i]->calcDists[k][j]);    }
                }
            }
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        return 0;
         
    }
    catch(exception& e) {
        m->errorOut(e, "SummarySharedCommand", "runCalcs");
        exit(1);
    }
}
/***********************************************************/
int SummarySharedCommand::process(SharedRAbundVectors* thisLookup, string sumFileName, string sumAllFileName, vector<string> currentLabels) {
	try {
        vector< vector< vector<seqDist> > > calcDistsTotals;  //each iter, one for each calc, then each groupCombos dists. this will be used to make .dist files
        vector< vector<seqDist>  > calcDists; calcDists.resize(numCalcs);
        
        SubSample sample;
        for (int thisIter = 0; thisIter < iters; thisIter++) {
            
            SharedRAbundVectors* thisItersLookup = new SharedRAbundVectors(*thisLookup);
            
            if (subsample) { //we want the summary results for the whole dataset, then the subsampling
                //make copy of lookup so we don't get access violations
                SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisItersLookup);
                
                if (withReplacement)    {  currentLabels = sample.getSampleWithReplacement(newLookup, subsampleSize);    }
                else                    {  currentLabels = sample.getSample(newLookup, subsampleSize);                   }
                
                delete thisItersLookup;
                thisItersLookup = newLookup;
            }
            
            runCalcs(thisItersLookup, sumFileName, sumAllFileName, calcDists);
            
            if (subsample) { //we want the summary results for the whole dataset, then the subsampling
                calcDistsTotals.push_back(calcDists);
                delete thisItersLookup;
            }else {
                if (createPhylip) {
                    for (int i = 0; i < calcDists.size(); i++) {
                        if (m->getControl_pressed()) { break; }
                        
                        //initialize matrix
                        vector< vector<double> > matrix; //square matrix to represent the distance
                        matrix.resize(thisLookup->size());
                        vector<string> GroupNames = thisLookup->getNamesGroups();
                        for (int k = 0; k < thisLookup->size(); k++) {  matrix[k].resize(thisLookup->size(), 0.0);  }
                        
                        for (int j = 0; j < calcDists[i].size(); j++) {
                            int row = calcDists[i][j].seq1;
                            int column = calcDists[i][j].seq2;
                            double dist = calcDists[i][j].dist;
                            
                            matrix[row][column] = dist;
                            matrix[column][row] = dist;
                        }
                        
                        map<string, string> variables; 
                        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
                        variables["[calc]"] = sumCalculatorsNames[i];
                        variables["[distance]"] = thisLookup->getLabel();
                        variables["[outputtag]"] = output;
                        variables["[tag2]"] = "";
                        string distFileName = getOutputFileName("phylip",variables);
                        outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                        ofstream outDist;
                        util.openOutputFile(distFileName, outDist);
                        outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
                        
                        printSims(outDist, matrix, GroupNames);
                        
                        outDist.close();
                    }
                }
            }
            for (int i = 0; i < calcDists.size(); i++) {  calcDists[i].clear(); }
		}

        if (subsample) {
            //we need to find the average distance and standard deviation for each groups distance
            vector< vector<seqDist>  > calcAverages = util.getAverages(calcDistsTotals);
            
            //find standard deviation
            vector< vector<seqDist>  > stdDev = util.getStandardDeviation(calcDistsTotals, calcAverages); 
            
            //print results
            for (int i = 0; i < calcDists.size(); i++) {
                vector< vector<double> > matrix; //square matrix to represent the distance
                matrix.resize(thisLookup->size());
                for (int k = 0; k < thisLookup->size(); k++) {  matrix[k].resize(thisLookup->size(), 0.0); }
                
                vector< vector<double> > stdmatrix; //square matrix to represent the stdDev
                stdmatrix.resize(thisLookup->size());
                vector<string> GroupNames = thisLookup->getNamesGroups();
                for (int k = 0; k < thisLookup->size(); k++) {  stdmatrix[k].resize(thisLookup->size(), 0.0);  }
                
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    int row = calcAverages[i][j].seq1;
                    int column = calcAverages[i][j].seq2;
                    float dist = calcAverages[i][j].dist;
                    float stdDist = stdDev[i][j].dist;
                    
                    matrix[row][column] = dist;
                    matrix[column][row] = dist;
                    stdmatrix[row][column] = stdDist;
                    stdmatrix[column][row] = stdDist;
                }
                
                map<string, string> variables; 
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
                variables["[calc]"] = sumCalculatorsNames[i];
                variables["[distance]"] = thisLookup->getLabel();
                variables["[outputtag]"] = output;
                variables["[tag2]"] = "ave";
                string distFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                ofstream outAve;
                util.openOutputFile(distFileName, outAve);
                outAve.setf(ios::fixed, ios::floatfield); outAve.setf(ios::showpoint);
                
                printSims(outAve, matrix, GroupNames); outAve.close();
                
                variables["[tag2]"] = "std";
                distFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                ofstream outSTD;
                util.openOutputFile(distFileName, outSTD);
                outSTD.setf(ios::fixed, ios::floatfield); outSTD.setf(ios::showpoint);
                
                printSims(outSTD, stdmatrix, GroupNames); outSTD.close();
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


