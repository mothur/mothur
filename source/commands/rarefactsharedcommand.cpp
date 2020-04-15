/*
 *  rarefactsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactsharedcommand.h"
#include "sharedsobs.h"
#include "sharednseqs.h"

#include "subsample.h"


//**********************************************************************************************************************
vector<string> RareFactSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pshared);
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pdesign);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pfreq);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pcalc("calc", "Multiple", "sharednseqs-sharedobserved", "sharedobserved", "", "", "","",true,false,true); parameters.push_back(pcalc);
        CommandParameter psubsampleiters("subsampleiters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(psubsampleiters);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
		CommandParameter pjumble("jumble", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pjumble);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
		CommandParameter pgroupmode("groupmode", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pgroupmode);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;    allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["sharedrarefaction"] = tempOutNames;
        outputTypes["sharedr_nseqs"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RareFactSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The rarefaction.shared command parameters are shared, design, label, iters, groups, sets, jumble, groupmode, processors and calc.  shared is required if there is no current sharedfile. \n";
        helpString += "The design parameter allows you to assign your groups to sets. If provided mothur will run rarefaction.shared on a per set basis. \n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
		helpString += "The rarefaction command should be in the following format: \n";
		helpString += "rarefaction.shared(label=yourLabel, iters=yourIters, calc=yourEstimators, jumble=yourJumble, groups=yourGroups).\n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += "Example rarefaction.shared(label=unique-0.01-0.03,  iters=10000, groups=B-C, jumble=T, calc=sharedobserved).\n";
		helpString += "The default values for iters is 1000, freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness.\n";
        helpString += "The subsampleiters parameter allows you to choose the number of times you would like to run the subsample.\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group.\n";
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
		helpString += "The default value for groups is all the groups in your groupfile, and jumble is true.\n";
		helpString += validCalculator.printCalc("sharedrarefaction");
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RareFactSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "sharedrarefaction") {  pattern = "[filename],shared.rarefaction"; }
        else if (type == "sharedr_nseqs") {  pattern = "[filename],shared.r_nseqs"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RareFactSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RareFactSharedCommand::RareFactSharedCommand(string option)  {
	try {
        if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required.\n");  abort = true; }
			}else { current->setSharedFile(sharedfile); }
            
            designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { abort = true; designfile = ""; }
			else if (designfile == "not found") {  	designfile = "";	}
			else { current->setDesignFile(designfile); }
			
			
			 
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
			if (calc == "not found") { calc = "sharedobserved";  }
			else { 
				 if (calc == "default")  {  calc = "sharedobserved";  }
			}
			util.splitAtDash(calc, Estimators);
			if (util.inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
            
            string sets = validParameter.valid(parameters, "sets");			
			if (sets == "not found") { sets = ""; }
			else { 
				util.splitAtDash(sets, Sets);
                if (Sets.size() != 0) { if (Sets[0] != "all") { Sets.clear(); } }
			}
			
			string temp;
			temp = validParameter.valid(parameters, "freq");			if (temp == "not found") { temp = "100"; }
			util.mothurConvert(temp, freq); 
			
			temp = validParameter.valid(parameters, "iters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, nIters); 
			
			temp = validParameter.valid(parameters, "jumble");			if (temp == "not found") { temp = "T"; }
			if (util.isTrue(temp)) { jumble = true; }
			else { jumble = false; }
            
            temp = validParameter.valid(parameters, "groupmode");		if (temp == "not found") { temp = "T"; }
			groupMode = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "subsampleiters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
            
            temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);

            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
            
            if (subsample == false) { iters = 1; }
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
				
		}

	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "RareFactSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RareFactSharedCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        DesignMap designMap;
        if (designfile == "") { //fake out designMap to run with process
            process(designMap, "");
        }else {
            designMap.read(designfile);
            
            if (Sets.size() == 0) {  Sets = designMap.getCategory(); }
        
            for (int i = 0; i < Sets.size(); i++) { process(designMap, Sets[i]); }
            
            if (groupMode) { outputNames = createGroupFile(outputNames); }
        }
                    
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int RareFactSharedCommand::process(DesignMap& designMap, string thisSet){
	try {
        Rarefact* rCurve;
        vector<Display*> rDisplays;
        
        InputData input(sharedfile, "sharedfile", Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
		SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        if (lookup->size() < 2) { m->mothurOut("[ERROR]: I cannot run the command without at least 2 valid groups.");  delete lookup; return 0; }
        
        string fileNameRoot = outputdir + util.getRootName(util.getSimpleName(sharedfile));
        
        vector<string> newGroups = lookup->getNamesGroups();
        if (thisSet != "") {  //make groups only filled with groups from this set so that's all inputdata will read
            vector<string> thisSets; thisSets.push_back(thisSet);
            newGroups = designMap.getNamesGroups(thisSets);
            fileNameRoot += thisSet + ".";
        }
        
        SharedRAbundVectors* subset = new SharedRAbundVectors();
        subset->setLabels(lookup->getLabel());
        subset->setOTUNames(lookup->getOTUNames());
        vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
        if (thisSet != "") {//remove unwanted groups
            for (int i = 0; i < data.size(); i++) { if (util.inUsersGroups(data[i]->getGroup(), newGroups)) { subset->push_back(data[i]); } }
            subset->eliminateZeroOTUS();
        }else { for (int i = 0; i < data.size(); i++) {  subset->push_back(data[i]); } }
        
        /******************************************************/
        if (subsample) {
            //user has not set size, set size = smallest samples size
            if (subsampleSize == -1) {
                subsampleSize = subset->getNumSeqsSmallestGroup();
                m->mothurOut("Setting subsample size to " + toString(subsampleSize) + ".\n\n");
            }
            
            subset->removeGroups(subsampleSize);
            newGroups = subset->getNamesGroups();
            
            if (subset->size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command.\n\n"); m->setControl_pressed(true); return 0; }
        }
        /******************************************************/
        
        map<string, string> variables; 
        variables["[filename]"] = fileNameRoot;
		ValidCalculators validCalculator;
		for (int i=0; i<Estimators.size(); i++) {
			if (validCalculator.isValidCalculator("sharedrarefaction", Estimators[i]) ) { 
				if (Estimators[i] == "sharedobserved") { 
					rDisplays.push_back(new RareDisplay(new SharedSobs(), new SharedThreeColumnFile(getOutputFileName("sharedrarefaction",variables), "")));
					outputNames.push_back(getOutputFileName("sharedrarefaction",variables)); outputTypes["sharedrarefaction"].push_back(getOutputFileName("sharedrarefaction",variables));
				}else if (Estimators[i] == "sharednseqs") { 
					rDisplays.push_back(new RareDisplay(new SharedNSeqs(), new SharedThreeColumnFile(getOutputFileName("sharedr_nseqs",variables), "")));
					outputNames.push_back(getOutputFileName("sharedr_nseqs",variables)); outputTypes["sharedr_nseqs"].push_back(getOutputFileName("sharedr_nseqs",variables));
				}
			}
            file2Group[outputNames.size()-1] = thisSet;
		}
		
		//if the users entered no valid calculators don't execute command
        if (rDisplays.size() == 0) { delete lookup;  delete subset; return 0; }
		
		if (m->getControl_pressed()) { 
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
			for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}
			delete lookup; delete subset; return 0;
		}
        
        while (subset != NULL) {
            
            if (m->getControl_pressed()) { delete subset; delete lookup; break; }
            
            rCurve = new Rarefact(subset, rDisplays, jumble, processors);
            rCurve->getSharedCurve(freq, nIters);
            delete rCurve;
            
            if (subsample) { subsampleLookup(subset, fileNameRoot);  }
            
           delete lookup; delete subset;
           lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
           
           if (lookup != NULL) {
               subset = new SharedRAbundVectors();
               data = lookup->getSharedRAbundVectors();
               if (thisSet != "") {//remove unwanted groups
                   for (int i = 0; i < data.size(); i++) { if (util.inUsersGroups(data[i]->getGroup(), newGroups)) { subset->push_back(data[i]); } }
                   subset->eliminateZeroOTUS();
               }else { for (int i = 0; i < data.size(); i++) {  subset->push_back(data[i]); } }
           }else {  subset = NULL; }
        }
        
		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
        
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {    util.mothurRemove(outputNames[i]);     } }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
int RareFactSharedCommand::subsampleLookup(SharedRAbundVectors*& thisLookup, string fileNameRoot) {
	try {
        
        map<string, vector<string> > filenames;
        SubSample sample;
        for (int thisIter = 0; thisIter < iters; thisIter++) {
            
            SharedRAbundVectors* thisItersLookup = new SharedRAbundVectors(*thisLookup);
            
            if (withReplacement)    {  sample.getSampleWithReplacement(thisItersLookup, subsampleSize);     }
            else                    {  sample.getSample(thisItersLookup, subsampleSize);                    }

            string thisfileNameRoot = fileNameRoot + toString(thisIter);
            map<string, string> variables; 
            variables["[filename]"] = thisfileNameRoot;
            
            vector<Display*> rDisplays;
            ValidCalculators validCalculator;
            for (int i=0; i<Estimators.size(); i++) {
                if (validCalculator.isValidCalculator("sharedrarefaction", Estimators[i]) ) { 
                    if (Estimators[i] == "sharedobserved") { 
                        rDisplays.push_back(new RareDisplay(new SharedSobs(), new SharedThreeColumnFile(getOutputFileName("sharedrarefaction",variables), "")));
                        filenames["sharedrarefaction"].push_back(getOutputFileName("sharedrarefaction",variables));
                    }else if (Estimators[i] == "sharednseqs") { 
                        rDisplays.push_back(new RareDisplay(new SharedNSeqs(), new SharedThreeColumnFile(getOutputFileName("sharedr_nseqs",variables), "")));
                        filenames["sharedr_nseqs"].push_back(getOutputFileName("sharedr_nseqs",variables));
                    }
                }
            }
            
            Rarefact rCurve(thisItersLookup, rDisplays, jumble, processors);
			rCurve.getSharedCurve(freq, nIters);
            
            //clean up memory
            for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
            delete thisItersLookup;
            
            if((thisIter+1) % 100 == 0)	{ m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 	}
        }
        
        //create std and ave outputs
        vector< vector< vector< double > > > results; //iter -> numSampled -> data
        for (map<string, vector<string> >::iterator it = filenames.begin(); it != filenames.end(); it++) {
            vector<string> thisTypesFiles = it->second;
            vector<string> columnHeaders;
            for (int i = 0; i < thisTypesFiles.size(); i++) {
                ifstream in;
                util.openInputFile(thisTypesFiles[i], in);
                
                string headers = util.getline(in); util.gobble(in);
                columnHeaders = util.splitWhiteSpace(headers);
                int numCols = columnHeaders.size();
                
                vector<vector<double> > thisFilesLines;
                while (!in.eof()) {
                    if (m->getControl_pressed()) { break; }
                    vector<double> data; data.resize(numCols, 0);
                    //read numSampled line
                    for (int j = 0; j < numCols; j++) { in >> data[j]; util.gobble(in); }
                    thisFilesLines.push_back(data);
                }
                in.close();
                results.push_back(thisFilesLines);
                util.mothurRemove(thisTypesFiles[i]);
            }
            
            if (!m->getControl_pressed()) {
                //process results
                map<string, string> variables; variables["[filename]"] = fileNameRoot + "ave-std." + thisLookup->getLabel() + ".";

                string outputFile = getOutputFileName(it->first,variables);
                ofstream out;
                util.openOutputFile(outputFile, out);
                outputNames.push_back(outputFile); outputTypes[it->first].push_back(outputFile);
                
                out << columnHeaders[0] << '\t' << "method";
                for (int i = 1; i < columnHeaders.size(); i++) { out  << '\t' << columnHeaders[i]; }
                out << endl;
            
                vector< vector<double> > aveResults; aveResults.resize(results[0].size());
                for (int i = 0; i < aveResults.size(); i++) { aveResults[i].resize(results[0][i].size(), 0.0); }
                
                for (int thisIter = 0; thisIter < iters; thisIter++) { //sum all groups dists for each calculator
                    for (int i = 0; i < aveResults.size(); i++) {  //initialize sums to zero.
                        aveResults[i][0] = results[thisIter][i][0];
                        for (int j = 1; j < aveResults[i].size(); j++) {
                            aveResults[i][j] += results[thisIter][i][j];
                        }
                    }
                }
                
                for (int i = 0; i < aveResults.size(); i++) {  //finds average.
                    for (int j = 1; j < aveResults[i].size(); j++) {
                        aveResults[i][j] /= (float) iters;
                    }
                }
                
                //standard deviation
                vector< vector<double> > stdResults; stdResults.resize(results[0].size());
                for (int i = 0; i < stdResults.size(); i++) { stdResults[i].resize(results[0][i].size(), 0.0); }
                
                for (int thisIter = 0; thisIter < iters; thisIter++) { //compute the difference of each dist from the mean, and square the result of each
                    for (int i = 0; i < stdResults.size(); i++) {  
                        stdResults[i][0] = aveResults[i][0];
                        for (int j = 1; j < stdResults[i].size(); j++) {
                            stdResults[i][j] += ((results[thisIter][i][j] - aveResults[i][j]) * (results[thisIter][i][j] - aveResults[i][j]));
                        }
                    }
                }
                
                for (int i = 0; i < stdResults.size(); i++) {  //finds average.
                    out << aveResults[i][0] << '\t' << "ave";
                    for (int j = 1; j < aveResults[i].size(); j++) { out  << '\t' << aveResults[i][j]; }
                    out << endl;
                    out << stdResults[i][0] << '\t' << "std";
                    for (int j = 1; j < stdResults[i].size(); j++) {
                        stdResults[i][j] /= (float) iters;
                        stdResults[i][j] = sqrt(stdResults[i][j]);
                        out << '\t' << stdResults[i][j];
                    }
                    out << endl;
                }
                out.close();
            }
        }
        
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "subsample");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactSharedCommand::createGroupFile(vector<string>& outputNames) {
	try {
		
		vector<string> newFileNames;
		
		//find different types of files
		map<string, map<string, string> > typesFiles;
        map<string, vector< vector<string> > > fileLabels; //combofile name to labels. each label is a vector because it may be unique lci hci.
        vector<string> groupNames;
		for (int i = 0; i < outputNames.size(); i++) {
            
			string extension = util.getExtension(outputNames[i]);
            string combineFileName = outputdir + util.getRootName(util.getSimpleName(sharedfile)) + "groups" + extension;
			util.mothurRemove(combineFileName); //remove old file
            
			ifstream in;
			util.openInputFile(outputNames[i], in);
			
            string labels = util.getline(in); util.gobble(in);
            vector<string> theseLabels = util.splitWhiteSpace(labels);
            
            vector< vector<string> > allLabels;
            vector<string> thisSet; thisSet.push_back(theseLabels[0]); allLabels.push_back(thisSet); thisSet.clear(); //makes "numSampled" its own grouping
            for (int j = 1; j < theseLabels.size()-1; j++) {
                if (theseLabels[j+1] == "lci") {
                    thisSet.push_back(theseLabels[j]); 
                    thisSet.push_back(theseLabels[j+1]); 
                    thisSet.push_back(theseLabels[j+2]);
                    j++; j++;
                }else{ //no lci or hci for this calc.
                    thisSet.push_back(theseLabels[j]); 
                }
                allLabels.push_back(thisSet); 
                thisSet.clear();
            }
            fileLabels[combineFileName] = allLabels;
            
            map<string, map<string, string> >::iterator itfind = typesFiles.find(extension);
            if (itfind != typesFiles.end()) {
                (itfind->second)[outputNames[i]] = file2Group[i];
            }else {
                map<string, string> temp;  
                temp[outputNames[i]] = file2Group[i];
                typesFiles[extension] = temp;
            }
            if (!(util.inUsersGroups(file2Group[i], groupNames))) {  groupNames.push_back(file2Group[i]); }
		}
		
		//for each type create a combo file
		
		for (map<string, map<string, string> >::iterator it = typesFiles.begin(); it != typesFiles.end(); it++) {
			
			ofstream out;
			string combineFileName = outputdir + util.getRootName(util.getSimpleName(sharedfile)) + "groups" + it->first;
			util.openOutputFileAppend(combineFileName, out);
			newFileNames.push_back(combineFileName);
			map<string, string> thisTypesFiles = it->second; //it->second maps filename to group
            set<int> numSampledSet;
            
			//open each type summary file
			map<string, map<int, vector< vector<string> > > > files; //maps file name to lines in file
			int maxLines = 0;
			for (map<string, string>::iterator itFileNameGroup = thisTypesFiles.begin(); itFileNameGroup != thisTypesFiles.end(); itFileNameGroup++) {
                
                string thisfilename = itFileNameGroup->first;
                string group = itFileNameGroup->second;
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + thisfilename + "\t" + group + "\n");  }
                
                ifstream temp;
                util.openInputFile(thisfilename, temp);
                
                //read through first line - labels
                string dummy = util.getline(temp);	util.gobble(temp);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + dummy + "\t" + toString(fileLabels[combineFileName].size()) + "\n");  }
				
				map<int, vector< vector<string> > > thisFilesLines;
				while (!temp.eof()){
                    float numSampled = 0;
                    string thisLineInfo = util.getline(temp); util.gobble(temp);
                    vector<string> parsedLine = util.splitWhiteSpace(thisLineInfo);
                    util.mothurConvert(parsedLine[0], numSampled);
                    
                    vector< vector<string> > theseReads;
                    vector<string> thisSet; thisSet.push_back(toString(numSampled)); theseReads.push_back(thisSet); thisSet.clear();
                    int columnIndex = 1; //0 -> numSampled, 1 -> 0.03, 2 -> 0.03lci, 3 -> 0.03hci, 4 -> 0.05, 5 -> 0.05lci, 6 -> 0.05hci
                    for (int k = 1; k < fileLabels[combineFileName].size(); k++) { //output thing like 0.03-A lci-A hci-A
                        vector<string> reads;
                        string next = "";
                        int numColumnsPerLabel = fileLabels[combineFileName][k].size();  // 0.03 lci hci  ... 0.05 lci hci -> 3 columns
                        for (int l = 0; l < numColumnsPerLabel; l++) {
                            reads.push_back(parsedLine[columnIndex]); columnIndex++;
                        }
                        theseReads.push_back(reads);
                        
                        if (m->getDebug()) { m->mothurOut("[DEBUG]: " + util.getStringFromVector(reads, " ") + "\n");  }
                    }
                    thisFilesLines[numSampled] = theseReads;
                    util.gobble(temp);
                    
                    numSampledSet.insert(numSampled);
				}
				
				files[group] = thisFilesLines;
				
				//save longest file for below
				if (maxLines < thisFilesLines.size()) { maxLines = thisFilesLines.size(); }
				
				temp.close();
				util.mothurRemove(thisfilename);
			}
			
            //output new labels line
            out << fileLabels[combineFileName][0][0];
            for (int k = 1; k < fileLabels[combineFileName].size(); k++) { //output thing like 0.03-A lci-A hci-A
                for (int n = 0; n < groupNames.size(); n++) { // for each group
                    for (int l = 0; l < fileLabels[combineFileName][k].size(); l++) { //output modified labels
                        out  << '\t' << fileLabels[combineFileName][k][l] << '-' << groupNames[n];
                    }
                }
            }
			out << endl;
            
			//for each label
			for (set<int>::iterator itNumSampled = numSampledSet.begin(); itNumSampled != numSampledSet.end(); itNumSampled++) {
				
                out << (*itNumSampled);
                
                if (m->getControl_pressed()) { break; }
                
                for (int k = 1; k < fileLabels[combineFileName].size(); k++) { //each chunk
				    //grab data for each group
                    for (map<string, map<int, vector< vector<string> > > >::iterator itFileNameGroup = files.begin(); itFileNameGroup != files.end(); itFileNameGroup++) {
                        
                        string group = itFileNameGroup->first;
                        
                        map<int, vector< vector<string> > >::iterator itLine = files[group].find(*itNumSampled);
                        if (itLine != files[group].end()) { 
                            for (int l = 0; l < (itLine->second)[k].size(); l++) { 
                                out  << '\t' << (itLine->second)[k][l];
                                
                            }                             
                        }else { 
                            for (int l = 0; l < fileLabels[combineFileName][k].size(); l++) { 
                                out << "\tNA";
                            } 
                        }
                    }
                }
                out << endl;
			}	
			out.close();
		}
		
		//return combine file name
		return newFileNames;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "createGroupFile");
		exit(1);
	}
}
//**********************************************************************************************************************
