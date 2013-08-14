//
//  lefsecommand.cpp
//  Mothur
//
//  Created by SarahsWork on 6/12/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "lefsecommand.h"
#include "linearalgebra.h"

//**********************************************************************************************************************
vector<string> LefseCommand::setParameters(){
	try {
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pdesign);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pshared);
        CommandParameter pclass("class", "String", "", "", "", "", "","",false,false); parameters.push_back(pclass);
        CommandParameter psubclass("subclass", "String", "", "", "", "", "","",false,false); parameters.push_back(psubclass);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pclasses("classes", "String", "", "", "", "", "","",false,false); parameters.push_back(pclasses);
        CommandParameter palpha("anova_alpha", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(palpha);
        CommandParameter pwalpha("wilcoxon_alpha", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(pwalpha);
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string LefseCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The lefse command allows you to ....\n";
		helpString += "The lefse command parameters are: shared, design, class, subclass, label, classes, wilcoxon_alpha, anova_alpha.\n";
		helpString += "The class parameter is used to indicate the which category you would like used for the Kruskal Wallis analysis. If none is provided first category is used.\n";
        helpString += "The subclass parameter is used to indicate the .....If none is provided second category is used, or if only one category subclass is ignored. \n";
        helpString += "The classes parameter is used to indicate the classes you would like to use. Clases should be inputted in the following format: classes=label<value1|value2|value3>-label<value1|value2>. For example to include groups from treatment with value early or late and age= young or old.  class=treatment<Early|Late>-age<young|old>.\n";
        helpString += "The label parameter is used to indicate which distances in the shared file you would like to use. labels are separated by dashes.\n";
		helpString += "The lefse command should be in the following format: lefse(shared=final.an.shared, design=final.design, class=treatment, subclass=age).\n";
        return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string LefseCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],[distance],lefse_summary"; }
        else if (type == "kruskall-wallis") {  pattern = "[filename],[distance],kruskall_wallis"; }
        else if (type == "wilcoxon") {  pattern = "[filename],[distance],wilcoxon"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
LefseCommand::LefseCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
        outputTypes["kruskall-wallis"] = tempOutNames;
        outputTypes["wilcoxon"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "LefseCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
LefseCommand::LefseCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
            outputTypes["kruskall-wallis"] = tempOutNames;
            outputTypes["wilcoxon"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["desing"] = inputDir + it->second;		}
				}
				
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
                    
            //get shared file, it is required
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile();
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
            
            //get shared file, it is required
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") {
				//if there is a current shared file, use it
				designfile = m->getDesignFile();
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current design file and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setDesignFile(designfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
				outputDir = m->hasPath(sharedfile); //if user entered a file with a path then preserve it
			}
            
            string label = validParameter.validFile(parameters, "label", false);
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
            
            mclass = validParameter.validFile(parameters, "class", false);
			if (mclass == "not found") { mclass = ""; }
			
            subclass = validParameter.validFile(parameters, "subclass", false);
			if (subclass == "not found") { subclass = ""; }
            
            classes = validParameter.validFile(parameters, "classes", false);
			if (classes == "not found") { classes = ""; }
            
            string temp = validParameter.validFile(parameters, "anova_alpha", false);
			if (temp == "not found") { temp = "0.05"; }
			m->mothurConvert(temp, anovaAlpha);
            
            temp = validParameter.validFile(parameters, "wilcoxon_alpha", false);
			if (temp == "not found") { temp = "0.05"; }
			m->mothurConvert(temp, wilcoxonAlpha);

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "LefseCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int LefseCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        DesignMap designMap(designfile);
        //if user set classes set groups=those classes
        if (classes != "") {
            map<string, vector<string> > thisClasses = m->parseClasses(classes);
            vector<string> groups = designMap.getNamesUnique(thisClasses);
            if (groups.size() != 0) { m->setGroups(groups); }
            else { m->mothurOut("[ERROR]: no groups meet your classes requirement, quitting.\n"); return 0; }
        }
        
        //if user did not select class use first column
        if (mclass == "") {  mclass = designMap.getDefaultClass(); m->mothurOut("\nYou did not provide a class, using " + mclass +".\n\n"); }
        
        InputData input(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
        string lastLabel = lookup[0]->getLabel();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
            
            if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
                
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup, designMap);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
            }
            
            if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup[0]->getLabel();
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                lookup = input.getSharedRAbundVectors(lastLabel);
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup, designMap);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
                
                //restore real lastlabel to save below
                lookup[0]->setLabel(saveLabel);
            }
            
            lastLabel = lookup[0]->getLabel();
            //prevent memory leak
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
            
            if (m->control_pressed) { return 0; }
            
            //get next line to process
            lookup = input.getSharedRAbundVectors();
        }
        
        if (m->control_pressed) {  return 0; }
        
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
            for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            process(lookup, designMap);
            
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
        }
                
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int LefseCommand::process(vector<SharedRAbundVector*>& lookup, DesignMap& designMap) {
	try {
        //run kruskal wallis on each otu
        vector<int> significantOtuLabels = runKruskalWallis(lookup, designMap);
        
        //check for subclass
        if (subclass != "") {  significantOtuLabels = runWilcoxon(lookup, designMap, significantOtuLabels);  }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************

vector<int> LefseCommand::runKruskalWallis(vector<SharedRAbundVector*>& lookup, DesignMap& designMap) {
	try {
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
		string outputFileName = getOutputFileName("kruskall-wallis",variables);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["kruskall-wallis"].push_back(outputFileName);
        out << "OTULabel\tKW\tPvalue\n";
        
        vector<int> significantOtuLabels;
        int numBins = lookup[0]->getNumBins();
        //sanity check to make sure each treatment has a group in the shared file
        set<string> treatments;
        for (int j = 0; j < lookup.size(); j++) {
            string group = lookup[j]->getGroup();
            string treatment = designMap.get(group, mclass); //get value for this group in this category
            treatments.insert(treatment);
        }
        if (treatments.size() < 2) { m->mothurOut("[ERROR]: need at least 2 things to classes to compare, quitting.\n"); m->control_pressed = true; }
        
        LinearAlgebra linear;
        for (int i = 0; i < numBins; i++) {
            if (m->control_pressed) { break; }
            
            vector<spearmanRank> values;
            for (int j = 0; j < lookup.size(); j++) {
                string group = lookup[j]->getGroup();
                string treatment = designMap.get(group, mclass); //get value for this group in this category
                spearmanRank temp(treatment, lookup[j]->getAbundance(i));
                values.push_back(temp);
            }
            
            double pValue = 0.0;
            double H = linear.calcKruskalWallis(values, pValue);
            
            //output H and signifigance
            out << m->currentBinLabels[i] << '\t' << H << '\t' << pValue << endl;
            
            if (pValue < anovaAlpha) {  significantOtuLabels.push_back(i);  }
        }
        out.close();
        
        return significantOtuLabels;
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "runKruskalWallis");
		exit(1);
	}
}
//**********************************************************************************************************************
//assumes not neccessarily paired
vector<int> LefseCommand::runWilcoxon(vector<SharedRAbundVector*>& lookup, DesignMap& designMap, vector<int> bins) {
    try {
        LinearAlgebra linear;
        vector<int> significantOtuLabels;
        //if it exists and meets the following requirements run Wilcoxon
        /*
         1. Subclass members all belong to same main class
         anything else
        */
        vector<string> subclasses;
        map<string, string> subclass2Class;
        map<string, int> subclassCounts;
        map<string, vector<int> > subClass2GroupIndex; //maps subclass name to vector of indexes in lookup from that subclass. old -> 1,2,3 means groups in location 1,2,3 of lookup are from old.  Saves time below.
        bool error = false;
        for (int j = 0; j < lookup.size(); j++) {
            string group = lookup[j]->getGroup();
            string treatment = designMap.get(group, mclass); //get value for this group in this category
            string thisSub = designMap.get(group, subclass);
            map<string, string>::iterator it = subclass2Class.find(thisSub);
            if (it == subclass2Class.end()) {
                subclass2Class[thisSub] = treatment;
                subclassCounts[thisSub] = 1;
                vector<int> temp; temp.push_back(j);
                subClass2GroupIndex[thisSub] = temp;
            }
            else {
                subclassCounts[thisSub]++;
                subClass2GroupIndex[thisSub].push_back(j);
                if (it->second != treatment) {
                    error = true;
                    m->mothurOut("[ERROR]: subclass " + thisSub + " has members in " + it->second + " and " + treatment + ". Subclass members must be from the same class. Ignoring wilcoxon.\n");
                }
            }
        }
        
        if (error) { return significantOtuLabels; }
        
        int numBins = lookup[0]->getNumBins();
        vector<compGroup> comp;
        //find comparisons and fill comp
        map<string, int>::iterator itB;
        for(map<string, int>::iterator it=subclassCounts.begin();it!=subclassCounts.end();it++){
            itB = it;itB++;
            for(itB;itB!=subclassCounts.end();itB++){
                compGroup temp(it->first,itB->first);
                comp.push_back(temp);
            }			
        }

        int numComp = comp.size();
        if (numComp < 2) {  m->mothurOut("[ERROR]: Need at least 2 subclasses, Ignoring Wilcoxon.\n");
            return significantOtuLabels;  }
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("wilcoxon",variables);
        
        ofstream out;
        m->openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["wilcoxon"].push_back(outputFileName);
        out << "OTULabel\tComparision\tWilcoxon\tPvalue\n";
        
        for (int i = 0; i < numBins; i++) {
            if (m->control_pressed) { break; }
            
            if (m->inUsersGroups(i, bins)) { //flagged in Kruskal Wallis
                
                bool sig = false;
                //for each subclass comparision
                for (int j = 0; j < numComp; j++) {
                    //fill x and y with this comparisons data
                    vector<double> x; vector<double> y; 
                    
                    cout << m->currentBinLabels[i] << '\t' << comp[j].getCombo() << " x <- (";
                    //fill x and y
                    vector<int> xIndexes = subClass2GroupIndex[comp[j].group1]; //indexes in lookup for this subclass
                    for (int k = 0; k < xIndexes.size(); k++) { x.push_back(lookup[xIndexes[k]]->getAbundance(i)); cout << lookup[xIndexes[k]]->getAbundance(i) << ", "; }
                    cout << ")\n";
                    
                    cout << m->currentBinLabels[i] << '\t' << comp[j].getCombo() << " y <- (";
                    
                    vector<int> yIndexes = subClass2GroupIndex[comp[j].group2]; //indexes in lookup for this subclass
                    for (int k = 0; k < yIndexes.size(); k++) { y.push_back(lookup[yIndexes[k]]->getAbundance(i)); cout << lookup[yIndexes[k]]->getAbundance(i) << ", ";}
                    cout << ")\n";
                    
                    double pValue = 0.0;
                    double H = linear.calcWilcoxon(x, y, pValue);
            
                    //output H and signifigance
                    if (!isnan(pValue)) { out << m->currentBinLabels[i] << '\t' << comp[j].getCombo() << '\t' << H << '\t' << pValue << endl; }
                    else { out << m->currentBinLabels[i] << '\t' << comp[j].getCombo() << '\t' << H << '\t' << "NA" << endl; }
                    
                    //set sig - not sure how yet
                }
                if (sig) {  significantOtuLabels.push_back(i);  }
            }
        }
        out.close();
        
        return significantOtuLabels;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "runWilcoxon");
        exit(1);
    }
}

//**********************************************************************************************************************


