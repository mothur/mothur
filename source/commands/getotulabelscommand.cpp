//
//  getotulabelscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/21/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "getotulabelscommand.h"

//**********************************************************************************************************************
vector<string> GetOtuLabelsCommand::setParameters(){	
	try {
        CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,true, true); parameters.push_back(paccnos);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "FNGLT", "none","constaxonomy",false,false, true); parameters.push_back(pconstaxonomy);
        CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false, true); parameters.push_back(plist);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "FNGLT", "none","shared",false,false, true); parameters.push_back(pshared);
		CommandParameter potucorr("otucorr", "InputTypes", "", "", "none", "FNGLT", "none","otucorr",false,false, true); parameters.push_back(potucorr);
        CommandParameter pcorraxes("corraxes", "InputTypes", "", "", "none", "FNGLT", "none","corraxes",false,false, true); parameters.push_back(pcorraxes);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetOtuLabelsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.otus command can be used to select specific otus with the output from classify.otu, otu.association, or corr.axes commands.  It can also be used to select a set of otus from a shared or list file.\n";
		helpString += "The get.otus parameters are: constaxonomy, otucorr, corraxes, shared, list, label and accnos.\n";
		helpString += "The constaxonomy parameter is used to input the results of the classify.otu command.\n";
        helpString += "The otucorr parameter is used to input the results of the otu.association command.\n";
        helpString += "The corraxes parameter is used to input the results of the corr.axes command.\n";
        helpString += "The label parameter is used to analyze specific labels in your input. \n";
		helpString += "The get.otus commmand should be in the following format: \n";
		helpString += "get.otus(accnos=yourListOfOTULabels, corraxes=yourCorrAxesFile)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetOtuLabelsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "constaxonomy")         {   pattern = "[filename],pick,[extension]";    }
        else if (type == "otucorr")         {   pattern = "[filename],pick,[extension]";    }
        else if (type == "corraxes")        {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")            {   pattern = "[filename],[distance],pick,[extension]";    }
        else if (type == "shared")          {   pattern = "[filename],[distance],pick,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetOtuLabelsCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
GetOtuLabelsCommand::GetOtuLabelsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["constaxonomy"] = tempOutNames; 
        outputTypes["otucorr"] = tempOutNames;
        outputTypes["corraxes"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "GetOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetOtuLabelsCommand::GetOtuLabelsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
        
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
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                //edit file types below to include only the types you added as parameters
                
				string path;
                it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("corraxes");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["corraxes"] = inputDir + it->second;		}
				}
                
                it = parameters.find("otucorr");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["otucorr"] = inputDir + it->second;		}
				}
                
                it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
            
            vector<string> tempOutNames;
            outputTypes["constaxonomy"] = tempOutNames; 
            outputTypes["otucorr"] = tempOutNames;
            outputTypes["corraxes"] = tempOutNames;
            outputTypes["shared"] = tempOutNames;
            outputTypes["list"] = tempOutNames;
            
 			//check for parameters
            accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = current->getAccnosFile(); 
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You have no valid accnos file and accnos is required."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else { current->setAccnosFile(accnosfile); }	
			
			constaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
			if (constaxonomyfile == "not open") { constaxonomyfile = ""; abort = true; }
			else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
            
            corraxesfile = validParameter.validFile(parameters, "corraxes");
			if (corraxesfile == "not open") { corraxesfile = ""; abort = true; }
			else if (corraxesfile == "not found") {  corraxesfile = "";  }
            
            otucorrfile = validParameter.validFile(parameters, "otucorr");
			if (otucorrfile == "not open") { otucorrfile = ""; abort = true; }
			else if (otucorrfile == "not found") {  otucorrfile = "";  }
            
            listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
            else { current->setListFile(listfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
            else { current->setSharedFile(sharedfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	 outputDir = ""; 	}
            
            if ((constaxonomyfile == "") && (corraxesfile == "") && (otucorrfile == "") && (sharedfile == "") && (listfile == ""))  { m->mothurOut("You must provide one of the following: constaxonomy, corraxes, otucorr, shared or list."); m->mothurOutEndLine(); abort = true; }
            
            if ((sharedfile != "") || (listfile != "")) {
                label = validParameter.valid(parameters, "label");			
                if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }
            }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "GetOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetOtuLabelsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        //get labels you want to keep
		labels = util.readAccnos(accnosfile);
        //simplfy labels
        set<string> newLabels;
        for (set<string>::iterator it = labels.begin(); it != labels.end(); it++) {  newLabels.insert(util.getSimpleLabel(*it)); }
        labels = newLabels;
        
		if (m->getControl_pressed()) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (constaxonomyfile != "")	{		readClassifyOtu();      }
		if (corraxesfile != "")		{		readCorrAxes();         }
		if (otucorrfile != "")		{		readOtuAssociation();	}
        if (listfile != "")         {		readList();             }
        if (sharedfile != "")		{		readShared();           }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
        string currentName = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }
        
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOtuLabelsCommand::readClassifyOtu(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(constaxonomyfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxonomyfile));
        variables["[extension]"] = util.getExtension(constaxonomyfile);
		string outputFileName = getOutputFileName("constaxonomy", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(constaxonomyfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        //read headers
        string headers = util.getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu >> size; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
            
            if (m->getDebug()) { m->mothurOut("Otu=" + otu + ", size=" + toString(size) + ", tax=" + tax + "\n"); }
            
            if (labels.count(util.getSimpleLabel(otu)) != 0) {
				wroteSomething = true;
				selectedCount++;
                
                out << otu << '\t' << size << '\t' << tax << endl;
            }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["constaxonomy"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " otus from your constaxonomy file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "readClassifyOtu");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOtuLabelsCommand::readOtuAssociation(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(otucorrfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(otucorrfile));
        variables["[extension]"] = util.getExtension(otucorrfile);
		string outputFileName = getOutputFileName("otucorr", variables);

		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(otucorrfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        //read headers
        string headers = util.getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string otu1 = ""; 
            string otu2 = ""; 
            in >> otu1 >> otu2;
            string line = util.getline(in); util.gobble(in);
            
            if ((labels.count(util.getSimpleLabel(otu1)) != 0) && (labels.count(util.getSimpleLabel(otu2)) != 0)){
				wroteSomething = true;
				selectedCount++;
                
                out << otu1 << '\t' << otu2 << '\t' << line << endl;
            }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["otucorr"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " lines from your otu.corr file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "readOtuAssociation");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOtuLabelsCommand::readCorrAxes(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(corraxesfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(corraxesfile));
        variables["[extension]"] = util.getExtension(corraxesfile);
		string outputFileName = getOutputFileName("corraxes", variables);

		ofstream out;
		util.openOutputFile(outputFileName, out);
		
        
		ifstream in;
		util.openInputFile(corraxesfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        //read headers
        string headers = util.getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string otu = ""; 
            in >> otu;
            string line = util.getline(in); util.gobble(in);
            
            if (labels.count(util.getSimpleLabel(otu)) != 0) {
				wroteSomething = true;
				selectedCount++;
                
                out << otu << '\t' << line << endl;
            }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["corraxes"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " lines from your corr.axes file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "readCorrAxes");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOtuLabelsCommand::readShared(){
	try {
        
        SharedRAbundVectors* lookup = getShared();
        
        if (m->getControl_pressed()) { delete lookup; return 0; }
          
        vector<string> newLabels;
        
        bool wroteSomething = false;
        int numSelected = 0;
        for (int i = 0; i < lookup->getNumBins();) {
            
            if (m->getControl_pressed()) { delete lookup; return 0; }
            
            //is this otu on the list
            if (labels.count(util.getSimpleLabel(lookup->getOTUNames()[i])) != 0) {
                numSelected++; wroteSomething = true;
                ++i;
            }else { lookup->removeOTU(i);  }
        }
    
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
        variables["[distance]"] = lookup->getLabel();
		string outputFileName = getOutputFileName("shared", variables); 
        ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		lookup->printHeaders(out);
        lookup->print(out);
		out.close();
        
        delete lookup;
        
        if (wroteSomething == false) { m->mothurOut("Your file does not contain any OTUs from the .accnos file."); m->mothurOutEndLine();  }

		m->mothurOut("Selected " + toString(numSelected) + " OTUs from your shared file."); m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "readShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOtuLabelsCommand::readList(){
	try {
        getListVector();
        
        if (m->getControl_pressed()) { delete list; return 0;}
        
        ListVector newList;
        newList.setLabel(list->getLabel());
        int selectedCount = 0;
        bool wroteSomething = false;
        
        vector<string> binLabels = list->getLabels();
        vector<string> newLabels;
        for (int i = 0; i < list->getNumBins(); i++) {
            
            if (m->getControl_pressed()) { delete list; return 0;}
            
            if (labels.count(util.getSimpleLabel(binLabels[i])) != 0) {
				selectedCount++;
                newList.push_back(list->get(i));
                newLabels.push_back(binLabels[i]);
            }
        }
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        variables["[distance]"] = list->getLabel();
		string outputFileName = getOutputFileName("list", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		delete list;
        //print new listvector
        if (newList.getNumBins() != 0) {
            wroteSomething = true;
            newList.setLabels(newLabels);
            newList.printHeaders(out);
            newList.print(out, false);
        }
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any OTUs from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " OTUs from your list file."); m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
            m->errorOut(e, "GetOtuLabelsCommand", "readList");
            exit(1);
        }
    }
//**********************************************************************************************************************
int GetOtuLabelsCommand::getListVector(){
	try {
		InputData input(listfile, "list", nullVector);
		list = input.getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {  return 0;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input.getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				break;
			}
			
			lastLabel = list->getLabel();			
			
			//get next line to process
			//prevent memory leak
			delete list; 
			list = input.getListVector();
		}
		
		
		if (m->getControl_pressed()) {  return 0;  }
		
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
		if (needToRun )  {
			delete list; 
			list = input.getListVector(lastLabel);
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "getListVector");	
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundVectors* GetOtuLabelsCommand::getShared(){
	try {
		InputData input(sharedfile, "sharedfile", nullVector);
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup->getLabel();
		
		if (label == "") { label = lastLabel;  return lookup; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {   delete lookup; return NULL;  }
			
			if(labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
                delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
            delete lookup;
			lookup = input.getSharedRAbundVectors();
		}
		
		
		if (m->getControl_pressed()) {  return 0;  }
		
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
		if (needToRun )  {
			delete lookup;
			lookup = input.getSharedRAbundVectors(lastLabel);
		}	
		
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "getShared");	
		exit(1);
	}
}
//**********************************************************************************************************************
