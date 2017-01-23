//
//  removeotulabels.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/21/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "removeotulabelscommand.h"

//**********************************************************************************************************************
vector<string> RemoveOtuLabelsCommand::setParameters(){	
	try {
        CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(paccnos);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "FNGLT", "none","constaxonomy",false,false); parameters.push_back(pconstaxonomy);
		CommandParameter potucorr("otucorr", "InputTypes", "", "", "none", "FNGLT", "none","otucorr",false,false); parameters.push_back(potucorr);
        CommandParameter pcorraxes("corraxes", "InputTypes", "", "", "none", "FNGLT", "none","corraxes",false,false); parameters.push_back(pcorraxes);
        CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false, true); parameters.push_back(plist);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "FNGLT", "none","shared",false,false, true); parameters.push_back(pshared);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveOtuLabelsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.otus command can be used to remove specific otus with the output from classify.otu, otu.association, or corr.axes. It can also be used to select a set of otus from a shared or list file.\n";
		helpString += "The remove.otus parameters are: constaxonomy, otucorr, corraxes, shared, list, label and accnos.\n";
		helpString += "The constaxonomy parameter is input the results of the classify.otu command.\n";
        helpString += "The otucorr parameter is input the results of the otu.association command.\n";
        helpString += "The corraxes parameter is input the results of the corr.axes command.\n";
        helpString += "The label parameter is used to analyze specific labels in your input. \n";
		helpString += "The remove.otus commmand should be in the following format: \n";
		helpString += "remove.otus(accnos=yourListOfOTULabels, corraxes=yourCorrAxesFile)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveOtuLabelsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "constaxonomy")         {   pattern = "[filename],pick,[extension]";                }
        else if (type == "otucorr")         {   pattern = "[filename],pick,[extension]";                }
        else if (type == "corraxes")        {   pattern = "[filename],pick,[extension]";                }
        else if (type == "list")            {   pattern = "[filename],[distance],pick,[extension]";     }
        else if (type == "shared")          {   pattern = "[filename],[distance],pick,[extension]";     }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveOtuLabelsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RemoveOtuLabelsCommand::RemoveOtuLabelsCommand(){	
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
		m->errorOut(e, "RemoveOtuLabelsCommand", "RemoveOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveOtuLabelsCommand::RemoveOtuLabelsCommand(string option)  {
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
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                //edit file types below to include only the types you added as parameters
                
				string path;
                it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("corraxes");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["corraxes"] = inputDir + it->second;		}
				}
                
                it = parameters.find("otucorr");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["otucorr"] = inputDir + it->second;		}
				}
                
                it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
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
            accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = m->getAccnosFile(); 
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You have no valid accnos file and accnos is required."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else { m->setAccnosFile(accnosfile); }	
			
			constaxonomyfile = validParameter.validFile(parameters, "constaxonomy", true);
			if (constaxonomyfile == "not open") { constaxonomyfile = ""; abort = true; }
			else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
            
            corraxesfile = validParameter.validFile(parameters, "corraxes", true);
			if (corraxesfile == "not open") { corraxesfile = ""; abort = true; }
			else if (corraxesfile == "not found") {  corraxesfile = "";  }
            
            otucorrfile = validParameter.validFile(parameters, "otucorr", true);
			if (otucorrfile == "not open") { otucorrfile = ""; abort = true; }
			else if (otucorrfile == "not found") {  otucorrfile = "";  }
            
            listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
            else { m->setListFile(listfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
            else { m->setSharedFile(sharedfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	 outputDir = ""; 	}
            
            if ((constaxonomyfile == "") && (corraxesfile == "") && (otucorrfile == "") && (sharedfile == "") && (listfile == ""))  { m->mothurOut("You must provide one of the following: constaxonomy, corraxes, otucorr, shared or list."); m->mothurOutEndLine(); abort = true; }
            
            if ((sharedfile != "") || (listfile != "")) {
                label = validParameter.validFile(parameters, "label", false);			
                if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }
            }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "RemoveOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveOtuLabelsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //get labels you want to keep
		otulabels = m->readAccnos(accnosfile);
        //simplfy labels
        set<string> newLabels;
        for (set<string>::iterator it = otulabels.begin(); it != otulabels.end(); it++) {  newLabels.insert(m->getSimpleLabel(*it)); }
        otulabels = newLabels;
        
        if (m->debug) { m->mothurOut("[DEBUG]: numlabels = " + toString(otulabels.size()) + "\n"); }
		
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (constaxonomyfile != "")	{		readClassifyOtu();      }
		if (corraxesfile != "")		{		readCorrAxes();         }
		if (otucorrfile != "")		{		readOtuAssociation();	}
        if (listfile != "")         {		readList();             }
        if (sharedfile != "")		{		readShared();           }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
        string current = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
        }
        
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setConsTaxonomyFile(current); }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetOtuLabelsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readClassifyOtu(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(constaxonomyfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(constaxonomyfile));
        variables["[extension]"] = m->getExtension(constaxonomyfile);
		string outputFileName = getOutputFileName("constaxonomy", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(constaxonomyfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in); m->gobble(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu >> size; m->gobble(in);
            tax = m->getline(in); m->gobble(in);
            
            if (m->debug) { m->mothurOut("[DEBUG]: " + otu + toString(size) + tax + "\n"); }
            
            if (otulabels.count(m->getSimpleLabel(otu)) == 0) {
				wroteSomething = true;
                out << otu << '\t' << size << '\t' << tax << endl;
            }else {  removedCount++;  }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file only contains labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["constaxonomy"].push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " otus from your constaxonomy file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readClassifyOtu");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readOtuAssociation(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(otucorrfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(otucorrfile));
        variables["[extension]"] = m->getExtension(otucorrfile);
		string outputFileName = getOutputFileName("otucorr", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(otucorrfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in); m->gobble(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu1 = ""; 
            string otu2 = ""; 
            in >> otu1 >> otu2;
            string line = m->getline(in); m->gobble(in);
            
            if ((otulabels.count(m->getSimpleLabel(otu1)) == 0) && (otulabels.count(m->getSimpleLabel(otu2)) == 0)){
				wroteSomething = true;
                
                out << otu1 << '\t' << otu2 << '\t' << line << endl;
            }else { removedCount++; }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file only contains labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["otucorr"].push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " lines from your otu.corr file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readOtuAssociation");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readCorrAxes(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(corraxesfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(corraxesfile));
        variables["[extension]"] = m->getExtension(corraxesfile);
		string outputFileName = getOutputFileName("corraxes", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
        
		ifstream in;
		m->openInputFile(corraxesfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in); m->gobble(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu = ""; 
            in >> otu;
            string line = m->getline(in); m->gobble(in);
            
            if (otulabels.count(m->getSimpleLabel(otu)) == 0) {
				wroteSomething = true;
                
                out << otu << '\t' << line << endl;
            }else { removedCount++; }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file only contains labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["corraxes"].push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " lines from your corr.axes file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readCorrAxes");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readShared(){
	try {
        
        getShared();
        
        if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } return 0; }
        
        vector<string> newLabels;
        
        //create new "filtered" lookup
        vector<SharedRAbundVector*> newLookup;
        for (int i = 0; i < lookup.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(lookup[i]->getLabel());
			temp->setGroup(lookup[i]->getGroup());
			newLookup.push_back(temp);
        }
        
        bool wroteSomething = false;
        int numRemoved = 0;
        for (int i = 0; i < lookup[0]->getNumBins(); i++) {
            
            if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) { delete newLookup[j]; } for (int j = 0; j < lookup.size(); j++) { delete lookup[j]; } return 0; }
            
            //is this otu on the list
            if (otulabels.count(m->getSimpleLabel(m->currentSharedBinLabels[i])) == 0) {
                wroteSomething = true;
                newLabels.push_back(m->currentSharedBinLabels[i]);
                for (int j = 0; j < newLookup.size(); j++) { //add this OTU to the new lookup
                    newLookup[j]->push_back(lookup[j]->getAbundance(i), lookup[j]->getGroup());
                }
            }else { numRemoved++; }
        }
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[extension]"] = m->getExtension(sharedfile);
        variables["[distance]"] = lookup[0]->getLabel();
		string outputFileName = getOutputFileName("shared", variables); 
        ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		for (int j = 0; j < lookup.size(); j++) { delete lookup[j]; }
        
        m->currentSharedBinLabels = newLabels;
        
		newLookup[0]->printHeaders(out);
		
		for (int i = 0; i < newLookup.size(); i++) {
			out << newLookup[i]->getLabel() << '\t' << newLookup[i]->getGroup() << '\t';
			newLookup[i]->print(out);
		}
		out.close();
        
        for (int j = 0; j < newLookup.size(); j++) { delete newLookup[j]; }
        
        if (wroteSomething == false) { m->mothurOut("Your file contains only OTUs from the .accnos file."); m->mothurOutEndLine();  }
        
		m->mothurOut("Removed " + toString(numRemoved) + " OTUs from your shared file."); m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readList(){
	try {
        getListVector();
        
        if (m->control_pressed) { delete list; return 0;}
        
        ListVector newList;
        newList.setLabel(list->getLabel());
        int removedCount = 0;
        bool wroteSomething = false;
        
        vector<string> binLabels = list->getLabels();
        vector<string> newLabels;
        for (int i = 0; i < list->getNumBins(); i++) {
            
            if (m->control_pressed) { delete list; return 0;}
            
            if (otulabels.count(m->getSimpleLabel(binLabels[i])) == 0) {
                newList.push_back(list->get(i));
                newLabels.push_back(binLabels[i]);
            }else { removedCount++; }
        }
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[extension]"] = m->getExtension(listfile);
        variables["[distance]"] = list->getLabel();
		string outputFileName = getOutputFileName("list", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		delete list;
        //print new listvector
        if (newList.getNumBins() != 0) {
            wroteSomething = true;
            newList.setLabels(newLabels);
            newList.printHeaders(out);
            newList.print(out, false);
        }
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file contains only OTUs from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " OTUs from your list file."); m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveOtuLabelsCommand", "readList");
        exit(1);
    }
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::getListVector(){
	try {
		InputData input(listfile, "list");
		list = input.getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->control_pressed) {  return 0;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
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
		
		
		if (m->control_pressed) {  return 0;  }
		
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
			delete list; 
			list = input.getListVector(lastLabel);
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "getListVector");	
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtuLabelsCommand::getShared(){
	try {
		InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && (userLabels.size() != 0)) {
			if (m->control_pressed) {   return 0;  }
			
			if(labels.count(lookup[0]->getLabel()) == 1){
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
				break;
			}
			
			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input.getSharedRAbundVectors();
		}
		
		
		if (m->control_pressed) {  return 0;  }
		
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
			lookup = input.getSharedRAbundVectors(lastLabel);
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "getShared");	
		exit(1);
	}
}
//**********************************************************************************************************************



