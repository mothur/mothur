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
        CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(paccnos);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pconstaxonomy);
		CommandParameter potucorr("otucorr", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(potucorr);
        CommandParameter pcorraxes("corraxes", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pcorraxes);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The remove.otulabels command can be used to remove specific otus with the output from classify.otu, otu.association, or corr.axes.\n";
		helpString += "The remove.otulabels parameters are: constaxonomy, otucorr, corraxes, and accnos.\n";
		helpString += "The constaxonomy parameter is input the results of the classify.otu command.\n";
        helpString += "The otucorr parameter is input the results of the otu.association command.\n";
        helpString += "The corraxes parameter is input the results of the corr.axes command.\n";
		helpString += "The remove.otulabels commmand should be in the following format: \n";
		helpString += "remove.otulabels(accnos=yourListOfOTULabels, corraxes=yourCorrAxesFile)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveOtuLabelsCommand::RemoveOtuLabelsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["contaxonomy"] = tempOutNames; 
        outputTypes["otu.corr"] = tempOutNames;
        outputTypes["corr.axes"] = tempOutNames;
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
            }
            
            vector<string> tempOutNames;
            outputTypes["contaxonomy"] = tempOutNames; 
            outputTypes["otu.corr"] = tempOutNames;
            outputTypes["corr.axes"] = tempOutNames;
            
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
            
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	 outputDir = ""; 	}
            
            if ((constaxonomyfile == "") && (corraxesfile == "") && (otucorrfile == ""))  { m->mothurOut("You must provide one of the following: constaxonomy, corraxes or otucorr."); m->mothurOutEndLine(); abort = true; }
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
		readAccnos();
		
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (constaxonomyfile != "")	{		readClassifyOtu();      }
		if (corraxesfile != "")		{		readCorrAxes();         }
		if (otucorrfile != "")		{		readOtuAssociation();	}
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(constaxonomyfile)) + "pick.taxonomy";
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(constaxonomyfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu >> size >> tax; m->gobble(in);
            
            if (labels.count(otu) == 0) {
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(otucorrfile)) + "pick.corr";
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(otucorrfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu1 = ""; 
            string otu2 = ""; 
            in >> otu1 >> otu2;
            string line = m->getline(in); m->gobble(in);
            
            if ((labels.count(otu1) == 0) && (labels.count(otu2) == 0)){
				wroteSomething = true;
                
                out << otu1 << '\t' << otu2 << '\t' << line << endl;
            }else { removedCount++; }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file only contains labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["otu.corr"].push_back(outputFileName);
		
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(corraxesfile)) + "pick.axes";
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
        
		ifstream in;
		m->openInputFile(corraxesfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        //read headers
        string headers = m->getline(in);
        out << headers << endl;
        
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string otu = ""; 
            in >> otu;
            string line = m->getline(in); m->gobble(in);
            
            if (labels.count(otu) == 0) {
				wroteSomething = true;
                
                out << otu << '\t' << line << endl;
            }else { removedCount++; }
        }
        in.close();
        out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file only contains labels from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["corr.axes"].push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " lines from your corr.axes file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readCorrAxes");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveOtuLabelsCommand::readAccnos(){
	try {
		
		ifstream in;
		m->openInputFile(accnosfile, in);
		string name;
		
		while(!in.eof()){
			in >> name;
            
			labels.insert(name);
			
			m->gobble(in);
		}
		in.close();	
		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtuLabelsCommand", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************



