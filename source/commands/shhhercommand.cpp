/*
 *  shhher.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "shhhercommand.h"

//**********************************************************************************************************************
vector<string> ShhherCommand::setParameters(){	
	try {
		CommandParameter pflow("flow", "InputTypes", "", "", "none", "fileflow", "none","fasta-name-group-counts-qfile",false,false,true); parameters.push_back(pflow);
		CommandParameter pfile("file", "InputTypes", "", "", "none", "fileflow", "none","fasta-name-group-counts-qfile",false,false,true); parameters.push_back(pfile);
		CommandParameter plookup("lookup", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(plookup);
		CommandParameter pcutoff("cutoff", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pmaxiter("maxiter", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(pmaxiter);
        CommandParameter plarge("large", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(plarge);
		CommandParameter psigma("sigma", "Number", "", "60", "", "", "","",false,false); parameters.push_back(psigma);
		CommandParameter pmindelta("mindelta", "Number", "", "0.000001", "", "", "","",false,false); parameters.push_back(pmindelta);
        CommandParameter porder("order", "Multiple", "A-B-I", "A", "", "", "","",false,false, true); parameters.push_back(porder);		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhherCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The shhh.flows command reads a file containing flowgrams and creates a file of corrected sequences.\n";
        helpString += "The shhh.flows command parameters are flow, file, lookup, cutoff, processors, large, maxiter, sigma, mindelta and order.\n";
        helpString += "The flow parameter is used to input your flow file.\n";
        helpString += "The file parameter is used to input the *flow.files file created by trim.flows.\n";
        helpString += "The lookup parameter is used specify the lookup file you would like to use. http://www.mothur.org/wiki/Lookup_files.\n";
        helpString += "The order parameter options are A, B or I.  Default=A. A = TACG and B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhherCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],shhh.fasta";   }
        else if (type == "name")    {   pattern = "[filename],shhh.names";   }
        else if (type == "group")        {   pattern = "[filename],shhh.groups";   }
        else if (type == "counts")        {   pattern = "[filename],shhh.counts";   }
        else if (type == "qfile")        {   pattern = "[filename],shhh.qual";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ShhherCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

ShhherCommand::ShhherCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["counts"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;

	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "ShhherCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

ShhherCommand::ShhherCommand(string option) {
	try {
        
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
            vector<string> tempOutNames;
            outputTypes["fasta"] = tempOutNames;
            outputTypes["name"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
            outputTypes["counts"] = tempOutNames;
            outputTypes["qfile"] = tempOutNames;

			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
				
				it = parameters.find("lookup");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["lookup"] = inputDir + it->second;		}
				}

				it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
			}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
            
			//check for required parameters
			flowFileName = validParameter.validFile(parameters, "flow");
			flowFilesFileName = validParameter.validFile(parameters, "file");
			if (flowFileName == "not found" && flowFilesFileName == "not found") {
				m->mothurOut("values for either flow or file must be provided for the shhh.flows command.");
				m->mothurOutEndLine();
				abort = true; 
			}
			else if (flowFileName == "not open" || flowFilesFileName == "not open") { abort = true; }
			
			if(flowFileName != "not found"){
				compositeFASTAFileName = "";	
				compositeNamesFileName = "";	
			}
			else{
				ofstream temp;
                
                string thisoutputDir = outputDir;
                if (outputDir == "") {  thisoutputDir =  util.hasPath(flowFilesFileName); } //if user entered a file with a path then preserve it
                
				//we want to rip off .files, and also .flow if its there
                string fileroot = util.getRootName(util.getSimpleName(flowFilesFileName));
                if (fileroot[fileroot.length()-1] == '.') {  fileroot = fileroot.substr(0, fileroot.length()-1); } //rip off dot
                string extension = util.getExtension(fileroot);
                if (extension == ".flow") { fileroot = util.getRootName(fileroot);  }
                else { fileroot += "."; } //add back if needed
                
				compositeFASTAFileName = thisoutputDir + fileroot + "shhh.fasta";
				util.openOutputFile(compositeFASTAFileName, temp);
				temp.close();
				
				compositeNamesFileName = thisoutputDir + fileroot + "shhh.names";
				util.openOutputFile(compositeNamesFileName, temp);
				temp.close();
			}
            
            if(flowFilesFileName != "not found"){
                string fName;
                
                ifstream flowFilesFile;
                util.openInputFile(flowFilesFileName, flowFilesFile);
                while(flowFilesFile){
                    fName = util.getline(flowFilesFile);
                    
                    //test if file is valid
                    ifstream in;
                    bool ableToOpen = util.openInputFile(fName, in, "noerror");
                    in.close();	
                    if (!ableToOpen) {
                        if (inputDir != "") { //default path is set
                            string tryPath = inputDir + fName;
                            m->mothurOut("Unable to open " + fName + ". Trying input directory " + tryPath); m->mothurOutEndLine();
                            ifstream in2;
                            ableToOpen = util.openInputFile(tryPath, in2, "noerror");
                            in2.close();
                            fName = tryPath;
                        }
                    }
                    
                    if (!ableToOpen) {
                        if (current->getDefaultPath() != "") { //default path is set
                            string tryPath = current->getDefaultPath() + util.getSimpleName(fName);
                            m->mothurOut("Unable to open " + fName + ". Trying default " + tryPath); m->mothurOutEndLine();
                            ifstream in2;
                            ableToOpen = util.openInputFile(tryPath, in2, "noerror");
                            in2.close();
                            fName = tryPath;
                        }
                    }
                    
                    //if you can't open it its not in current working directory or inputDir, try mothur excutable location
                    if (!ableToOpen) {
                        string exepath = current->getProgramPath();
                        //string tempPath = exepath;
                        //for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
                        //exepath = exepath.substr(0, (tempPath.find_last_of('m')));
                        
                        string tryPath = util.getFullPathName(exepath) + util.getSimpleName(fName);
                        m->mothurOut("Unable to open " + fName + ". Trying mothur's executable location " + tryPath); m->mothurOutEndLine();
                        ifstream in2;
                        ableToOpen = util.openInputFile(tryPath, in2, "noerror");
                        in2.close();
                        fName = tryPath;
                    }
                    
                    if (!ableToOpen) {  m->mothurOut("Unable to open " + fName + ". Disregarding. "); m->mothurOutEndLine();  }
                    else { flowFileVector.push_back(fName); }
                    util.gobble(flowFilesFile);
                }
                flowFilesFile.close();
                if (flowFileVector.size() == 0) {  m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
            }
            else{
                if (outputDir == "") { outputDir = util.hasPath(flowFileName); }
                flowFileVector.push_back(flowFileName);
            }
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "lookup");
			if (temp == "not found")	{	
				string path = current->getProgramPath();
                //string tempPath = path;
                //for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
                //path = path.substr(0, (tempPath.find_last_of('m')));
                
#if defined NON_WINDOWS
                path += "lookupFiles/";
#else
                path += "lookupFiles\\";
#endif
				lookupFileName = util.getFullPathName(path) + "LookUp_Titanium.pat";
				
				bool ableToOpen;
				ifstream in;
				ableToOpen = util.openInputFile(lookupFileName, in, "noerror");
				in.close();	
				
				//if you can't open it, try input location
				if (!ableToOpen) {
					if (inputDir != "") { //default path is set
						string tryPath = inputDir + util.getSimpleName(lookupFileName);
						m->mothurOut("Unable to open " + lookupFileName + ". Trying input directory " + tryPath); m->mothurOutEndLine();
						ifstream in2;
						ableToOpen = util.openInputFile(tryPath, in2, "noerror");
						in2.close();
						lookupFileName = tryPath;
					}
				}
				
				//if you can't open it, try default location
				if (!ableToOpen) {
					if (current->getDefaultPath() != "") { //default path is set
						string tryPath = current->getDefaultPath() + util.getSimpleName(lookupFileName);
						m->mothurOut("Unable to open " + lookupFileName + ". Trying default " + tryPath); m->mothurOutEndLine();
						ifstream in2;
						ableToOpen = util.openInputFile(tryPath, in2, "noerror");
						in2.close();
						lookupFileName = tryPath;
					}
				}
				
				//if you can't open it its not in current working directory or inputDir, try mothur excutable location
				if (!ableToOpen) {
					string exepath = current->getProgramPath();
					//string tempPath = exepath;
					//for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
					//exepath = exepath.substr(0, (tempPath.find_last_of('m')));
					
					string tryPath = util.getFullPathName(exepath) + util.getSimpleName(lookupFileName);
					m->mothurOut("Unable to open " + lookupFileName + ". Trying mothur's executable location " + tryPath); m->mothurOutEndLine();
					ifstream in2;
					ableToOpen = util.openInputFile(tryPath, in2, "noerror");
					in2.close();
					lookupFileName = tryPath;
				}
				
				if (!ableToOpen) {  m->mothurOut("Unable to open " + lookupFileName + "."); m->mothurOutEndLine(); abort=true;  }
			}
			else if(temp == "not open")	{	
				
				lookupFileName = validParameter.valid(parameters, "lookup");
				
				//if you can't open it its not inputDir, try mothur excutable location
				string exepath = current->getProgramPath();
				//string tempPath = exepath;
				//for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
				//exepath = exepath.substr(0, (tempPath.find_last_of('m')));
					
				string tryPath = util.getFullPathName(exepath) + util.getSimpleName(lookupFileName);
				m->mothurOut("Unable to open " + lookupFileName + ". Trying mothur's executable location " + tryPath); m->mothurOutEndLine();
				ifstream in2;
				bool ableToOpen = util.openInputFile(tryPath, in2, "noerror");
				in2.close();
				lookupFileName = tryPath;
				
				if (!ableToOpen) {  m->mothurOut("Unable to open " + lookupFileName + "."); m->mothurOutEndLine(); abort=true;  }
			}else						{	lookupFileName = temp;	}
			
			temp = validParameter.valid(parameters, "cutoff");	if (temp == "not found"){	temp = "0.01";		}
			util.mothurConvert(temp, cutoff); 
			
			temp = validParameter.valid(parameters, "mindelta");	if (temp == "not found"){	temp = "0.000001";	}
			util.mothurConvert(temp, minDelta); 

			temp = validParameter.valid(parameters, "maxiter");	if (temp == "not found"){	temp = "1000";		}
			util.mothurConvert(temp, maxIters); 
            
            temp = validParameter.valid(parameters, "large");	if (temp == "not found"){	temp = "0";		}
			util.mothurConvert(temp, largeSize); 
            if (largeSize != 0) { large = true; }
            else { large = false;  }
            if (largeSize < 0) {  m->mothurOut("The value of the large cannot be negative.\n"); }

			temp = validParameter.valid(parameters, "sigma"); if (temp == "not found")	{	temp = "60";		}
			util.mothurConvert(temp, sigma); 
			
			temp = validParameter.valid(parameters, "order");  if (temp == "not found"){ 	temp = "A";	}
            if (temp.length() > 1) {  m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
            }
            else {
                if (toupper(temp[0]) == 'A') {  flowOrder = "TACG";   }
                else if(toupper(temp[0]) == 'B'){
                    flowOrder = "TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC";   }
                else if(toupper(temp[0]) == 'I'){
                    flowOrder = "TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC";   }
                else {
                    m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
                }
            }

			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "ShhherCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ShhherCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		getSingleLookUp();	if (m->getControl_pressed()) { return 0; }
		getJointLookUp();	if (m->getControl_pressed()) { return 0; }
		
        driver(flowFileVector, compositeFASTAFileName, compositeNamesFileName);

		if(compositeFASTAFileName != ""){
			outputNames.push_back(compositeFASTAFileName); outputTypes["fasta"].push_back(compositeFASTAFileName);
			outputNames.push_back(compositeNamesFileName); outputTypes["name"].push_back(compositeNamesFileName);
		}

		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "execute");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts biggest to smallest
inline bool compareFileSizes(string left, string right){
    
    FILE * pFile;
    long leftsize = 0;
    
    //get num bytes in file
    string filename = left;
    pFile = fopen (filename.c_str(),"rb");
    string error = "Error opening " + filename;
    if (pFile==NULL) perror (error.c_str());
    else{
        fseek (pFile, 0, SEEK_END);
        leftsize=ftell (pFile);
        fclose (pFile);
    }
    
    FILE * pFile2;
    long rightsize = 0;
    
    //get num bytes in file
    filename = right;
    pFile2 = fopen (filename.c_str(),"rb");
    error = "Error opening " + filename;
    if (pFile2==NULL) perror (error.c_str());
    else{
        fseek (pFile2, 0, SEEK_END);
        rightsize=ftell (pFile2);
        fclose (pFile2);
    }
    
    return (leftsize > rightsize);	
} 
/**************************************************************************************************/

vector<string> ShhherCommand::parseFlowFiles(string filename){
    try {
        vector<string> files;
        int count = 0;
        
        ifstream in;
        util.openInputFile(filename, in);
        
        int thisNumFLows = 0;
        in >> thisNumFLows; util.gobble(in);
        
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            ofstream out;
            string outputFileName = filename + toString(count) + ".temp";
            util.openOutputFile(outputFileName, out);
            out << thisNumFLows << endl;
            files.push_back(outputFileName);
            
            int numLinesWrote = 0;
            for (int i = 0; i < largeSize; i++) {
                if (in.eof()) { break; }
                string line = util.getline(in); util.gobble(in);
                out << line << endl;
                numLinesWrote++;
            }
            out.close();
            
            if (numLinesWrote == 0) {  util.mothurRemove(outputFileName); files.pop_back();  }
            count++;
        }
        in.close();
        
        if (m->getControl_pressed()) { for (int i = 0; i < files.size(); i++) { util.mothurRemove(files[i]); }  files.clear(); }
        
        m->mothurOut("\nDivided " + filename + " into " + toString(files.size()) + " files.\n\n"); 
        
        return files;
    }
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "parseFlowFiles");
		exit(1);
	}
}
/**************************************************************************************************/

int ShhherCommand::driver(vector<string> filenames, string thisCompositeFASTAFileName, string thisCompositeNamesFileName){
    try {
        
        int numCompleted = 0;
        
        for(int i=0;i<filenames.size();i++){
			
			if (m->getControl_pressed()) { break; }
			
            vector<string> theseFlowFileNames; theseFlowFileNames.push_back(filenames[i]);
            if (large) {  theseFlowFileNames = parseFlowFiles(filenames[i]);  }
            
            if (m->getControl_pressed()) { break; }
            
            double begClock = clock();
            unsigned long long begTime;
            
            string fileNameForOutput = filenames[i];
            
            for (int g = 0; g < theseFlowFileNames.size(); g++) {
                
                string flowFileName = theseFlowFileNames[g];
                m->mothurOut("\n>>>>>\tProcessing " + flowFileName + " (file " + toString(i+1) + " of " + toString(filenames.size()) + ")\t<<<<<\n");
                m->mothurOut("Reading flowgrams...\n");
                
                vector<string> seqNameVector;
                vector<int> lengths;
                vector<short> flowDataIntI;
                vector<double> flowDataPrI;
                map<string, int> nameMap;
                vector<short> uniqueFlowgrams;
                vector<int> uniqueCount;
                vector<int> mapSeqToUnique;
                vector<int> mapUniqueToSeq;
                vector<int> uniqueLengths;
                int numFlowCells;
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: About to read flowgrams.\n"); }
                int numSeqs = getFlowData(flowFileName, seqNameVector, lengths, flowDataIntI, nameMap, numFlowCells);
                
                if (m->getControl_pressed()) { break; }
                
                m->mothurOut("Identifying unique flowgrams...\n");
                int numUniques = getUniques(numSeqs, numFlowCells, uniqueFlowgrams, uniqueCount, uniqueLengths, mapSeqToUnique, mapUniqueToSeq, lengths, flowDataPrI, flowDataIntI);
                
                if (m->getControl_pressed()) { break; }
                
                m->mothurOut("Calculating distances between flowgrams...\n");
                string distFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.dist";
                begTime = time(NULL);
               
                
                flowDistParentFork(numFlowCells, distFileName, numUniques, mapUniqueToSeq, mapSeqToUnique, lengths, flowDataPrI, flowDataIntI);
                
                m->mothurOutEndLine();
                m->mothurOut("Total time: " + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/CLOCKS_PER_SEC) + '\n');
                
                
                string namesFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.names";
                createNamesFile(numSeqs, numUniques, namesFileName, seqNameVector, mapSeqToUnique, mapUniqueToSeq);
                
                if (m->getControl_pressed()) { break; }
                
                m->mothurOut("\nClustering flowgrams...\n");
                string listFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.list";
                cluster(listFileName, distFileName, namesFileName);
                
                if (m->getControl_pressed()) { break; }
                
                vector<int> otuData;
                vector<int> cumNumSeqs;
                vector<int> nSeqsPerOTU;
                vector<vector<int> > aaP;	//tMaster->aanP:	each row is a different otu / each col contains the sequence indices
                vector<vector<int> > aaI;	//tMaster->aanI:	that are in each otu - can't differentiate between aaP and aaI 
                vector<int> seqNumber;		//tMaster->anP:		the sequence id number sorted by OTU
                vector<int> seqIndex;		//tMaster->anI;		the index that corresponds to seqNumber
                
                
                int numOTUs = getOTUData(numSeqs, listFileName, otuData, cumNumSeqs, nSeqsPerOTU, aaP, aaI, seqNumber, seqIndex, nameMap);
                
                if (m->getControl_pressed()) { break; }
                
                util.mothurRemove(distFileName);
                util.mothurRemove(namesFileName);
                util.mothurRemove(listFileName);
                
                vector<double> dist;		//adDist - distance of sequences to centroids
                vector<short> change;		//did the centroid sequence change? 0 = no; 1 = yes
                vector<int> centroids;		//the representative flowgram for each cluster m
                vector<double> weight;
                vector<double> singleTau;	//tMaster->adTau:	1-D Tau vector (1xnumSeqs)
                vector<int> nSeqsBreaks;
                vector<int> nOTUsBreaks;
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: numSeqs = " + toString(numSeqs) + " numOTUS = " + toString(numOTUs) + " about to alloc a dist vector with size = " + toString((numSeqs * numOTUs)) + ".\n"); }
                
                dist.assign(numSeqs * numOTUs, 0);
                change.assign(numOTUs, 1);
                centroids.assign(numOTUs, -1);
                weight.assign(numOTUs, 0);
                singleTau.assign(numSeqs, 1.0);
                
                nSeqsBreaks.assign(2, 0);
                nOTUsBreaks.assign(2, 0);
                
                nSeqsBreaks[0] = 0;
                nSeqsBreaks[1] = numSeqs;
                nOTUsBreaks[1] = numOTUs;
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: done allocating memory, about to denoise.\n"); }
                
                if (m->getControl_pressed()) { break; }
                
                double maxDelta = 0;
                int iter = 0;
                
                begClock = clock();
                begTime = time(NULL);
                
                m->mothurOut("\nDenoising flowgrams...\n");
                m->mothurOut("iter\tmaxDelta\tnLL\t\tcycletime\n");
                
                while((maxIters == 0 && maxDelta > minDelta) || iter < MIN_ITER || (maxDelta > minDelta && iter < maxIters)){
                    
                    if (m->getControl_pressed()) { break; }
                    
                    double cycClock = clock();
                    unsigned long long cycTime = time(NULL);
                    fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);
                    
                    if (m->getControl_pressed()) { break; }
                    
                    calcCentroidsDriver(numOTUs, cumNumSeqs, nSeqsPerOTU, seqIndex, change, centroids, singleTau, mapSeqToUnique, uniqueFlowgrams, flowDataIntI, lengths, numFlowCells, seqNumber);
                    
                    if (m->getControl_pressed()) { break; }
                    
                    maxDelta = getNewWeights(numOTUs, cumNumSeqs, nSeqsPerOTU, singleTau, seqNumber, weight);  
                    
                    if (m->getControl_pressed()) { break; }
                    
                    double nLL = getLikelihood(numSeqs, numOTUs, nSeqsPerOTU, seqNumber, cumNumSeqs, seqIndex, dist, weight); 
                    
                    if (m->getControl_pressed()) { break; }
                    
                    checkCentroids(numOTUs, centroids, weight);
                    
                    if (m->getControl_pressed()) { break; }
                    
                    calcNewDistances(numSeqs, numOTUs, nSeqsPerOTU,  dist, weight, change, centroids, aaP, singleTau, aaI, seqNumber, seqIndex, uniqueFlowgrams, flowDataIntI, numFlowCells, lengths);
                    
                    if (m->getControl_pressed()) { break; }
                    
                    iter++;
                    
                    m->mothurOut(toString(iter) + '\t' + toString(maxDelta) + '\t' + toString(nLL) + '\t' + toString(time(NULL) - cycTime) + '\t' + toString((clock() - cycClock)/(double)CLOCKS_PER_SEC) + '\n');
                    
                }	
                
                if (m->getControl_pressed()) { break; }
                
                m->mothurOut("\nFinalizing...\n");
                fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: done fill().\n"); }

                if (m->getControl_pressed()) { break; }
                
                setOTUs(numOTUs, numSeqs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, otuData, singleTau, dist, aaP, aaI);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: done setOTUs().\n"); }
                
                if (m->getControl_pressed()) { break; }
                
                vector<int> otuCounts(numOTUs, 0);
                for(int j=0;j<numSeqs;j++)	{	otuCounts[otuData[j]]++;	}
                
                calcCentroidsDriver(numOTUs, cumNumSeqs, nSeqsPerOTU, seqIndex, change, centroids, singleTau, mapSeqToUnique, uniqueFlowgrams, flowDataIntI, lengths, numFlowCells, seqNumber);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: done calcCentroidsDriver().\n"); }
                
                if (m->getControl_pressed()) { break; }
                
                if ((large) && (g == 0)) {  flowFileName = filenames[i]; theseFlowFileNames[0] = filenames[i]; }
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir = util.hasPath(flowFileName);  }
                map<string, string> variables; 
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(flowFileName));
                string qualityFileName = getOutputFileName("qfile",variables);
                string fastaFileName = getOutputFileName("fasta",variables);
                string nameFileName = getOutputFileName("name",variables);
                string otuCountsFileName = getOutputFileName("counts",variables);
                string fileRoot = util.getRootName(util.getSimpleName(flowFileName));
                int pos = fileRoot.find_first_of('.');
                string fileGroup = fileRoot;
                if (pos != string::npos) {  fileGroup = fileRoot.substr(pos+1, (fileRoot.length()-1-(pos+1)));  }
                string groupFileName = getOutputFileName("group",variables);

                
                writeQualities(numOTUs, numFlowCells, qualityFileName, otuCounts, nSeqsPerOTU, seqNumber, singleTau, flowDataIntI, uniqueFlowgrams, cumNumSeqs, mapUniqueToSeq, seqNameVector, centroids, aaI); if (m->getControl_pressed()) { break; }
                writeSequences(thisCompositeFASTAFileName, numOTUs, numFlowCells, fastaFileName, otuCounts, uniqueFlowgrams, seqNameVector, aaI, centroids);if (m->getControl_pressed()) { break; }
                writeNames(thisCompositeNamesFileName, numOTUs, nameFileName, otuCounts, seqNameVector, aaI, nSeqsPerOTU);				if (m->getControl_pressed()) { break; }
                writeClusters(otuCountsFileName, numOTUs, numFlowCells,otuCounts, centroids, uniqueFlowgrams, seqNameVector, aaI, nSeqsPerOTU, lengths, flowDataIntI);			if (m->getControl_pressed()) { break; }
                writeGroups(groupFileName, fileGroup, numSeqs, seqNameVector);						if (m->getControl_pressed()) { break; }
                
                if (large) {
                    if (g > 0) {
                        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(theseFlowFileNames[0]));
                        util.appendFiles(qualityFileName, getOutputFileName("qfile",variables));
                        util.mothurRemove(qualityFileName);
                        util.appendFiles(fastaFileName, getOutputFileName("fasta",variables));
                        util.mothurRemove(fastaFileName);
                        util.appendFiles(nameFileName, getOutputFileName("name",variables));
                        util.mothurRemove(nameFileName);
                        util.appendFiles(otuCountsFileName, getOutputFileName("counts",variables));
                        util.mothurRemove(otuCountsFileName);
                        util.appendFiles(groupFileName, getOutputFileName("group",variables));
                        util.mothurRemove(groupFileName);
                    }
                    util.mothurRemove(theseFlowFileNames[g]);
                }
			}
            
            numCompleted++;
			m->mothurOut("Total time to process " + fileNameForOutput + ":\t" + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/(double)CLOCKS_PER_SEC) + '\n');
		}
		
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
        
        return numCompleted;
        
    }catch(exception& e) {
            m->errorOut(e, "ShhherCommand", "driver");
            exit(1);
    }
}

/**************************************************************************************************/
int ShhherCommand::getFlowData(string filename, vector<string>& thisSeqNameVector, vector<int>& thisLengths, vector<short>& thisFlowDataIntI, map<string, int>& thisNameMap, int& numFlowCells){
	try{
       
		ifstream flowFile;
       
		util.openInputFile(filename, flowFile);
		
		string seqName;
		int currentNumFlowCells;
		float intensity;
        thisSeqNameVector.clear();
		thisLengths.clear();
		thisFlowDataIntI.clear();
		thisNameMap.clear();
		
		string numFlowTest;
        flowFile >> numFlowTest;
        
        if (!util.isContainingOnlyDigits(numFlowTest)) { m->mothurOut("[ERROR]: expected a number and got " + numFlowTest + ", quitting. Did you use the flow parameter instead of the file parameter?"); m->mothurOutEndLine(); exit(1); }
        else { convert(numFlowTest, numFlowCells); }
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: numFlowCells = " + toString(numFlowCells) + ".\n"); }
		int index = 0;//pcluster
		while(!flowFile.eof()){
			
			if (m->getControl_pressed()) { break; }
			
			flowFile >> seqName >> currentNumFlowCells;
            
			thisLengths.push_back(currentNumFlowCells);
           
			thisSeqNameVector.push_back(seqName);
			thisNameMap[seqName] = index++;//pcluster
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: seqName = " + seqName + " length = " + toString(currentNumFlowCells) + " index = " + toString(index) + "\n"); }
            
			for(int i=0;i<numFlowCells;i++){
				flowFile >> intensity;
				if(intensity > 9.99)	{	intensity = 9.99;	}
				int intI = int(100 * intensity + 0.0001);
				thisFlowDataIntI.push_back(intI);
			}
			util.gobble(flowFile);
		}
		flowFile.close();
		
		int numSeqs = thisSeqNameVector.size();		
		
		for(int i=0;i<numSeqs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int iNumFlowCells = i * numFlowCells;
			for(int j=thisLengths[i];j<numFlowCells;j++){
				thisFlowDataIntI[iNumFlowCells + j] = 0;
			}
		}
        
        return numSeqs;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getFlowData");
		exit(1);
	}
}
/**************************************************************************************************/

int ShhherCommand::flowDistParentFork(int numFlowCells, string distFileName, int stopSeq, vector<int>& mapUniqueToSeq, vector<int>& mapSeqToUnique, vector<int>& lengths, vector<double>& flowDataPrI, vector<short>& flowDataIntI){
	try{		
        
		ostringstream outStream;
		outStream.setf(ios::fixed, ios::floatfield);
		outStream.setf(ios::dec, ios::basefield);
		outStream.setf(ios::showpoint);
		outStream.precision(6);
		
		int begTime = time(NULL);
		double begClock = clock();
        
		for(int i=0;i<stopSeq;i++){
			
			if (m->getControl_pressed()) { break; }
			
			for(int j=0;j<i;j++){
				float flowDistance = calcPairwiseDist(numFlowCells, mapUniqueToSeq[i], mapUniqueToSeq[j], mapSeqToUnique, lengths, flowDataPrI, flowDataIntI);
                
				if(flowDistance < 1e-6){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << 0.000000 << endl;
				}
				else if(flowDistance <= cutoff){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << flowDistance << endl;
				}
			}
			if(i % 100 == 0){
				m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - begTime));
				m->mothurOutJustToScreen("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC)+"\n");
			}
		}
		
		ofstream distFile(distFileName.c_str());
		distFile << outStream.str();		
		distFile.close();
		
		if (m->getControl_pressed()) {}
		else {
			m->mothurOutJustToScreen(toString(stopSeq-1) + "\t" + toString(time(NULL) - begTime));
			m->mothurOutJustToScreen("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC)+"\n");
		}
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "flowDistParentFork");
		exit(1);
	}
}
/**************************************************************************************************/

float ShhherCommand::calcPairwiseDist(int numFlowCells, int seqA, int seqB, vector<int>& mapSeqToUnique, vector<int>& lengths, vector<double>& flowDataPrI, vector<short>& flowDataIntI){
	try{
		int minLength = lengths[mapSeqToUnique[seqA]];
		if(lengths[seqB] < minLength){	minLength = lengths[mapSeqToUnique[seqB]];	}
		
		int ANumFlowCells = seqA * numFlowCells;
		int BNumFlowCells = seqB * numFlowCells;
		
		float dist = 0;
		
		for(int i=0;i<minLength;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int flowAIntI = flowDataIntI[ANumFlowCells + i];
			float flowAPrI = flowDataPrI[ANumFlowCells + i];
			
			int flowBIntI = flowDataIntI[BNumFlowCells + i];
			float flowBPrI = flowDataPrI[BNumFlowCells + i];
			dist += jointLookUp[flowAIntI * NUMBINS + flowBIntI] - flowAPrI - flowBPrI;
		}
		
		dist /= (float) minLength;
		return dist;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcPairwiseDist");
		exit(1);
	}
}

/**************************************************************************************************/

int ShhherCommand::getUniques(int numSeqs, int numFlowCells, vector<short>& uniqueFlowgrams, vector<int>& uniqueCount, vector<int>& uniqueLengths, vector<int>& mapSeqToUnique, vector<int>& mapUniqueToSeq, vector<int>& lengths, vector<double>& flowDataPrI, vector<short>& flowDataIntI){
	try{
		int numUniques = 0;
		uniqueFlowgrams.assign(numFlowCells * numSeqs, -1);
		uniqueCount.assign(numSeqs, 0);							//	anWeights
		uniqueLengths.assign(numSeqs, 0);
		mapSeqToUnique.assign(numSeqs, -1);
		mapUniqueToSeq.assign(numSeqs, -1);
		
		vector<short> uniqueFlowDataIntI(numFlowCells * numSeqs, -1);
		
		for(int i=0;i<numSeqs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int index = 0;
			
			vector<short> current(numFlowCells);
			for(int j=0;j<numFlowCells;j++){
				current[j] = short(((flowDataIntI[i * numFlowCells + j] + 50.0)/100.0));
			}
            
			for(int j=0;j<numUniques;j++){
				int offset = j * numFlowCells;
				bool toEnd = 1;
				
				int shorterLength;
				if(lengths[i] < uniqueLengths[j])	{	shorterLength = lengths[i];			}
				else								{	shorterLength = uniqueLengths[j];	}
                
				for(int k=0;k<shorterLength;k++){
					if(current[k] != uniqueFlowgrams[offset + k]){
						toEnd = 0;
						break;
					}
				}
				
				if(toEnd){
					mapSeqToUnique[i] = j;
					uniqueCount[j]++;
					index = j;
					if(lengths[i] > uniqueLengths[j])	{	uniqueLengths[j] = lengths[i];	}
					break;
				}
				index++;
			}
			
			if(index == numUniques){
				uniqueLengths[numUniques] = lengths[i];
				uniqueCount[numUniques] = 1;
				mapSeqToUnique[i] = numUniques;//anMap
				mapUniqueToSeq[numUniques] = i;//anF
				
				for(int k=0;k<numFlowCells;k++){
					uniqueFlowgrams[numUniques * numFlowCells + k] = current[k];
					uniqueFlowDataIntI[numUniques * numFlowCells + k] = flowDataIntI[i * numFlowCells + k];
				}
				
				numUniques++;
			}
		}
		uniqueFlowDataIntI.resize(numFlowCells * numUniques);
		uniqueLengths.resize(numUniques);	
		
		flowDataPrI.resize(numSeqs * numFlowCells, 0);
		for(int i=0;i<flowDataPrI.size();i++)	{	if (m->getControl_pressed()) { break; } flowDataPrI[i] = getProbIntensity(flowDataIntI[i]);		}
        
        return numUniques;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getUniques");
		exit(1);
	}
}
/**************************************************************************************************/
int ShhherCommand::createNamesFile(int numSeqs, int numUniques, string filename, vector<string>& seqNameVector, vector<int>& mapSeqToUnique, vector<int>& mapUniqueToSeq){
	try{
		
		vector<string> duplicateNames(numUniques, "");
		for(int i=0;i<numSeqs;i++){
			duplicateNames[mapSeqToUnique[i]] += seqNameVector[i] + ',';
		}
		
		ofstream nameFile;
		util.openOutputFile(filename, nameFile);
		
		for(int i=0;i<numUniques;i++){
			
			if (m->getControl_pressed()) { break; }
			
            //			nameFile << seqNameVector[mapUniqueToSeq[i]] << '\t' << duplicateNames[i].substr(0, duplicateNames[i].find_last_of(',')) << endl;
			nameFile << mapUniqueToSeq[i] << '\t' << duplicateNames[i].substr(0, duplicateNames[i].find_last_of(',')) << endl;
		}
		
		nameFile.close();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************

int ShhherCommand::cluster(string filename, string distFileName, string namesFileName){
	try {
		
		ReadMatrix* read = new ReadColumnMatrix(distFileName); 	
		read->setCutoff(cutoff);
		
		NameAssignment* clusterNameMap = new NameAssignment(namesFileName);
		clusterNameMap->readMap();
		read->read(clusterNameMap);
        
		ListVector* list = read->getListVector();
		SparseDistanceMatrix* matrix = read->getDMatrix();
		
		delete read; 
		delete clusterNameMap; 
        
		RAbundVector* rabund = new RAbundVector(list->getRAbundVector());
		
        float adjust = -1.0;
		Cluster* cluster = new CompleteLinkage(rabund, list, matrix, cutoff, "furthest", adjust);
		string tag = cluster->getTag();
		
		double clusterCutoff = cutoff;
		while (matrix->getSmallDist() <= clusterCutoff && matrix->getNNodes() > 0){
			
			if (m->getControl_pressed()) { break; }
			
			cluster->update(clusterCutoff);
		}
		
		list->setLabel(toString(cutoff));
		
		ofstream listFile;
		util.openOutputFile(filename, listFile);
		list->print(listFile, true);
		listFile.close();
		
		delete matrix;	delete cluster;	delete rabund; delete list;
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "cluster");
		exit(1);	
	}		
}
/**************************************************************************************************/

int ShhherCommand::getOTUData(int numSeqs, string fileName,  vector<int>& otuData,
                               vector<int>& cumNumSeqs,
                               vector<int>& nSeqsPerOTU,
                               vector<vector<int> >& aaP,	//tMaster->aanP:	each row is a different otu / each col contains the sequence indices
                               vector<vector<int> >& aaI,	//tMaster->aanI:	that are in each otu - can't differentiate between aaP and aaI 
                               vector<int>& seqNumber,		//tMaster->anP:		the sequence id number sorted by OTU
                               vector<int>& seqIndex,
                               map<string, int>& nameMap){
	try {
        InputData input(fileName, "list", nullVector);
        ListVector* list = input.getListVector();
        
        string label = list->getLabel();
        int numOTUs = list->getNumBins();
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: Getting OTU Data...\n"); }
        
		otuData.assign(numSeqs, 0);
		cumNumSeqs.assign(numOTUs, 0);
		nSeqsPerOTU.assign(numOTUs, 0);
		aaP.clear();aaP.resize(numOTUs);
		
		seqNumber.clear();
		aaI.clear();
		seqIndex.clear();
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
            if (m->getDebug()) { m->mothurOut("[DEBUG]: processing OTU " + toString(i) + ".\n"); }
            
			string singleOTU = list->get(i);
			
            vector<string> otuSeqs; util.splitAtComma(singleOTU, otuSeqs);
            
			for(int j=0;j<otuSeqs.size();j++){
				
                string seqName = otuSeqs[j];
                map<string,int>::iterator nmIt = nameMap.find(seqName);
                int index = nmIt->second;
						
                nameMap.erase(nmIt);
                otuData[index] = i;
                nSeqsPerOTU[i]++;
                aaP[i].push_back(index);
            }
			
			sort(aaP[i].begin(), aaP[i].end());
			for(int j=0;j<nSeqsPerOTU[i];j++)       { seqNumber.push_back(aaP[i][j]);   }
			for(int j=nSeqsPerOTU[i];j<numSeqs;j++) { aaP[i].push_back(0);              }
		}
		
		for(int i=1;i<numOTUs;i++){ cumNumSeqs[i] = cumNumSeqs[i-1] + nSeqsPerOTU[i-1]; }
		aaI = aaP;
		seqIndex = seqNumber;
        delete list;
      
        return numOTUs;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getOTUData");
		exit(1);	
	}		
}
/**************************************************************************************************/

int ShhherCommand::calcCentroidsDriver(int numOTUs, 
                                          vector<int>& cumNumSeqs,
                                          vector<int>& nSeqsPerOTU,
                                          vector<int>& seqIndex,
                                          vector<short>& change,		//did the centroid sequence change? 0 = no; 1 = yes
                                          vector<int>& centroids,		//the representative flowgram for each cluster m
                                          vector<double>& singleTau,	//tMaster->adTau:	1-D Tau vector (1xnumSeqs)
                                          vector<int>& mapSeqToUnique,
                                          vector<short>& uniqueFlowgrams,
                                          vector<short>& flowDataIntI,
                                          vector<int>& lengths,
                                          int numFlowCells,
                                          vector<int>& seqNumber){                          
	
	//this function gets the most likely homopolymer length at a flow position for a group of sequences
	//within an otu
	
	try{
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			double count = 0;
			int position = 0;
			int minFlowGram = 100000000;
			double minFlowValue = 1e8;
			change[i] = 0; //FALSE
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				count += singleTau[seqNumber[cumNumSeqs[i] + j]];
			}
            
			if(nSeqsPerOTU[i] > 0 && count > MIN_COUNT){
				vector<double> adF(nSeqsPerOTU[i]);
				vector<int> anL(nSeqsPerOTU[i]);
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int index = cumNumSeqs[i] + j;
					int nI = seqIndex[index];
					int nIU = mapSeqToUnique[nI];
					
					int k;
					for(k=0;k<position;k++){
						if(nIU == anL[k]){
							break;
						}
					}
					if(k == position){
						anL[position] = nIU;
						adF[position] = 0.0000;
						position++;
					}						
				}
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int index = cumNumSeqs[i] + j;
					int nI = seqIndex[index];
					
					double tauValue = singleTau[seqNumber[index]];
					
					for(int k=0;k<position;k++){
						double dist = getDistToCentroid(anL[k], nI, lengths[nI], uniqueFlowgrams, flowDataIntI, numFlowCells);
						adF[k] += dist * tauValue;
					}
				}
				
				for(int j=0;j<position;j++){
					if(adF[j] < minFlowValue){
						minFlowGram = j;
						minFlowValue = adF[j];
					}
				}
				
				if(centroids[i] != anL[minFlowGram]){
					change[i] = 1;
					centroids[i] = anL[minFlowGram];
				}
			}
			else if(centroids[i] != -1){
				change[i] = 1;
				centroids[i] = -1;			
			}
		}
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcCentroidsDriver");
		exit(1);	
	}		
}
/**************************************************************************************************/

double ShhherCommand::getDistToCentroid(int cent, int flow, int length, vector<short>& uniqueFlowgrams,
                                        vector<short>& flowDataIntI, int numFlowCells){
	try{
		
		int flowAValue = cent * numFlowCells;
		int flowBValue = flow * numFlowCells;
		
		double dist = 0;
        
		for(int i=0;i<length;i++){
			dist += singleLookUp[uniqueFlowgrams[flowAValue] * NUMBINS + flowDataIntI[flowBValue]];
			flowAValue++;
			flowBValue++;
		}
		
		return dist / (double)length;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getDistToCentroid");
		exit(1);	
	}		
}
/**************************************************************************************************/

double ShhherCommand::getNewWeights(int numOTUs, vector<int>& cumNumSeqs, vector<int>& nSeqsPerOTU, vector<double>& singleTau, vector<int>& seqNumber, vector<double>& weight){
	try{
		
		double maxChange = 0;
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			double difference = weight[i];
			weight[i] = 0;
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				double tauValue = singleTau[seqNumber[index]];
				weight[i] += tauValue;
			}
			
			difference = fabs(weight[i] - difference);
			if(difference > maxChange){	maxChange = difference;	}
		}
		return maxChange;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getNewWeights");
		exit(1);	
	}		
}

/**************************************************************************************************/

double ShhherCommand::getLikelihood(int numSeqs, int numOTUs, vector<int>& nSeqsPerOTU, vector<int>& seqNumber, vector<int>& cumNumSeqs, vector<int>& seqIndex, vector<double>& dist, vector<double>& weight){
	
	try{
		
		vector<long double> P(numSeqs, 0);
		int effNumOTUs = 0;
		
		for(int i=0;i<numOTUs;i++){
			if(weight[i] > MIN_WEIGHT){
				effNumOTUs++;
			}
		}
		
		string hold;
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				int nI = seqIndex[index];
				double singleDist = dist[seqNumber[index]];
				
				P[nI] += weight[i] * exp(-singleDist * sigma);
			}
		}
		double nLL = 0.00;
		for(int i=0;i<numSeqs;i++){
			if(P[i] == 0){	P[i] = DBL_EPSILON;	}
            
			nLL += -log(P[i]);
		}
		
		nLL = nLL -(double)numSeqs * log(sigma);
        
		return nLL; 
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getNewWeights");
		exit(1);	
	}		
}

/**************************************************************************************************/

int ShhherCommand::checkCentroids(int numOTUs, vector<int>& centroids, vector<double>& weight){
	try{
		vector<int> unique(numOTUs, 1);
		
		for(int i=0;i<numOTUs;i++){
			if(centroids[i] == -1 || weight[i] < MIN_WEIGHT){
				unique[i] = -1;
			}
		}
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			if(unique[i] == 1){
				for(int j=i+1;j<numOTUs;j++){
					if(unique[j] == 1){
						
						if(centroids[j] == centroids[i]){
							unique[j] = 0;
							centroids[j] = -1;
							
							weight[i] += weight[j];
							weight[j] = 0.0;
						}
					}
				}
			}
		}
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "checkCentroids");
		exit(1);	
	}		
}
/**************************************************************************************************/

void ShhherCommand::calcNewDistances(int numSeqs, int numOTUs, vector<int>& nSeqsPerOTU, vector<double>& dist, 
                                     vector<double>& weight, vector<short>& change, vector<int>& centroids,
                                     vector<vector<int> >& aaP,	vector<double>& singleTau, vector<vector<int> >& aaI,	
                                     vector<int>& seqNumber, vector<int>& seqIndex,
                                     vector<short>& uniqueFlowgrams,
                                     vector<short>& flowDataIntI, int numFlowCells, vector<int>& lengths){
	
	try{
		
		int total = 0;
		vector<double> newTau(numOTUs,0);
		vector<double> norms(numSeqs, 0);
		nSeqsPerOTU.assign(numOTUs, 0);
        
		for(int i=0;i<numSeqs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int indexOffset = i * numOTUs;
            
			double offset = 1e8;
			
			for(int j=0;j<numOTUs;j++){
                
				if(weight[j] > MIN_WEIGHT && change[j] == 1){
					dist[indexOffset + j] = getDistToCentroid(centroids[j], i, lengths[i], uniqueFlowgrams, flowDataIntI, numFlowCells);
				}
                
				if(weight[j] > MIN_WEIGHT && dist[indexOffset + j] < offset){
					offset = dist[indexOffset + j];
				}
			}
            
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT){
					newTau[j] = exp(sigma * (-dist[indexOffset + j] + offset)) * weight[j];
					norms[i] += newTau[j];
				}
				else{
					newTau[j] = 0.0;
				}
			}
            
			for(int j=0;j<numOTUs;j++){
				newTau[j] /= norms[i];
			}
            
			for(int j=0;j<numOTUs;j++){
				if(newTau[j] > MIN_TAU){
					
					int oldTotal = total;
					
					total++;
					
					singleTau.resize(total, 0);
					seqNumber.resize(total, 0);
					seqIndex.resize(total, 0);
					
					singleTau[oldTotal] = newTau[j];
					
					aaP[j][nSeqsPerOTU[j]] = oldTotal;
					aaI[j][nSeqsPerOTU[j]] = i;
					nSeqsPerOTU[j]++;
				}
			}
            
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistances");
		exit(1);	
	}		
}
/**************************************************************************************************/

int ShhherCommand::fill(int numOTUs, vector<int>& seqNumber, vector<int>& seqIndex, vector<int>& cumNumSeqs, vector<int>& nSeqsPerOTU, vector<vector<int> >& aaP, vector<vector<int> >& aaI){
	try {
		int index = 0;
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			cumNumSeqs[i] = index;
			for(int j=0;j<nSeqsPerOTU[i];j++){
				seqNumber[index] = aaP[i][j];
				seqIndex[index] = aaI[i][j];
				
				index++;
			}
		}
        
        return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "fill");
		exit(1);	
	}		
}
/**************************************************************************************************/

void ShhherCommand::setOTUs(int numOTUs, int numSeqs, vector<int>& seqNumber, vector<int>& seqIndex, vector<int>& cumNumSeqs, vector<int>& nSeqsPerOTU,
                            vector<int>& otuData, vector<double>& singleTau, vector<double>& dist, vector<vector<int> >& aaP, vector<vector<int> >& aaI){
	
	try {
		vector<double> bigTauMatrix(numOTUs * numSeqs, 0.0000);
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				double tauValue = singleTau[seqNumber[index]];
				int sIndex = seqIndex[index];
				bigTauMatrix[sIndex * numOTUs + i] = tauValue;				
			}
		}
		
		for(int i=0;i<numSeqs;i++){
			double maxTau = -1.0000;
			int maxOTU = -1;
			for(int j=0;j<numOTUs;j++){
				if(bigTauMatrix[i * numOTUs + j] > maxTau){
					maxTau = bigTauMatrix[i * numOTUs + j];
					maxOTU = j;
				}
			}
			
			otuData[i] = maxOTU;
		}
		
		nSeqsPerOTU.assign(numOTUs, 0);		
		
		for(int i=0;i<numSeqs;i++){
			int index = otuData[i];
			
			singleTau[i] = 1.0000;
			dist[i] = 0.0000;
			
			aaP[index][nSeqsPerOTU[index]] = i;
			aaI[index][nSeqsPerOTU[index]] = i;
			
			nSeqsPerOTU[index]++;
		}
        
		fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);	
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "setOTUs");
		exit(1);	
	}		
}
/**************************************************************************************************/

void ShhherCommand::writeQualities(int numOTUs, int numFlowCells, string qualityFileName, vector<int> otuCounts, vector<int>& nSeqsPerOTU, vector<int>& seqNumber,
                                   vector<double>& singleTau, vector<short>& flowDataIntI, vector<short>& uniqueFlowgrams, vector<int>& cumNumSeqs,
                                   vector<int>& mapUniqueToSeq, vector<string>& seqNameVector, vector<int>& centroids, vector<vector<int> >& aaI){
	
	try {
        
		ofstream qualityFile;
		util.openOutputFile(qualityFileName, qualityFile);
        
		qualityFile.setf(ios::fixed, ios::floatfield);
		qualityFile.setf(ios::showpoint);
		qualityFile << setprecision(6);
		
		vector<vector<int> > qualities(numOTUs);
		vector<double> pr(HOMOPS, 0);
		
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int index = 0;
            
			if(nSeqsPerOTU[i] > 0){
				
				while(index < numFlowCells){
                    
					double maxPrValue = 1e8;
					short maxPrIndex = -1;
					double count = 0.0000;
					
					pr.assign(HOMOPS, 0);
					
					for(int j=0;j<nSeqsPerOTU[i];j++){
						int lIndex = cumNumSeqs[i] + j;
						double tauValue = singleTau[seqNumber[lIndex]];
						int sequenceIndex = aaI[i][j];
						short intensity = flowDataIntI[sequenceIndex * numFlowCells + index];
						
						count += tauValue;
						
						for(int s=0;s<HOMOPS;s++){
							pr[s] += tauValue * singleLookUp[s * NUMBINS + intensity];
						}
					}
                    
					maxPrIndex = uniqueFlowgrams[centroids[i] * numFlowCells + index];
					maxPrValue = pr[maxPrIndex];
					
					if(count > MIN_COUNT){
						double U = 0.0000;
						double norm = 0.0000;
						
						for(int s=0;s<HOMOPS;s++){
							norm += exp(-(pr[s] - maxPrValue));
						}
                        
						for(int s=1;s<=maxPrIndex;s++){
							int value = 0;
							double temp = 0.0000;
							
							U += exp(-(pr[s-1]-maxPrValue))/norm;
							
							if(U>0.00){
								temp = log10(U);
							}
							else{
								temp = -10.1;
							}
							temp = floor(-10 * temp);
							value = (int)floor(temp);
							if(value > 100){	value = 100;	}
							
                            qualities[i].push_back((int)value);
						}
					}//end if
					
					index++;
                    
				}//end while
                
			}//end if
			
            
			if(otuCounts[i] > 0){
				qualityFile << '>' << seqNameVector[mapUniqueToSeq[i]] << endl;
                //need to get past the first four bases
                for (int j = 4; j < qualities[i].size(); j++) { qualityFile << qualities[i][j] << ' '; }
				qualityFile << endl;
			}
		}//end for
		qualityFile.close();
		outputNames.push_back(qualityFileName); outputTypes["qfile"].push_back(qualityFileName);
        
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeQualities");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeSequences(string thisCompositeFASTAFileName, int numOTUs, int numFlowCells, string fastaFileName, vector<int> otuCounts, vector<short>& uniqueFlowgrams, vector<string>& seqNameVector, vector<vector<int> >& aaI, vector<int>& centroids){
	try {
		
		ofstream fastaFile;
		util.openOutputFile(fastaFileName, fastaFile);
		
		vector<string> names(numOTUs, "");
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			int index = centroids[i];
			
			if(otuCounts[i] > 0){
				fastaFile << '>' << seqNameVector[aaI[i][0]] << endl;
				
				string newSeq = "";
				
				for(int j=0;j<numFlowCells;j++){
					
					char base = flowOrder[j % flowOrder.length()];
					for(int k=0;k<uniqueFlowgrams[index * numFlowCells + j];k++){
						newSeq += base;
					}
				}
				
				if (newSeq.length() >= 4) {  fastaFile << newSeq.substr(4) << endl;  }
                else {  fastaFile << "NNNN" << endl;  }
			}
		}
		fastaFile.close();
        
		outputNames.push_back(fastaFileName); outputTypes["fasta"].push_back(fastaFileName);
        
		if(thisCompositeFASTAFileName != ""){
			util.appendFiles(fastaFileName, thisCompositeFASTAFileName);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeSequences");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeNames(string thisCompositeNamesFileName, int numOTUs, string nameFileName, vector<int> otuCounts, vector<string>& seqNameVector, vector<vector<int> >& aaI, vector<int>& nSeqsPerOTU){
	try {
		
		ofstream nameFile;
		util.openOutputFile(nameFileName, nameFile);
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { break; }
			
			if(otuCounts[i] > 0){
				nameFile << seqNameVector[aaI[i][0]] << '\t' << seqNameVector[aaI[i][0]];
				
				for(int j=1;j<nSeqsPerOTU[i];j++){
					nameFile << ',' << seqNameVector[aaI[i][j]];
				}
				
				nameFile << endl;
			}
		}
		nameFile.close();
		outputNames.push_back(nameFileName); outputTypes["name"].push_back(nameFileName);
		
		
		if(thisCompositeNamesFileName != ""){
			util.appendFiles(nameFileName, thisCompositeNamesFileName);
		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeNames");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeGroups(string groupFileName, string fileRoot, int numSeqs, vector<string>& seqNameVector){
	try {
        ofstream groupFile;
		util.openOutputFile(groupFileName, groupFile);
		
		for(int i=0;i<numSeqs;i++){
			if (m->getControl_pressed()) { break; }
			groupFile << seqNameVector[i] << '\t' << fileRoot << endl;
		}
		groupFile.close();
		outputNames.push_back(groupFileName); outputTypes["group"].push_back(groupFileName);
        
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeGroups");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeClusters(string otuCountsFileName, int numOTUs, int numFlowCells, vector<int> otuCounts, vector<int>& centroids, vector<short>& uniqueFlowgrams, vector<string>& seqNameVector, vector<vector<int> >& aaI, vector<int>& nSeqsPerOTU, vector<int>& lengths, vector<short>& flowDataIntI){
	try {
		ofstream otuCountsFile;
		util.openOutputFile(otuCountsFileName, otuCountsFile);
		
		string bases = flowOrder;
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) {
				break;
			}
			//output the translated version of the centroid sequence for the otu
			if(otuCounts[i] > 0){
				int index = centroids[i];
				
				otuCountsFile << "ideal\t";
				for(int j=8;j<numFlowCells;j++){
					char base = bases[j % bases.length()];
					for(int s=0;s<uniqueFlowgrams[index * numFlowCells + j];s++){
						otuCountsFile << base;
					}
				}
				otuCountsFile << endl;
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int sequence = aaI[i][j];
					otuCountsFile << seqNameVector[sequence] << '\t';
					
					string newSeq = "";
					
					for(int k=0;k<lengths[sequence];k++){
						char base = bases[k % bases.length()];
						int freq = int(0.01 * (double)flowDataIntI[sequence * numFlowCells + k] + 0.5);
                        
						for(int s=0;s<freq;s++){
							newSeq += base;
							//otuCountsFile << base;
						}
					}
					
                    if (newSeq.length() >= 4) {  otuCountsFile << newSeq.substr(4) << endl;  }
                    else {  otuCountsFile << "NNNN" << endl;  }
				}
				otuCountsFile << endl;
			}
		}
		otuCountsFile.close();
		outputNames.push_back(otuCountsFileName); outputTypes["counts"].push_back(otuCountsFileName);
        
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeClusters");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::getSingleLookUp(){
	try{
		//	these are the -log probabilities that a signal corresponds to a particular homopolymer length
		singleLookUp.assign(HOMOPS * NUMBINS, 0);
		
		int index = 0;
		ifstream lookUpFile;
		util.openInputFile(lookupFileName, lookUpFile);
		
		for(int i=0;i<HOMOPS;i++){
			
			if (m->getControl_pressed()) { break; }
			
			float logFracFreq;
			lookUpFile >> logFracFreq;
			
			for(int j=0;j<NUMBINS;j++)	{
				lookUpFile >> singleLookUp[index];
				index++;
			}
		}	
		lookUpFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getSingleLookUp");
		exit(1);
	}
}

/**************************************************************************************************/

void ShhherCommand::getJointLookUp(){
	try{
		
		//	the most likely joint probability (-log) that two intenities have the same polymer length
		jointLookUp.resize(NUMBINS * NUMBINS, 0);
		
		for(int i=0;i<NUMBINS;i++){
			
			if (m->getControl_pressed()) { break; }
			
			for(int j=0;j<NUMBINS;j++){		
				
				double minSum = 100000000;
				
				for(int k=0;k<HOMOPS;k++){
					double sum = singleLookUp[k * NUMBINS + i] + singleLookUp[k * NUMBINS + j];
					
					if(sum < minSum)	{	minSum = sum;		}
				}	
				jointLookUp[i * NUMBINS + j] = minSum;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getJointLookUp");
		exit(1);
	}
}

/**************************************************************************************************/

double ShhherCommand::getProbIntensity(int intIntensity){                          
	try{

		double minNegLogProb = 100000000; 

		
		for(int i=0;i<HOMOPS;i++){//loop signal strength
			
			if (m->getControl_pressed()) { break; }
			
			float negLogProb = singleLookUp[i * NUMBINS + intIntensity];
			if(negLogProb < minNegLogProb)	{	minNegLogProb = negLogProb; }
		}
		
		return minNegLogProb;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getProbIntensity");
		exit(1);
	}
}




