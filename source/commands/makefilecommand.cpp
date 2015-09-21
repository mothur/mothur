//
//  makefilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/24/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "makefilecommand.h"

//**********************************************************************************************************************
vector<string> MakeFileCommand::setParameters(){
    try {
        CommandParameter ptype("type", "Multiple", "fastq-gz", "fastq", "", "", "","",false,false); parameters.push_back(ptype);
        CommandParameter pnumcols("numcols", "Multiple", "2-3", "3", "", "", "","",false,false, true); parameters.push_back(pnumcols);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MakeFileCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The make.file command takes a input directory and creates a file file containing the fastq or gz files in the directory.\n";
        helpString += "The make.fastq command parameters are inputdir, numcols and type.  inputdir is required.\n";
        helpString += "May create more than one file. Mothur will attempt to match paired files. \n";
        helpString += "The type parameter allows you to set the type of files to look for. Options are fastq or gz.  Default=fastq. \n";
        helpString += "The numcols parameter allows you to set number of columns you mothur to make in the file.  Default=3, meaning groupName forwardFastq reverseFastq. The groupName is made from the beginning part of the forwardFastq file. Everything up to the first '_' or if no '_' is found then the root of the forwardFastq filename.\n";
        helpString += "The make.file command should be in the following format: \n";
        helpString += "make.file(inputdir=yourInputDirectory). \n";
        helpString += "Example make.group(inputdir=fastqFiles)\n";
        helpString += "Note: No spaces between parameter labels (i.e. inputdir), '=' and parameters (i.e. yourInputDirectory).\n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string MakeFileCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "file") {  pattern = "[filename],[tag],file"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeFileCommand::MakeFileCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["file"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "MakeFileCommand");
        exit(1);
    }
}

//**********************************************************************************************************************

MakeFileCommand::MakeFileCommand(string option)  {
    try {
        
        abort = false; calledHelp = false;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            map<string, string>::iterator it;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
            }
            
            //initialize outputTypes
            vector<string> tempOutNames;
            outputTypes["file"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            inputDir = validParameter.validFile(parameters, "inputdir", false);
            if (inputDir == "not found"){	inputDir = "";	m->mothurOut("[ERROR]: The inputdir parameter is required, aborting."); m->mothurOutEndLine(); abort = true;	}
            else {
                if (m->dirCheck(inputDir)) {} // all set
                else { abort = true; }
            }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = inputDir;		}
            
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            typeFile = validParameter.validFile(parameters, "type", false);
            if (typeFile == "not found"){	typeFile = "fastq";		}
            
            if ((typeFile != "fastq") && (typeFile != "gz")) { m->mothurOut(typeFile + " is not a valid type. Options are fastq or gz. I will use fastq."); m->mothurOutEndLine(); typeFile = "fastq"; }
            
            string temp = validParameter.validFile(parameters, "numcols", false);		if(temp == "not found"){	temp = "3"; }
            if ((temp != "2") && (temp != "3")) { m->mothurOut(temp + " is not a valid numcols. Options are 2 or 3. I will use 3."); m->mothurOutEndLine(); temp = "3";  }
            m->mothurConvert(temp, numCols);
            
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "MakeFileCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int MakeFileCommand::execute(){
    try {
        if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //find all .fastq files
        string tempFile = inputDir + "fileList.temp";
        string findCommand = "find \"" + inputDir.substr(0, inputDir.length()-1) + "\" -maxdepth 1 -name \"*." + typeFile + "\"  > \"" + tempFile + "\"";
        if (m->debug) { m->mothurOut(findCommand + "\n"); }
        system(findCommand.c_str());
        
        //read in list of files
        vector<string> fastqFiles;
        m->readAccnos(tempFile, fastqFiles, "no error");
        m->mothurRemove(tempFile);
        
        if (fastqFiles.size() == 0) { m->mothurOut("[WARNING]: Unable to find any " + typeFile + " files in your directory.\n"); }
        else {
            
            //sort into alpha order to put pairs togther if they exist
            sort(fastqFiles.begin(), fastqFiles.end());
            
            vector< vector<string> > paired;
            vector<string> singles;
            string lastFile = "";
            for (int i = 0; i < fastqFiles.size()-1; i++) {
                
                if (m->debug) { m->mothurOut("[DEBUG]: File " + toString(i) + " = " + fastqFiles[i] + ".\n"); }
                
                if (m->control_pressed) { break; }
                
                string simpleName1 = m->getRootName(m->getSimpleName(fastqFiles[i]));
                string simpleName2 = m->getRootName(m->getSimpleName(fastqFiles[i+1]));
                
                //possible pair
                if (simpleName1.length() == simpleName2.length()) {
                    int numDiffs = 0;
                    for (int j = 0; j < simpleName1.length(); j++) {
                        if (numDiffs > 1) { break; }
                        else if (simpleName1[j] != simpleName2[j]) { numDiffs++; }
                    }
                    if (numDiffs > 1) { singles.push_back(fastqFiles[i]); lastFile = fastqFiles[i]; }
                    else { //only one diff = paired files
                        int pos = simpleName1.find("R1");
                        int pos2 = simpleName2.find("R2");
                        if ((pos != string::npos) && (pos2 != string::npos)){
                            vector<string> temp;
                            if (numCols == 3) {
                                string groupName = "noGroup"+toString(i);
                                int posUnderscore = fastqFiles[i].find_first_of('_');
                                if (posUnderscore == string::npos) {   groupName = m->getSimpleName(m->getRootName(fastqFiles[i]));  }
                                else{  groupName = m->getSimpleName(fastqFiles[i].substr(0, posUnderscore));  }
                                temp.push_back(groupName);
                            }
                            temp.push_back(fastqFiles[i]); temp.push_back(fastqFiles[i+1]); lastFile = fastqFiles[i+1];
                            paired.push_back(temp);
                            i++;
                        }else {
                            singles.push_back(fastqFiles[i]); lastFile = fastqFiles[i];
                        }
                    }
                }else{
                    singles.push_back(fastqFiles[i]); lastFile = fastqFiles[i];
                }
            }
            if (lastFile != fastqFiles[fastqFiles.size()-1]) { singles.push_back(fastqFiles[fastqFiles.size()-1]); }
            
            if (singles.size() != 0) {
                map<string, string> variables;
                variables["[filename]"] = outputDir + "fileList.";
                variables["[tag]"] = "single";
                string filename = getOutputFileName("file",variables);
                ofstream out;
                m->openOutputFile(filename, out);
                outputNames.push_back(filename); outputTypes["file"].push_back(filename);
                m->setFileFile(filename);
                
                for (int i = 0; i < singles.size(); i++) {
                    out << singles[i] << endl;
                }
                out.close();
            }
            
            if (paired.size() != 0) {
                map<string, string> variables;
                variables["[filename]"] = outputDir + "fileList.";
                variables["[tag]"] = "paired";
                string filename = getOutputFileName("file",variables);
                ofstream out;
                m->openOutputFile(filename, out);
                outputNames.push_back(filename); outputTypes["file"].push_back(filename);
                m->setFileFile(filename);
                
                for (int i = 0; i < paired.size(); i++) {
                    for (int j = 0; j < paired[i].size(); j++) {
                        out << paired[i][j] << '\t';
                    }
                    out << endl;
                }
                out.close();
            }
            
        }
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************


