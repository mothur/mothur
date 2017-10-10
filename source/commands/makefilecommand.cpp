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
        CommandParameter pprefix("prefix", "String", "", "", "", "", "","",false,false); parameters.push_back(pprefix);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pdelim("delim", "String", "", "_", "", "", "","",false,false); parameters.push_back(pdelim);
        
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
        helpString += "The make.file command parameters are inputdir, numcols, type and prefix.  inputdir is required.\n";
        helpString += "May create more than one file. Mothur will attempt to match paired files. \n";
        helpString += "The type parameter allows you to set the type of files to look for. Options are fastq or gz.  Default=fastq. \n";
        helpString += "The numcols parameter allows you to set number of columns you mothur to make in the file.  Default=3, meaning groupName forwardFastq reverseFastq. The groupName is made from the beginning part of the forwardFastq file. Everything up to the first '_' or if no '_' is found then the root of the forwardFastq filename.\n";
        helpString += "The prefix parameter allows you to enter your own prefix for the output filename. Default=stability.";
        helpString += "The delim parameter allow you to enter the character you would like to use to create the sample name. Default='_'. For example, M6D7_S163_L001_R2_001.fastq.gz would produce the sample name M6D7. Set delim=* to indicate you want mothur to create unique names for each file pair. (no pooling)\n";
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
        
        if (type == "file") {  pattern = "[filename],[tag],files-[filename],files"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
            
            prefix = validParameter.validFile(parameters, "prefix", false);		if (prefix == "not found") { prefix = "stability"; }
            
            delim = validParameter.validFile(parameters, "delim", false);			if (delim == "not found") { delim = "_"; }
            
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
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        //find all .fastq files
        string tempFile = inputDir + "fileList.temp";
        fillAccnosFile(tempFile);
        
        //read in list of files
        vector<string> fastqFiles;
        m->readAccnos(tempFile, fastqFiles, "no error");
        m->mothurRemove(tempFile);
        
        if (m->getDebug()) {
            m->mothurOut("[DEBUG]: Found " + toString(fastqFiles.size()) + " files of type " + typeFile + ".\n");
            for (int i = 0; i < fastqFiles.size(); i++) { m->mothurOut("[DEBUG]: " + toString(i) + " = " + fastqFiles[i] + "\n");}
        }
        
        if (fastqFiles.size() == 0) { m->mothurOut("[WARNING]: Unable to find any " + typeFile + " files in your directory.\n"); }
        else {
            
            //sort into alpha order to put pairs togther if they exist
            sort(fastqFiles.begin(), fastqFiles.end());
            
            vector< vector<string> > paired;
            vector<string> singles;
            set<string> groups;
            string lastFile = "";
            for (int i = 0; i < fastqFiles.size()-1; i++) {
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(i) + ".\n"); }
                if (m->getDebug()) { m->mothurOut("[DEBUG]: File " + toString(i) + " = " + fastqFiles[i] + ".\n"); }
                
                if (m->getControl_pressed()) { break; }
                
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
                        vector<string> temp;
                        temp.push_back(fastqFiles[i]); temp.push_back(fastqFiles[i+1]); lastFile = fastqFiles[i+1];
                        if (m->getDebug()) { m->mothurOut("[DEBUG]: Pairing " + fastqFiles[i] + " with " + fastqFiles[i+1] + ".\n"); }
                        paired.push_back(temp);
                        i++;
                    }
                }else{
                    singles.push_back(fastqFiles[i]); lastFile = fastqFiles[i];
                }
            }
            if (lastFile != fastqFiles[fastqFiles.size()-1]) { singles.push_back(fastqFiles[fastqFiles.size()-1]); }
            
            if (singles.size() != 0) {
                map<string, string> variables;
                variables["[filename]"] = outputDir + prefix + ".";
                if (paired.size() != 0) { variables["[tag]"] = "single"; }
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
            
            //generates unique group names
            if (numCols == 3) { paired = findGroupNames(paired); }

            if (paired.size() != 0) {
                map<string, string> variables;
                variables["[filename]"] = outputDir + prefix + ".";
                if (singles.size() != 0) { variables["[tag]"] = "paired"; }
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
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
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
//groupName defaults to "noGroup"+toString(i);
vector< vector<string> > MakeFileCommand::findGroupNames(vector< vector<string> > paired){
    try {
        vector< vector<string> > results; results.resize(paired.size());
        
        if (delim == "*") {
            //remove any "words" in filenames that is the same in all filenames separated by delim(_ .)
            //MI.M00833_0261.001.FLD0207.TRIN-META_16S_R2.fastq
            //MI.M00833_0261.001.FLD0223.ERIFF-META_16S_R2.fastq
            //would become...
            //FLD0207.TRIN-META
            //FLD0223.ERIFF-META
            
            //split all forward names into pieces
            vector<vector<string> > words; words.resize(paired.size());
            map<int, set<string> > posToWord;
            
            for (int i = 0; i < paired.size(); i++) {
                if (m->getControl_pressed()) { break; }
                
                string filename = m->getRootName(m->getSimpleName(paired[i][0]));
                
                int pos = 0;
                string individual = "";
                for(int j=0;j<filename.length();j++){
                    if((filename[j] == '.') || (filename[j] == '_')){
                        words[i].push_back(individual);
                        
                        map<int, set<string> >::iterator it = posToWord.find(pos);
                        if (it != posToWord.end()) { posToWord[pos].insert(individual); }
                        else {
                            set<string> temp; temp.insert(individual);
                            posToWord[pos] = temp;
                        }
                        individual = "";
                        pos++;
                    }
                    else{
                        individual += filename[j];
                    }
                }
                if (!m->allSpaces(individual)) {
                    words[i].push_back(individual);
                    map<int, set<string> >::iterator it = posToWord.find(pos);
                    if (it != posToWord.end()) { posToWord[pos].insert(individual); }
                    else {
                        set<string> temp; temp.insert(individual);
                        posToWord[pos] = temp;
                    }
                }
            }
            
            //remove duplicate pieces
            set<int> goodIndexes;
            for (map<int, set<string> >::iterator it = posToWord.begin(); it != posToWord.end(); it++) {
                set<string> w = it->second;;
                if (w.size() != 1) { goodIndexes.insert(it->first);  }
            }
            
            set<string> groups;
            for (int i = 0; i < words.size(); i++) {
                
                //assemble groupNames
                string groupName = "";
                for (int j = 0; j < words[i].size(); j++) {
                    //include word
                    if (goodIndexes.count(j) != 0) { groupName += words[i][j] + "_"; }
                }
                
                if (groupName != "") { groupName = groupName.substr(0, groupName.length()-1); }
                
                //is this name unique
                if (groups.count(groupName) == 0) {  groups.insert(groupName);  }
                else { groupName = "Group_"+ toString(i); groups.insert(groupName); }
                
                results[i].push_back(groupName); results[i].push_back(paired[i][0]); results[i].push_back(paired[i][1]);
            }
            
        }else { //separate by the user selected deliminator. default='_'
            for (int i = 0; i < paired.size(); i++) {
                
                string groupName = "Group_" + toString(i);
                string filename = m->getSimpleName(paired[i][0]);
                int pos = filename.find(delim);
                
                if (pos != string::npos) { groupName = filename.substr(0, pos); }
                
                results[i].push_back(groupName); results[i].push_back(paired[i][0]); results[i].push_back(paired[i][1]);
            }
        }
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "findGroupName");
        exit(1);
    }
}

//**********************************************************************************************************************

int MakeFileCommand::fillAccnosFile(string tempFile){
    try {
        
        string findCommand = "";
        string tempOut = tempFile;
        tempFile = "\"" + tempFile + "\"";
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)

        findCommand = "find \"" + inputDir.substr(0, inputDir.length()-1) + "\" -maxdepth 1 -name \"*." + typeFile + "\" > " + tempFile;
        if (m->getDebug()) { m->mothurOut(findCommand + "\n"); }
        system(findCommand.c_str());
        
#else
        //use ls command
        findCommand = "find \"\" "  + inputDir.substr(0, inputDir.length()-1) + "\\*." + typeFile + " > " + tempFile;
        if (m->getDebug()) { m->mothurOut(findCommand + "\n"); }
        system(findCommand.c_str());

        ifstream in;
        ofstream out;
        tempOut += ".temp";
        
        tempFile = tempFile.substr(1, tempFile.length()-2); //remove ""
        
        m->openOutputFile(tempOut, out);
        
        m->openInputFile(tempFile, in);
        
        string junk, filename;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            in >> junk; m->gobble(in);
            in >> filename; m->gobble(in);
            
            out << filename << endl;
        }
        in.close();
        out.close();
        
        m->mothurRemove(tempFile);
        m->renameFile(tempOut, tempFile);
        
#endif
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "fillAccnosFile");
        exit(1);
    }
}
//**********************************************************************************************************************


