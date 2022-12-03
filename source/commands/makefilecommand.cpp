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
        CommandParameter ptype("type", "Multiple", "fasta-fastq-gz", "fastq", "", "", "","",false,false); parameters.push_back(ptype);
        CommandParameter pnumcols("numcols", "Multiple", "2-3", "3", "", "", "","",false,false, true); parameters.push_back(pnumcols);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pprefix("prefix", "String", "", "", "", "", "","",false,false); parameters.push_back(pprefix);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pdelim("delim", "String", "", "_", "", "", "","",false,false); parameters.push_back(pdelim);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["file"] = tempOutNames;
        
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
        helpString += "The make.file command takes a input directory and creates a file file containing the fasta, fastq or gz files in the directory.\n";
        helpString += "The make.file command parameters are inputdir, numcols, type and prefix.  inputdir is required.\n";
        helpString += "May create more than one file. Mothur will attempt to match paired files. \n";
        helpString += "The type parameter allows you to set the type of files to look for. Options are fasta, fastq or gz.  Default=fastq. \n";
        helpString += "The numcols parameter allows you to set number of columns you mothur to make in the file.  Default=3 for fastq and gz, meaning groupName forwardFastq reverseFastq. The groupName is made from the beginning part of the forwardFastq file. Everything up to the first '_' or if no '_' is found then the root of the forwardFastq filename. Default=2 for fasta files, meaning groupName fastaFile.\n";
        helpString += "The prefix parameter allows you to enter your own prefix for the output filename. Default=stability.";
        helpString += "The delim parameter allow you to enter the character you would like to use to create the sample name. Default='_'. For example, M6D7_S163_L001_R2_001.fastq.gz would produce the sample name M6D7. Set delim=* to indicate you want mothur to create unique names for each file pair. (no pooling)\n";
        helpString += "The make.file command should be in the following format: \n";
        helpString += "make.file(inputdir=yourInputDirectory). \n";
        helpString += "Example make.file(inputdir=fastqFiles)\n";
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
MakeFileCommand::MakeFileCommand(string option) : Command()  {
    try {
         if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            inputDir = validParameter.validPath(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";	m->mothurOut("[ERROR]: The inputdir parameter is required, aborting.\n");  abort = true;	}
            else {
                if (util.dirCheckExists(inputDir)) {} // all set
                else { abort = true; }
            }
            
            if (outputdir == ""){    outputdir = inputDir;		}
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            typeFile = validParameter.valid(parameters, "type");
            if (typeFile == "not found"){	typeFile = "fastq";		}
            
            if ((typeFile != "fastq") && (typeFile != "gz") && (typeFile != "fasta")) { m->mothurOut(typeFile + " is not a valid type. Options are fasta, fastq or gz. I will use fastq.\n");  typeFile = "fastq"; }
            
            string temp = validParameter.valid(parameters, "numcols");		if(temp == "not found"){ temp = "3"; if (typeFile == "fasta") {	temp = "2"; } }
            if ((temp != "2") && (temp != "3")) { m->mothurOut(temp + " is not a valid numcols. Options are 2 or 3. I will use 3.\n");  temp = "3";  }
            util.mothurConvert(temp, numCols);
            cout << numCols << endl;
            prefix = validParameter.valid(parameters, "prefix");		if (prefix == "not found") { prefix = "stability"; }
            
            delim = validParameter.valid(parameters, "delim");			if (delim == "not found") { delim = "_"; }
            
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
        vector<string> files;
        util.readAccnos(tempFile, files, "no error");
        util.mothurRemove(tempFile);
        
        if (m->getDebug()) {
            m->mothurOut("[DEBUG]: Found " + toString(files.size()) + " files of type " + typeFile + ".\n");
            for (int i = 0; i < files.size(); i++) { m->mothurOut("[DEBUG]: " + toString(i) + " = " + files[i] + "\n");}
        }
        
        if (files.size() == 0) { m->mothurOut("[WARNING]: Unable to find any " + typeFile + " files in your directory.\n"); }
        else {
            
            //sort into alpha order to put pairs togther if they exist
            sort(files.begin(), files.end());
            
            if (typeFile == "fasta")    { createFastaFile(files);   }
            else                        { createFastqGzFile(files); }
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]+"\n"); 	}
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

void MakeFileCommand::createFastaFile(vector<string> files){
    try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + prefix + ".";
        string filename = getOutputFileName("file",variables);
        
        ofstream out; util.openOutputFile(filename, out);
        outputNames.push_back(filename); outputTypes["file"].push_back(filename);
        current->setFileFile(filename);
        
        vector< vector<string> > outputs;
        if (numCols == 3) {
            vector<string> singles;
            vector< vector<string> > paired = pairFiles(files, singles);
            outputs = prepareOutputs(paired);
        }else {
            outputs = prepareOutputs(files);
        }
        
        for (int i = 0; i < outputs.size(); i++) {
            for (int j = 0; j < outputs[i].size(); j++) { out << outputs[i][j] << '\t'; } out << endl;
        }
        out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "createFastaFile");
        exit(1);
    }
}
//**********************************************************************************************************************

void MakeFileCommand::createFastqGzFile(vector<string> files){
    try {
        vector<string> singles;
        vector< vector<string> > paired = pairFiles(files, singles);
        
        if (singles.size() != 0) {
            map<string, string> variables;
            variables["[filename]"] = outputdir + prefix + ".";
            if (paired.size() != 0) { variables["[tag]"] = "single"; }
            string filename = getOutputFileName("file",variables);
            
            ofstream out; util.openOutputFile(filename, out);
            
            for (int i = 0; i < singles.size(); i++) { out << singles[i] << endl; } out.close();
            
            if (util.isBlank(filename)) {  util.mothurRemove(filename); }
            else {
                outputNames.push_back(filename); outputTypes["file"].push_back(filename);
                m->mothurOut("\n[WARNNG]: mothur found unpaired files in your input directory. Outputting list of filenames to " + filename + " for your review.\n\n");
            }
        }
        
        //generates unique group names
        if (numCols == 3) { paired = prepareOutputs(paired); }

        if (paired.size() != 0) {
            map<string, string> variables;
            variables["[filename]"] = outputdir + prefix + ".";
            string filename = getOutputFileName("file",variables);
            
            ofstream out; util.openOutputFile(filename, out);
            outputNames.push_back(filename); outputTypes["file"].push_back(filename);
            current->setFileFile(filename);
            
            for (int i = 0; i < paired.size(); i++) {
                for (int j = 0; j < paired[i].size(); j++) { out << paired[i][j] << '\t'; } out << endl;
            }
            out.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "createFastqGzFile");
        exit(1);
    }
}
//**********************************************************************************************************************
vector< vector<string> > MakeFileCommand::pairFiles(vector<string> files, vector<string>& singles){
    try {
        vector< vector<string> > paired;
        string lastFile = "";
        
        for (int i = 0; i < files.size()-1; i++) {
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: File " + toString(i) + " = " + files[i] + ".\n"); }
            
            if (m->getControl_pressed()) { break; }
            
            string simpleName1 = util.getRootName(util.getSimpleName(files[i]));
            string simpleName2 = util.getRootName(util.getSimpleName(files[i+1]));
            
            //possible pair
            if (simpleName1.length() == simpleName2.length()) {
                int numDiffs = 0;
                for (int j = 0; j < simpleName1.length(); j++) {
                    if (numDiffs > 1) { break; }
                    else if (simpleName1[j] != simpleName2[j]) { numDiffs++; }
                }
                if (numDiffs > 1) { singles.push_back(util.getSimpleName(files[i])); lastFile = files[i]; }
                else { //only one diff = paired files
                    vector<string> temp;
                    temp.push_back(util.getSimpleName(files[i])); temp.push_back(util.getSimpleName(files[i+1])); lastFile = files[i+1];
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: Pairing " + files[i] + " with " + files[i+1] + ".\n"); }
                    paired.push_back(temp);
                    i++;
                }
            }else{
                if (m->getDebug()) { m->mothurOut("[DEBUG]: Adding single " + files[i] + ".\n"); }
                singles.push_back(util.getSimpleName(files[i])); lastFile = files[i];
            }
        }
        
        if (lastFile != files[files.size()-1]) {
            if (m->getDebug()) { m->mothurOut("[DEBUG]: Adding single " + files[files.size()-1] + ".\n"); }
            singles.push_back(util.getSimpleName(files[files.size()-1]));
        }
        
        return paired;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "pairFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
//groupName defaults to "noGroup"+toString(i);
vector< vector<string> > MakeFileCommand::prepareOutputs(vector<string> files){
    try {
        vector< vector<string> > results; results.resize(files.size());
        
        vector<string> groupNames;
        if (delim == "*") {
            groupNames = denoiseGroupNames(files);
        }else {
            groupNames = createDefaultGroupNames(files);
        }
        
        for (int i = 0; i < files.size(); i++) {
            results[i].push_back(groupNames[i]); results[i].push_back(util.getSimpleName(files[i]));
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "prepareOutputs");
        exit(1);
    }
}
//**********************************************************************************************************************
//groupName defaults to "noGroup"+toString(i);
vector< vector<string> > MakeFileCommand::prepareOutputs(vector< vector<string> > paired){
    try {
        vector< vector<string> > results; results.resize(paired.size());
        
        vector<string> forwardNames;
        for (int i = 0; i < paired.size(); i++) {
            forwardNames.push_back(util.getRootName(util.getSimpleName(paired[i][0])));
        }
        
        vector<string> groupNames;
        if (delim == "*") {
            groupNames = denoiseGroupNames(forwardNames);
        }else {
            groupNames = createDefaultGroupNames(forwardNames);
        }
        
        for (int i = 0; i < paired.size(); i++) {
            results[i].push_back(groupNames[i]); results[i].push_back(paired[i][0]); results[i].push_back(paired[i][1]);
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "prepareOutputs");
        exit(1);
    }
}
//**********************************************************************************************************************

vector<string> MakeFileCommand::createDefaultGroupNames(vector<string> files){
    try {
        set<string> groups;
        vector<string> groupNames;
        for (int i = 0; i < files.size(); i++) {
            
            string groupName = "Group_" + toString(i);
            string filename = util.getSimpleName(files[i]);
            int pos = filename.find(delim);
            
            if (pos != string::npos) { groupName = filename.substr(0, pos); }
            
            if (groups.count(groupName) == 0) { util.checkName(groupName); groups.insert(groupName);  }
            else {
                //look for another delim
                string tempFilename = filename.substr(pos+1); //grab rest of name
                pos = tempFilename.find(delim);
                if (pos != string::npos) {
                    groupName += "_" + tempFilename.substr(0, pos);
                    if (groups.count(groupName) != 0) {  groupName += "_"+ toString(i);  } //already have this name
                }
                else { groupName += "_"+ toString(i); }
                util.checkName(groupName);
                groups.insert(groupName);
            }
            
            groupNames.push_back(groupName);
        }
        
        return groupNames;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "createGroupNames");
        exit(1);
    }
}
//**********************************************************************************************************************

vector<string> MakeFileCommand::denoiseGroupNames(vector<string> files){
    try {
        //remove any "words" in filenames that is the same in all filenames separated by delim(_ .)
        //MI.M00833_0261.001.FLD0207.TRIN-META_16S_R2.fastq
        //MI.M00833_0261.001.FLD0223.ERIFF-META_16S_R2.fastq
        //would become...
        //FLD0207.TRIN-META
        //FLD0223.ERIFF-META
                     
        //split all names into pieces
        vector<vector<string> > words; words.resize(files.size());
        map<int, set<string> > posToWord;
        
        for (int i = 0; i < files.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            string filename = util.getRootName(util.getSimpleName(files[i]));
            
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
            if (!util.allSpaces(individual)) {
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
        vector<string> orderedGroups;
        for (int i = 0; i < words.size(); i++) {
            
            //assemble groupNames
            string groupName = "";
            for (int j = 0; j < words[i].size(); j++) {
                //include word
                if (goodIndexes.count(j) != 0) { groupName += words[i][j] + "_"; }
            }
            
            if (groupName != "") { groupName = groupName.substr(0, groupName.length()-1); }
           
            //is this name unique
            if (groups.count(groupName) == 0) { util.checkName(groupName); groups.insert(groupName);  }
            else { groupName = "Group_"+ toString(i); util.checkName(groupName); groups.insert(groupName); }
            
            orderedGroups.push_back(groupName);
        }
        
        return orderedGroups;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "denoiseGroupNames");
        exit(1);
    }
}
//**********************************************************************************************************************

int MakeFileCommand::fillAccnosFile(string tempFile){
    try {
        string findCommand = "";
        string tempOut = tempFile;
        tempFile = "\"" + tempFile + "\"";
        string wrappedInput = "\"" + inputDir + "\"";
        
#if defined NON_WINDOWS

        findCommand = "ls " + wrappedInput + "*." + typeFile + " > " + tempFile;
        //findCommand = "find \"" + inputDir.substr(0, inputDir.length()-1) + "\" -maxdepth 1 -name \"*." + typeFile + "\" > " + tempFile;
        if (m->getDebug()) { m->mothurOut(findCommand + "\n"); }
        system(findCommand.c_str());
        
#else
        //use ls command
        findCommand = "dir /B "  + wrappedInput + "*." + typeFile + " > " + tempFile;
        //findCommand = "dir /B \""  + inputDir.substr(0, inputDir.length()-1) + "\\*.\"" + typeFile + " > " + tempFile + "\"";
        if (m->getDebug()) { m->mothurOut(findCommand + "\n"); }
        system(findCommand.c_str());
        
        tempOut += ".temp";
        tempFile = tempFile.substr(1, tempFile.length()-2); //remove ""
        
        ofstream out; util.openOutputFile(tempOut, out);
        ifstream in; util.openInputFile(tempFile, in);
        
        string junk, filename;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
           
            in >> filename; gobble(in);

            if (m->getDebug()) { m->mothurOut("[DEBUG]: Found file " + filename + ".\n"); }
            //ignore hidden files
            if (filename[0] != '.') { out << filename << endl; }
        }
        in.close();
        out.close();
        
        util.mothurRemove(tempFile);
        util.renameFile(tempOut, tempFile);
        
#endif
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFileCommand", "fillAccnosFile");
        exit(1);
    }
}
//**********************************************************************************************************************


