/*
 *  mergefilecommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/14/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mergefilecommand.h"

//**********************************************************************************************************************
vector<string> MergeFileCommand::setParameters(){	
	try {
		CommandParameter pinput("input", "String", "", "", "", "", "","",false,true,true); parameters.push_back(pinput);
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,true,true); parameters.push_back(poutput);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter ptaxonomy("taxonomy", "", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptaxonomy);
        CommandParameter pfasta("fasta", "", "", "", "none", "none", "none","taxonomy",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MergeFileCommand::getHelpString(){	
	try {
		string helpString = "";
        helpString += "The merge.file command takes a list of files separated by dashes and appends them into one file. Altternatively, the merge file command can combine the data of several files. For example, you can combine a fasta, taxonomy and name or count field to achieve outputs like: GQY1XT001C44N8 3677 Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Porphyromonadaceae_unclassified; C-G--T-T--GA-A-A-C-T-G-G--CG-T-T-C--T-T-G-A-G-T-G-G-GC-GA-G-A-A-G-T-A--TG-C-GG-A-ATG-C-G-T-G-GT-GT-A-G-CGGT-G-AAA--...";
		helpString += "The merge.file command parameters are input and output or fasta, taxonomy, name and count.";
		helpString += "Example merge.file(input=small.fasta-large.fasta, output=all.fasta).";
        helpString += "Example merge.file(fasta=final.fasta, name=final.names, taxonomy=final.taxonomy).";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MergeFileCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],merged,[extension]"; }
        else if (type == "merge") { pattern = ""; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeFileCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MergeFileCommand::MergeFileCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["merge"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "MergeFileCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MergeFileCommand::MergeFileCommand(string option)  {
	try {
        abort = false; calledHelp = false;   appendMode = true;
		
		if(option == "help") {
			help();
			abort = true; calledHelp = true;
		}else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["merge"] = tempOutNames;
            outputTypes["fasta"] = tempOutNames;
			
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found")	{	outputDir = "";		}

            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){    inputDir = "";        }
            
			string fileList = validParameter.valid(parameters, "input");
            if(fileList == "not found") { appendMode = false; fileList = "";  }
			else{ 	util.splitAtDash(fileList, fileNames);	}
            
            outputFileName = validParameter.valid(parameters, "output");
            if (outputFileName == "not found") { appendMode = false; outputFileName = "";   }

            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { fastafile = ""; abort = true; }
            else if (fastafile == "not found") {  fastafile = "";  }
            else { current->setFastaFile(fastafile); appendMode = false; }
            
            namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { namefile = ""; abort = true; }
            else if (namefile == "not found") {  namefile = "";  }
            else { current->setNameFile(namefile); appendMode = false; }
            
            taxfile = validParameter.validFile(parameters, "taxonomy");
            if (taxfile == "not open") { taxfile = ""; abort = true; }
            else if (taxfile == "not found") {  taxfile = "";  }
            else { current->setTaxonomyFile(taxfile); appendMode = false; }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }	
            else { current->setCountFile(countfile); appendMode = false; }
            
            if (!appendMode) { //if you are not appending, fasta is required as well as at least one of taxonomy, name or count
                if (fastafile == "") { //look for current
                    fastafile = current->getFastaFile();
                    if (fastafile != "") {  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                    else { 	m->mothurOut("[ERROR]: You have no current fastafile and the fasta parameter is required.\n"); abort = true; }
                }
                
                if ((namefile == "") && (countfile == "") && (taxfile == "")) {
                    taxfile = current->getTaxonomyFile();
                    if (taxfile != "") {  m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter.\n");  }
                    else {
                        countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                        else {
                            namefile = current->getNameFile();
                            if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
                            else { m->mothurOut("[ERROR]: You have no current taxonomy, name or count files. At least one is required. \n"); abort = true; }
                        }
                    }
                }
                
                if ((namefile != "") && (countfile != "")) {
                    m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n"); abort = true;
                }
            }else {
                numInputFiles = fileNames.size();
                ifstream testFile;
                if(numInputFiles == 0){  m->mothurOut("you must enter two or more file names and you entered " + toString(fileNames.size()) +  " file names\n"); abort=true; }
                else{
                    for(int i=0;i<numInputFiles;i++){
                        string path = util.hasPath(fileNames[i]);
                        if (inputDir != "") {
                            //if the user has not given a path then, add inputdir. else leave path alone.
                            if (path == "") { fileNames[i] = inputDir + fileNames[i]; }
                        }
                        
                        if (util.checkLocations(fileNames[i], current->getLocations())) { }
                        else { fileNames.erase(fileNames.begin()+i); i--; } //erase from file list
                        
                        path = util.hasPath(fileNames[i]);
                        if (path != "") { if (outputDir == "") { outputDir = path; } }
                    }
                    if (outputDir != "") { outputFileName = outputDir + util.getSimpleName(outputFileName);  }
                }
            }
			
        }
			
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "MergeFileCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeFileCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        if (appendMode) {
            util.mothurRemove(outputFileName);
            for(int i=0;i<numInputFiles;i++){  util.appendFiles(fileNames[i], outputFileName);  }
            outputNames.push_back(outputFileName); outputTypes["merge"].push_back(outputFileName);
        }else { outputFileName = mergeFileData(); }
		
        if (m->getControl_pressed()) {  util.mothurRemove(outputFileName); return 0;  }
        
        //set taxonomy file as new current taxonomyfile
        string currentName = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); } }

        
		m->mothurOut("\nOutput File Names: \n"); 
		m->mothurOut(outputFileName); m->mothurOutEndLine();
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

string MergeFileCommand::mergeFileData(){
    try {
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
        outputFileName = getOutputFileName("fasta", variables);
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        //extract seq counts from name or count file
        map<string, int> nameMap; map<string, int>::iterator itCount;
        bool useNameMap = true;
        if (countfile != "") {
            CountTable ct; ct.readTable(countfile, false, false);
            nameMap = ct.getNameMap();
        }else if (namefile != "") {
            nameMap = util.readNames(namefile);
        }else { useNameMap = false;  }
        
        map<string, string> taxMap; map<string, string>::iterator itTax;
        bool useTax = false;
        if (taxfile != "") {  util.readTax(taxfile, taxMap, false);  useTax = true;  }
        
        ifstream in;
        util.openInputFile(fastafile, in);

        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            Sequence currSeq(in); util.gobble(in);
            
            string comment = " ";
            
            if (useNameMap) {
                itCount = nameMap.find(currSeq.getName());
                if (itCount != nameMap.end()) {
                    comment += toString(itCount->second) + " ";
                    nameMap.erase(itCount);
                }else { m->mothurOut("[ERROR]: Missing count data for " + currSeq.getName() + ", please correct.\n"); m->setControl_pressed(true); }
            }
            
            
            if (useTax) {
                itTax = taxMap.find(currSeq.getName());
                if (itTax != taxMap.end()) {
                    comment += itTax->second;
                    taxMap.erase(itTax);
                }else { m->mothurOut("[ERROR]: Missing taxonomy for " + currSeq.getName() + ", please correct.\n"); m->setControl_pressed(true); }
            }
            
            currSeq.setComment(comment);
            currSeq.printSequence(out);
        }
        in.close(); out.close();
        
        return outputFileName;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeFileCommand", "mergeFileData");
        exit(1);
    }
}
//**********************************************************************************************************************
