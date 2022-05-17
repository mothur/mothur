//
//  mergesfffilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/31/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "mergesfffilecommand.h"
#include "endiannessmacros.h"

//********************************************************************************
MergeSfffilesCommand::~MergeSfffilesCommand(){
    for (int i = 0; i < commonHeaders.size(); i++) { delete commonHeaders[i]; }
    commonHeaders.clear();
}
//********************************************************************************
vector<string> MergeSfffilesCommand::setParameters(){
	try {
		CommandParameter psff("sff", "InputTypes", "", "", "sffFile", "sffFile", "none","sff",false,false); parameters.push_back(psff);
        CommandParameter pfile("file", "InputTypes", "", "", "sffFile", "sffFile", "none","sff",false,false); parameters.push_back(pfile);
        CommandParameter pkeytrim("keytrim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeytrim);
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,true,true); parameters.push_back(poutput);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["sff"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************
string MergeSfffilesCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The merge.sfffiles command reads a sff file or a file containing a list of sff files and merges the individual files into a single sff file. \n";
		helpString += "The merge.sfffiles command parameters are sff, file and output. sff or file is required. \n";
		helpString += "The sff parameter allows you to enter the sff list of sff files separated by -'s.\n";
		helpString += "The file parameter allows you to provide a file containing a list of sff files to merge.  \n";
        helpString += "The keytrim parameter allows you to mergesff files with different keysequence by trimming them to the first 4 characters. Provided the first 4 match.  \n";
        helpString += "The output parameter allows you to provide an output filename.  \n";
		helpString += "Example sffinfo(sff=mySffFile.sff-mySecond.sff).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "getHelpString");
		exit(1);
	}
}
//*******************************************************************************
string MergeSfffilesCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "sff")            {   pattern =  "[filename],";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeSfffilesCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************
MergeSfffilesCommand::MergeSfffilesCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
 			
			ValidParameters validParameter;
			
            string inputDir = validParameter.validPath(parameters, "inputdir");
            if (inputDir == "not found"){    inputDir = "";        }
			
			sffFilename = validParameter.validPath(parameters, "sff");
			if (sffFilename == "not found") { sffFilename = "";  }
			else {
				util.splitAtDash(sffFilename, filenames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < filenames.size(); i++) {
					bool ignore = false;
					if (filenames[i] == "current") {
						filenames[i] = current->getSFFFile();
						if (filenames[i] != "") {  m->mothurOut("Using " + filenames[i] + " as input file for the sff parameter where you had given current.\n");  }
						else {
							m->mothurOut("You have no current sfffile, ignoring current.\n");  ignore=true;
							//erase from file list
							filenames.erase(filenames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						if (inputDir != "") {
							string path = util.hasPath(filenames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	filenames[i] = inputDir + filenames[i];		}
						}
                        
                        bool ableToOpen = util.checkLocations(filenames[i], current->getLocations());
                        
						if (!ableToOpen) {
							m->mothurOut("Unable to open " + filenames[i] + ". It will be disregarded.\n");
							filenames.erase(filenames.begin()+i); //erase from file list
							i--;
						}else { current->setSFFFile(filenames[i]); }
					}
				}
			}
			
			file = validParameter.validFile(parameters, "file");
			if (file == "not open") {  abort = true; }
			else if (file == "not found") { file = "";  }
            
            if ((file == "") && (filenames.size() == 0)) { m->mothurOut("[ERROR]: no valid files.\n");  abort = true; }
            
            if ((file != "") && (filenames.size() != 0)) { //both are given
                m->mothurOut("[ERROR]: cannot use file option and sff option at the same time, choose one.\n");  abort = true;
            }
            
            outputFile = validParameter.validPath(parameters, "output");
			if (outputFile == "not found") { m->mothurOut("you must enter an output file name\n");   abort=true;  }
			if (outputdir != "") { outputFile = outputdir + util.getSimpleName(outputFile);  }
            
            string temp = validParameter.valid(parameters, "keytrim");				if (temp == "not found") { temp = "F"; }
            keyTrim = util.isTrue(temp);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "MergeSfffilesCommand");
		exit(1);
	}
}
//*****************************************************************************
int MergeSfffilesCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (file != "") {
            readFile();
            if (outputdir == "") { outputdir = util.hasPath(file); }
        }
        ofstream out;
        map<string, string> variables;
        string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(outputFile);  }
        variables["[filename]"] = thisOutputDir + util.getSimpleName(outputFile);
		outputFile = getOutputFileName("sff",variables);
        util.openOutputFileBinary(outputFile, out);
        outputNames.push_back(outputFile); outputTypes["sff"].push_back(outputFile);
        outputFileHeader = outputFile + ".headers";
        numTotalReads = 0;
        
		for (int s = 0; s < filenames.size(); s++) {
			
			if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
			
			long start = time(nullptr);
			
            filenames[s] = util.getFullPathName(filenames[s]);
			m->mothurOut("\nMerging info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
            
			int numReads = mergeSffInfo(filenames[s], out);
            
			m->mothurOut("It took " + toString(time(nullptr) - start) + " secs to merge " + toString(numReads) + ".\n");
		}
        out.close();
        
        //create new common header and add to merged file
        adjustCommonHeader();

		if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
		
		//set sff file as new current sff file
		string currentName = "";
		itTypes = outputTypes.find("sff");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSFFFile(currentName); }
		}
		
		//report output filenames
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "execute");
		exit(1);
	}
}
//*****************************************************************************
int MergeSfffilesCommand::mergeSffInfo(string input, ofstream& out){
	try {
		currentFileName = input;
        
		ifstream in; util.openInputFileBinary(input, in);
		
		SffCommonHeader* header = new SffCommonHeader();
        bool goodHeader = header->read(in);
        
		if (!goodHeader) {  return 0; }
    
        commonHeaders.push_back(header); //save for adjustHeader sanity check
        
		//read through the sff file
        int count = 0; int numFlows = header->getNumFlows();
		while (!in.eof()) {
            
			//read data
			SffRead* read = new SffRead(numFlows);
            
            bool okay = read->readSff(in);
            
            if (!okay) { break; }
			
            read->printSff(out); numTotalReads++; count++;
            delete read;
            
			//report progress
			if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
            
			if (m->getControl_pressed()) { count = 0; break;   }
			
			if (count >= header->getNumReads()) { break; }
		}
        
		//report progress
		if (!m->getControl_pressed()) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "mergeSffInfo");
		exit(1);
	}
}
//****************************************************************************
void MergeSfffilesCommand::adjustCommonHeader(){
	try {
        //sanity check
        bool okayMagic = true;
        bool okayVersion = true;
        bool okayHeader = true;
        bool okayKeyLength = true;
        bool okayNumFlows = true;
        bool okayformatCode = true;
        bool okayflowChar = true;
        bool okayKeySequence = true;
        if (commonHeaders.size() != 0) {
            unsigned int magicN = commonHeaders[0]->getMagicNumber();
            string version = commonHeaders[0]->getVersion();
            unsigned short headerLength = commonHeaders[0]->getHeaderLength();
            unsigned short keyLength = commonHeaders[0]->getKeyLength();
            unsigned short numFlows = commonHeaders[0]->getNumFlows();
            int flowCode = commonHeaders[0]->getFlowgramFormat();
            string flowChars = commonHeaders[0]->getFlows();
            string keySeq = commonHeaders[0]->getKeySequence();
            
            for (int i = 1; i < commonHeaders.size(); i++) {
                if (commonHeaders[i]->getMagicNumber() != magicN)             { okayMagic = false;  m->mothurOut("[ERROR]: merge issue with common headers. Magic numbers do not match. " + filenames[0] + " magic number is " + toString(magicN) + ", but " + filenames[i] + " magic number is " + toString(commonHeaders[i]->getMagicNumber()) + ".\n");  }
                if (commonHeaders[i]->getVersion() != version)                { okayVersion = false;   m->mothurOut("[ERROR]: merge issue with common headers. Versions do not match. " + filenames[0] + " version is " + version + ", but " + filenames[i] + " version is " + commonHeaders[i]->getVersion() + ".\n");     }
                if (commonHeaders[i]->getHeaderLength() != headerLength)      { okayHeader = false;    m->mothurOut("[ERROR]: merge issue with common headers. Header lengths do not match. " + filenames[0] + " header length is " + toString(headerLength) + ", but " + filenames[i] + " header length is " + toString(commonHeaders[i]->getHeaderLength()) + ".\n");    }
                if (commonHeaders[i]->getKeyLength() != keyLength)            { okayKeyLength = false;  m->mothurOut("[ERROR]: merge issue with common headers. Key Lengths do not match. " + filenames[0] + " Key length is " + toString(keyLength) + ", but " + filenames[i] + " key length is " + toString(commonHeaders[i]->getKeyLength()) + ".\n");    }
                if (commonHeaders[i]->getNumFlows() != numFlows)       { okayNumFlows = false;   m->mothurOut("[ERROR]: merge issue with common headers. Number of flows per read do not match. " + filenames[0] + " number of flows is " + toString(numFlows) + ", but " + filenames[i] + " number of flows is " + toString(commonHeaders[i]->getNumFlows()) + ".\n");     }
                if (commonHeaders[i]->getFlowgramFormat() != flowCode)     { okayformatCode = false;    m->mothurOut("[ERROR]: merge issue with common headers. Flow format codes do not match. " + filenames[0] + " Flow format code is " + toString(flowCode) + ", but " + filenames[i] + " flow format code is " + toString(commonHeaders[i]->getFlowgramFormat()) + ".\n");    }
                if (commonHeaders[i]->getFlows() != flowChars)            { okayflowChar = false;   m->mothurOut("[ERROR]: merge issue with common headers. Flow characters do not match. " + filenames[0] + " Flow characters are " + flowChars + ", but " + filenames[i] + " flow characters are " + commonHeaders[i]->getFlows() + ".\n");    }
                if (commonHeaders[i]->getKeySequence() != keySeq)             { okayKeySequence = false;
                    if (keyTrim) {
                        m->mothurOut("[WARNING]: merge issue with common headers. Key sequences do not match. " + filenames[0] + " Key sequence is " + keySeq + ", but " + filenames[i] + " key sequence is " + commonHeaders[i]->getKeySequence() + ". We will attempt to trim them.\n");
                    }else { m->mothurOut("[ERROR]: merge issue with common headers. Key sequences do not match. " + filenames[0] + " Key sequence is " + keySeq + ", but " + filenames[i] + " key sequence is " + commonHeaders[i]->getKeySequence() + ".\n");
                    }
                }
            }
        }else { m->setControl_pressed(true); return; } //should never get here
        
        bool modify = false;
        if (!okayMagic || !okayVersion || !okayHeader || !okayKeyLength || !okayNumFlows || !okayformatCode || !okayflowChar) { m->setControl_pressed(true); return; }
        if (!okayKeySequence) {
            bool okayKeySequence2 = true;
            string keySeq = commonHeaders[0]->getKeySequence().substr(0,4);
            for (int i = 1; i < commonHeaders.size(); i++) {
                if ((commonHeaders[i]->getKeySequence().substr(0,4)) != keySeq)          { okayKeySequence2 = false;   }
            }
            if (okayKeySequence2 && keyTrim) {  modify = true;
                m->mothurOut("We are able to trim the key sequences. Merged key seqeunce will be " + keySeq + ".\n");
            }
        }
        
        ofstream out;
        util.openOutputFileBinaryAppend(outputFileHeader, out);
        commonHeaders[0]->printSampleCommonHeader(out, numTotalReads);
        out.close();
        
        util.appendSFFFiles(outputFile, outputFileHeader);
        util.renameFile(outputFileHeader, outputFile);
        util.mothurRemove(outputFileHeader);
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "adjustCommonHeader");
		exit(1);
	}
}
//*************************************************************************************
void MergeSfffilesCommand::readFile(){
	try {
        ifstream in; util.openInputFile(file, in);
        
        string filename;
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { return; }
            
            in >> filename; gobble(in);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: filename = " + filename + ".\n"); }
            
            bool ableToOpen = util.checkLocations(filename, current->getLocations());
            
            if (!ableToOpen) { //can't find it
                m->mothurOut("[WARNING]: can't find " + filename + ", ignoring.\n");
            }else{  filenames.push_back(filename); }
            
        }
        in.close();
    }
    catch(exception& e) {
        m->errorOut(e, "MergeSfffilesCommand", "readFile");
        exit(1);
    }
}
//******************************************************************************************
