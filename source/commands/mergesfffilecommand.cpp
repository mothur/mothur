//
//  mergesfffilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/31/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "mergesfffilecommand.h"
#include "endiannessmacros.h"

//**********************************************************************************************************************
vector<string> MergeSfffilesCommand::setParameters(){
	try {
		CommandParameter psff("sff", "InputTypes", "", "", "sffFile", "sffFile", "none","sff",false,false); parameters.push_back(psff);
        CommandParameter pfile("file", "InputTypes", "", "", "sffFile", "sffFile", "none","sff",false,false); parameters.push_back(pfile);
        CommandParameter pkeytrim("keytrim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeytrim);
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,true,true); parameters.push_back(poutput);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
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

//**********************************************************************************************************************
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
//**********************************************************************************************************************
MergeSfffilesCommand::MergeSfffilesCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
        outputTypes["sff"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "MergeSfffilesCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MergeSfffilesCommand::MergeSfffilesCommand(string option)  {
	try {
		abort = false; calledHelp = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
            map<string,string>::iterator it;
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
            outputTypes["sff"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.valid(parameters, "inputdir");	  if (inputDir == "not found"){	inputDir = "";		}
            else {
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					string path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
            }
            
			sffFilename = validParameter.valid(parameters, "sff");
			if (sffFilename == "not found") { sffFilename = "";  }
			else {
				util.splitAtDash(sffFilename, filenames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < filenames.size(); i++) {
					bool ignore = false;
					if (filenames[i] == "current") {
						filenames[i] = current->getSFFFile();
						if (filenames[i] != "") {  m->mothurOut("Using " + filenames[i] + " as input file for the sff parameter where you had given current."); m->mothurOutEndLine(); }
						else {
							m->mothurOut("You have no current sfffile, ignoring current."); m->mothurOutEndLine(); ignore=true;
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
                        
						ifstream in;
						bool ableToOpen = util.openInputFile(filenames[i], in, "noerror");
                        
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (current->getDefaultPath() != "") { //default path is set
								string tryPath = current->getDefaultPath() + util.getSimpleName(filenames[i]);
								m->mothurOut("Unable to open " + filenames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = util.openInputFile(tryPath, in2, "noerror");
								in2.close();
								filenames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (current->getOutputDir() != "") { //default path is set
								string tryPath = current->getOutputDir() + util.getSimpleName(filenames[i]);
								m->mothurOut("Unable to open " + filenames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = util.openInputFile(tryPath, in2, "noerror");
								in2.close();
								filenames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (!ableToOpen) {
							m->mothurOut("Unable to open " + filenames[i] + ". It will be disregarded."); m->mothurOutEndLine();
							//erase from file list
							filenames.erase(filenames.begin()+i);
							i--;
						}else { current->setSFFFile(filenames[i]); }
					}
				}
			}
			
			file = validParameter.validFile(parameters, "file");
			if (file == "not open") {  abort = true; }
			else if (file == "not found") { file = "";  }
            
            if ((file == "") && (filenames.size() == 0)) {
                m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((file != "") && (filenames.size() != 0)) { //both are given
                m->mothurOut("[ERROR]: cannot use file option and sff option at the same time, choose one."); m->mothurOutEndLine(); abort = true;
            }
            
            outputFile = validParameter.valid(parameters, "output");
			if (outputFile == "not found") { m->mothurOut("you must enter an output file name"); m->mothurOutEndLine();  abort=true;  }
			if (outputDir != "") { outputFile = outputDir + util.getSimpleName(outputFile);  }
            
            string temp = validParameter.valid(parameters, "keytrim");				if (temp == "not found") { temp = "F"; }
            keyTrim = util.isTrue(temp);
            
		}
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "MergeSfffilesCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int MergeSfffilesCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (file != "") {
            readFile();
            if (outputDir == "") { outputDir = util.hasPath(file); }
        }
        ofstream out;
        map<string, string> variables;
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(outputFile);  }
        variables["[filename]"] = thisOutputDir + util.getSimpleName(outputFile);
		outputFile = getOutputFileName("sff",variables);
        util.openOutputFileBinary(outputFile, out);
        outputNames.push_back(outputFile); outputTypes["sff"].push_back(outputFile);
        outputFileHeader = outputFile + ".headers";
        numTotalReads = 0;
        
		for (int s = 0; s < filenames.size(); s++) {
			
			if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
			
			long start = time(NULL);
			
            filenames[s] = util.getFullPathName(filenames[s]);
			m->mothurOut("\nMerging info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
            
			int numReads = mergeSffInfo(filenames[s], out);
            
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to merge " + toString(numReads) + ".\n");
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
//**********************************************************************************************************************
int MergeSfffilesCommand::mergeSffInfo(string input, ofstream& out){
	try {
		currentFileName = input;
        
		ifstream in;
		util.openInputFileBinary(input, in);
		
		CommonHeader header;
		readCommonHeader(in, header);
        
		int count = 0;
		
		//check magic number and version
		if (header.magicNumber != 779314790) { m->mothurOut("Magic Number is not correct, not a valid .sff file"); m->mothurOutEndLine(); return count; }
		if (header.version != "0001") { m->mothurOut("Version is not supported, only support version 0001."); m->mothurOutEndLine(); return count; }
		
        //save for adjustHeader sanity check
        commonHeaders.push_back(header);
        
		//read through the sff file
		while (!in.eof()) {
            
			//read data
			seqRead read;  Header readheader;
            readSeqData(in, read, header.numFlowsPerRead, readheader, out);
            
            bool okay = sanityCheck(readheader, read);
            if (!okay) { break; }
			
			count++;
            
			//report progress
			if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
            
			if (m->getControl_pressed()) { count = 0; break;   }
			
			if (count >= header.numReads) { break; }
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
//**********************************************************************************************************************
int MergeSfffilesCommand::readCommonHeader(ifstream& in, CommonHeader& header){
	try {
        
        if (!in.eof()) {
            
            //read magic number
            char buffer[4];
            in.read(buffer, 4);
            header.magicNumber = be_int4(*(unsigned int *)(&buffer));
            
            //read version
            char buffer9[4];
            in.read(buffer9, 4);
            header.version = "";
            for (int i = 0; i < 4; i++) {  header.version += toString((int)(buffer9[i]));  }
            
            //read offset
            char buffer2 [8];
            in.read(buffer2, 8);
            header.indexOffset =  be_int8(*(unsigned long long *)(&buffer2));
            
            
            //read index length
            char buffer3 [4];
            in.read(buffer3, 4);
            header.indexLength =  be_int4(*(unsigned int *)(&buffer3));
            
            //read num reads
            char buffer4 [4];
            in.read(buffer4, 4);
            header.numReads =  be_int4(*(unsigned int *)(&buffer4));
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: numReads = " + toString(header.numReads) + "\n"); }
            
            //read header length
            char buffer5 [2];
            in.read(buffer5, 2);
            header.headerLength =  be_int2(*(unsigned short *)(&buffer5));
            
            //read key length
            char buffer6 [2];
            in.read(buffer6, 2);
            header.keyLength = be_int2(*(unsigned short *)(&buffer6));
            //cout << "header key length = " << header.keyLength << endl;
            
            //read number of flow reads
            char buffer7 [2];
            in.read(buffer7, 2);
            header.numFlowsPerRead =  be_int2(*(unsigned short *)(&buffer7));
            
            //read format code
            char buffer8 [1];
            in.read(buffer8, 1);
            header.flogramFormatCode = (int)(buffer8[0]);
            
            //read flow chars
            char* tempBuffer = new char[header.numFlowsPerRead];
            in.read(&(*tempBuffer), header.numFlowsPerRead);
            header.flowChars = tempBuffer;
            if (header.flowChars.length() > header.numFlowsPerRead) { header.flowChars = header.flowChars.substr(0, header.numFlowsPerRead);  }
            delete[] tempBuffer;
            
            //read key
            char* tempBuffer2 = new char[header.keyLength];
            in.read(&(*tempBuffer2), header.keyLength);
            header.keySequence = tempBuffer2;
            //cout << "key sequence " <<header.keySequence << endl;
            if (header.keySequence.length() > header.keyLength) { header.keySequence = header.keySequence.substr(0, header.keyLength);  }
            delete[] tempBuffer2;
            //cout << "key sequence " <<header.keySequence << endl;
            
            /* Pad to 8 chars */
            unsigned long long spotInFile = in.tellg();
            unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
            in.seekg(spot);
            
        }else{
            m->mothurOut("Error reading sff common header."); m->mothurOutEndLine();
        }
        
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "readCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int MergeSfffilesCommand::adjustCommonHeader(){
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
            unsigned int magicN = commonHeaders[0].magicNumber;
            string version = commonHeaders[0].version;
            unsigned short headerLength = commonHeaders[0].headerLength;
            unsigned short keyLength = commonHeaders[0].keyLength;
            unsigned short numFlows = commonHeaders[0].numFlowsPerRead;
            int flowCode = commonHeaders[0].flogramFormatCode;
            string flowChars = commonHeaders[0].flowChars;
            string keySeq = commonHeaders[0].keySequence;
            
            for (int i = 1; i < commonHeaders.size(); i++) {
                if (commonHeaders[i].magicNumber != magicN)             { okayMagic = false;  m->mothurOut("[ERROR]: merge issue with common headers. Magic numbers do not match. " + filenames[0] + " magic number is " + toString(commonHeaders[0].magicNumber) + ", but " + filenames[i] + " magic number is " + toString(commonHeaders[i].magicNumber) + ".\n");  }
                if (commonHeaders[i].version != version)                { okayVersion = false;   m->mothurOut("[ERROR]: merge issue with common headers. Versions do not match. " + filenames[0] + " version is " + commonHeaders[0].version + ", but " + filenames[i] + " version is " + commonHeaders[i].version + ".\n");     }
                if (commonHeaders[i].headerLength != headerLength)      { okayHeader = false;    m->mothurOut("[ERROR]: merge issue with common headers. Header lengths do not match. " + filenames[0] + " header length is " + toString(commonHeaders[0].headerLength) + ", but " + filenames[i] + " header length is " + toString(commonHeaders[i].headerLength) + ".\n");    }
                if (commonHeaders[i].keyLength != keyLength)            { okayKeyLength = false;  m->mothurOut("[ERROR]: merge issue with common headers. Key Lengths do not match. " + filenames[0] + " Key length is " + toString(commonHeaders[0].keyLength) + ", but " + filenames[i] + " key length is " + toString(commonHeaders[i].keyLength) + ".\n");    }
                if (commonHeaders[i].numFlowsPerRead != numFlows)       { okayNumFlows = false;   m->mothurOut("[ERROR]: merge issue with common headers. Number of flows per read do not match. " + filenames[0] + " number of flows is " + toString(commonHeaders[0].numFlowsPerRead) + ", but " + filenames[i] + " number of flows is " + toString(commonHeaders[i].numFlowsPerRead) + ".\n");     }
                if (commonHeaders[i].flogramFormatCode != flowCode)     { okayformatCode = false;    m->mothurOut("[ERROR]: merge issue with common headers. Flow format codes do not match. " + filenames[0] + " Flow format code is " + toString(commonHeaders[0].flogramFormatCode) + ", but " + filenames[i] + " flow format code is " + toString(commonHeaders[i].flogramFormatCode) + ".\n");    }
                if (commonHeaders[i].flowChars != flowChars)            { okayflowChar = false;   m->mothurOut("[ERROR]: merge issue with common headers. Flow characters do not match. " + filenames[0] + " Flow characters are " + commonHeaders[0].flowChars + ", but " + filenames[i] + " flow characters are " + commonHeaders[i].flowChars + ".\n");    }
                if (commonHeaders[i].keySequence != keySeq)             { okayKeySequence = false;
                    if (keyTrim) {
                        m->mothurOut("[WARNING]: merge issue with common headers. Key sequences do not match. " + filenames[0] + " Key sequence is " + commonHeaders[0].keySequence + ", but " + filenames[i] + " key sequence is " + commonHeaders[i].keySequence + ". We will attempt to trim them.\n");
                    }else { m->mothurOut("[ERROR]: merge issue with common headers. Key sequences do not match. " + filenames[0] + " Key sequence is " + commonHeaders[0].keySequence + ", but " + filenames[i] + " key sequence is " + commonHeaders[i].keySequence + ".\n");
                    }
                }
            }
        }else { m->setControl_pressed(true); return 0; } //should never get here
        
        bool modify = false;
        if (!okayMagic || !okayVersion || !okayHeader || !okayKeyLength || !okayNumFlows || !okayformatCode || !okayflowChar) { m->setControl_pressed(true); return 0; }
        if (!okayKeySequence) {
            bool okayKeySequence2 = true;
            string keySeq = commonHeaders[0].keySequence.substr(0,4);
            for (int i = 1; i < commonHeaders.size(); i++) {
                if ((commonHeaders[i].keySequence.substr(0,4)) != keySeq)          { okayKeySequence2 = false;   }
            }
            if (okayKeySequence2 && keyTrim) {  modify = true;
                m->mothurOut("We are able to trim the key sequences. Merged key seqeunce will be " + keySeq + ".\n");
            }
        }
        
        string endian = util.findEdianness();
        
        char* mybuffer = new char[4];
        ifstream in;
        util.openInputFileBinary(currentFileName, in);
        
        //magic number
        in.read(mybuffer,4);
        ofstream out;
        util.openOutputFileBinaryAppend(outputFileHeader, out);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
        
        //version
        mybuffer = new char[4];
        in.read(mybuffer,4);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
        
        //offset
        mybuffer = new char[8];
        in.read(mybuffer,8);
        unsigned long long offset = 0;
        char* thisbuffer = new char[8];
        thisbuffer[0] = (offset >> 56) & 0xFF;
        thisbuffer[1] = (offset >> 48) & 0xFF;
        thisbuffer[2] = (offset >> 40) & 0xFF;
        thisbuffer[3] = (offset >> 32) & 0xFF;
        thisbuffer[4] = (offset >> 24) & 0xFF;
        thisbuffer[5] = (offset >> 16) & 0xFF;
        thisbuffer[6] = (offset >> 8) & 0xFF;
        thisbuffer[7] = offset & 0xFF;
        out.write(thisbuffer, 8);
        delete[] thisbuffer;
        delete[] mybuffer;
        
        //read index length
		mybuffer = new char[4];
        in.read(mybuffer,4);
        offset = 0;
        char* thisbuffer2 = new char[4];
        thisbuffer2[0] = (offset >> 24) & 0xFF;
        thisbuffer2[1] = (offset >> 16) & 0xFF;
        thisbuffer2[2] = (offset >> 8) & 0xFF;
        thisbuffer2[3] = offset & 0xFF;
        out.write(thisbuffer2, 4);
        delete[] thisbuffer2;
        delete[] mybuffer;
		
        //change num reads
        mybuffer = new char[4];
        in.read(mybuffer,4);
        delete[] mybuffer;
        thisbuffer2 = new char[4];
        thisbuffer2[0] = (numTotalReads >> 24) & 0xFF;
        thisbuffer2[1] = (numTotalReads >> 16) & 0xFF;
        thisbuffer2[2] = (numTotalReads >> 8) & 0xFF;
        thisbuffer2[3] = numTotalReads & 0xFF;
        
        out.write(thisbuffer2, 4);
        delete[] thisbuffer2;
        
        //read header length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
       
        
        //read key length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        if (modify) {
            unsigned short fourL = 4;
            thisbuffer2 = new char[2];
            thisbuffer2[0] = (fourL >> 8) & 0xFF;
            thisbuffer2[1] = fourL & 0xFF;
            out.write(thisbuffer2, in.gcount());
            delete[] thisbuffer2;
        }else {
            out.write(mybuffer, in.gcount());
        }
        delete[] mybuffer;
        
        //read number of flow reads
        mybuffer = new char[2];
        in.read(mybuffer,2);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
        
        //read format code
        mybuffer = new char[1];
        in.read(mybuffer,1);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
        
        //read flow chars
        mybuffer = new char[commonHeaders[0].numFlowsPerRead];
        in.read(mybuffer,commonHeaders[0].numFlowsPerRead);
        out.write(mybuffer, in.gcount());
        delete[] mybuffer;
        
        //read key
        mybuffer = new char[commonHeaders[0].keyLength];
        in.read(mybuffer,commonHeaders[0].keyLength);
        if (modify) {
            out.write(mybuffer, 4);
        }else {
            out.write(mybuffer, in.gcount());
        }
        delete[] mybuffer;
        
        /* Pad to 8 chars */
        unsigned long long spotInFile = in.tellg();
        if (modify) {  spotInFile -= commonHeaders[0].keyLength - 4;  }
        unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
        in.seekg(spot);
        
        mybuffer = new char[spot-spotInFile];
        out.write(mybuffer, spot-spotInFile);
        delete[] mybuffer;
        in.close();
        out.close();
        
        util.appendSFFFiles(outputFile, outputFileHeader);
        util.renameFile(outputFileHeader, outputFile);
        util.mothurRemove(outputFileHeader);
        
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "adjustCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
bool MergeSfffilesCommand::readSeqData(ifstream& in, seqRead& read, int numFlowReads, Header& header, ofstream& out){
	try {
        unsigned long long startSpotInFile = in.tellg();
		if (!in.eof()) {
            
            /*****************************************/
            //read header
            
            //read header length
			char buffer [2];
			in.read(buffer, 2);
			header.headerLength = be_int2(*(unsigned short *)(&buffer));
            
			//read name length
			char buffer2 [2];
			in.read(buffer2, 2);
			header.nameLength = be_int2(*(unsigned short *)(&buffer2));
            
			//read num bases
			char buffer3 [4];
			in.read(buffer3, 4);
			header.numBases =  be_int4(*(unsigned int *)(&buffer3));
            
			
			//read clip qual left
			char buffer4 [2];
			in.read(buffer4, 2);
			header.clipQualLeft =  be_int2(*(unsigned short *)(&buffer4));
			header.clipQualLeft = 5;
            
			
			//read clip qual right
			char buffer5 [2];
			in.read(buffer5, 2);
			header.clipQualRight =  be_int2(*(unsigned short *)(&buffer5));
            
            
			//read clipAdapterLeft
			char buffer6 [2];
			in.read(buffer6, 2);
			header.clipAdapterLeft = be_int2(*(unsigned short *)(&buffer6));
            
            
			//read clipAdapterRight
			char buffer7 [2];
			in.read(buffer7, 2);
			header.clipAdapterRight = be_int2(*(unsigned short *)(&buffer7));
            
            
			//read name
			char* tempBuffer = new char[header.nameLength];
			in.read(&(*tempBuffer), header.nameLength);
			header.name = tempBuffer;
			if (header.name.length() > header.nameLength) { header.name = header.name.substr(0, header.nameLength);  }
            
			delete[] tempBuffer;
						
			/* Pad to 8 chars */
			unsigned long long spotInFile = in.tellg();
			unsigned long long spot = (spotInFile + 7)& ~7;
			in.seekg(spot);
            
            /*****************************************/
            //sequence read
            
			//read flowgram
			read.flowgram.resize(numFlowReads);
			for (int i = 0; i < numFlowReads; i++) {
				char buffer [2];
				in.read(buffer, 2);
				read.flowgram[i] = be_int2(*(unsigned short *)(&buffer));
			}
            
			//read flowIndex
			read.flowIndex.resize(header.numBases);
			for (int i = 0; i < header.numBases; i++) {
				char temp[1];
				in.read(temp, 1);
				read.flowIndex[i] = be_int1(*(unsigned char *)(&temp));
			}
            
			//read bases
			char* tempBuffer6 = new char[header.numBases];
			in.read(&(*tempBuffer6), header.numBases);
			read.bases = tempBuffer6;
			if (read.bases.length() > header.numBases) { read.bases = read.bases.substr(0, header.numBases);  }
			delete[] tempBuffer6;
            
			//read qual scores
			read.qualScores.resize(header.numBases);
			for (int i = 0; i < header.numBases; i++) {
				char temp[1];
				in.read(temp, 1);
				read.qualScores[i] = be_int1(*(unsigned char *)(&temp));
			}
            
			/* Pad to 8 chars */
			spotInFile = in.tellg();
			spot = (spotInFile + 7)& ~7;
			in.seekg(spot);
            
            char * mybuffer;
            mybuffer = new char [spot-startSpotInFile];
            
            ifstream in2;
            util.openInputFileBinary(currentFileName, in2);
            in2.seekg(startSpotInFile);
            in2.read(mybuffer,spot-startSpotInFile);
            
            out.write(mybuffer, in2.gcount());
            numTotalReads++;
            
            delete[] mybuffer;
            in2.close();
            
		}else{
			m->mothurOut("Error reading."); m->mothurOutEndLine();
		}
        
        if (in.eof()) {  return true; }
        
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "readSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
bool MergeSfffilesCommand::sanityCheck(Header& header, seqRead& read) {
	try {
        bool okay = true;
        string message = "[WARNING]: Your sff file may be corrupted! Sequence: " + header.name + "\n";
        
        if (header.clipQualLeft > read.bases.length()) {
            okay = false; message += "Clip Qual Left = " + toString(header.clipQualLeft) + ", but we only read " + toString(read.bases.length()) + " bases.\n";
        }
        if (header.clipQualRight > read.bases.length()) {
            okay = false; message += "Clip Qual Right = " + toString(header.clipQualRight) + ", but we only read " + toString(read.bases.length()) + " bases.\n";
        }
        if (header.clipQualLeft > read.qualScores.size()) {
            okay = false; message += "Clip Qual Left = " + toString(header.clipQualLeft) + ", but we only read " + toString(read.qualScores.size()) + " quality scores.\n";
        }
        if (header.clipQualRight > read.qualScores.size()) {
            okay = false; message += "Clip Qual Right = " + toString(header.clipQualRight) + ", but we only read " + toString(read.qualScores.size()) + " quality scores.\n";
        }
        
        if (okay == false) {
            m->mothurOut(message); m->mothurOutEndLine();
        }
        
		return okay;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeSfffilesCommand", "sanityCheck");
		exit(1);
	}
}
//**********************************************************************************************************************
int MergeSfffilesCommand::readFile(){
	try {
        
        string filename;
        
        ifstream in;
        util.openInputFile(file, in);
        
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { return 0; }
            
            in >> filename; util.gobble(in);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: filename = " + filename + ".\n"); }
            
            //check to make sure both are able to be opened
            ifstream in2;
            bool openForward = util.openInputFile(filename, in2, "noerror");
            
            //if you can't open it, try default location
            if (!openForward) {
                if (current->getDefaultPath() != "") { //default path is set
                    string tryPath = current->getDefaultPath() + util.getSimpleName(filename);
                    m->mothurOut("Unable to open " + filename + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openForward = util.openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    filename = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (!openForward) {
                if (current->getOutputDir() != "") { //default path is set
                    string tryPath = current->getOutputDir() + util.getSimpleName(filename);
                    m->mothurOut("Unable to open " + filename + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = util.openInputFile(tryPath, in4, "noerror");
                    filename = tryPath;
                    in4.close();
                }
            }
            
            if (!openForward) { //can't find it
                m->mothurOut("[WARNING]: can't find " + filename + ", ignoring.\n");
            }else{  filenames.push_back(filename); }
            
        }
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeSfffilesCommand", "readFileNames");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeSfffilesCommand::printCommonHeaderForDebug(CommonHeader& header, ofstream& out, int numReads){
    try {
        string endian = util.findEdianness();
        
        ifstream in;
        util.openInputFileBinary(currentFileName, in);
        
        //magic number
        char* mybuffer = new char[4];
        in.read(mybuffer,4);
        out.write(mybuffer, in.gcount());
        string contents = mybuffer;
        m->mothurOut("magicNumber = " + contents + "\n");
        delete[] mybuffer;
        
        //version
        char* mybuffer1 = new char[4];
        in.read(mybuffer1,4);
        out.write(mybuffer1, in.gcount());
        contents = mybuffer1;
        m->mothurOut("version = " + contents + "\n");
        m->mothurOut("version = " + header.version + "\n");
        delete[] mybuffer1;
        
        //offset
        char* mybuffer2 = new char[8];
        in.read(mybuffer2,8);
        unsigned long long offset = 0;
        char* thisbuffer = new char[8];
        thisbuffer[0] = (offset >> 56) & 0xFF;
        thisbuffer[1] = (offset >> 48) & 0xFF;
        thisbuffer[2] = (offset >> 40) & 0xFF;
        thisbuffer[3] = (offset >> 32) & 0xFF;
        thisbuffer[4] = (offset >> 24) & 0xFF;
        thisbuffer[5] = (offset >> 16) & 0xFF;
        thisbuffer[6] = (offset >> 8) & 0xFF;
        thisbuffer[7] = offset & 0xFF;
        out.write(thisbuffer, 8);
        delete[] thisbuffer;
        delete[] mybuffer2;
        m->mothurOut("index offset = " + toString(header.indexOffset) + "\n");
        
        //read index length
        char* mybuffer3 = new char[4];
        in.read(mybuffer3,4);
        offset = 0;
        char* thisbuffer2 = new char[4];
        thisbuffer2[0] = (offset >> 24) & 0xFF;
        thisbuffer2[1] = (offset >> 16) & 0xFF;
        thisbuffer2[2] = (offset >> 8) & 0xFF;
        thisbuffer2[3] = offset & 0xFF;
        out.write(thisbuffer2, 4);
        delete[] thisbuffer2;
        delete[] mybuffer3;
        m->mothurOut("index read length = " + toString(header.indexLength) + "\n");
        
        //change num reads
        char* mybuffer4 = new char[4];
        in.read(mybuffer4,4);
        
        char* thisbuffer3 = new char[4];
        if (endian == "BIG_ENDIAN") {
            thisbuffer3[0] = (numReads >> 24) & 0xFF;
            thisbuffer3[1] = (numReads >> 16) & 0xFF;
            thisbuffer3[2] = (numReads >> 8) & 0xFF;
            thisbuffer3[3] = numReads & 0xFF;
        }else {
            thisbuffer3[0] = numReads & 0xFF;
            thisbuffer3[1] = (numReads >> 8) & 0xFF;
            thisbuffer3[2] = (numReads >> 16) & 0xFF;
            thisbuffer3[3] = (numReads >> 24) & 0xFF;
        }
        out.write(thisbuffer3, 4);
        contents = mybuffer4;
        m->mothurOut("numReads = " + contents + "\n");
        unsigned int numTReads = be_int4(*(unsigned int *)(mybuffer4));
        m->mothurOut("numReads = " + toString(numTReads) + "\n");
        m->mothurOut("numReads = " + toString(header.numReads) + "\n");
        delete[] thisbuffer3;
        delete[] mybuffer4;
        
        //read header length
        char* mybuffer5 = new char[2];
        in.read(mybuffer5,2);
        out.write(mybuffer5, in.gcount());
        contents = mybuffer5;
        m->mothurOut("readLength = " + contents + "\n");
        m->mothurOut("readLength = " + toString(header.headerLength) + "\n");
        delete[] mybuffer5;
        
        //read key length
        char* mybuffer6 = new char[2];
        in.read(mybuffer6,2);
        out.write(mybuffer6, in.gcount());
        contents = mybuffer6;
        m->mothurOut("key length = " + contents + "\n");
        m->mothurOut("key length = " + toString(header.keyLength) + "\n");
        delete[] mybuffer6;
        
        //read number of flow reads
        char* mybuffer7 = new char[2];
        in.read(mybuffer7,2);
        out.write(mybuffer7, in.gcount());
        contents = mybuffer7;
        m->mothurOut("num flow reads = " + contents + "\n");
        int numFlowReads = be_int2(*(unsigned short *)(mybuffer7));
        m->mothurOut("num flow Reads = " + toString(numFlowReads) + "\n");
        delete[] mybuffer7;
        
        //read format code
        char* mybuffer8 = new char[1];
        in.read(mybuffer8,1);
        out.write(mybuffer8, in.gcount());
        contents = mybuffer8;
        m->mothurOut("read format code = " + contents + "\n");
        m->mothurOut("read format code = " + toString(header.flogramFormatCode) + "\n");
        delete[] mybuffer8;
        
        //read flow chars
        char* mybuffer9 = new char[header.numFlowsPerRead];
        in.read(mybuffer9,header.numFlowsPerRead);
        out.write(mybuffer9, in.gcount());
        contents = mybuffer9;
        m->mothurOut("flow chars = " + contents + "\n");
        m->mothurOut("flow chars = " + header.flowChars + "\n");
        delete[] mybuffer9;
        
        //read key
        char* mybuffer10 = new char[header.keyLength];
        in.read(mybuffer10,header.keyLength);
        out.write(mybuffer10, in.gcount());
        contents = mybuffer10;
        m->mothurOut("key = " + contents + "\n");
        m->mothurOut("key = " + header.keySequence + "\n");
        delete[] mybuffer10;
        
        
        /* Pad to 8 chars */
        unsigned long long spotInFile = in.tellg();
        unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
        in.seekg(spot);
        
        char* mybuffer11 = new char[spot-spotInFile];
        out.write(mybuffer11, spot-spotInFile);
        delete[] mybuffer11;
        in.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeSfffilesCommand", "printCommonHeaderForDebug");
        exit(1);
    }
}
//**********************************************************************************************************************




