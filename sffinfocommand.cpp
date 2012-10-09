/*
 *  sffinfocommand.cpp
 *  Mothur
 *
 *  Created by westcott on 7/7/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "sffinfocommand.h"
#include "endiannessmacros.h"
#include "trimoligos.h"
#include "sequence.hpp"
#include "qualityscores.h"

//**********************************************************************************************************************
vector<string> SffInfoCommand::setParameters(){	
	try {		
		CommandParameter psff("sff", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(psff);
        CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(poligos);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(paccnos);
		CommandParameter psfftxt("sfftxt", "String", "", "", "", "", "",false,false); parameters.push_back(psfftxt);
		CommandParameter pflow("flow", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pflow);
		CommandParameter ptrim("trim", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(ptrim);
		CommandParameter pfasta("fasta", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pfasta);
		CommandParameter pqfile("name", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pqfile);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ptdiffs);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffInfoCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sffinfo command reads a sff file and extracts the sequence data, or you can use it to parse a sfftxt file.\n";
		helpString += "The sffinfo command parameters are sff, fasta, qfile, accnos, flow, sfftxt, oligos, bdiffs, tdiffs, ldiffs, sdiffs, pdiffs and trim. sff is required. \n";
		helpString += "The sff parameter allows you to enter the sff file you would like to extract data from.  You may enter multiple files by separating them by -'s.\n";
		helpString += "The fasta parameter allows you to indicate if you would like a fasta formatted file generated.  Default=True. \n";
		helpString += "The qfile parameter allows you to indicate if you would like a quality file generated.  Default=True. \n";
        helpString += "The oligos parameter allows you to provide an oligos file to split your sff file into separate sff files by barcode. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The flow parameter allows you to indicate if you would like a flowgram file generated.  Default=True. \n";
		helpString += "The sfftxt parameter allows you to indicate if you would like a sff.txt file generated.  Default=False. \n";
		helpString += "If you want to parse an existing sfftxt file into flow, fasta and quality file, enter the file name using the sfftxt parameter. \n";
		helpString += "The trim parameter allows you to indicate if you would like a sequences and quality scores trimmed to the clipQualLeft and clipQualRight values.  Default=True. \n";
		helpString += "The accnos parameter allows you to provide a accnos file containing the names of the sequences you would like extracted. You may enter multiple files by separating them by -'s. \n";
		helpString += "Example sffinfo(sff=mySffFile.sff, trim=F).\n";
		helpString += "Note: No spaces between parameter labels (i.e. sff), '=' and parameters (i.e.yourSffFileName).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffInfoCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "fasta")            {   outputFileName =  "fasta";   }
            else if (type == "flow")    {   outputFileName =  "flow";   }
            else if (type == "sfftxt")        {   outputFileName =  "sff.txt";   }
            else if (type == "sff")        {   outputFileName =  "sff";   }
            else if (type == "qfile")       {   outputFileName =  "qual";   }
             else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "getOutputFileNameTag");
		exit(1);
	}
}


//**********************************************************************************************************************
SffInfoCommand::SffInfoCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["flow"] = tempOutNames;
		outputTypes["sfftxt"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["sff"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "SffInfoCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

SffInfoCommand::SffInfoCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		hasAccnos = false; hasOligos = false;
        split = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["flow"] = tempOutNames;
			outputTypes["sfftxt"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["sff"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);	  if (inputDir == "not found"){	inputDir = "";		}

			sffFilename = validParameter.validFile(parameters, "sff", false);
			if (sffFilename == "not found") { sffFilename = "";  }
			else { 
				m->splitAtDash(sffFilename, filenames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < filenames.size(); i++) {
					bool ignore = false;
					if (filenames[i] == "current") { 
						filenames[i] = m->getSFFFile(); 
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
							string path = m->hasPath(filenames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	filenames[i] = inputDir + filenames[i];		}
						}
		
						ifstream in;
						int ableToOpen = m->openInputFile(filenames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(filenames[i]);
								m->mothurOut("Unable to open " + filenames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								filenames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(filenames[i]);
								m->mothurOut("Unable to open " + filenames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								filenames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + filenames[i] + ". It will be disregarded."); m->mothurOutEndLine();
							//erase from file list
							filenames.erase(filenames.begin()+i);
							i--;
						}else { m->setSFFFile(filenames[i]); }
					}
				}
				
				//make sure there is at least one valid file left
				if (filenames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			accnosName = validParameter.validFile(parameters, "accnos", false);
			if (accnosName == "not found") { accnosName = "";  }
			else { 
				hasAccnos = true;
				m->splitAtDash(accnosName, accnosFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < accnosFileNames.size(); i++) {
					bool ignore = false;
					if (accnosFileNames[i] == "current") { 
						accnosFileNames[i] = m->getAccnosFile(); 
						if (accnosFileNames[i] != "") {  m->mothurOut("Using " + accnosFileNames[i] + " as input file for the accnos parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current accnosfile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							accnosFileNames.erase(accnosFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
					
						if (inputDir != "") {
							string path = m->hasPath(accnosFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	accnosFileNames[i] = inputDir + accnosFileNames[i];		}
						}
		
						ifstream in;
						int ableToOpen = m->openInputFile(accnosFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(accnosFileNames[i]);
								m->mothurOut("Unable to open " + accnosFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								accnosFileNames[i] = tryPath;
							}
						}
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(accnosFileNames[i]);
								m->mothurOut("Unable to open " + accnosFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								accnosFileNames[i] = tryPath;
							}
						}
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + accnosFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();
							//erase from file list
							accnosFileNames.erase(accnosFileNames.begin()+i);
							i--;
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (accnosFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
            
            oligosfile = validParameter.validFile(parameters, "oligos", false);
			if (oligosfile == "not found") { oligosfile = "";  }
			else { 
				hasOligos = true;
				m->splitAtDash(oligosfile, oligosFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < oligosFileNames.size(); i++) {
					bool ignore = false;
					if (oligosFileNames[i] == "current") { 
						oligosFileNames[i] = m->getOligosFile(); 
						if (oligosFileNames[i] != "") {  m->mothurOut("Using " + oligosFileNames[i] + " as input file for the accnos parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current oligosfile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							oligosFileNames.erase(oligosFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
                        
						if (inputDir != "") {
							string path = m->hasPath(oligosFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	oligosFileNames[i] = inputDir + oligosFileNames[i];		}
						}
                        
						ifstream in;
						int ableToOpen = m->openInputFile(oligosFileNames[i], in, "noerror");
                        
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(oligosFileNames[i]);
								m->mothurOut("Unable to open " + oligosFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								oligosFileNames[i] = tryPath;
							}
						}
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(oligosFileNames[i]);
								m->mothurOut("Unable to open " + oligosFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								oligosFileNames[i] = tryPath;
							}
						}
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + oligosFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();
							//erase from file list
							oligosFileNames.erase(oligosFileNames.begin()+i);
							i--;
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (oligosFileNames.size() == 0) { m->mothurOut("no valid oligos files."); m->mothurOutEndLine(); abort = true; }
			}

			if (hasOligos) {
                split = 2;
				if (oligosFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide a oligos file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
            
			if (hasAccnos) {
				if (accnosFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide a accnos file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
			
			string temp = validParameter.validFile(parameters, "qfile", false);			if (temp == "not found"){	temp = "T";				}
			qual = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "fasta", false);				if (temp == "not found"){	temp = "T";				}
			fasta = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "flow", false);					if (temp == "not found"){	temp = "T";				}
			flow = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "trim", false);					if (temp == "not found"){	temp = "T";				}
			trim = m->isTrue(temp); 
            
            temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, pdiffs);
            
            temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, ldiffs);
            
            temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, sdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
            
			temp = validParameter.validFile(parameters, "sfftxt", false);				
			if (temp == "not found")	{	temp = "F";	 sfftxt = false; sfftxtFilename = "";		}
			else if (m->isTrue(temp))	{	sfftxt = true;		sfftxtFilename = "";				}
			else {
				//you are a filename
				if (inputDir != "") {
					map<string,string>::iterator it = parameters.find("sfftxt");
					//user has given a template file
					if(it != parameters.end()){ 
						string path = m->hasPath(it->second);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	parameters["sfftxt"] = inputDir + it->second;		}
					}
				}
				
				sfftxtFilename = validParameter.validFile(parameters, "sfftxt", true);
				if (sfftxtFilename == "not found") { sfftxtFilename = "";  }
				else if (sfftxtFilename == "not open") { sfftxtFilename = "";  }
			}
			
			if ((sfftxtFilename == "") && (filenames.size() == 0)) {  
				//if there is a current sff file, use it
				string filename = m->getSFFFile(); 
				if (filename != "") { filenames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the sff parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("[ERROR]: you must provide a valid sff or sfftxt file."); m->mothurOutEndLine(); abort=true;  }
			}
            
            
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "SffInfoCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < filenames.size(); s++) {
			
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
			
			int start = time(NULL);
			
            filenames[s] = m->getFullPathName(filenames[s]);
			m->mothurOut("Extracting info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
			
			string accnos = "";
			if (hasAccnos) { accnos = accnosFileNames[s]; }
            
            string oligos = "";
            if (hasOligos) { oligos = oligosFileNames[s]; }
			
			int numReads = extractSffInfo(filenames[s], accnos, oligos);

			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to extract " + toString(numReads) + ".");
		}
		
		if (sfftxtFilename != "") {  parseSffTxt(); }
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
		}
		
		itTypes = outputTypes.find("flow");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFlowFile(current); }
		}
		
		//report output filenames
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::extractSffInfo(string input, string accnos, string oligos){
	try {
		currentFileName = input;
		if (outputDir == "") {  outputDir += m->hasPath(input); }
		
		if (accnos != "")	{  readAccnosFile(accnos);  }
		else				{	seqNames.clear();		}
         
        if (oligos != "")   {   readOligos(oligos);  split = 2;   }

		ofstream outSfftxt, outFasta, outQual, outFlow;
		string outFastaFileName, outQualFileName;
        string rootName = outputDir + m->getRootName(m->getSimpleName(input));
        if(rootName.find_last_of(".") == rootName.npos){ rootName += "."; }
        
		string sfftxtFileName = outputDir + m->getRootName(m->getSimpleName(input)) + getOutputFileNameTag("sfftxt");
		string outFlowFileName = outputDir + m->getRootName(m->getSimpleName(input)) + getOutputFileNameTag("flow");
		if (trim) {
			outFastaFileName = outputDir + m->getRootName(m->getSimpleName(input)) + getOutputFileNameTag("fasta");
			outQualFileName = outputDir + m->getRootName(m->getSimpleName(input)) + getOutputFileNameTag("qfile");
		}else{
			outFastaFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "raw." + getOutputFileNameTag("fasta");
			outQualFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "raw." + getOutputFileNameTag("qfile");
		}
		
		if (sfftxt) { m->openOutputFile(sfftxtFileName, outSfftxt); outSfftxt.setf(ios::fixed, ios::floatfield); outSfftxt.setf(ios::showpoint);  outputNames.push_back(sfftxtFileName);  outputTypes["sfftxt"].push_back(sfftxtFileName); }
		if (fasta)	{ m->openOutputFile(outFastaFileName, outFasta);	outputNames.push_back(outFastaFileName); outputTypes["fasta"].push_back(outFastaFileName); }
		if (qual)	{ m->openOutputFile(outQualFileName, outQual);		outputNames.push_back(outQualFileName); outputTypes["qfile"].push_back(outQualFileName);  }
		if (flow)	{ m->openOutputFile(outFlowFileName, outFlow);		outputNames.push_back(outFlowFileName);  outFlow.setf(ios::fixed, ios::floatfield); outFlow.setf(ios::showpoint); outputTypes["flow"].push_back(outFlowFileName);  }
		
		ifstream in;
		in.open(input.c_str(), ios::binary);
		
		CommonHeader header; 
		readCommonHeader(in, header);
	
		int count = 0;
		mycount = 0;
		
		//check magic number and version
		if (header.magicNumber != 779314790) { m->mothurOut("Magic Number is not correct, not a valid .sff file"); m->mothurOutEndLine(); return count; }
		if (header.version != "0001") { m->mothurOut("Version is not supported, only support version 0001."); m->mothurOutEndLine(); return count; }
	
		//print common header
		if (sfftxt) {	printCommonHeader(outSfftxt, header);		}
		if (flow)	{	outFlow << header.numFlowsPerRead << endl;	}
			
		//read through the sff file
		while (!in.eof()) {
			
			bool print = true;
						
			//read data
			seqRead read;  Header readheader;
			readSeqData(in, read, header.numFlowsPerRead, readheader);
            bool okay = sanityCheck(readheader, read);
            if (!okay) { break; }
            
			//if you have provided an accosfile and this seq is not in it, then dont print
			if (seqNames.size() != 0) {   if (seqNames.count(readheader.name) == 0) { print = false; }  }
			
			//print 
			if (print) {
				if (sfftxt) { printHeader(outSfftxt, readheader); printSffTxtSeqData(outSfftxt, read, readheader); }
				if (fasta)	{	printFastaSeqData(outFasta, read, readheader);	}
				if (qual)	{	printQualSeqData(outQual, read, readheader);	}
				if (flow)	{	printFlowSeqData(outFlow, read, readheader);	}
			}
			
			count++;
			mycount++;
        
			//report progress
			if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
		
			if (m->control_pressed) { count = 0; break;   }
			
			if (count >= header.numReads) { break; }
		}
		
		//report progress
		if (!m->control_pressed) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
		
		in.close();
		
		if (sfftxt) {  outSfftxt.close();	}
		if (fasta)	{  outFasta.close();	}
		if (qual)	{  outQual.close();		}
		if (flow)	{  outFlow.close();		}
		
        if (split > 1) {
            //create new common headers for each file with the correct number of reads
            adjustCommonHeader(header);
            
			map<string, string>::iterator it;
			set<string> namesToRemove;
			for(int i=0;i<filehandles.size();i++){
				for(int j=0;j<filehandles[0].size();j++){
					if (filehandles[i][j] != "") {
						if (namesToRemove.count(filehandles[i][j]) == 0) {
							if(m->isBlank(filehandles[i][j])){
								m->mothurRemove(filehandles[i][j]);
                                m->mothurRemove(filehandlesHeaders[i][j]);
								namesToRemove.insert(filehandles[i][j]);
                            }
						}
					}
				}
			}
            
            //append new header to reads
            for (int i = 0; i < filehandles.size(); i++) {
                for (int j = 0; j < filehandles[i].size(); j++) {
                    m->appendFiles(filehandles[i][j], filehandlesHeaders[i][j]);
                    m->renameFile(filehandlesHeaders[i][j], filehandles[i][j]);
                    m->mothurRemove(filehandlesHeaders[i][j]);
                    if (numSplitReads[i][j] == 0) { m->mothurRemove(filehandles[i][j]); }
                }
            }
			
			//remove names for outputFileNames, just cleans up the output
			for(int i = 0; i < outputNames.size(); i++) { 
                if (namesToRemove.count(outputNames[i]) != 0) { 
                    outputNames.erase(outputNames.begin()+i);
                    i--;
                } 
            }
            
            if(m->isBlank(noMatchFile)){  m->mothurRemove(noMatchFile); }
            else { outputNames.push_back(noMatchFile); outputTypes["sff"].push_back(noMatchFile); }
        }
        
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "extractSffInfo");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::readCommonHeader(ifstream& in, CommonHeader& header){
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
				
			//read header length
			char buffer5 [2];
			in.read(buffer5, 2);
			header.headerLength =  be_int2(*(unsigned short *)(&buffer5));
					
			//read key length
			char buffer6 [2];
			in.read(buffer6, 2);
			header.keyLength = be_int2(*(unsigned short *)(&buffer6));
			
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
			if (header.keySequence.length() > header.keyLength) { header.keySequence = header.keySequence.substr(0, header.keyLength);  }
			delete[] tempBuffer2;
			
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
		m->errorOut(e, "SffInfoCommand", "readCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::adjustCommonHeader(CommonHeader header){
	try {

        char* mybuffer = new char[4];
        ifstream in;
        in.open(currentFileName.c_str(), ios::binary);
        
        //magic number
        in.read(mybuffer,4);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
        
        //version
        mybuffer = new char[4];
        in.read(mybuffer,4);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
        
        //offset
        mybuffer = new char[8];
        in.read(mybuffer,8);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
            
			
        //read index length
		mybuffer = new char[4];
        in.read(mybuffer,4);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
		
        //change num reads
        mybuffer = new char[4];
        in.read(mybuffer,4);
        delete[] mybuffer;
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                //convert number of reads to 4 byte char*
                char* thisbuffer = new char[4];
                thisbuffer[0] = (numSplitReads[i][j] >> 24) & 0xFF;
                thisbuffer[1] = (numSplitReads[i][j] >> 16) & 0xFF;
                thisbuffer[2] = (numSplitReads[i][j] >> 8) & 0xFF;
                thisbuffer[3] = numSplitReads[i][j] & 0xFF;
                out.write(thisbuffer, 4);
                out.close();
                delete[] thisbuffer;
            }
        }
            
        //read header length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
            
        //read key length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
			
        //read number of flow reads
        mybuffer = new char[2];
        in.read(mybuffer,2);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
            
        //read format code
        mybuffer = new char[1];
        in.read(mybuffer,1);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
			
        //read flow chars
        mybuffer = new char[header.numFlowsPerRead];
        in.read(mybuffer,header.numFlowsPerRead);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
			
        //read key
        mybuffer = new char[header.keyLength];
        in.read(mybuffer,header.keyLength);
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, in.gcount()); 
                out.close();
            }
        }
        delete[] mybuffer;
        
			
        /* Pad to 8 chars */
        unsigned long long spotInFile = in.tellg();
        unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
        in.seekg(spot);
        
        mybuffer = new char[spot-spotInFile];
        for (int i = 0; i < filehandlesHeaders.size(); i++) { 
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                m->openOutputFileAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, spot-spotInFile); 
                out.close();
            }
        }
        delete[] mybuffer;
        in.close();
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "adjustCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::readSeqData(ifstream& in, seqRead& read, int numFlowReads, Header& header){
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
			
			//extract info from name
			decodeName(header.timestamp, header.region, header.xy, header.name);
			
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
            
            if (split > 1) {
                char * mybuffer;
                mybuffer = new char [spot-startSpotInFile];
                ifstream in2;
                m->openInputFile(currentFileName, in2);
                in2.seekg(startSpotInFile);
                in2.read(mybuffer,spot-startSpotInFile);
                in2.close();
                
                int barcodeIndex, primerIndex;
                int trashCodeLength = findGroup(header, read, barcodeIndex, primerIndex);
                                
                if(trashCodeLength == 0){
                    ofstream out;
                    m->openOutputFileAppend(filehandles[barcodeIndex][primerIndex], out);
                    out.write(mybuffer, in2.gcount()); 
                    out.close();
                    delete[] mybuffer;
                    numSplitReads[barcodeIndex][primerIndex]++;
				}
				else{
					ofstream out;
                    m->openOutputFileAppend(noMatchFile, out);
                    out.write(mybuffer, in2.gcount()); 
                    out.close();
                    delete[] mybuffer;
				}
				
			}
		}else{
			m->mothurOut("Error reading."); m->mothurOutEndLine();
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::findGroup(Header header, seqRead read, int& barcode, int& primer) {
	try {
        //find group read belongs to
        TrimOligos trimOligos(pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimer, linker, spacer);
        
        int success = 1;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        string seq = read.bases;
        
        if (trim) {
            if(header.clipQualRight < header.clipQualLeft){
                seq = "NNNN";
            }
            else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
                seq = seq.substr((header.clipQualLeft-1), (header.clipQualRight-header.clipQualLeft));
            }
            else {
                seq = seq.substr(header.clipQualLeft-1);
            }
        }else{
            //if you wanted the sfftxt then you already converted the bases to the right case
            if (!sfftxt) {
                //make the bases you want to clip lowercase and the bases you want to keep upper case
                if(header.clipQualRight == 0){	header.clipQualRight = seq.length();	}
                for (int i = 0; i < (header.clipQualLeft-1); i++) { seq[i] = tolower(seq[i]);  }
                for (int i = (header.clipQualLeft-1); i < (header.clipQualRight-1); i++)  {   seq[i] = toupper(seq[i]);  }
                for (int i = (header.clipQualRight-1); i < seq.length(); i++) {   seq[i] = tolower(seq[i]);  }
            }
        }
        
        Sequence currSeq(header.name, seq);
        QualityScores currQual;
        
        if(numLinkers != 0){
            success = trimOligos.stripLinker(currSeq, currQual);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(barcodes.size() != 0){
            success = trimOligos.stripBarcode(currSeq, currQual, barcode);
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numSpacers != 0){
            success = trimOligos.stripSpacer(currSeq, currQual);
            if(success > sdiffs)		{	trashCode += 's';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numFPrimers != 0){
            success = trimOligos.stripForward(currSeq, currQual, primer, true);
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
        
        if(revPrimer.size() != 0){
            success = trimOligos.stripReverse(currSeq, currQual);
            if(!success)				{	trashCode += 'r';	}
        }

        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "findGroup");
		exit(1);
	}
}     
//**********************************************************************************************************************
int SffInfoCommand::decodeName(string& timestamp, string& region, string& xy, string name) {
	try {
		
		if (name.length() >= 6) {
			string time = name.substr(0, 6);
			unsigned int timeNum = m->fromBase36(time);
			
			int q1 = timeNum / 60;
			int sec = timeNum - 60 * q1;
			int q2 = q1 / 60;
			int minute = q1 - 60 * q2;
			int q3 = q2 / 24;
			int hr = q2 - 24 * q3;
			int q4 = q3 / 32;
			int day = q3 - 32 * q4;
			int q5 = q4 / 13;
			int mon = q4 - 13 * q5;
			int year = 2000 + q5;
		
			timestamp = toString(year) + "_" + toString(mon) + "_" + toString(day) + "_" + toString(hr) + "_" + toString(minute) + "_" + toString(sec);
		}
		
		if (name.length() >= 9) {
			region = name.substr(7, 2);
		
			string xyNum = name.substr(9);
			unsigned int myXy = m->fromBase36(xyNum);
			int x = myXy >> 12;
			int y = myXy & 4095;
		
			xy = toString(x) + "_" + toString(y);
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "decodeName");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printCommonHeader(ofstream& out, CommonHeader& header) {
	try {
	
		out << "Common Header:\nMagic Number: " << header.magicNumber << endl;
		out << "Version: " << header.version << endl;
		out << "Index Offset: " << header.indexOffset << endl;
		out << "Index Length: " << header.indexLength << endl;
		out << "Number of Reads: " << header.numReads << endl;
		out << "Header Length: " << header.headerLength << endl;
		out << "Key Length: " << header.keyLength << endl;
		out << "Number of Flows: " << header.numFlowsPerRead << endl;
		out << "Format Code: " << header.flogramFormatCode << endl;
		out << "Flow Chars: " << header.flowChars << endl;
		out << "Key Sequence: " << header.keySequence << endl << endl;
			
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printHeader(ofstream& out, Header& header) {
	try {
		
		out << ">" << header.name << endl;
		out << "Run Prefix: " << header.timestamp << endl;
		out << "Region #:  " << header.region << endl;
		out << "XY Location: " << header.xy << endl << endl;
		
		out << "Run Name:  " << endl;
		out << "Analysis Name:  " << endl;
		out << "Full Path: " << endl << endl;
		
		out << "Read Header Len: " << header.headerLength << endl;
		out << "Name Length: " << header.nameLength << endl;
		out << "# of Bases: " << header.numBases << endl;
		out << "Clip Qual Left: " << header.clipQualLeft << endl;
		out << "Clip Qual Right: " << header.clipQualRight << endl;
		out << "Clip Adap Left: " << header.clipAdapterLeft << endl;
		out << "Clip Adap Right: " << header.clipAdapterRight << endl << endl;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
bool SffInfoCommand::sanityCheck(Header& header, seqRead& read) {
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
		m->errorOut(e, "SffInfoCommand", "sanityCheck");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printSffTxtSeqData(ofstream& out, seqRead& read, Header& header) {
	try {
		out << "Flowgram: ";
		for (int i = 0; i < read.flowgram.size(); i++) { out << setprecision(2) << (read.flowgram[i]/(float)100) << '\t';  }
		
		out << endl <<  "Flow Indexes: ";
		int sum = 0;
		for (int i = 0; i < read.flowIndex.size(); i++) {  sum +=  read.flowIndex[i];  out << sum << '\t'; }
		
		//make the bases you want to clip lowercase and the bases you want to keep upper case
		if(header.clipQualRight == 0){	header.clipQualRight = read.bases.length();	}
		for (int i = 0; i < (header.clipQualLeft-1); i++) { read.bases[i] = tolower(read.bases[i]); }
		for (int i = (header.clipQualLeft-1); i < (header.clipQualRight-1); i++) {   read.bases[i] = toupper(read.bases[i]);  }
		for (int i = (header.clipQualRight-1); i < read.bases.length(); i++) {   read.bases[i] = tolower(read.bases[i]);  }
		
		out << endl <<  "Bases: " << read.bases << endl << "Quality Scores: ";
		for (int i = 0; i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';  }
	
		
		out << endl << endl;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printSffTxtSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printFastaSeqData(ofstream& out, seqRead& read, Header& header) {
	try {
		string seq = read.bases;
		
        if (trim) {
			if(header.clipQualRight < header.clipQualLeft){
				seq = "NNNN";
			}
			else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
				seq = seq.substr((header.clipQualLeft-1), (header.clipQualRight-header.clipQualLeft));
			}
			else {
				seq = seq.substr(header.clipQualLeft-1);
			}
		}else{
			//if you wanted the sfftxt then you already converted the bases to the right case
			if (!sfftxt) {
				//make the bases you want to clip lowercase and the bases you want to keep upper case
				if(header.clipQualRight == 0){	header.clipQualRight = seq.length();	}
				for (int i = 0; i < (header.clipQualLeft-1); i++) { seq[i] = tolower(seq[i]);  }
				for (int i = (header.clipQualLeft-1); i < (header.clipQualRight-1); i++)  {   seq[i] = toupper(seq[i]);  }
				for (int i = (header.clipQualRight-1); i < seq.length(); i++) {   seq[i] = tolower(seq[i]);  }
			}
		}
		
		out << ">" << header.name  << " xy=" << header.xy << endl;
		out << seq << endl;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printFastaSeqData");
		exit(1);
	}
}

//**********************************************************************************************************************
int SffInfoCommand::printQualSeqData(ofstream& out, seqRead& read, Header& header) {
	try {
		
		if (trim) {
			if(header.clipQualRight < header.clipQualLeft){
				out << ">" << header.name << " xy=" << header.xy << endl;
				out << "0\t0\t0\t0";
			}
			else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
				out << ">" << header.name << " xy=" << header.xy << " length=" << (header.clipQualRight-header.clipQualLeft) << endl;
				for (int i = (header.clipQualLeft-1); i < (header.clipQualRight-1); i++) {   out << read.qualScores[i] << '\t';	}
			}
			else{
				out << ">" << header.name << " xy=" << header.xy << " length=" << (header.clipQualRight-header.clipQualLeft) << endl;
				for (int i = (header.clipQualLeft-1); i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';	}			
			}
		}else{
			out << ">" << header.name << " xy=" << header.xy << " length=" << read.qualScores.size() << endl;
			for (int i = 0; i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';  }
		}
		
		out << endl;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printQualSeqData");
		exit(1);
	}
}

//**********************************************************************************************************************
int SffInfoCommand::printFlowSeqData(ofstream& out, seqRead& read, Header& header) {
	try {
		if(header.clipQualRight > header.clipQualLeft){
			
			int rightIndex = 0;
			for (int i = 0; i < header.clipQualRight; i++) {  rightIndex +=  read.flowIndex[i];	}

			out << header.name << ' ' << rightIndex;
			for (int i = 0; i < read.flowgram.size(); i++) { out << setprecision(2) << ' ' << (read.flowgram[i]/(float)100);  }
			out << endl;
		}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printFlowSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::readAccnosFile(string filename) {
	try {
		//remove old names
		seqNames.clear();
		
		ifstream in;
		m->openInputFile(filename, in);
		string name;
		
		while(!in.eof()){
			in >> name; m->gobble(in);
						
			seqNames.insert(name);
			
			if (m->control_pressed) { seqNames.clear(); break; }
		}
		in.close();		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readAccnosFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::parseSffTxt() {
	try {
		
		ifstream inSFF;
		m->openInputFile(sfftxtFilename, inSFF);
		
		if (outputDir == "") {  outputDir += m->hasPath(sfftxtFilename); }
		
		//output file names
		ofstream outFasta, outQual, outFlow;
		string outFastaFileName, outQualFileName;
		string fileRoot = m->getRootName(m->getSimpleName(sfftxtFilename));
		if (fileRoot.length() > 0) {
			//rip off last .
			fileRoot = fileRoot.substr(0, fileRoot.length()-1);
			fileRoot = m->getRootName(fileRoot);
		}
		
		string outFlowFileName = outputDir + fileRoot + getOutputFileNameTag("flow");
		if (trim) {
			outFastaFileName = outputDir + fileRoot + getOutputFileNameTag("fasta");
			outQualFileName = outputDir + fileRoot + getOutputFileNameTag("qfile");
		}else{
			outFastaFileName = outputDir + fileRoot + "raw." + getOutputFileNameTag("fasta");
			outQualFileName = outputDir + fileRoot + "raw." + getOutputFileNameTag("qfile");
		}
		
		if (fasta)	{ m->openOutputFile(outFastaFileName, outFasta);	outputNames.push_back(outFastaFileName); outputTypes["fasta"].push_back(outFastaFileName); }
		if (qual)	{ m->openOutputFile(outQualFileName, outQual);		outputNames.push_back(outQualFileName); outputTypes["qfile"].push_back(outQualFileName);  }
		if (flow)	{ m->openOutputFile(outFlowFileName, outFlow);		outputNames.push_back(outFlowFileName);  outFlow.setf(ios::fixed, ios::floatfield); outFlow.setf(ios::showpoint); outputTypes["flow"].push_back(outFlowFileName);  }
		
		//read common header
		string commonHeader = m->getline(inSFF);
		string magicNumber = m->getline(inSFF);	
		string version = m->getline(inSFF);
		string indexOffset = m->getline(inSFF);
		string indexLength = m->getline(inSFF);
		int numReads = parseHeaderLineToInt(inSFF);
		string headerLength = m->getline(inSFF);
		string keyLength = m->getline(inSFF);
		int numFlows = parseHeaderLineToInt(inSFF);
		string flowgramCode = m->getline(inSFF);
		string flowChars = m->getline(inSFF);
		string keySequence = m->getline(inSFF);
		m->gobble(inSFF);
		
		string seqName;
		
		if (flow)	{	outFlow << numFlows << endl;	}
		
		for(int i=0;i<numReads;i++){
			
			//sanity check
			if (inSFF.eof()) { m->mothurOut("[ERROR]: Expected " + toString(numReads) + " but reached end of file at " + toString(i+1) + "."); m->mothurOutEndLine(); break; }
			
			Header header;
			
			//parse read header
			inSFF >> seqName;
			seqName = seqName.substr(1);
			m->gobble(inSFF);
			header.name = seqName;
			
			string runPrefix = parseHeaderLineToString(inSFF);		header.timestamp = runPrefix;
			string regionNumber = parseHeaderLineToString(inSFF);	header.region = regionNumber;
			string xyLocation = parseHeaderLineToString(inSFF);		header.xy = xyLocation;
			m->gobble(inSFF);
				
			string runName = parseHeaderLineToString(inSFF);
			string analysisName = parseHeaderLineToString(inSFF);
			string fullPath = parseHeaderLineToString(inSFF);
			m->gobble(inSFF);
			
			string readHeaderLen = parseHeaderLineToString(inSFF);  convert(readHeaderLen, header.headerLength);
			string nameLength = parseHeaderLineToString(inSFF);		convert(nameLength, header.nameLength);
			int numBases = parseHeaderLineToInt(inSFF);				header.numBases = numBases;
			string clipQualLeft = parseHeaderLineToString(inSFF);	convert(clipQualLeft, header.clipQualLeft);
			int clipQualRight = parseHeaderLineToInt(inSFF);		header.clipQualRight = clipQualRight;
			string clipAdapLeft = parseHeaderLineToString(inSFF);	convert(clipAdapLeft, header.clipAdapterLeft);
			string clipAdapRight = parseHeaderLineToString(inSFF);	convert(clipAdapRight, header.clipAdapterRight);
			m->gobble(inSFF);
				
			seqRead read;
			
			//parse read
			vector<unsigned short> flowVector = parseHeaderLineToFloatVector(inSFF, numFlows);	read.flowgram = flowVector;
			vector<unsigned int> flowIndices = parseHeaderLineToIntVector(inSFF, numBases);	
			
			//adjust for print
			vector<unsigned int> flowIndicesAdjusted; flowIndicesAdjusted.push_back(flowIndices[0]);
			for (int j = 1; j < flowIndices.size(); j++) {   flowIndicesAdjusted.push_back(flowIndices[j] - flowIndices[j-1]);   }
			read.flowIndex = flowIndicesAdjusted;
			
			string bases = parseHeaderLineToString(inSFF);										read.bases = bases;
			vector<unsigned int> qualityScores = parseHeaderLineToIntVector(inSFF, numBases);	read.qualScores = qualityScores;
			m->gobble(inSFF);
					
			//if you have provided an accosfile and this seq is not in it, then dont print
			bool print = true;
			if (seqNames.size() != 0) {   if (seqNames.count(header.name) == 0) { print = false; }  }
			
			//print 
			if (print) {
				if (fasta)	{	printFastaSeqData(outFasta, read, header);	}
				if (qual)	{	printQualSeqData(outQual, read, header);	}
				if (flow)	{	printFlowSeqData(outFlow, read, header);	}
			}
			
			//report progress
			if((i+1) % 10000 == 0){	m->mothurOut(toString(i+1)); m->mothurOutEndLine();		}
			
			if (m->control_pressed) {  break;  }
		}
		
		//report progress
		if (!m->control_pressed) {   if((numReads) % 10000 != 0){	m->mothurOut(toString(numReads)); m->mothurOutEndLine();		}  }
		
		inSFF.close();
		
		if (fasta)	{  outFasta.close();	}
		if (qual)	{  outQual.close();		}
		if (flow)	{  outFlow.close();		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseSffTxt");
		exit(1);
	}
}
//**********************************************************************************************************************

int SffInfoCommand::parseHeaderLineToInt(ifstream& file){
	try {
		int number;
		
		while (!file.eof())	{
			
			char c = file.get(); 
			if (c == ':'){
				file >> number;
				break;
			}
			
		}
		m->gobble(file);
		return number;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToInt");
		exit(1);
	}
	
}

//**********************************************************************************************************************

string SffInfoCommand::parseHeaderLineToString(ifstream& file){
	try {
		string text;
		
		while (!file.eof())	{
			char c = file.get(); 
			
			if (c == ':'){
				//m->gobble(file);
				//text = m->getline(file);	
				file >> text;
				break;
			}
		}
		m->gobble(file);
		
		return text;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToString");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<unsigned short> SffInfoCommand::parseHeaderLineToFloatVector(ifstream& file, int length){
	try {
		vector<unsigned short> floatVector(length);
		
		while (!file.eof())	{
			char c = file.get(); 
			if (c == ':'){
				float temp;
				for(int i=0;i<length;i++){
					file >> temp;
					floatVector[i] = temp * 100;
				}
				break;
			}
		}
		m->gobble(file);	
		return floatVector;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToFloatVector");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<unsigned int> SffInfoCommand::parseHeaderLineToIntVector(ifstream& file, int length){
	try {
		vector<unsigned int> intVector(length);
		
		while (!file.eof())	{
			char c = file.get(); 
			if (c == ':'){
				for(int i=0;i<length;i++){
					file >> intVector[i];
				}
				break;
			}
		}
		m->gobble(file);	
		return intVector;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToIntVector");
		exit(1);
	}
}
//***************************************************************************************************************

bool SffInfoCommand::readOligos(string oligoFile){
	try {
        filehandles.clear();
        numSplitReads.clear();
        filehandlesHeaders.clear();
        
		ifstream inOligos;
		m->openInputFile(oligoFile, inOligos);
		
		string type, oligo, group;
        
		int indexPrimer = 0;
		int indexBarcode = 0;
		
		while(!inOligos.eof()){
            
			inOligos >> type; 
            
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				m->gobble(inOligos);
			}
			else{
				m->gobble(inOligos);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					group = "";
					
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{	
						char c = inOligos.get(); 
						if (c == 10 || c == 13){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
					
					//check for repeat barcodes
					map<string, int>::iterator itPrime = primers.find(oligo);
					if (itPrime != primers.end()) { m->mothurOut("primer " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
					
					primers[oligo]=indexPrimer; indexPrimer++;		
					primerNameVector.push_back(group);
				}else if(type == "REVERSE"){
					//Sequence oligoRC("reverse", oligo);
					//oligoRC.reverseComplement();
                    string oligoRC = reverseOligo(oligo);
					revPrimer.push_back(oligoRC);
				}
				else if(type == "BARCODE"){
					inOligos >> group;
					
					//check for repeat barcodes
					map<string, int>::iterator itBar = barcodes.find(oligo);
					if (itBar != barcodes.end()) { m->mothurOut("barcode " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
                    
					barcodes[oligo]=indexBarcode; indexBarcode++;
					barcodeNameVector.push_back(group);
				}else if(type == "LINKER"){
					linker.push_back(oligo);
				}else if(type == "SPACER"){
					spacer.push_back(oligo);
				}
				else{	m->mothurOut("[WARNING]: " + type + " is not recognized as a valid type. Choices are forward, reverse, and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine(); }
			}
			m->gobble(inOligos);
		}	
		inOligos.close();
		
		if(barcodeNameVector.size() == 0 && primerNameVector[0] == ""){	split = 1;	}
		
		//add in potential combos
		if(barcodeNameVector.size() == 0){
			barcodes[""] = 0;
			barcodeNameVector.push_back("");			
		}
		
		if(primerNameVector.size() == 0){
			primers[""] = 0;
			primerNameVector.push_back("");			
		}
		
		filehandles.resize(barcodeNameVector.size());
		for(int i=0;i<filehandles.size();i++){
			filehandles[i].assign(primerNameVector.size(), "");
		}
			
		if(split > 1){
			set<string> uniqueNames; //used to cleanup outputFileNames
			for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
				for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
					
					string primerName = primerNameVector[itPrimer->second];
					string barcodeName = barcodeNameVector[itBar->second];
					
					string comboGroupName = "";
					string fastaFileName = "";
					string qualFileName = "";
					string nameFileName = "";
					
					if(primerName == ""){
						comboGroupName = barcodeNameVector[itBar->second];
					}
					else{
						if(barcodeName == ""){
							comboGroupName = primerNameVector[itPrimer->second];
						}
						else{
							comboGroupName = barcodeNameVector[itBar->second] + "." + primerNameVector[itPrimer->second];
						}
					}
					
					ofstream temp;
					string thisFilename = outputDir + m->getRootName(m->getSimpleName(currentFileName)) + comboGroupName + "." + getOutputFileNameTag("sff");
					if (uniqueNames.count(thisFilename) == 0) {
						outputNames.push_back(thisFilename);
						outputTypes["sff"].push_back(thisFilename);
						uniqueNames.insert(thisFilename);
					}
					
					filehandles[itBar->second][itPrimer->second] = thisFilename;
					m->openOutputFile(thisFilename, temp);		temp.close();
				}
			}
		}
		numFPrimers = primers.size();
        numLinkers = linker.size();
        numSpacers = spacer.size();
		noMatchFile = outputDir + m->getRootName(m->getSimpleName(currentFileName)) + "scrap." + getOutputFileNameTag("sff");
        m->mothurRemove(noMatchFile);
        
		bool allBlank = true;
		for (int i = 0; i < barcodeNameVector.size(); i++) {
			if (barcodeNameVector[i] != "") {
				allBlank = false;
				break;
			}
		}
		for (int i = 0; i < primerNameVector.size(); i++) {
			if (primerNameVector[i] != "") {
				allBlank = false;
				break;
			}
		}
		
        filehandlesHeaders.resize(filehandles.size());
        numSplitReads.resize(filehandles.size());
        for (int i = 0; i < filehandles.size(); i++) { 
            numSplitReads[i].resize(filehandles[i].size(), 0); 
            for (int j = 0; j < filehandles[i].size(); j++) {
                filehandlesHeaders[i].push_back(filehandles[i][j]+"headers");
            }
        }
                             
		if (allBlank) {
			m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a split the sff file."); m->mothurOutEndLine();
			split = 1;
			return false;
		}
		
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readOligos");
		exit(1);
	}
}
//********************************************************************/
string SffInfoCommand::reverseOligo(string oligo){
	try {
        string reverse = "";
        
        for(int i=oligo.length()-1;i>=0;i--){
            
            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}
            
            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}
            
            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}
            
            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}
            
            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}
            
            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}
            
            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}
            
            else						{	reverse += 'N';	}
        }
        
        
        return reverse;
    }
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "reverseOligo");
		exit(1);
	}
}

//**********************************************************************************************************************


				
				
