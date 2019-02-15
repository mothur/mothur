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
		CommandParameter psff("sff", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(psff);
        CommandParameter poligos("oligos", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
        CommandParameter pgroup("group", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(pgroup);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);
		CommandParameter psfftxt("sfftxt", "String", "", "", "", "", "","",false,false); parameters.push_back(psfftxt);
		CommandParameter pflow("flow", "Boolean", "", "T", "", "", "","flow",false,false); parameters.push_back(pflow);
		CommandParameter ptrim("trim", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(ptrim);
		CommandParameter pfasta("fasta", "Boolean", "", "T", "", "", "","fasta",false,false); parameters.push_back(pfasta);
		CommandParameter pqfile("qfile", "Boolean", "", "T", "", "", "","qfile",false,false); parameters.push_back(pqfile);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The sffinfo command reads a sff file and extracts the sequence data, or you can use it to parse a sff file.\n";
		helpString += "The sffinfo command parameters are sff, fasta, qfile, accnos, flow, sfftxt, oligos, group, bdiffs, tdiffs, ldiffs, sdiffs, pdiffs, checkorient and trim. sff is required. \n";
		helpString += "The sff parameter allows you to enter the sff file you would like to extract data from.\n";
		helpString += "The fasta parameter allows you to indicate if you would like a fasta formatted file generated.  Default=True. \n";
		helpString += "The qfile parameter allows you to indicate if you would like a quality file generated.  Default=True. \n";
        helpString += "The oligos parameter allows you to provide an oligos file to split your sff file into separate sff files by barcode. \n";
        helpString += "The group parameter allows you to provide a group file to split your sff file into separate sff files by group. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. The default is false.\n";
		helpString += "The flow parameter allows you to indicate if you would like a flowgram file generated.  Default=True. \n";
		helpString += "The sfftxt parameter allows you to indicate if you would like a sff.txt file generated.  Default=False. \n";
		helpString += "If you want to parse an existing sfftxt file into flow, fasta and quality file, enter the file name using the sfftxt parameter. \n";
		helpString += "The trim parameter allows you to indicate if you would like a sequences and quality scores trimmed to the clipQualLeft and clipQualRight values.  Default=True. \n";
		helpString += "The accnos parameter allows you to provide a accnos file containing the names of the sequences you would like extracted.\n";
		helpString += "Example sffinfo(sff=mySffFile.sff, trim=F).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
string SffInfoCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern =  "[filename],fasta-[filename],[tag],fasta";   }
        else if (type == "flow")    {   pattern =  "[filename],flow";   }
        else if (type == "sfftxt")        {   pattern =  "[filename],sff.txt";   }
        else if (type == "sff")        {   pattern =  "[filename],[group],sff";   }
        else if (type == "qfile")       {   pattern =  "[filename],qual-[filename],[tag],qual";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SffInfoCommand", "getOutputPattern");
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
		hasAccnos = false; hasOligos = false; hasGroup = false;
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["flow"] = tempOutNames;
			outputTypes["sfftxt"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["sff"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");	  if (inputDir == "not found"){	inputDir = "";		}

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
							m->mothurOut("You have no current sfffile, ignoring current.\n"); ignore=true; 
							//erase from file list
							filenames.erase(filenames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(filenames[i], current->getLocations())) { current->setSFFFile(filenames[i]); }
                        else { m->mothurOut("Unable to open " + filenames[i] + ". It will be disregarded.\n"); filenames.erase(filenames.begin()+i); i--; } //erase from file list
                    }
				}
				
				//make sure there is at least one valid file left
				if (filenames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			accnosName = validParameter.valid(parameters, "accnos");
			if (accnosName == "not found") { accnosName = "";  }
			else { 
				hasAccnos = true;
				util.splitAtDash(accnosName, accnosFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < accnosFileNames.size(); i++) {
					bool ignore = false;
					if (accnosFileNames[i] == "current") { 
						accnosFileNames[i] = current->getAccnosFile(); 
						if (accnosFileNames[i] != "") {  m->mothurOut("Using " + accnosFileNames[i] + " as input file for the accnos parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current accnosfile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							accnosFileNames.erase(accnosFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(accnosFileNames[i], current->getLocations())) { current->setAccnosFile(accnosFileNames[i]); }
                        else { m->mothurOut("Unable to open " + accnosFileNames[i] + ". It will be disregarded.\n"); accnosFileNames.erase(accnosFileNames.begin()+i); i--; } //erase from file list
                    }
				}
				
				//make sure there is at least one valid file left
				if (accnosFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
            
            oligosfile = validParameter.valid(parameters, "oligos");
			if (oligosfile == "not found") { oligosfile = "";  }
			else { 
				hasOligos = true;
				util.splitAtDash(oligosfile, oligosFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < oligosFileNames.size(); i++) {
					bool ignore = false;
					if (oligosFileNames[i] == "current") { 
						oligosFileNames[i] = current->getOligosFile(); 
						if (oligosFileNames[i] != "") {  m->mothurOut("Using " + oligosFileNames[i] + " as input file for the oligos parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current oligosfile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							oligosFileNames.erase(oligosFileNames.begin()+i);
							i--;
						}
					}
                    
                    if (!ignore) {
                        if (util.checkLocations(oligosFileNames[i], current->getLocations())) { current->setOligosFile(oligosFileNames[i]); }
                        else { m->mothurOut("Unable to open " + oligosFileNames[i] + ". It will be disregarded.\n"); oligosFileNames.erase(oligosFileNames.begin()+i); i--; } //erase from file list
                    }
				}
				
				//make sure there is at least one valid file left
				if (oligosFileNames.size() == 0) { m->mothurOut("no valid oligos files."); m->mothurOutEndLine(); abort = true; }
			}
            
            groupfile = validParameter.valid(parameters, "group");
			if (groupfile == "not found") { groupfile = "";  }
			else {
				hasGroup = true;
				util.splitAtDash(groupfile, groupFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupFileNames.size(); i++) {
					bool ignore = false;
					if (groupFileNames[i] == "current") {
						groupFileNames[i] = current->getGroupFile();
						if (groupFileNames[i] != "") {  m->mothurOut("Using " + groupFileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else {
							m->mothurOut("You have no current group file, ignoring current."); m->mothurOutEndLine(); ignore=true;
							//erase from file list
							groupFileNames.erase(groupFileNames.begin()+i);
							i--;
						}
					}
                    
                    if (!ignore) {
                        if (util.checkLocations(groupFileNames[i], current->getLocations())) { current->setGroupFile(groupFileNames[i]); }
                        else { m->mothurOut("Unable to open " + groupFileNames[i] + ". It will be disregarded.\n"); groupFileNames.erase(groupFileNames.begin()+i); i--; } //erase from file list
                    }
				}
				
				//make sure there is at least one valid file left
				if (groupFileNames.size() == 0) { m->mothurOut("no valid group files."); m->mothurOutEndLine(); abort = true; }
			}

			if (hasGroup) {
                split = 2;
				if (groupFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide a group file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
            
            if (hasOligos) {
                split = 2;
				if (oligosFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide an oligos file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
            
            if (hasGroup && hasOligos) { m->mothurOut("You must enter ONLY ONE of the following: oligos or group."); m->mothurOutEndLine(); abort = true;}
            
			if (hasAccnos) {
				if (accnosFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide a accnos file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
			
			string temp = validParameter.valid(parameters, "qfile");			if (temp == "not found"){	temp = "T";				}
			qual = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "fasta");				if (temp == "not found"){	temp = "T";				}
			fasta = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "flow");					if (temp == "not found"){	temp = "T";				}
			flow = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "trim");					if (temp == "not found"){	temp = "T";				}
			trim = util.isTrue(temp); 
            
            temp = validParameter.valid(parameters, "bdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, bdiffs);
			
			temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs);
            
            temp = validParameter.valid(parameters, "ldiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, ldiffs);
            
            temp = validParameter.valid(parameters, "sdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, sdiffs);
			
			temp = validParameter.valid(parameters, "tdiffs");		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			util.mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
            
			temp = validParameter.valid(parameters, "sfftxt");
			if (temp == "not found")	{	temp = "F";	 sfftxt = false; sfftxtFilename = "";		}
			else if (util.isTrue(temp))	{	sfftxt = true;		sfftxtFilename = "";				}
			else {
				//you are a filename
				if (inputDir != "") {
					map<string,string>::iterator it = parameters.find("sfftxt");
					//user has given a template file
					if(it != parameters.end()){ 
						string path = util.hasPath(it->second);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	parameters["sfftxt"] = inputDir + it->second;		}
					}
				}
				
				sfftxtFilename = validParameter.validFile(parameters, "sfftxt");
				if (sfftxtFilename == "not found") { sfftxtFilename = "";  }
				else if (sfftxtFilename == "not open") { sfftxtFilename = "";  }
			}
			
			if ((sfftxtFilename == "") && (filenames.size() == 0)) {  
				//if there is a current sff file, use it
				string filename = current->getSFFFile(); 
				if (filename != "") { filenames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the sff parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("[ERROR]: you must provide a valid sff or sfftxt file."); m->mothurOutEndLine(); abort=true;  }
			}
            
            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "F"; }
			reorient = util.isTrue(temp);
            
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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
     
		for (int s = 0; s < filenames.size(); s++) {
			
			if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
			
			long start = time(NULL);
			
            filenames[s] = util.getFullPathName(filenames[s]);
			m->mothurOut("Extracting info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
			
			string accnos = "";
			if (hasAccnos) { accnos = accnosFileNames[s]; }
            
            string oligos = "";
            if (hasOligos) { oligos = oligosFileNames[s]; }
            if (hasGroup) { oligos = groupFileNames[s]; }
            
			int numReads = extractSffInfo(filenames[s], accnos, oligos);

			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to extract " + toString(numReads) + ".\n");
		}
		
		if (sfftxtFilename != "") {  parseSffTxt(); }
		
		if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
		}
		
		itTypes = outputTypes.find("flow");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFlowFile(currentName); }
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
        oligosObject = new Oligos();
		currentFileName = input;
		if (outputDir == "") {  outputDir += util.hasPath(input); }
		
		if (accnos != "")	{  readAccnosFile(accnos);  }
		else				{	seqNames.clear();		}
        
        TrimOligos* trimOligos = NULL; TrimOligos* rtrimOligos = NULL;
        if (hasOligos)   {
            readOligos(oligos);    split = 2;
            if (m->getControl_pressed()) { delete oligosObject; return 0; }
            trimOligos = new TrimOligos(pdiffs, bdiffs, ldiffs, sdiffs, oligosObject->getPrimers(), oligosObject->getBarcodes(), oligosObject->getReversePrimers(), oligosObject->getLinkers(), oligosObject->getSpacers());  numFPrimers = oligosObject->getPrimers().size(); numBarcodes = oligosObject->getBarcodes().size();
            if (reorient) {
                rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligosObject->getReorientedPairedPrimers(), oligosObject->getReorientedPairedBarcodes(), false); numBarcodes = oligosObject->getReorientedPairedBarcodes().size();
            }
        }
        
        if (hasGroup)    {   readGroup(oligos);     split = 2;      }
        
		ofstream outSfftxt, outFasta, outQual, outFlow;
		string outFastaFileName, outQualFileName;
        string rootName = outputDir + util.getRootName(util.getSimpleName(input));
        if(rootName.find_last_of(".") == rootName.npos){ rootName += "."; }
        
        map<string, string> variables; 
		variables["[filename]"] = rootName;
		string sfftxtFileName = getOutputFileName("sfftxt",variables);
		string outFlowFileName = getOutputFileName("flow",variables);
		if (!trim) { variables["[tag]"] = "raw"; }
		outFastaFileName = getOutputFileName("fasta",variables);
        outQualFileName = getOutputFileName("qfile",variables);
        
		if (sfftxt) { util.openOutputFile(sfftxtFileName, outSfftxt); outSfftxt.setf(ios::fixed, ios::floatfield); outSfftxt.setf(ios::showpoint);  outputNames.push_back(sfftxtFileName);  outputTypes["sfftxt"].push_back(sfftxtFileName); }
		if (fasta)	{ util.openOutputFile(outFastaFileName, outFasta);	outputNames.push_back(outFastaFileName); outputTypes["fasta"].push_back(outFastaFileName); }
		if (qual)	{ util.openOutputFile(outQualFileName, outQual);		outputNames.push_back(outQualFileName); outputTypes["qfile"].push_back(outQualFileName);  }
		if (flow)	{ util.openOutputFile(outFlowFileName, outFlow);		outputNames.push_back(outFlowFileName);  outFlow.setf(ios::fixed, ios::floatfield); outFlow.setf(ios::showpoint); outputTypes["flow"].push_back(outFlowFileName);  }
        
		ifstream in;
		util.openInputFileBinary(input, in);
		
		CommonHeader header;
		readCommonHeader(in, header);
        
		int count = 0;
		
		//check magic number and version
		if (header.magicNumber != 779314790) { m->mothurOut("Magic Number is not correct, not a valid .sff file"); m->mothurOutEndLine(); delete oligosObject; if (hasOligos)   { delete trimOligos; if (reorient) { delete  rtrimOligos; } } return count; }
		if (header.version != "0001") { m->mothurOut("Version is not supported, only support version 0001."); m->mothurOutEndLine(); delete oligosObject; if (hasOligos)   { delete trimOligos; if (reorient) { delete  rtrimOligos; } } return count; }
	
		//print common header
		if (sfftxt) {	printCommonHeader(outSfftxt, header);		}
		if (flow)	{	outFlow << header.numFlowsPerRead << endl;	}
        
		//read through the sff file
		while (!in.eof()) {
			
			bool print = true;
			
			//read data
			seqRead read;  Header readheader;
            readSeqData(in, read, header.numFlowsPerRead, readheader, trimOligos, rtrimOligos);
            
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
            
			//report progress
			if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
		
			if (m->getControl_pressed()) { count = 0; break;   }
			
			if (count >= header.numReads) { break; }
		}
		
		//report progress
		if (!m->getControl_pressed()) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
		
		in.close();
		
		if (sfftxt) {  outSfftxt.close();	}
		if (fasta)	{  outFasta.close();	}
		if (qual)	{  outQual.close();		}
		if (flow)	{  outFlow.close();		}
        
        if (split > 1) {
            //create new common headers for each file with the correct number of reads
            adjustCommonHeader(header);
            
            if (hasGroup) { delete groupMap; }
            
            
			map<string, string>::iterator it;
			set<string> namesToRemove;
			for(int i=0;i<filehandles.size();i++){
				for(int j=0;j<filehandles[0].size();j++){
                    
					if (filehandles[i][j] != "") {
						if (namesToRemove.count(filehandles[i][j]) == 0) {
							if(util.isBlank(filehandles[i][j])){
                               
								util.mothurRemove(filehandles[i][j]);
                                util.mothurRemove(filehandlesHeaders[i][j]);
								namesToRemove.insert(filehandles[i][j]);
                            }
						}
					}
				}
			}
            
            //append new header to reads
            for (int i = 0; i < filehandles.size(); i++) {
                for (int j = 0; j < filehandles[i].size(); j++) {
                    if (filehandles[i][j] != "") {
                        util.appendSFFFiles(filehandles[i][j], filehandlesHeaders[i][j]);
                        util.renameFile(filehandlesHeaders[i][j], filehandles[i][j]);
                        util.mothurRemove(filehandlesHeaders[i][j]);
                        
                        if (numSplitReads[i][j] == 0) { util.mothurRemove(filehandles[i][j]); }
                    }
                }
            }
			
			//remove names for outputFileNames, just cleans up the output
			for(int i = 0; i < outputNames.size(); i++) { 
                if (namesToRemove.count(outputNames[i]) != 0) {
                    
                    outputNames.erase(outputNames.begin()+i);
                    i--;
                }else { outputTypes["sff"].push_back(outputNames[i]); }
            }
            
            if(util.isBlank(noMatchFile)){  util.mothurRemove(noMatchFile); }
            else { outputNames.push_back(noMatchFile); outputTypes["sff"].push_back(noMatchFile); }
        }
        
        
        delete oligosObject;
        if (hasOligos)   { delete trimOligos; if (reorient) { delete  rtrimOligos; } }
        
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
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: numReads = " + toString(header.numReads) + "\n"); }
				
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
        string endian = util.findEdianness();
        
        char* mybuffer = new char[4];
        ifstream in;
        util.openInputFileBinary(currentFileName, in);
        
        ofstream outNoMatchHeader;
        string tempNoHeader = "tempNoMatchHeader";
        util.openOutputFileBinary(tempNoHeader, outNoMatchHeader);
        
        //magic number
        in.read(mybuffer,4);
        streamsize lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
        
        //version
        mybuffer = new char[4];
        in.read(mybuffer,4);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
        
        //offset
        mybuffer = new char[8];
        in.read(mybuffer,8);
        lengthRead = in.gcount();
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
        for (int i = 0; i < filehandlesHeaders.size(); i++) {
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(thisbuffer, 8);
                out.close();
            }
        }
        outNoMatchHeader.write(thisbuffer, 8);
        delete[] thisbuffer;
        delete[] mybuffer;
            
			
        //read index length
		mybuffer = new char[4];
        in.read(mybuffer,4);
        lengthRead = in.gcount();
        offset = 0;
        char* thisbuffer2 = new char[4];
        thisbuffer2[0] = (offset >> 24) & 0xFF;
        thisbuffer2[1] = (offset >> 16) & 0xFF;
        thisbuffer2[2] = (offset >> 8) & 0xFF;
        thisbuffer2[3] = offset & 0xFF;
        for (int i = 0; i < filehandlesHeaders.size(); i++) {
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(thisbuffer2, 4);
                out.close();
            }
        }
        outNoMatchHeader.write(thisbuffer2, 4);
        delete[] thisbuffer2;
        delete[] mybuffer;
		
        //change num reads
        mybuffer = new char[4];
        in.read(mybuffer,4);
        lengthRead = in.gcount();
        delete[] mybuffer;
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                
                char* thisbuffer = new char[4];
                thisbuffer[0] = (numSplitReads[i][j] >> 24) & 0xFF;
                thisbuffer[1] = (numSplitReads[i][j] >> 16) & 0xFF;
                thisbuffer[2] = (numSplitReads[i][j] >> 8) & 0xFF;
                thisbuffer[3] = numSplitReads[i][j] & 0xFF;
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(thisbuffer, 4);
                out.close();
                delete[] thisbuffer;
            }
        }
        char* thisbuffer3 = new char[4];
        thisbuffer3[0] = (numNoMatch >> 24) & 0xFF;
        thisbuffer3[1] = (numNoMatch >> 16) & 0xFF;
        thisbuffer3[2] = (numNoMatch >> 8) & 0xFF;
        thisbuffer3[3] = numNoMatch & 0xFF;
        outNoMatchHeader.write(thisbuffer3, 4);
        delete[] thisbuffer3;
        
        
        //read header length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
            
        //read key length
        mybuffer = new char[2];
        in.read(mybuffer,2);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
			
        //read number of flow reads
        mybuffer = new char[2];
        in.read(mybuffer,2);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
            
        //read format code
        mybuffer = new char[1];
        in.read(mybuffer,1);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
			
        //read flow chars
        mybuffer = new char[header.numFlowsPerRead];
        in.read(mybuffer,header.numFlowsPerRead);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
			
        //read key
        mybuffer = new char[header.keyLength];
        in.read(mybuffer,header.keyLength);
        lengthRead = in.gcount();
        for (int i = 0; i < filehandlesHeaders.size(); i++) {  
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                
                out.write(mybuffer, lengthRead);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, lengthRead);
        delete[] mybuffer;
        
			
        /* Pad to 8 chars */
        unsigned long long spotInFile = in.tellg();
        unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
        in.seekg(spot);
        
        mybuffer = new char[spot-spotInFile];
        for (int i = 0; i < filehandlesHeaders.size(); i++) { 
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                ofstream out;
                util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                out.write(mybuffer, spot-spotInFile);
                out.close();
            }
        }
        outNoMatchHeader.write(mybuffer, spot-spotInFile);
        outNoMatchHeader.close();
        delete[] mybuffer;
        in.close();
        
        util.appendSFFFiles(noMatchFile, tempNoHeader);
        util.renameFile(tempNoHeader, noMatchFile);
        util.mothurRemove(tempNoHeader);
        
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "adjustCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printCommonHeaderForDebug(CommonHeader& header, ofstream& out, int numReads){
    try {
        string endian = util.findEdianness();
      
        ifstream in;
        util.openInputFileBinary(currentFileName, in);
        
        //magic number
        char* mybuffer = new char[4];
        in.read(mybuffer,4);
        out.write(mybuffer, in.gcount());
        unsigned int magic = be_int4(*(unsigned int *)(mybuffer));
        string contents = toString(magic);
        m->mothurOut("magicNumber = " + contents + "\n");
        delete[] mybuffer;
        
        //version
        char* mybuffer1 = new char[4];
        in.read(mybuffer1,4);
        out.write(mybuffer1, in.gcount());
        contents = "";
        for (int i = 0; i < 4; i++) {  contents += toString((int)(mybuffer1[i]));  }
        m->mothurOut("version = " + contents + "\n");
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
        //m->mothurOut("numReads = " + contents + "\n");
        unsigned int numTReads = be_int4(*(unsigned int *)(mybuffer4));
        m->mothurOut("numReads = " + toString(numTReads) + "\n");
        //m->mothurOut("numReads = " + toString(header.numReads) + "\n");
        delete[] thisbuffer3;
        delete[] mybuffer4;
        
        //read header length
        char* mybuffer5 = new char[2];
        in.read(mybuffer5,2);
        out.write(mybuffer5, in.gcount());
        unsigned short hl = be_int2(*(unsigned short *)(mybuffer5));
        contents = toString(hl);
        m->mothurOut("readLength = " + contents + "\n");
        delete[] mybuffer5;
        
        //read key length
        char* mybuffer6 = new char[2];
        in.read(mybuffer6,2);
        out.write(mybuffer6, in.gcount());
        unsigned short kl = be_int2(*(unsigned short *)(mybuffer6));
        contents = toString(kl);
        m->mothurOut("key length = " + contents + "\n");
        delete[] mybuffer6;
        
        //read number of flow reads
        char* mybuffer7 = new char[2];
        in.read(mybuffer7,2);
        out.write(mybuffer7, in.gcount());
        contents = mybuffer7;
        //m->mothurOut("num flow reads = " + contents + "\n");
        int numFlowReads = be_int2(*(unsigned short *)(mybuffer7));
        m->mothurOut("num flow Reads = " + toString(numFlowReads) + "\n");
        delete[] mybuffer7;
        
        //read format code
        char* mybuffer8 = new char[1];
        in.read(mybuffer8,1);
        out.write(mybuffer8, in.gcount());
        int fc = (int)(mybuffer8[0]);
        contents = toString(fc);
        m->mothurOut("read format code = " + contents + "\n");
        delete[] mybuffer8;
        
        //read flow chars
        char* mybuffer9 = new char[header.numFlowsPerRead];
        in.read(mybuffer9,header.numFlowsPerRead);
        out.write(mybuffer9, in.gcount());
        contents = mybuffer9;
        m->mothurOut("flow chars = " + contents + "\n");
        delete[] mybuffer9;
        
        //read key
        char* mybuffer10 = new char[header.keyLength];
        in.read(mybuffer10,header.keyLength);
        out.write(mybuffer10, in.gcount());
        contents = mybuffer10;
        m->mothurOut("key = " + contents + "\n");
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
        m->errorOut(e, "SffInfoCommand", "printCommonHeaderForDebug");
        exit(1);
    }
}
//**********************************************************************************************************************
bool SffInfoCommand::readSeqData(ifstream& in, seqRead& read, int numFlowReads, Header& header, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos){
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
			//header.clipQualLeft = 5;
			
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
			read.qualScores.resize(header.numBases, 0);
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
               
                int barcodeIndex, primerIndex, trashCodeLength;
                
                if (hasOligos)      {  trashCodeLength = findGroup(header, read, barcodeIndex, primerIndex, trimOligos, rtrimOligos);                }
                else if (hasGroup)  {  trashCodeLength = findGroup(header, read, barcodeIndex, primerIndex, "groupMode");   }
                else {  m->mothurOut("[ERROR]: uh oh, we shouldn't be here...\n"); }

                char * mybuffer;
                mybuffer = new char [spot-startSpotInFile];
                
                ifstream in2;
                util.openInputFileBinary(currentFileName, in2);
                in2.seekg(startSpotInFile);
                in2.read(mybuffer,spot-startSpotInFile);
                
                if(trashCodeLength == 0){
                    ofstream out;
                    util.openOutputFileBinaryAppend(filehandles[barcodeIndex][primerIndex], out);
                    out.write(mybuffer, in2.gcount());
                    out.close();
                    numSplitReads[barcodeIndex][primerIndex]++;
				}
				else{
					ofstream out;
                    util.openOutputFileBinaryAppend(noMatchFile, out);
                    out.write(mybuffer, in2.gcount());
                    out.close();
                    numNoMatch++;
				}
				delete[] mybuffer;
                in2.close();
        }    
            
		}else{ m->mothurOut("Error reading.\n"); }
        
        if (in.eof()) {  return true; }
        
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::findGroup(Header header, seqRead read, int& barcode, int& primer, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos) {
	try {
        
        int success = 1;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        string seq = read.bases;
        vector<int> quals;
        for (int i = 0; i < read.qualScores.size(); i++) {   quals.push_back((int)read.qualScores[i]); }

        for (int i = 0; i < seq.length(); i++) { seq[i] = toupper(seq[i]); }
        
        if (trim) {
            if(header.clipQualRight < header.clipQualLeft){
                if (header.clipQualRight == 0) { //don't trim right
                    
                    vector<int>::iterator first = quals.begin() + (header.clipQualLeft-1);
                    vector<int>::iterator last = quals.begin() + (quals.size()-1);
                    vector<int> newVec(first, last);
                    quals = newVec;
                    seq = seq.substr(header.clipQualLeft-1);
                }else {
                    seq = "NNNN";
                    quals.clear(); quals.resize(4, 0);
                }
            }
            else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
                seq = seq.substr((header.clipQualLeft-1), (header.clipQualRight-header.clipQualLeft+1));
                vector<int>::iterator first = quals.begin() + (header.clipQualLeft-1);
                vector<int>::iterator last = quals.begin() + (header.clipQualRight-header.clipQualLeft+1);
                vector<int> newVec(first, last);
                quals = newVec;
            }
            else {
                vector<int>::iterator first = quals.begin() + (header.clipQualLeft-1);
                vector<int>::iterator last = quals.begin() + (quals.size()-1);
                vector<int> newVec(first, last);
                quals = newVec;
                seq = seq.substr(header.clipQualLeft-1);
            }
        }
        
        Sequence currSeq(header.name, seq);
        QualityScores currQual(header.name, quals);
        
        //for reorient
        Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
        QualityScores savedQual(currQual.getName(), currQual.getScores());
        
        if(numLinkers != 0){
            success = trimOligos->stripLinker(currSeq, currQual);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numBarcodes != 0){
            vector<int> results = trimOligos->stripBarcode(currSeq, currQual, barcode);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numSpacers != 0){
            success = trimOligos->stripSpacer(currSeq, currQual);
            if(success > sdiffs)		{	trashCode += 's';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numFPrimers != 0){
            vector<int> results = trimOligos->stripForward(currSeq, currQual, primer, true);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numRPrimers != 0){
            vector<int> results = trimOligos->stripReverse(currSeq, currQual);
            success = results[0];
            if(success > pdiffs)		{	trashCode += 'r';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
        

        if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
            int thisSuccess = 0;
            string thisTrashCode = "";
            int thisCurrentSeqsDiffs = 0;
            
            int thisBarcodeIndex = 0;
            int thisPrimerIndex = 0;
            
            if(numBarcodes != 0){
                vector<int> results = rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                }
                if(thisSuccess > bdiffs)		{ thisTrashCode += "b"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if(numFPrimers != 0){
                vector<int> results = rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, true);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                }
                if(thisSuccess > pdiffs)		{ thisTrashCode += "f"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if (thisCurrentSeqsDiffs > tdiffs)	{	thisTrashCode += 't';   }
            
            if (thisTrashCode == "") {
                trashCode = thisTrashCode;
                success = thisSuccess;
                currentSeqsDiffs = thisCurrentSeqsDiffs;
                barcode = thisBarcodeIndex;
                primer = thisPrimerIndex;
                savedSeq.reverseComplement();
                currSeq.setAligned(savedSeq.getAligned());
                savedQual.flipQScores();
                currQual.setScores(savedQual.getScores());
            }else { trashCode += "(" + thisTrashCode + ")";  }
        }

        if (trashCode.length() == 0) { //is this sequence in the ignore group
            string thisGroup = oligosObject->getGroupName(barcode, primer);
            
            int pos = thisGroup.find("ignore");
            if (pos != string::npos) {  trashCode += "i"; }
        }
        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::findGroup(Header header, seqRead read, int& barcode, int& primer, string groupMode) {
	try {
        string trashCode = "";
        primer = 0;
        
        string group = groupMap->getGroup(header.name);
        if (group == "not found") {     trashCode += "g";   } //scrap for group
        else {  barcode = GroupToFile[group]; }
        
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
			unsigned int timeNum = util.fromBase36(time);
			
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
			unsigned int myXy = util.fromBase36(xyNum);
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
        
        if (!okay) { m->mothurOut(message); m->mothurOutEndLine(); }
        
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
        int endValue = header.clipQualRight;
		if(endValue == 0){	endValue = read.bases.length();	}
		for (int i = 0; i < (header.clipQualLeft-1); i++) { read.bases[i] = tolower(read.bases[i]); }
		for (int i = (header.clipQualLeft-1); i < (endValue-1); i++) {   read.bases[i] = toupper(read.bases[i]);  }
		for (int i = (endValue-1); i < read.bases.length(); i++) {   read.bases[i] = tolower(read.bases[i]);  }
		
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
				if (header.clipQualRight == 0) { //don't trim right
                    seq = seq.substr(header.clipQualLeft-1);
                }else {
                    seq = "NNNN";
                }
			}
			else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
				seq = seq.substr((header.clipQualLeft-1), (header.clipQualRight-header.clipQualLeft+1));
			}
			else {
				seq = seq.substr(header.clipQualLeft-1);
			}
		}else{
			//if you wanted the sfftxt then you already converted the bases to the right case
			if (!sfftxt) {
                int endValue = header.clipQualRight;
				//make the bases you want to clip lowercase and the bases you want to keep upper case
				if(endValue == 0){	endValue = seq.length();	}
				for (int i = 0; i < (header.clipQualLeft-1); i++) { seq[i] = tolower(seq[i]);  }
				for (int i = (header.clipQualLeft-1); i < (endValue-1); i++)  {   seq[i] = toupper(seq[i]);  }
				for (int i = (endValue-1); i < seq.length(); i++) {   seq[i] = tolower(seq[i]);  }
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
                if (header.clipQualRight == 0) { //don't trim right
                    out << ">" << header.name << " xy=" << header.xy << " length=" << (read.qualScores.size()-header.clipQualLeft) << endl;
                    for (int i = (header.clipQualLeft-1); i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';	}	
                }else {
                    out << ">" << header.name << " xy=" << header.xy << endl;
                    out << "0\t0\t0\t0";
                }
			}
			else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
				out << ">" << header.name << " xy=" << header.xy << " length=" << (header.clipQualRight-header.clipQualLeft+1) << endl;
				for (int i = (header.clipQualLeft-1); i < (header.clipQualRight); i++) {   out << read.qualScores[i] << '\t';	}
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
        
        int endValue = header.clipQualRight;
        if (header.clipQualRight == 0) {
            endValue = read.flowIndex.size();
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + header.name + " has clipQualRight=0.\n"); }
        }
        if(endValue > header.clipQualLeft){
            
            int rightIndex = 0;
            for (int i = 0; i < endValue; i++) {  rightIndex +=  read.flowIndex[i];	 }
            
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
		util.openInputFile(filename, in);
		string name;
		
		while(!in.eof()){
			in >> name; util.gobble(in);
						
			seqNames.insert(name);
			
			if (m->getControl_pressed()) { seqNames.clear(); break; }
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
		util.openInputFile(sfftxtFilename, inSFF);
		
		if (outputDir == "") {  outputDir += util.hasPath(sfftxtFilename); }
		
		//output file names
		ofstream outFasta, outQual, outFlow;
		string outFastaFileName, outQualFileName;
		string fileRoot = util.getRootName(util.getSimpleName(sfftxtFilename));
		if (fileRoot.length() > 0) {
			//rip off last .
			fileRoot = fileRoot.substr(0, fileRoot.length()-1);
			fileRoot = util.getRootName(fileRoot);
		}
		
        map<string, string> variables; 
		variables["[filename]"] = fileRoot;
		string sfftxtFileName = getOutputFileName("sfftxt",variables);
		string outFlowFileName = getOutputFileName("flow",variables);
		if (!trim) { variables["[tag]"] = "raw"; }
		outFastaFileName = getOutputFileName("fasta",variables);
        outQualFileName = getOutputFileName("qfile",variables);
		
		if (fasta)	{ util.openOutputFile(outFastaFileName, outFasta);	outputNames.push_back(outFastaFileName); outputTypes["fasta"].push_back(outFastaFileName); }
		if (qual)	{ util.openOutputFile(outQualFileName, outQual);		outputNames.push_back(outQualFileName); outputTypes["qfile"].push_back(outQualFileName);  }
		if (flow)	{ util.openOutputFile(outFlowFileName, outFlow);		outputNames.push_back(outFlowFileName);  outFlow.setf(ios::fixed, ios::floatfield); outFlow.setf(ios::showpoint); outputTypes["flow"].push_back(outFlowFileName);  }
		
		//read common header
		string commonHeader = util.getline(inSFF);
		string magicNumber = util.getline(inSFF);	
		string version = util.getline(inSFF);
		string indexOffset = util.getline(inSFF);
		string indexLength = util.getline(inSFF);
		int numReads = parseHeaderLineToInt(inSFF);
		string headerLength = util.getline(inSFF);
		string keyLength = util.getline(inSFF);
		int numFlows = parseHeaderLineToInt(inSFF);
		string flowgramCode = util.getline(inSFF);
		string flowChars = util.getline(inSFF);
		string keySequence = util.getline(inSFF);
		util.gobble(inSFF);
		
		string seqName;
		
		if (flow)	{	outFlow << numFlows << endl;	}
		
		for(int i=0;i<numReads;i++){
			
			//sanity check
			if (inSFF.eof()) { m->mothurOut("[ERROR]: Expected " + toString(numReads) + " but reached end of file at " + toString(i+1) + "."); m->mothurOutEndLine(); break; }
			
			Header header;
			
			//parse read header
			inSFF >> seqName;
			seqName = seqName.substr(1);
			util.gobble(inSFF);
			header.name = seqName;
			
			string runPrefix = parseHeaderLineToString(inSFF);		header.timestamp = runPrefix;
			string regionNumber = parseHeaderLineToString(inSFF);	header.region = regionNumber;
			string xyLocation = parseHeaderLineToString(inSFF);		header.xy = xyLocation;
			util.gobble(inSFF);
				
			string runName = parseHeaderLineToString(inSFF);
			string analysisName = parseHeaderLineToString(inSFF);
			string fullPath = parseHeaderLineToString(inSFF);
			util.gobble(inSFF);
			
			string readHeaderLen = parseHeaderLineToString(inSFF);  convert(readHeaderLen, header.headerLength);
			string nameLength = parseHeaderLineToString(inSFF);		convert(nameLength, header.nameLength);
			int numBases = parseHeaderLineToInt(inSFF);				header.numBases = numBases;
			string clipQualLeft = parseHeaderLineToString(inSFF);	convert(clipQualLeft, header.clipQualLeft);
			int clipQualRight = parseHeaderLineToInt(inSFF);		header.clipQualRight = clipQualRight;
			string clipAdapLeft = parseHeaderLineToString(inSFF);	convert(clipAdapLeft, header.clipAdapterLeft);
			string clipAdapRight = parseHeaderLineToString(inSFF);	convert(clipAdapRight, header.clipAdapterRight);
			util.gobble(inSFF);
				
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
			util.gobble(inSFF);
					
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
			
			if (m->getControl_pressed()) {  break;  }
		}
		
		//report progress
		if (!m->getControl_pressed()) {   if((numReads) % 10000 != 0){	m->mothurOut(toString(numReads)); m->mothurOutEndLine();		}  }
		
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
		util.gobble(file);
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
				//util.gobble(file);
				//text = util.getline(file);	
				file >> text;
				break;
			}
		}
		util.gobble(file);
		
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
		util.gobble(file);	
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
		util.gobble(file);	
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
        
		bool allBlank = false;
        oligosObject->read(oligoFile);
        
        if (m->getControl_pressed()) { return false; } //error in reading oligos
        
        if (oligosObject->hasPairedPrimers() || oligosObject->hasPairedBarcodes()) {
            pairedOligos = true;
            m->mothurOut("[ERROR]: sffinfo does not support paired barcodes and primers, aborting.\n"); m->setControl_pressed(true); return true;
        }else {
            pairedOligos = false;
            numFPrimers = oligosObject->getPrimers().size();
            numBarcodes = oligosObject->getBarcodes().size();
        }
        
        numLinkers = oligosObject->getLinkers().size();
        numSpacers = oligosObject->getSpacers().size();
        numRPrimers = oligosObject->getReversePrimers().size();
        
        vector<string> groupNames = oligosObject->getGroupNames();
        if (groupNames.size() == 0) { allBlank = true;  }
        
        filehandles.resize(oligosObject->getBarcodeNames().size());
		for(int i=0;i<filehandles.size();i++){
            for(int j=0;j<oligosObject->getPrimerNames().size();j++){  filehandles[i].push_back(""); }
		}
        
        if(split > 1){
            set<string> uniqueNames; //used to cleanup outputFileNames
            map<string, int> barcodes = oligosObject->getBarcodes() ;
            map<string, int> primers = oligosObject->getPrimers();
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligosObject->getPrimerName(itPrimer->second);
                    string barcodeName = oligosObject->getBarcodeName(itBar->second);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        if(itPrimer->first == ""){
                            comboName = itBar->first;
                        }else{
                            if(itBar->first == ""){
                                comboName = itPrimer->first;
                            }
                            else{
                                comboName = itBar->first + "." + itPrimer->first;
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        
                        ofstream temp;
                        map<string, string> variables;
                        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(currentFileName));
                        variables["[group]"] = comboGroupName;
                        string thisFilename = getOutputFileName("sff",variables);
                        if (uniqueNames.count(thisFilename) == 0) {
                            outputNames.push_back(thisFilename);
                            uniqueNames.insert(thisFilename);
                        }
                        
                        filehandles[itBar->second][itPrimer->second] = thisFilename;
                        util.openOutputFileBinary(thisFilename, temp);		temp.close();
                    }
                }
            }
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(currentFileName));
        variables["[group]"] = "scrap";
		noMatchFile = getOutputFileName("sff",variables);
        util.mothurRemove(noMatchFile);
        numNoMatch = 0;
		
        filehandlesHeaders.resize(filehandles.size());
        numSplitReads.resize(filehandles.size());
        for (int i = 0; i < filehandles.size(); i++) {
            numSplitReads[i].resize(filehandles[i].size(), 0);
            for (int j = 0; j < filehandles[i].size(); j++) {
                filehandlesHeaders[i].push_back(filehandles[i][j]+"headers");
                ofstream temp;
                util.openOutputFileBinary(filehandles[i][j]+"headers", temp); temp.close();
            }
        }
        
		if (allBlank) {
			m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a split the sff file.\n"); 
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
//***************************************************************************************************************

bool SffInfoCommand::readGroup(string oligoFile){
	try {
        filehandles.clear();
        numSplitReads.clear();
        filehandlesHeaders.clear();
        
        groupMap = new GroupMap();
        groupMap->readMap(oligoFile);
    
        //like barcodeNameVector - no primer names
        vector<string> groups = groupMap->getNamesOfGroups();
		
		filehandles.resize(groups.size());
        for (int i = 0; i < filehandles.size(); i++) {
                map<string, string> variables;
                variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(currentFileName));
                variables["[group]"] = groups[i];
                string thisFilename = getOutputFileName("sff",variables);
                outputNames.push_back(thisFilename);
               
                ofstream temp;
                util.openOutputFileBinary(thisFilename, temp); temp.close();
                filehandles[i].push_back(thisFilename);
                GroupToFile[groups[i]] = i;
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(currentFileName));
        variables["[group]"] = "scrap";
		noMatchFile = getOutputFileName("sff",variables);
        util.mothurRemove(noMatchFile);
        numNoMatch = 0;
        
		
        filehandlesHeaders.resize(groups.size());
        numSplitReads.resize(filehandles.size());
        for (int i = 0; i < filehandles.size(); i++) {
            numSplitReads[i].resize(filehandles[i].size(), 0);
            for (int j = 0; j < filehandles[i].size(); j++) {
                string thisHeader = filehandles[i][j]+"headers";
                filehandlesHeaders[i].push_back(thisHeader);
                ofstream temp;
                util.openOutputFileBinary(thisHeader, temp); temp.close();
            }
        }
		
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readGroup");
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


				
				
