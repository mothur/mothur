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


//**************************************************************************************
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
        
        abort = false; calledHelp = false;
        hasAccnos = false; hasOligos = false; hasGroup = false;
        split = 1;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["flow"] = tempOutNames;
        outputTypes["sfftxt"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
        outputTypes["sff"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "setParameters");
		exit(1);
	}
}
//****************************************************************************************
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

//************************************************************************************
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
//*******************************************************************************
SffInfoCommand::SffInfoCommand(string option)  {
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

            sffFilename = validParameter.validFile(parameters, "sff");
            if (sffFilename == "not found") {
                sffFilename = current->getSFFFile();
                if (sffFilename != "") { m->mothurOut("Using " + sffFilename + " as input file for the sff parameter.\n"); }
                else { 	m->mothurOut("[ERROR]: You have no current sff file and the sff parameter is required.\n");  abort = true; }
            }
            else if (sffFilename == "not open") { abort = true; }
            else { current->setSFFFile(sffFilename); }
			
            accnosName = validParameter.validFile(parameters, "accnos");
            if (accnosName == "not found") { accnosName = ""; }
            else if (accnosName == "not open") { accnosName = ""; abort = true; }
            else { current->setAccnosFile(accnosName);  hasAccnos = true; }
            
            oligosfile = validParameter.validFile(parameters, "oligos");
			if (oligosfile == "not found") { oligosfile = "";  }
            else if (oligosfile == "not open") { oligosfile = ""; abort = true; }
            else { current->setOligosFile(oligosfile); hasOligos = true;  }
            
            groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") { groupfile = "";  }
            else if (groupfile == "not open") { groupfile = ""; abort = true; }
            else { current->setGroupFile(groupfile);  hasGroup = true; }
            
            sfftxtFilename = validParameter.valid(parameters, "sfftxt");
            if (sfftxtFilename == "not found")      { sfftxt = false; sfftxtFilename = "";          }
            else if (util.isTrue(sfftxtFilename))	{	sfftxt = true;		sfftxtFilename = "";	}
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
                else if (sfftxtFilename == "not open") { sfftxtFilename = "";  abort = true; }
            }

			if ((hasGroup) || (hasOligos)) { split = 2; }
            
            if (hasGroup && hasOligos) { m->mothurOut("[ERROR]: You may enter ONLY ONE of the following: oligos or group.\n"); abort = true; }
			
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
            
            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "F"; }
			reorient = util.isTrue(temp);
            
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "SffInfoCommand");
		exit(1);
	}
}
//********************************************************************************
int SffInfoCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
	
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
        
        long start = time(NULL);
        
        sffFilename = util.getFullPathName(sffFilename);
        m->mothurOut("Extracting info from " + sffFilename + " ...\n" );
        
        string oligos = "";
        if (hasOligos) { oligos = oligosfile; }
        if (hasGroup) { oligos = groupfile; }
        
        int numReads = extractSffInfo(sffFilename, accnosName, oligos);
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to extract " + toString(numReads) + ".\n");
		
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
		m->mothurOut("\nOutput File Names:\n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]+"\n"); 	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************
int SffInfoCommand::extractSffInfo(string input, string accnos, string oligos){
	try {
        oligosObject = new Oligos();
		currentFileName = input;
		if (outputdir == "") {  outputdir += util.hasPath(input); }
		
        if (accnos != "")	{  seqNames.clear(); seqNames = util.readAccnos(accnos);    }
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
        string rootName = outputdir + util.getRootName(util.getSimpleName(input));
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
        
		SffCommonHeader* header = new SffCommonHeader();
        bool goodHeader = header->read(in);
		
		if (!goodHeader) { delete oligosObject; if (hasOligos)   { delete trimOligos; if (reorient) { delete  rtrimOligos; } } return 0; }
	
		//print common header
        if (sfftxt) {	header->printSFFTxt(outSfftxt);		        }
		if (flow)	{	outFlow << header->getNumFlows() << endl;	}
        
		//read through the sff file
        int count = 0; int numFlows = header->getNumFlows();
		while (!in.eof()) {
			
			bool print = true;

			SffRead* read = new SffRead(in, numFlows);
            
            if (!read->isOkay()) { break; }
            
            if (split > 1) { assignToSample(read, trimOligos, rtrimOligos); }
            
			//if you have provided an accosfile and this seq is not in it, then dont print
            if (seqNames.size() != 0) {   if (seqNames.count(read->getName()) == 0) { print = false; }  }
			
			//print 
			if (print) {
                if (sfftxt) {   read->printSffTxt(outSfftxt);        }
                if (fasta)	{	read->printFasta(outFasta, trim);	}
                if (qual)	{   read->printQuality(outQual, trim);	}
                if (flow)	{	read->printFlow(outFlow);	        }
			}
			count++;
            
            delete read;
            
			//report progress
			if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)+"\n"); 		}
		
			if (m->getControl_pressed()) { count = 0; break;   }
			
			if (count >= header->getNumReads()) { break; }
            
            //if (count >= 10000) { break; } //debug
		}
		
		//report progress
		if (!m->getControl_pressed()) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)+"\n"); 	}  }
		
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
                            if (numSplitReads[i][j] == 0) {
                                util.mothurRemove(filehandles[i][j]);
                                util.mothurRemove(filehandlesHeaders[i][j]);
                                namesToRemove.insert(filehandles[i][j]);
                            }else {
                                if(util.isBlank(filehandles[i][j])){
                                    
                                    util.mothurRemove(filehandles[i][j]);
                                    util.mothurRemove(filehandlesHeaders[i][j]);
                                    namesToRemove.insert(filehandles[i][j]);
                                }
                            }
						}
					}
				}
			}
            
            //append new header to reads
            for (int i = 0; i < filehandles.size(); i++) {
                for (int j = 0; j < filehandles[i].size(); j++) {
                    if (filehandles[i][j] != "") {
                        
                        if (namesToRemove.count(filehandles[i][j]) == 0) {
                            util.appendSFFFiles(filehandles[i][j], filehandlesHeaders[i][j]);
                            util.renameFile(filehandlesHeaders[i][j], filehandles[i][j]);
                            util.mothurRemove(filehandlesHeaders[i][j]);
                        }
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
        
        delete header;
        delete oligosObject;
        if (hasOligos)   { delete trimOligos; if (reorient) { delete  rtrimOligos; } }
        
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "extractSffInfo");
		exit(1);
	}
}
//****************************************************************************************
int SffInfoCommand::adjustCommonHeader(SffCommonHeader*& header){
	try {
        string endian = util.findEdianness();
        
        for (int i = 0; i < filehandlesHeaders.size(); i++) {
            for (int j = 0; j < filehandlesHeaders[i].size(); j++) {
                
                ofstream out; util.openOutputFileBinaryAppend(filehandlesHeaders[i][j], out);
                                
                header->printSampleCommonHeader(out, numSplitReads[i][j]);
                out.close();
            }
        }
        
        ofstream outNoMatchHeader;
        string tempNoHeader = "tempNoMatchHeader";
        util.openOutputFileBinary(tempNoHeader, outNoMatchHeader);
        
        header->printSampleCommonHeader(outNoMatchHeader, numNoMatch); outNoMatchHeader.close();
        
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
//***********************************************************************************************
void SffInfoCommand::assignToSample(SffRead*& read, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos){
    try {
        
        int barcodeIndex, primerIndex, trashCodeLength;
        
        if (hasOligos)      {  trashCodeLength = findGroup(read, barcodeIndex, primerIndex, trimOligos, rtrimOligos);                }
        else if (hasGroup)  {  trashCodeLength = findGroup(read, barcodeIndex, primerIndex, "groupMode");   }
        else {  m->mothurOut("[ERROR]: uh oh, we shouldn't be here...\n"); }
        
        if(trashCodeLength == 0){
            ofstream out; util.openOutputFileBinaryAppend(filehandles[barcodeIndex][primerIndex], out);
            read->printSff(out); out.close();
            numSplitReads[barcodeIndex][primerIndex]++;
        }
        else{
            ofstream out; util.openOutputFileBinaryAppend(noMatchFile, out);
            read->printSff(out); out.close();
            numNoMatch++;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SffInfoCommand", "assignToSample");
        exit(1);
    }
}
//***************************************************************************************************
int SffInfoCommand::findGroup(SffRead*& read, int& barcode, int& primer, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos) {
	try {
        
        int success = 1;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        string seq = read->getBases();
        int readLength = read->getBases().length();
        
        unsigned short clipLeft = read->getClipQualLeft();
        unsigned short clipRight = read->getClipQualRight();
        
        for (int i = 0; i < readLength; i++) {   seq[i] = toupper(seq[i]);                      }
        
        if (trim) {
            if(clipRight < clipLeft){
                //don't trim right
                if (clipRight == 0) { seq = seq.substr(clipLeft-1); }
                else { seq = "NNNN"; }
            }
            else if((clipRight != 0) && ((clipRight-clipLeft) >= 0)){  seq = seq.substr((clipLeft-1), (clipRight-clipLeft+1)); }
            else { seq = seq.substr(clipLeft-1); }
        }
        
        Sequence currSeq(read->getName(), seq);
        Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
        
        if(numLinkers != 0){
            success = trimOligos->stripLinker(currSeq);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numBarcodes != 0){
            vector<int> results = trimOligos->stripBarcode(currSeq, barcode);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numSpacers != 0){
            success = trimOligos->stripSpacer(currSeq);
            if(success > sdiffs)		{	trashCode += 's';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numFPrimers != 0){
            vector<int> results = trimOligos->stripForward(currSeq, primer);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numRPrimers != 0){
            vector<int> results = trimOligos->stripReverse(currSeq);
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
                vector<int> results = rtrimOligos->stripBarcode(savedSeq, thisBarcodeIndex);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                }
                if(thisSuccess > bdiffs)		{ thisTrashCode += "b"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if(numFPrimers != 0){
                vector<int> results = rtrimOligos->stripForward(savedSeq, thisPrimerIndex);
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
//***************************************************************************************
int SffInfoCommand::findGroup(SffRead*& read, int& barcode, int& primer, string groupMode) {
	try {
        string trashCode = "";
        primer = 0;
        
        string group = groupMap->getGroup(read->getName());
        if (group == "not found") {     trashCode += "g";   } //scrap for group
        else {  barcode = GroupToFile[group]; }
        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "findGroup");
		exit(1);
	}
}
//***********************************************************************************
int SffInfoCommand::parseSffTxt() {
	try {
		
		ifstream inSFF;
		util.openInputFile(sfftxtFilename, inSFF);
		
		if (outputdir == "") {  outputdir += util.hasPath(sfftxtFilename); }
		
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
			if (inSFF.eof()) { m->mothurOut("[ERROR]: Expected " + toString(numReads) + " but reached end of file at " + toString(i+1) + ".\n");  break; }
			
			SffRead read(numFlows);
			
			//parse read header
			inSFF >> seqName;
			seqName = seqName.substr(1);
			util.gobble(inSFF);
			read.setName(seqName);
			
			string runPrefix = parseHeaderLineToString(inSFF);		read.setTimeStamp(runPrefix);
			string regionNumber = parseHeaderLineToString(inSFF);	read.setRegion(regionNumber);
			string xyLocation = parseHeaderLineToString(inSFF);		read.setXY(xyLocation);
			util.gobble(inSFF);
				
			string runName = parseHeaderLineToString(inSFF);
			string analysisName = parseHeaderLineToString(inSFF);
			string fullPath = parseHeaderLineToString(inSFF);
			util.gobble(inSFF);
			
            unsigned short readHeaderLen = parseHeaderLineToShort(inSFF); read.setHeaderLength(readHeaderLen);
            unsigned short nameLength = parseHeaderLineToShort(inSFF); read.setNameLength(nameLength);
			int numBases = parseHeaderLineToInt(inSFF);	read.setNumBases(numBases);
            unsigned short clipQualLeft = parseHeaderLineToShort(inSFF); read.setClipQualLeft(clipQualLeft);
			unsigned short clipQualRight = parseHeaderLineToShort(inSFF); read.setClipQualRight(clipQualRight);
            unsigned short clipAdapLeft = parseHeaderLineToShort(inSFF); read.setClipAdapterLeft(clipAdapLeft);
            unsigned short clipAdapRight = parseHeaderLineToShort(inSFF); read.setClipAdapterRight(clipAdapRight);
			util.gobble(inSFF);
				
			//parse read
            vector<unsigned short> flowVector = parseHeaderLineToFloatVector(inSFF, numFlows);	read.setFlowgrams(flowVector);
			vector<unsigned int> flowIndices = parseHeaderLineToIntVector(inSFF, numBases);	
			
			//adjust for print
			vector<unsigned int> flowIndicesAdjusted; flowIndicesAdjusted.push_back(flowIndices[0]);
			for (int j = 1; j < flowIndices.size(); j++) {   flowIndicesAdjusted.push_back(flowIndices[j] - flowIndices[j-1]);   }
			read.setFlowIndex(flowIndicesAdjusted);
			
			string bases = parseHeaderLineToString(inSFF); read.setBases(bases);
			vector<unsigned int> qualityScores = parseHeaderLineToIntVector(inSFF, numBases); read.setQualScores(qualityScores);
			util.gobble(inSFF);
					
			//if you have provided an accosfile and this seq is not in it, then dont print
			bool print = true;
			if (seqNames.size() != 0) {   if (seqNames.count(read.getName()) == 0) { print = false; }  }
			
			//print 
			if (print) {
				if (fasta)	{	read.printFasta(outFasta, trim);	}
				if (qual)	{	read.printQuality(outQual, trim);	}
				if (flow)	{	read.printFlow(outFlow);	        }
			}
			
			//report progress
			if((i+1) % 10000 == 0){	m->mothurOut(toString(i+1)+"\n"); 		}
			
			if (m->getControl_pressed()) {  break;  }
		}
		
		//report progress
		if (!m->getControl_pressed()) {   if((numReads) % 10000 != 0){	m->mothurOut(toString(numReads)+"\n"); 		}  }
		
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
//*****************************************************************************
int SffInfoCommand::parseHeaderLineToInt(ifstream& file){
	try {
		int number;
		
		while (!file.eof())	{
			char c = file.get(); 
			if (c == ':'){ file >> number; break; }
		}
		util.gobble(file);
		return number;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToInt");
		exit(1);
	}
	
}
//*****************************************************************************
unsigned short SffInfoCommand::parseHeaderLineToShort(ifstream& file){
    try {
        string text;
        
        while (!file.eof())    {
            char c = file.get();
            
            if (c == ':'){ file >> text; break; }
        }
        util.gobble(file);
        
        unsigned short value;
        util.mothurConvert(text, value);
        return value;
    }
    catch(exception& e) {
        m->errorOut(e, "SffInfoCommand", "parseHeaderLineToShort");
        exit(1);
    }
}
//*****************************************************************************
string SffInfoCommand::parseHeaderLineToString(ifstream& file){
	try {
		string text;
		
		while (!file.eof())	{
			char c = file.get(); 
			
			if (c == ':'){ file >> text; break; }
		}
		util.gobble(file);
		
		return text;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "parseHeaderLineToString");
		exit(1);
	}
}
//********************************************************************************************
vector<unsigned short> SffInfoCommand::parseHeaderLineToFloatVector(ifstream& file, int length){
	try {
		vector<unsigned short> floatVector(length);
		
		while (!file.eof())	{
			char c = file.get(); 
			if (c == ':'){
				float temp;
				for(int i=0;i<length;i++){ file >> temp; floatVector[i] = temp * 100; }
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

//******************************************************************************************
vector<unsigned int> SffInfoCommand::parseHeaderLineToIntVector(ifstream& file, int length){
	try {
		vector<unsigned int> intVector(length);
		
		while (!file.eof())	{
			char c = file.get(); 
			if (c == ':'){
				for(int i=0;i<length;i++){ file >> intVector[i]; }
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
//***********************************************************************************************
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
                        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(currentFileName));
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
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(currentFileName));
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
//***************************************************************************************************
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
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(currentFileName));
                variables["[group]"] = groups[i];
                string thisFilename = getOutputFileName("sff",variables);
                outputNames.push_back(thisFilename);
               
                ofstream temp;
                util.openOutputFileBinary(thisFilename, temp); temp.close();
                filehandles[i].push_back(thisFilename);
                GroupToFile[groups[i]] = i;
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(currentFileName));
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

