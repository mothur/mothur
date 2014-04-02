/*
 *  parsefastaqcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "parsefastaqcommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> ParseFastaQCommand::setParameters(){	
	try {
		CommandParameter pfastq("fastq", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pfastq);
        CommandParameter poligos("oligos", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter pgroup("group", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter pfasta("fasta", "Boolean", "", "T", "", "", "","fasta",false,false); parameters.push_back(pfasta);
		CommandParameter pqual("qfile", "Boolean", "", "T", "", "", "","qfile",false,false); parameters.push_back(pqual);
        CommandParameter ppacbio("pacbio", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ppacbio);
 		CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "sanger", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The fastq.info command reads a fastq file and creates a fasta and quality file.\n";
		helpString += "The fastq.info command parameters are fastq, fasta, qfile, oligos, group and format; fastq is required.\n";
        helpString += "The fastq.info command should be in the following format: fastq.info(fastaq=yourFastaQFile).\n";
        helpString += "The oligos parameter allows you to provide an oligos file to split your fastq file into separate fastq files by barcode and primers. \n";
        helpString += "The group parameter allows you to provide a group file to split your fastq file into separate fastq files by group. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the reads. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=sanger.\n";
        helpString += "The fasta parameter allows you to indicate whether you want a fasta file generated. Default=T.\n";
        helpString += "The qfile parameter allows you to indicate whether you want a quality file generated. Default=T.\n";
        helpString += "The pacbio parameter allows you to indicate .... When set to true, quality scores of 0 will results in a corresponding base of N. Default=F.\n";
		helpString += "Example fastq.info(fastaq=test.fastaq).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fastq), '=' and yourFastQFile.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],fasta"; } 
        else if (type == "qfile") {  pattern = "[filename],qual"; }
        else if (type == "fastq") {  pattern = "[filename],[group],fastq"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ParseFastaQCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["fastq"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(string option){
	try {
		abort = false; calledHelp = false;
        split = 1;
		
		if(option == "help") {	help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;

			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["fastq"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
				}
                
                it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fastaQFile = validParameter.validFile(parameters, "fastq", true);
			if (fastaQFile == "not found") {	m->mothurOut("fastq is a required parameter for the fastq.info command.");	m->mothurOutEndLine();	abort = true;	}
			else if (fastaQFile == "not open")	{	fastaQFile = ""; abort = true;	}
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found") {	oligosfile = "";	}
			else if (oligosfile == "not open")	{	oligosfile = ""; abort = true;	}
            else { m->setOligosFile(oligosfile); split = 2; }
            
            groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not found") {	groupfile = "";	}
			else if (groupfile == "not open")	{	groupfile = ""; abort = true;	}
            else { m->setGroupFile(groupfile); split = 2; }
            
            if ((groupfile != "") && (oligosfile != "")) { m->mothurOut("You must enter ONLY ONE of the following: oligos or group."); m->mothurOutEndLine(); abort = true;  }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);	if (outputDir == "not found"){	outputDir = m->hasPath(fastaQFile); 	}
			
			string temp;
			temp = validParameter.validFile(parameters, "fasta", false);	if(temp == "not found"){	temp = "T";	}
			fasta = m->isTrue(temp); 

			temp = validParameter.validFile(parameters, "qfile", false);	if(temp == "not found"){	temp = "T";	}
			qual = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "pacbio", false);	if(temp == "not found"){	temp = "F";	}
			pacbio = m->isTrue(temp);

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

			
            format = validParameter.validFile(parameters, "format", false);		if (format == "not found"){	format = "sanger";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  { 
				m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
				abort=true;
			}

			if ((!fasta) && (!qual)) { m->mothurOut("[ERROR]: no outputs selected. Aborting."); m->mothurOutEndLine(); abort=true; }

		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ParseFastaQCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//open Output Files
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
		string fastaFile = getOutputFileName("fasta",variables);
		string qualFile = getOutputFileName("qfile",variables);
		ofstream outFasta, outQual;
		
		if (fasta) { m->openOutputFile(fastaFile, outFasta);  outputNames.push_back(fastaFile); outputTypes["fasta"].push_back(fastaFile);	}
		if (qual) { m->openOutputFile(qualFile, outQual);	outputNames.push_back(qualFile);  outputTypes["qfile"].push_back(qualFile);		}
        
        TrimOligos* trimOligos = NULL;
        int numBarcodes, numPrimers; numBarcodes = 0; numPrimers = 0;
        if (oligosfile != "")       {
            readOligos(oligosfile);
            numPrimers = primers.size(); numBarcodes = barcodes.size();
            //find group read belongs to
            if (pairedOligos)   {   trimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, pairedPrimers, pairedBarcodes); numBarcodes = pairedBarcodes.size(); numPrimers = pairedPrimers.size(); }
            else                {   trimOligos = new TrimOligos(pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimer, linker, spacer);  }

        }
        else if (groupfile != "")   { readGroup(groupfile);     }
		
		ifstream in;
		m->openInputFile(fastaQFile, in);
        
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) { 
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
		
        
        int count = 0;
		while (!in.eof()) {
			
			if (m->control_pressed) { break; }
            
            bool ignore;
            fastqRead2 thisRead = readFastq(in, ignore);
            
            if (!ignore) {
                vector<int> qualScores;
                if (qual) {
                    qualScores = convertQual(thisRead.quality);
                    outQual << ">" << thisRead.seq.getName() << endl;
                    for (int i = 0; i < qualScores.size(); i++) { outQual << qualScores[i] << " "; }
                    outQual << endl;
                }
                
                if (m->control_pressed) { break; }
                
                if (pacbio) {
                    if (!qual) { qualScores = convertQual(thisRead.quality); } //convert if not done
                    string sequence = thisRead.seq.getAligned();
                    for (int i = 0; i < qualScores.size(); i++) {
                        if (qualScores[i] == 0){ sequence[i] = 'N'; }
                    }
                    thisRead.seq.setAligned(sequence);
                }
                
                //print sequence info to files
                if (fasta) { thisRead.seq.printSequence(outFasta); }
                
                if (split > 1) {
                    int barcodeIndex, primerIndex, trashCodeLength;
                    if (oligosfile != "")      {  trashCodeLength = findGroup(thisRead, barcodeIndex, primerIndex, trimOligos, numBarcodes, numPrimers);    }
                    else if (groupfile != "")  {  trashCodeLength = findGroup(thisRead, barcodeIndex, primerIndex, "groupMode");   }
                    else {  m->mothurOut("[ERROR]: uh oh, we shouldn't be here...\n"); }
                    
                    if(trashCodeLength == 0){
                        ofstream out;
                        m->openOutputFileAppend(fastqFileNames[barcodeIndex][primerIndex], out);
                        out << thisRead.wholeRead;
                        out.close();
                    }else{
                        ofstream out;
                        m->openOutputFileAppend(noMatchFile, out);
                        out << thisRead.wholeRead;
                        out.close();
                    }
                }
                //report progress
                if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
                count++;
			}
		}
		
		in.close();
		if (fasta)	{ outFasta.close();	}
		if (qual)	{ outQual.close();	}
        
        //report progress
		if (!m->control_pressed) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
        
        if (split > 1) {
            
            if (groupfile != "")        { delete groupMap;      }
            else if (oligosfile != "")  { delete trimOligos;    }
           
			map<string, string>::iterator it;
			set<string> namesToRemove;
			for(int i=0;i<fastqFileNames.size();i++){
				for(int j=0;j<fastqFileNames[0].size();j++){
					if (fastqFileNames[i][j] != "") {
						if (namesToRemove.count(fastqFileNames[i][j]) == 0) {
							if(m->isBlank(fastqFileNames[i][j])){
								m->mothurRemove(fastqFileNames[i][j]);
								namesToRemove.insert(fastqFileNames[i][j]);
                            }
						}
					}
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
            else { outputNames.push_back(noMatchFile); outputTypes["fastq"].push_back(noMatchFile); }
        }
		
		if (m->control_pressed) { outputTypes.clear(); outputNames.clear(); m->mothurRemove(fastaFile); m->mothurRemove(qualFile); return 0; }
		
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
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
fastqRead2 ParseFastaQCommand::readFastq(ifstream& in, bool& ignore){
    try {
        ignore = false;
        string wholeRead = "";
        
        //read sequence name
        string line = m->getline(in); m->gobble(in); if (split > 1) { wholeRead += line + "\n"; }
        vector<string> pieces = m->splitWhiteSpace(line);
        string name = "";  if (pieces.size() != 0) { name = pieces[0]; }
        if (name == "") {  m->mothurOut("[WARNING]: Blank fasta name, ignoring read."); m->mothurOutEndLine(); ignore=true;  }
        else if (name[0] != '@') { m->mothurOut("[WARNING]: reading " + name + " expected a name with @ as a leading character, ignoring read."); m->mothurOutEndLine(); ignore=true; }
        else { name = name.substr(1); }
        
        //read sequence
        string sequence = m->getline(in); m->gobble(in); if (split > 1) { wholeRead += sequence + "\n"; }
        if (sequence == "") {  m->mothurOut("[WARNING]: missing sequence for " + name + ", ignoring."); ignore=true; }
        
        //read sequence name
        line = m->getline(in); m->gobble(in); if (split > 1) { wholeRead += line + "\n"; }
        pieces = m->splitWhiteSpace(line);
        string name2 = "";  if (pieces.size() != 0) { name2 = pieces[0]; }
        if (name2 == "") {  m->mothurOut("[WARNING]: expected a name with + as a leading character, ignoring."); ignore=true; }
        else if (name2[0] != '+') { m->mothurOut("[WARNING]: reading " + name2 + " expected a name with + as a leading character, ignoring."); ignore=true; }
        else { name2 = name2.substr(1); if (name2 == "") { name2 = name; } }
        
                
        //read quality scores
        string quality = m->getline(in); m->gobble(in); if (split > 1) { wholeRead += quality + "\n"; }
        if (quality == "") {  m->mothurOut("[WARNING]: missing quality for " + name2 + ", ignoring."); ignore=true; }
        
        //sanity check sequence length and number of quality scores match
        if (name2 != "") { if (name != name2) { m->mothurOut("[WARNING]: names do not match. read " + name + " for fasta and " + name2 + " for quality, ignoring."); ignore=true; } }
        if (quality.length() != sequence.length()) { m->mothurOut("[WARNING]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores, ignoring read."); ignore=true; }
        
        m->checkName(name);
        Sequence seq(name, sequence);
        fastqRead2 read(seq, quality, wholeRead);
            
        if (m->debug) { m->mothurOut("[DEBUG]: " + read.seq.getName() + " " + read.seq.getAligned() + " " + quality + "\n"); }
        
        return read;
    }
    catch(exception& e) {
        m->errorOut(e, "ParseFastaQCommand", "readFastq");
        exit(1);
    }
}

//**********************************************************************************************************************
vector<int> ParseFastaQCommand::convertQual(string qual) {
	try {
		vector<int> qualScores;
		
        bool negativeScores = false;
        
		for (int i = 0; i < qual.length(); i++) { 
            
            int temp = 0;
            temp = int(qual[i]);
            if (format == "illumina") {
                temp -= 64; //char '@'
            }else if (format == "illumina1.8+") {
                temp -= int('!'); //char '!'
            }else if (format == "solexa") {
                temp = int(convertTable[temp]); //convert to sanger
                temp -= int('!'); //char '!'
            }else {
                temp -= int('!'); //char '!'
            }
            if (temp < -5) { negativeScores = true; }
			qualScores.push_back(temp);
		}
		
        if (negativeScores) { m->mothurOut("[ERROR]: finding negative quality scores, do you have the right format selected? http://en.wikipedia.org/wiki/FASTQ_format#Encoding \n");  m->control_pressed = true;  }
        
		return qualScores;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "convertQual");
		exit(1);
	}
}
//**********************************************************************************************************************
int ParseFastaQCommand::findGroup(fastqRead2 thisRead, int& barcode, int& primer, TrimOligos*& trimOligos, int numBarcodes, int numPrimers) {
	try {
        int success = 1;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        Sequence currSeq(thisRead.seq.getName(), thisRead.seq.getAligned());
        QualityScores currQual; currQual.setScores(convertQual(thisRead.quality));
        
        if(linker.size() != 0){
            success = trimOligos->stripLinker(currSeq, currQual);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numBarcodes != 0){
            success = trimOligos->stripBarcode(currSeq, currQual, barcode);
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(spacer.size() != 0){
            success = trimOligos->stripSpacer(currSeq, currQual);
            if(success > sdiffs)		{	trashCode += 's';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numPrimers != 0){
            success = trimOligos->stripForward(currSeq, currQual, primer, true);
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
        
        if(revPrimer.size() != 0){
            success = trimOligos->stripReverse(currSeq, currQual);
            if(!success)				{	trashCode += 'r';	}
        }
        
        if (trashCode.length() == 0) { //is this sequence in the ignore group
            string thisGroup = "";
            
            if(barcodes.size() != 0){
                thisGroup = barcodeNameVector[barcode];
                if (numPrimers != 0) {
                    if (primerNameVector[primer] != "") {
                        if(thisGroup != "") {
                            thisGroup += "." + primerNameVector[primer];
                        }else {
                            thisGroup = primerNameVector[primer];
                        }
                    }
                }
            }
            
            int pos = thisGroup.find("ignore");
            if (pos != string::npos) {  trashCode += "i"; }
        }

        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int ParseFastaQCommand::findGroup(fastqRead2 thisRead, int& barcode, int& primer, string groupMode) {
	try {
        string trashCode = "";
        primer = 0;
        
        string group = groupMap->getGroup(thisRead.seq.getName());
        if (group == "not found") {     trashCode += "g";   } //scrap for group
        else { //find file group
            map<string, int>::iterator it = barcodes.find(group);
            if (it != barcodes.end()) {
                barcode = it->second;
            }else { trashCode += "g"; }
        }
        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//***************************************************************************************************************

bool ParseFastaQCommand::readOligos(string oligoFile){
	try {
		ifstream inOligos;
		m->openInputFile(oligoFile, inOligos);
		
		string type, oligo, roligo, group;
        bool hasPrimer = false; bool hasPairedBarcodes = false; pairedOligos = false;
        
		int indexPrimer = 0;
		int indexBarcode = 0;
        int indexPairedPrimer = 0;
		int indexPairedBarcode = 0;
        set<string> uniquePrimers;
        set<string> uniqueBarcodes;
		
		while(!inOligos.eof()){
            
			inOligos >> type;
            
		 	if (m->debug) { m->mothurOut("[DEBUG]: reading type - " + type + ".\n"); }
            
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				m->gobble(inOligos);
			}
			else{
				m->gobble(inOligos);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
                
                if (m->debug) { m->mothurOut("[DEBUG]: reading - " + oligo + ".\n"); }
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					group = "";
					
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					}
					
					//check for repeat barcodes
					map<string, int>::iterator itPrime = primers.find(oligo);
					if (itPrime != primers.end()) { m->mothurOut("primer " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
					
                    if (m->debug) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer " + oligo + ".\n"); }  }
                    
					primers[oligo]=indexPrimer; indexPrimer++;
					primerNameVector.push_back(group);
				}
                else if (type == "PRIMER"){
                    m->gobble(inOligos);
					
                    inOligos >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    roligo = reverseOligo(roligo);
                    
                    group = "";
                    
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					}
                    
                    oligosPair newPrimer(oligo, roligo);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: primer pair " + newPrimer.forward + " " + newPrimer.reverse + ", and group = " + group + ".\n"); }
					
					//check for repeat barcodes
                    string tempPair = oligo+roligo;
                    if (uniquePrimers.count(tempPair) != 0) { m->mothurOut("primer pair " + newPrimer.forward + " " + newPrimer.reverse + " is in your oligos file already."); m->mothurOutEndLine();  }
                    else { uniquePrimers.insert(tempPair); }
					
                    if (m->debug) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer pair " + newPrimer.forward + " " + newPrimer.reverse + ".\n"); }  }
                    
					pairedPrimers[indexPairedPrimer]=newPrimer; indexPairedPrimer++;
					primerNameVector.push_back(group);
                    hasPrimer = true;
                }
				else if(type == "REVERSE"){
					//Sequence oligoRC("reverse", oligo);
					//oligoRC.reverseComplement();
                    string oligoRC = reverseOligo(oligo);
					revPrimer.push_back(oligoRC);
				}
				else if(type == "BARCODE"){
					inOligos >> group;
                    
                    //barcode lines can look like   BARCODE   atgcatgc   groupName  - for 454 seqs
                    //or                            BARCODE   atgcatgc   atgcatgc    groupName  - for illumina data that has forward and reverse info
                    
                    string temp = "";
                    while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	temp += c;  }
					}
					
                    //then this is illumina data with 4 columns
                    if (temp != "") {
                        hasPairedBarcodes = true;
                        string reverseBarcode = group; //reverseOligo(group); //reverse barcode
                        group = temp;
                        
                        for(int i=0;i<reverseBarcode.length();i++){
                            reverseBarcode[i] = toupper(reverseBarcode[i]);
                            if(reverseBarcode[i] == 'U')	{	reverseBarcode[i] = 'T';	}
                        }
                        
                        reverseBarcode = reverseOligo(reverseBarcode);
                        oligosPair newPair(oligo, reverseBarcode);
                        
                        if (m->debug) { m->mothurOut("[DEBUG]: barcode pair " + newPair.forward + " " + newPair.reverse + ", and group = " + group + ".\n"); }
                        //check for repeat barcodes
                        string tempPair = oligo+reverseBarcode;
                        if (uniqueBarcodes.count(tempPair) != 0) { m->mothurOut("barcode pair " + newPair.forward + " " + newPair.reverse +  " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                        else { uniqueBarcodes.insert(tempPair); }
                        
                        pairedBarcodes[indexPairedBarcode]=newPair; indexPairedBarcode++;
                        barcodeNameVector.push_back(group);
                    }else {
                        //check for repeat barcodes
                        map<string, int>::iterator itBar = barcodes.find(oligo);
                        if (itBar != barcodes.end()) { m->mothurOut("barcode " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
                        
                        barcodes[oligo]=indexBarcode; indexBarcode++;
                        barcodeNameVector.push_back(group);
                    }
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
		
        if (hasPairedBarcodes || hasPrimer) {
            pairedOligos = true;
            if ((primers.size() != 0) || (barcodes.size() != 0) || (linker.size() != 0) || (spacer.size() != 0) || (revPrimer.size() != 0)) { m->control_pressed = true;  m->mothurOut("[ERROR]: cannot mix paired primers and barcodes with non paired or linkers and spacers, quitting."); m->mothurOutEndLine();  return 0; }
        }
        
		//add in potential combos
		if(barcodeNameVector.size() == 0){
			barcodes[""] = 0;
			barcodeNameVector.push_back("");
		}
		
		if(primerNameVector.size() == 0){
			primers[""] = 0;
			primerNameVector.push_back("");
		}
		
		fastqFileNames.resize(barcodeNameVector.size());
		for(int i=0;i<fastqFileNames.size();i++){
			fastqFileNames[i].assign(primerNameVector.size(), "");
		}
		
		
			set<string> uniqueNames; //used to cleanup outputFileNames
            if (pairedOligos) {
                for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                    for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                        
                        string primerName = primerNameVector[itPrimer->first];
                        string barcodeName = barcodeNameVector[itBar->first];
                        
                        if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                        else {
                            string comboGroupName = "";
                            string fastqFileName = "";
                            
                            if(primerName == ""){
                                comboGroupName = barcodeNameVector[itBar->first];
                            }
                            else{
                                if(barcodeName == ""){
                                    comboGroupName = primerNameVector[itPrimer->first];
                                }
                                else{
                                    comboGroupName = barcodeNameVector[itBar->first] + "." + primerNameVector[itPrimer->first];
                                }
                            }
                            
                            
                            ofstream temp;
                            map<string, string> variables;
                            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
                            variables["[group]"] = comboGroupName;
                            fastqFileName = getOutputFileName("fastq", variables);
                            if (uniqueNames.count(fastqFileName) == 0) {
                                outputNames.push_back(fastqFileName);
                                outputTypes["fastq"].push_back(fastqFileName);
                                uniqueNames.insert(fastqFileName);
                            }
                            
                            fastqFileNames[itBar->first][itPrimer->first] = fastqFileName;
                            m->openOutputFile(fastqFileName, temp);		temp.close();
                            
                        }
                    }
                }
            }else {
                for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                    for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                        
                        string primerName = primerNameVector[itPrimer->second];
                        string barcodeName = barcodeNameVector[itBar->second];
                        
                        if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                        else {
                            string comboGroupName = "";
                            string fastqFileName = "";
                            
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
                            map<string, string> variables;
                            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
                            variables["[group]"] = comboGroupName;
                            fastqFileName = getOutputFileName("fastq", variables);
                            if (uniqueNames.count(fastqFileName) == 0) {
                                outputNames.push_back(fastqFileName);
                                outputTypes["fastq"].push_back(fastqFileName);
                                uniqueNames.insert(fastqFileName);
                            }
                            
                            fastqFileNames[itBar->second][itPrimer->second] = fastqFileName;
                            m->openOutputFile(fastqFileName, temp);		temp.close();
                            
                        }
                    }
                }
            }
		
        ofstream temp;
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
        variables["[group]"] = "scrap";
        noMatchFile = getOutputFileName("fastq", variables);
        m->openOutputFile(noMatchFile, temp);		temp.close();

		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getOligos");
		exit(1);
	}
}
//***************************************************************************************************************
bool ParseFastaQCommand::readGroup(string groupfile){
	try {
        fastqFileNames.clear();
        
        groupMap = new GroupMap();
        groupMap->readMap(groupfile);
        
        //like barcodeNameVector - no primer names
        vector<string> groups = groupMap->getNamesOfGroups();
		
		fastqFileNames.resize(groups.size());
        for (int i = 0; i < fastqFileNames.size(); i++) {
            for (int j = 0; j < 1; j++) {
                
                map<string, string> variables;
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
                variables["[group]"] = groups[i];
                string thisFilename = getOutputFileName("fastq",variables);
                outputNames.push_back(thisFilename);
                outputTypes["fastq"].push_back(thisFilename);
                
                ofstream temp;
                m->openOutputFileBinary(thisFilename, temp); temp.close();
                fastqFileNames[i].push_back(thisFilename);
                barcodes[groups[i]] = i;
            }
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaQFile));
        variables["[group]"] = "scrap";
		noMatchFile = getOutputFileName("fastq",variables);
        m->mothurRemove(noMatchFile);
		
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "readGroup");
		exit(1);
	}
}
//********************************************************************/
string ParseFastaQCommand::reverseOligo(string oligo){
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
		m->errorOut(e, "ParseFastaQCommand", "reverseOligo");
		exit(1);
	}
}


//**********************************************************************************************************************



