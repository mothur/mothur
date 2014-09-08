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
        CommandParameter pfile("file", "InputTypes", "", "", "fastqFile", "fastqFile", "none","",false,false,true); parameters.push_back(pfile);
		CommandParameter pfastq("fastq", "InputTypes", "", "", "fastqFile", "fastqFile", "none","",false,false,true); parameters.push_back(pfastq);
        CommandParameter poligos("oligos", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter pgroup("group", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
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
		helpString += "The fastq.info command parameters are file, fastq, fasta, qfile, oligos, group and format; file or fastq is required.\n";
        helpString += "The fastq.info command should be in the following format: fastq.info(fastaq=yourFastaQFile).\n";
        helpString += "The oligos parameter allows you to provide an oligos file to split your fastq file into separate fastq files by barcode and primers. \n";
        helpString += "The group parameter allows you to provide a group file to split your fastq file into separate fastq files by group. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the reads. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. If found the sequence is flipped. The default is false.\n";
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
			inputDir = validParameter.validFile(parameters, "inputdir", false);		
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
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
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
			if (fastaQFile == "not found") {	fastaQFile= "";	}
			else if (fastaQFile == "not open")	{	fastaQFile = ""; abort = true;	}
            else { inputfile = fastaQFile; }
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not found") {	file = "";	}
			else if (file == "not open")	{	file = ""; abort = true;	}
            else { inputfile = file; }
            
            if ((file == "") && (fastaQFile == "")) {  m->mothurOut("You must provide a file or fastq option."); m->mothurOutEndLine(); abort = true;  }

            
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
            
            temp = validParameter.validFile(parameters, "checkorient", false);		if (temp == "not found") { temp = "F"; }
			reorient = m->isTrue(temp);

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
		
        if (file != "") {
            //read file
            vector< vector<string> > files = readFile();
            
            if (m->control_pressed) { return 0; }
            
            for (int i = 0; i < files.size(); i++) { //process each pair
                
                if (m->control_pressed) { break; }
                
            }
           
            //if oligos, then parse, otherwise just create fasta and qual for each.
            
            //process forward and reverse
            
            //parse each set of file with oligos file
            
            
            
        }else {
            
            //open Output Files
            map<string, string> variables;
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
            string fastaFile = getOutputFileName("fasta",variables);
            string qualFile = getOutputFileName("qfile",variables);
            ofstream outFasta, outQual;
            
            if (fasta) { m->openOutputFile(fastaFile, outFasta);  outputNames.push_back(fastaFile); outputTypes["fasta"].push_back(fastaFile);	}
            if (qual) { m->openOutputFile(qualFile, outQual);	outputNames.push_back(qualFile);  outputTypes["qfile"].push_back(qualFile);		}
            
            TrimOligos* trimOligos = NULL; TrimOligos* rtrimOligos = NULL;
            pairedOligos = false; numBarcodes = 0; numPrimers= 0; numLinkers= 0; numSpacers = 0; numRPrimers = 0;
            if (oligosfile != "")       {
                readOligos(oligosfile);
                //find group read belongs to
                if (pairedOligos)   {   trimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos.getPairedPrimers(), oligos.getPairedBarcodes()); numBarcodes = oligos.getPairedBarcodes().size(); numPrimers = oligos.getPairedPrimers().size(); }
                else                {   trimOligos = new TrimOligos(pdiffs, bdiffs, ldiffs, sdiffs, oligos.getPrimers(), oligos.getBarcodes(), oligos.getReversePrimers(), oligos.getLinkers(), oligos.getSpacers());  numPrimers = oligos.getPrimers().size(); numBarcodes = oligos.getBarcodes().size();  }
                
                if (reorient) {
                    rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos.getReorientedPairedPrimers(), oligos.getReorientedPairedBarcodes()); numBarcodes = oligos.getReorientedPairedBarcodes().size();
                }
                
            }else if (groupfile != "")   { readGroup(groupfile);     }

            
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
                        if (oligosfile != "")      {  trashCodeLength = findGroup(thisRead, barcodeIndex, primerIndex, trimOligos, rtrimOligos, numBarcodes, numPrimers);    }
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
                else if (oligosfile != "")  { delete trimOligos; if (reorient) { delete rtrimOligos; }   }
                
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
		}
		
		
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
int ParseFastaQCommand::findGroup(fastqRead2 thisRead, int& barcode, int& primer, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos, int numBarcodes, int numPrimers) {
	try {
        int success = 1;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        Sequence currSeq(thisRead.seq.getName(), thisRead.seq.getAligned());
        QualityScores currQual; currQual.setScores(convertQual(thisRead.quality));
        
        //for reorient
        Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
        QualityScores savedQual(currQual.getName(), currQual.getScores());
        
        if(numLinkers != 0){
            success = trimOligos->stripLinker(currSeq, currQual);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numBarcodes != 0){
            success = trimOligos->stripBarcode(currSeq, currQual, barcode);
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numSpacers != 0){
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
        
        if(numRPrimers != 0){
            success = trimOligos->stripReverse(currSeq, currQual);
            if(!success)				{	trashCode += 'r';	}
        }
        
        if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
            int thisSuccess = 0;
            string thisTrashCode = "";
            int thisCurrentSeqsDiffs = 0;
            
            int thisBarcodeIndex = 0;
            int thisPrimerIndex = 0;
            //cout << currSeq.getName() << '\t' << savedSeq.getUnaligned() << endl;
            if(numBarcodes != 0){
                thisSuccess = rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                if(thisSuccess > bdiffs)		{ thisTrashCode += "b"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            //cout << currSeq.getName() << '\t' << savedSeq.getUnaligned() << endl;
            if(numPrimers != 0){
                thisSuccess = rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, true);
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
            string thisGroup = oligos.getGroupName(barcode, primer);
            
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
                
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************

/*
 
 file option 1
 
 ffastqfile1 rfastqfile1
 ffastqfile2 rfastqfile2
 ...
 
 file option 2
 
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 ...
 
 */

vector< vector<string> > ParseFastaQCommand::readFile(){
	try {
        
        vector< vector<string> > files;
        
        ifstream in;
        m->openInputFile(file, in);
        
        while(!in.eof()) {
            
            if (m->control_pressed) { return files; }
            
            string line = m->getline(in);  m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            string group = "";
            string thisFileName1, thisFileName2; thisFileName1 = ""; thisFileName2 = "";
            if (pieces.size() == 2) {
                thisFileName1 = pieces[0];
                thisFileName2 = pieces[1];
            }else if (pieces.size() == 3) {
                thisFileName1 = pieces[1];
                thisFileName2 = pieces[2];
                string group = pieces[0];
            }else {
                m->mothurOut("[ERROR]: file lines can be 2 or 3 columns. The 2 column files are forwardFastqFile and reverseFastqFile. The 3 column are groupName, forwardFastqFile reverseFastqFile. You may have multiple lines in the file. \n"); m->control_pressed = true;
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + ", thisFileName1 = " + thisFileName1 + ", thisFileName2 = " + thisFileName2  + ".\n"); }
            
            if (inputDir != "") {
                string path = m->hasPath(thisFileName1);
                if (path == "") {  thisFileName1 = inputDir + thisFileName1;  }
                
                path = m->hasPath(thisFileName2);
                if (path == "") {  thisFileName2 = inputDir + thisFileName2;  }
            }
            
            //check to make sure both are able to be opened
            ifstream in2;
            int openForward = m->openInputFile(thisFileName1, in2, "noerror");
            
            //if you can't open it, try default location
            if (openForward == 1) {
                
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openForward = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    thisFileName1 = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openForward == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = m->openInputFile(tryPath, in4, "noerror");
                    thisFileName1 = tryPath;
                    in4.close();
                }
            }
            
            if (openForward == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + thisFileName1 + ", ignoring.\n");
            }else{  in2.close();  }
            
            int openReverse = 1;
            
            ifstream in3;
            openReverse = m->openInputFile(thisFileName2, in3, "noerror");
            
            //if you can't open it, try default location
            if (openReverse == 1) {
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName2);
                    m->mothurOut("Unable to open " + thisFileName2 + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openReverse = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    thisFileName2 = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openReverse == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName2);
                    m->mothurOut("Unable to open " + thisFileName2 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openReverse = m->openInputFile(tryPath, in4, "noerror");
                    thisFileName2 = tryPath;
                    in4.close();
                }
            }
            
            if (openReverse == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + thisFileName2 + ", ignoring pair.\n");
            }else{  in3.close();  }
            
            
            if ((pieces.size() == 2) && (openForward != 1) && (openReverse != 1)) { //good pair a
                vector<string> temp;
                temp.push_back(thisFileName1); temp.push_back(thisFileName2);
                files.push_back(temp);
            }else if((pieces.size() == 3) && (openForward != 1) && (openReverse != 1)) { //good pair
                vector<string> temp;
                temp.push_back(thisFileName1); temp.push_back(thisFileName2);
                files.push_back(temp);
                Groups.push_back(group);
            }
        }
        in.close();
        
        return files;
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "readFile");
		exit(1);
	}
}

//***************************************************************************************************************

bool ParseFastaQCommand::readOligos(string oligoFile){
	try {
        bool allBlank = false;
        oligos.read(oligosfile);
        
        if (m->control_pressed) { return false; } //error in reading oligos
        
        if (oligos.hasPairedBarcodes()) {
            pairedOligos = true;
            numPrimers = oligos.getPairedPrimers().size();
            numBarcodes = oligos.getPairedBarcodes().size();
        }else {
            pairedOligos = false;
            numPrimers = oligos.getPrimers().size();
            numBarcodes = oligos.getBarcodes().size();
        }
        
        numLinkers = oligos.getLinkers().size();
        numSpacers = oligos.getSpacers().size();
        numRPrimers = oligos.getReversePrimers().size();
        
        vector<string> groupNames = oligos.getGroupNames();
        if (groupNames.size() == 0) { allBlank = true;  }
        
        
        fastqFileNames.resize(oligos.getBarcodeNames().size());
		for(int i=0;i<fastqFileNames.size();i++){
            for(int j=0;j<oligos.getPrimerNames().size();j++){  fastqFileNames[i].push_back(""); }
		}
        
        set<string> uniqueNames; //used to cleanup outputFileNames
        if (pairedOligos) {
            map<int, oligosPair> barcodes = oligos.getPairedBarcodes();
            map<int, oligosPair> primers = oligos.getPairedPrimers();
            for(map<int, oligosPair>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->first);
                    string barcodeName = oligos.getBarcodeName(itBar->first);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        
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
                        
                        ofstream temp;
                        map<string, string> variables;
                        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
                        variables["[group]"] = comboGroupName;
                        string fastqFileName = getOutputFileName("fastq", variables);
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
            map<string, int> barcodes = oligos.getBarcodes() ;
            map<string, int> primers = oligos.getPrimers();
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->second);
                    string barcodeName = oligos.getBarcodeName(itBar->second);
                   
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        
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
                        
                        ofstream temp;
                        map<string, string> variables;
                        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
                        variables["[group]"] = comboGroupName;
                        string fastqFileName = getOutputFileName("fastq", variables);
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
        
        if (allBlank) {
            m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a groupfile."); m->mothurOutEndLine();
            return false;
        }
       
        ofstream temp;
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
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
//**********************************************************************************************************************



