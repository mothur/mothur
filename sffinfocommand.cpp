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

//**********************************************************************************************************************

SffInfoCommand::SffInfoCommand(string option)  {
	try {
		abort = false;
		hasAccnos = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"sff","qfile","fasta","flow","trim","accnos","sfftxt","outputdir","inputdir", "outputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);	  if (inputDir == "not found"){	inputDir = "";		}

			sffFilename = validParameter.validFile(parameters, "sff", false);
			if (sffFilename == "not found") { m->mothurOut("sff is a required parameter for the sffinfo command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(sffFilename, filenames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < filenames.size(); i++) {
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
							ableToOpen = m->openInputFile(tryPath, in, "noerror");
							filenames[i] = tryPath;
						}
					}
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + filenames[i] + ". It will be disregarded."); m->mothurOutEndLine();
						//erase from file list
						filenames.erase(filenames.begin()+i);
						i--;
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
							ableToOpen = m->openInputFile(tryPath, in, "noerror");
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
				
				//make sure there is at least one valid file left
				if (accnosFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasAccnos) {
				if (accnosFileNames.size() != filenames.size()) { abort = true; m->mothurOut("If you provide a accnos file, you must have one for each sff file."); m->mothurOutEndLine(); }
			}
			
			string temp = validParameter.validFile(parameters, "qfile", false);			if (temp == "not found"){	temp = "T";				}
			qual = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "fasta", false);				if (temp == "not found"){	temp = "T";				}
			fasta = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "flow", false);					if (temp == "not found"){	temp = "F";				}
			flow = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "trim", false);					if (temp == "not found"){	temp = "T";				}
			trim = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "sfftxt", false);				if (temp == "not found"){	temp = "F";				}
			sfftxt = m->isTrue(temp); 
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "SffInfoCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SffInfoCommand::help(){
	try {
		m->mothurOut("The sffinfo command reads a sff file and extracts the sequence data.\n");
		m->mothurOut("The sffinfo command parameters are sff, fasta, qfile, accnos, flow, sfftxt, and trim. sff is required. \n");
		m->mothurOut("The sff parameter allows you to enter the sff file you would like to extract data from.  You may enter multiple files by separating them by -'s.\n");
		m->mothurOut("The fasta parameter allows you to indicate if you would like a fasta formatted file generated.  Default=True. \n");
		m->mothurOut("The qfile parameter allows you to indicate if you would like a quality file generated.  Default=True. \n");
		m->mothurOut("The flow parameter allows you to indicate if you would like a flowgram file generated.  Default=False. \n");
		m->mothurOut("The sfftxt parameter allows you to indicate if you would like a sff.txt file generated.  Default=False. \n");
		m->mothurOut("The trim parameter allows you to indicate if you would like a sequences and quality scores trimmed to the clipQualLeft and clipQualRight values.  Default=True. \n");
		m->mothurOut("The accnos parameter allows you to provide a accnos file containing the names of the sequences you would like extracted. You may enter multiple files by separating them by -'s. \n");
		m->mothurOut("Example sffinfo(sff=mySffFile.sff, trim=F).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. sff), '=' and parameters (i.e.yourSffFileName).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************

SffInfoCommand::~SffInfoCommand(){}

//**********************************************************************************************************************
int SffInfoCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		for (int s = 0; s < filenames.size(); s++) {
			
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }
			
			int start = time(NULL);
			
			m->mothurOut("Extracting info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
			
			string accnos = "";
			if (hasAccnos) { accnos = accnosFileNames[s]; }
			
			int numReads = extractSffInfo(filenames[s], accnos);

			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to extract " + toString(numReads) + ".");
		}
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }
		
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
int SffInfoCommand::extractSffInfo(string input, string accnos){
	try {
		
		if (outputDir == "") {  outputDir += m->hasPath(input); }
		
		if (accnos != "")	{  readAccnosFile(accnos);  }
		else				{	seqNames.clear();		}

		ofstream outSfftxt, outFasta, outQual, outFlow;
		string outFastaFileName, outQualFileName;
		string sfftxtFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "sff.txt";
		string outFlowFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "flow";
		if (trim) {
			outFastaFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "fasta";
			outQualFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "qual";
		}else{
			outFastaFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "raw.fasta";
			outQualFileName = outputDir + m->getRootName(m->getSimpleName(input)) + "raw.qual";
		}
		
		if (sfftxt) { m->openOutputFile(sfftxtFileName, outSfftxt); outSfftxt.setf(ios::fixed, ios::floatfield); outSfftxt.setf(ios::showpoint);  outputNames.push_back(sfftxtFileName); }
		if (fasta)	{ m->openOutputFile(outFastaFileName, outFasta);	outputNames.push_back(outFastaFileName); }
		if (qual)	{ m->openOutputFile(outQualFileName, outQual);		outputNames.push_back(outQualFileName);  }
		if (flow)	{ m->openOutputFile(outFlowFileName, outFlow);		outputNames.push_back(outFlowFileName);  }
		
		ifstream in;
		in.open(input.c_str(), ios::binary);
		
		CommonHeader header; 
		readCommonHeader(in, header);
		
		int count = 0;
		
		//check magic number and version
		if (header.magicNumber != 779314790) { m->mothurOut("Magic Number is not correct, not a valid .sff file"); m->mothurOutEndLine(); return count; }
		if (header.version != "0001") { m->mothurOut("Version is not supported, only support version 0001."); m->mothurOutEndLine(); return count; }
	
		//print common header
		if (sfftxt) { printCommonHeader(outSfftxt, header); }
	
		//read through the sff file
		while (!in.eof()) {
			
			bool print = true;
			
			//read header
			Header readheader;
			readHeader(in, readheader);
			
			//read data
			seqRead read; 
			readSeqData(in, read, header.numFlowsPerRead, readheader.numBases);
				
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
			for (int i = 0; i < 4; i++) {  header.version += toString((int)(buffer9[i])); }
				
			//read offset
			char buffer2 [8];
			in.read(buffer2, 8);
			header.indexOffset =  be_int8(*(unsigned long int *)(&buffer2));
			
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
			unsigned long int spotInFile = in.tellg();
			unsigned long int spot = (spotInFile + 7)& ~7;  // ~ inverts
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
int SffInfoCommand::readHeader(ifstream& in, Header& header){
	try {
	
		if (!in.eof()) {
			
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
			unsigned long int spotInFile = in.tellg();
			unsigned long int spot = (spotInFile + 7)& ~7;
			in.seekg(spot);
			
		}else{
			m->mothurOut("Error reading sff header info."); m->mothurOutEndLine();
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::readSeqData(ifstream& in, seqRead& read, int numFlowReads, int numBases){
	try {
	
		if (!in.eof()) {
	
			//read flowgram
			read.flowgram.resize(numFlowReads);
			for (int i = 0; i < numFlowReads; i++) {  
				char buffer [2];
				in.read(buffer, 2);
				read.flowgram[i] = be_int2(*(unsigned short *)(&buffer));
			}
	
			//read flowIndex
			read.flowIndex.resize(numBases);
			for (int i = 0; i < numBases; i++) {  
				char temp[1];
				in.read(temp, 1);
				read.flowIndex[i] = be_int1(*(unsigned char *)(&temp));
			}
	
			//read bases
			char* tempBuffer = new char[numBases];
			in.read(&(*tempBuffer), numBases);
			read.bases = tempBuffer;
			if (read.bases.length() > numBases) { read.bases = read.bases.substr(0, numBases);  }
			delete[] tempBuffer;

			//read qual scores
			read.qualScores.resize(numBases);
			for (int i = 0; i < numBases; i++) {  
				char temp[1];
				in.read(temp, 1);
				read.qualScores[i] = be_int1(*(unsigned char *)(&temp));
			}
	
			/* Pad to 8 chars */
			unsigned long int spotInFile = in.tellg();
			unsigned long int spot = (spotInFile + 7)& ~7;
			in.seekg(spot);
			
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
		out << "Run Prefix: " << endl;
		out << "Region #:  " << endl;
		out << "XY Location: " << endl << endl;
		
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
		
		out << ">" << header.name << endl;
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
				out << "0\t0\t0\t0";
			}
			else if((header.clipQualRight != 0) && ((header.clipQualRight-header.clipQualLeft) >= 0)){
				out << ">" << header.name << " length=" << (header.clipQualRight-header.clipQualLeft) << endl;
				for (int i = (header.clipQualLeft-1); i < (header.clipQualRight-1); i++) {   out << read.qualScores[i] << '\t';	}
			}
			else{
				out << ">" << header.name << " length=" << (header.clipQualRight-header.clipQualLeft) << endl;
				for (int i = (header.clipQualLeft-1); i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';	}			
			}
		}else{
			out << ">" << header.name << " length=" << read.qualScores.size() << endl;
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
		
		out << ">" << header.name << endl;
		for (int i = 0; i < read.flowgram.size(); i++) { out << setprecision(2) << (read.flowgram[i]/(float)100) << '\t';  }
		out << endl;
		
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
//**********************************************************************************************************************/
