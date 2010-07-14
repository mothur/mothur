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
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"sff","outputdir","inputdir", "outputdir"};
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
				splitAtDash(sffFilename, filenames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < filenames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(filenames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	filenames[i] = inputDir + filenames[i];		}
					}
	
					ifstream in;
					int ableToOpen = openInputFile(filenames[i], in);
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut(filenames[i] + " will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						filenames.erase(filenames.begin()+i);
						i--;
					}
				}
				
				//make sure there is at least one valid file left
				if (filenames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
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
		m->mothurOut("The sffinfo command reads a sff file and outputs a .sff.txt file.\n");
		
		m->mothurOut("Example sffinfo(sff=...).\n");
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
			
			m->mothurOut("Extracting info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
			
			if (outputDir == "") {  outputDir += hasPath(filenames[s]); }
			string outputFileName = outputDir + getRootName(getSimpleName(filenames[s])) + "sff.txt";
						
			extractSffInfo(filenames[s], outputFileName);
			
			outputNames.push_back(outputFileName);

			m->mothurOut("Done."); m->mothurOutEndLine();
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
int SffInfoCommand::extractSffInfo(string input, string output){
	try {
		ofstream out;
		openOutputFile(output, out);
		
		ifstream in;
		in.open(input.c_str(), ios::binary);
		
		CommonHeader* header = readCommonHeader(in);
		
		//cout << "magic = " << header->magicNumber << endl << "version = " << header->version << endl << "index offset = " << header->indexOffset << endl << "index length = "<< header->indexLength << endl << "numreads = " << header->numReads << endl << "header length = " << header->headerLength << endl << "key length = " << header->keyLength << endl;
//cout << "numflowreads = "<< header->numFlowsPerRead << endl << "flow format code = "<< header->flogramFormatCode << endl << "flow chars = " << header->flowChars << endl << "key sequence = " << header->keySequence << endl << endl;		
		cout << in.tellg() << endl;
		//read through the sff file
		while (!in.eof()) {
		
			//read header
			Header* readheader = readHeader(in);
			
			//read data
			seqRead* read = readSeqData(in, header->numFlowsPerRead, readheader->numBases);
			
			cout << in.tellg() << endl;
			
			//print common header
			printCommonHeader(out, header, true);
			
			//print header
			printHeader(out, readheader, true);
			
			//print data
			printSeqData(out, read, true);
	
		}
		
		
		in.close();
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "extractSffInfo");
		exit(1);
	}
}
//**********************************************************************************************************************
CommonHeader* SffInfoCommand::readCommonHeader(ifstream& in){
	try {
		CommonHeader* header = new CommonHeader();
	
		if (!in.eof()) {
			string tempBuf = "";
			
			//read magic number
			char* buffer = new char(sizeof(header->magicNumber));
			in.read(buffer, sizeof(header->magicNumber));
			header->magicNumber = be_int4(*(uint32_t *)(buffer));
			delete[] buffer;
			
			//read version
			header->version = new char(4);
			in.read(header->version, 4);
			tempBuf = buffer;
			if (tempBuf.length() > 4) { tempBuf = tempBuf.substr(0, 4); strcpy(header->version, tempBuf.c_str());  }
			
			//read offset
			buffer = new char(sizeof(header->indexOffset));
			in.read(buffer, sizeof(header->indexOffset));
			header->indexOffset =  be_int8(*(uint64_t *)(buffer));
			delete[] buffer;
			
			//read index length
			buffer = new char(sizeof(header->indexLength));
			in.read(buffer, sizeof(header->indexLength));
			header->indexLength =  be_int4(*(uint32_t *)(buffer));
			delete[] buffer;
			
			//read num reads
			buffer = new char(sizeof(header->numReads));
			in.read(buffer, sizeof(header->numReads));
			header->numReads =  be_int4(*(uint32_t *)(buffer));
			delete[] buffer;	
			
			//read header length
			buffer = new char(sizeof(header->headerLength));
			in.read(buffer, sizeof(header->headerLength));
			header->headerLength =  be_int2(*(uint16_t *)(buffer));
			delete[] buffer;
			
			//read key length
			buffer = new char(sizeof(header->keyLength));
			in.read(buffer, sizeof(header->keyLength));
			header->keyLength = be_int2(*(uint16_t *)(buffer));
			delete[] buffer;

			//read number of flow reads
			buffer = new char(sizeof(header->numFlowsPerRead));
			in.read(buffer, sizeof(header->numFlowsPerRead));
			header->numFlowsPerRead =  be_int2(*(uint16_t *)(buffer));
			delete[] buffer;
			
			//read format code
			buffer = new char(sizeof(header->flogramFormatCode));
			in.read(buffer, sizeof(header->flogramFormatCode));
			header->flogramFormatCode = be_int1(*(uint8_t *)(buffer));
			delete[] buffer;
			
			//read flow chars
			header->flowChars = new char(header->numFlowsPerRead);
			in.read(header->flowChars, header->numFlowsPerRead);
			tempBuf = buffer;
			if (tempBuf.length() > header->numFlowsPerRead) { tempBuf = tempBuf.substr(0, header->numFlowsPerRead); strcpy(header->flowChars, tempBuf.c_str());  }
			
			//read key
			header->keySequence = new char(header->keyLength);
			in.read(header->keySequence, header->keyLength);
			tempBuf = header->keySequence;
			if (tempBuf.length() > header->keyLength) { tempBuf = tempBuf.substr(0, header->keyLength); strcpy(header->keySequence, tempBuf.c_str());  }
			
		}else{
			m->mothurOut("Error reading sff common header."); m->mothurOutEndLine();
		}

		return header;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
Header* SffInfoCommand::readHeader(ifstream& in){
	try {
		Header* header = new Header();
	
		if (!in.eof()) {
			string tempBuf = "";
			
			//read header length
			char* buffer = new char(sizeof(header->headerLength));
			in.read(buffer, sizeof(header->headerLength));
			header->headerLength = be_int2(*(unsigned short *)(buffer));
			delete[] buffer;
						
			//read name length
			buffer = new char(sizeof(header->nameLength));
			in.read(buffer, sizeof(header->nameLength));
			header->nameLength = be_int2(*(unsigned short *)(buffer));
			delete[] buffer;

			//read num bases
			buffer = new char(sizeof(header->numBases));
			in.read(buffer, sizeof(header->numBases));
			header->numBases =  be_int4(*(unsigned int *)(buffer));
			delete[] buffer;
			
			//read clip qual left
			buffer = new char(sizeof(header->clipQualLeft));
			in.read(buffer, sizeof(header->clipQualLeft));
			header->clipQualLeft =  be_int2(*(unsigned short *)(buffer));
			delete[] buffer;
			
			//read clip qual right
			buffer = new char(sizeof(header->clipQualRight));
			in.read(buffer, sizeof(header->clipQualRight));
			header->clipQualRight =  be_int2(*(unsigned short *)(buffer));
			delete[] buffer;
			
			//read clipAdapterLeft
			buffer = new char(sizeof(header->clipAdapterLeft));
			in.read(buffer, sizeof(header->clipAdapterLeft));
			header->clipAdapterLeft = be_int2(*(unsigned short *)(buffer));
			delete[] buffer;

			//read clipAdapterRight
			buffer = new char(sizeof(header->clipAdapterRight));
			in.read(buffer, sizeof(header->clipAdapterRight));
			header->clipAdapterRight =  be_int2(*(unsigned short *)(buffer));
			delete[] buffer;
			
			//read name
			header->name = new char(header->nameLength);
			in.read(header->name, header->nameLength);
			tempBuf = header->name;
			if (tempBuf.length() > header->nameLength) { tempBuf = tempBuf.substr(0, header->nameLength); strcpy(header->name, tempBuf.c_str());  }
	
		}else{
			m->mothurOut("Error reading sff header info."); m->mothurOutEndLine();
		}

		return header;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
seqRead* SffInfoCommand::readSeqData(ifstream& in, int numFlowReads, int numBases){
	try {
		seqRead* read = new seqRead();
	
		if (!in.eof()) {
						
			string tempBuf = "";
			char* buffer;
			
			//read flowgram
			read->flowgram.resize(numFlowReads);
			for (int i = 0; i < numFlowReads; i++) {  
				buffer = new char((sizeof(unsigned short)));
				in.read(buffer, (sizeof(unsigned short)));
				read->flowgram[i] = be_int2(*(unsigned short *)(buffer));
				delete[] buffer;
			}
			
			//read flowgram
			read->flowIndex.resize(numBases);
			for (int i = 0; i < numBases; i++) {  
				buffer = new char(1);
				in.read(buffer, 1);
				read->flowgram[i] = be_int1(*(unsigned int *)(buffer));
				delete[] buffer;
			}
			
			//read bases
			read->bases = new char(numBases);
			in.read(read->bases, numBases);
			tempBuf = buffer;
			if (tempBuf.length() > numBases) { tempBuf = tempBuf.substr(0, numBases); strcpy(read->bases, tempBuf.c_str());  }
			

			//read flowgram
			read->qualScores.resize(numBases);
			for (int i = 0; i < numBases; i++) {  
				buffer = new char(1);
				in.read(buffer, 1);
				read->qualScores[i] = be_int1(*(unsigned int *)(buffer));
				delete[] buffer;
			}
			
			/* Pad to 8 chars */
			int spotInFile = in.tellg();
	cout << spotInFile << endl;
			int spot = floor((spotInFile + 7) /(float) 8) * 8;
	cout << spot << endl;
			in.seekg(spot);
			
		}else{
			m->mothurOut("Error reading."); m->mothurOutEndLine();
		}

		return read;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "readSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printCommonHeader(ofstream& out, CommonHeader* header, bool debug) {
	try {
	
		string output = "Common Header:\nMagic Number: ";
		output += toString(header->magicNumber) + '\n';
		output += "Version: " + toString(header->version) + '\n';
		output += "Index Offset: " + toString(header->indexOffset) + '\n';
		output += "Index Length: " + toString(header->indexLength) + '\n';	
		output += "Number of Reads: " + toString(header->numReads) + '\n';
		output += "Header Length: " + toString(header->headerLength) + '\n';
		output += "Key Length: " + toString(header->keyLength) + '\n';
		output += "Number of Flows: " + toString(header->numFlowsPerRead) + '\n';
		output += "Format Code: " + toString(header->flogramFormatCode) + '\n';
		output += "Flow Chars: " + toString(header->flowChars) + '\n';
		output += "Key Sequence: " + toString(header->keySequence) + '\n';
		
		out << output << endl;
		
		if (debug) { cout << output << endl; }
			
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printCommonHeader");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffInfoCommand::printHeader(ofstream& out, Header* header, bool debug) {
	try {
		string name = header->name;
		string output = ">" + name + "\nRead Header Length: " + toString(header->headerLength) + '\n';
		output += "Name Length: " + toString(header->nameLength) + '\n';
		output += "Number of Bases: " + toString(header->numBases) + '\n';
		output += "Clip Qual Left: " + toString(header->clipQualLeft) + '\n';
		output += "Clip Qual Right: " + toString(header->clipQualLeft) + '\n';
		output += "Clip Adap Left: " + toString(header->clipQualLeft) + '\n';
		output += "Clip Adap Right: " + toString(header->clipQualLeft) + '\n';
		
		out << output << endl;
		
		if (debug) { cout << output << endl; }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printHeader");
		exit(1);
	}
}

//**********************************************************************************************************************
int SffInfoCommand::printSeqData(ofstream& out, seqRead* read, bool debug) {
	try {

		string output = "FlowGram: ";
		for (int i = 0; i < read->flowgram.size(); i++) {  output += toString(read->flowgram[i]) +'\t';  }
		output += "\nFlow Indexes: ";
		for (int i = 0; i < read->flowIndex.size(); i++) {  output += toString(read->flowIndex[i]) +'\t';  }
		string bases = read->bases;
		output += "\nBases: " + bases + '\n';
		for (int i = 0; i < read->qualScores.size(); i++) {  output += toString(read->qualScores[i]) +'\t';  }
		output += '\n';
		
		out << output << endl;
		
		if (debug) { cout << output << endl; }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************/
