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
			
			int start = time(NULL);
			
			m->mothurOut("Extracting info from " + filenames[s] + " ..." ); m->mothurOutEndLine();
			
			if (outputDir == "") {  outputDir += hasPath(filenames[s]); }
			string outputFileName = outputDir + getRootName(getSimpleName(filenames[s])) + "sff.txt";
						
			int numReads = extractSffInfo(filenames[s], outputFileName);
			
			outputNames.push_back(outputFileName);

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
int SffInfoCommand::extractSffInfo(string input, string output){
	try {
		
		ofstream out;
		openOutputFile(output, out);
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		ifstream in;
		in.open(input.c_str(), ios::binary);
		
		CommonHeader header; 
		readCommonHeader(in, header);
		
		//print common header
		printCommonHeader(out, header);
		
		int count = 0;
		
		//check magic number and version
		if (header.magicNumber != 779314790) { m->mothurOut("Magic Number is not correct, not a valid .sff file"); m->mothurOutEndLine(); return count; }
		if (header.version != "0001") { m->mothurOut("Version is not supported, only support version 0001."); m->mothurOutEndLine(); return count; }
		
		//read through the sff file
		while (!in.eof()) {
			
			//read header
			Header readheader;
			readHeader(in, readheader);
		
			//print header
			printHeader(out, readheader);
			
			//read data
			seqRead read; 
			readSeqData(in, read, header.numFlowsPerRead, readheader.numBases);
		
			//print data
			printSeqData(out, read);
			
			count++;
			
			//report progress
			if((count+1) % 500 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
		
			if (m->control_pressed) { count = 0; break;   }
			
			if (count >= header.numReads) { break; }
		}
		
		//report progress
		if (!m->control_pressed) {   if((count) % 500 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
		
		in.close();
		out.close();
		
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
			char buffer[sizeof(header.magicNumber)];
			in.read(buffer, sizeof(header.magicNumber));
			header.magicNumber = be_int4(*(unsigned int *)(&buffer));
			
			//read version
			char buffer9[4];
			in.read(buffer9, 4);
			header.version = "";
			for (int i = 0; i < 4; i++) {  header.version += toString((int)(buffer9[i])); }
				
			//read offset
			char buffer2 [sizeof(header.indexOffset)];
			in.read(buffer2, sizeof(header.indexOffset));
			header.indexOffset =  be_int8(*(unsigned long int *)(&buffer2));
			
			//read index length
			char buffer3 [sizeof(header.indexLength)];
			in.read(buffer3, sizeof(header.indexLength));
			header.indexLength =  be_int4(*(unsigned int *)(&buffer3));
			
			//read num reads
			char buffer4 [sizeof(header.numReads)];
			in.read(buffer4, sizeof(header.numReads));
			header.numReads =  be_int4(*(unsigned int *)(&buffer4));
				
			//read header length
			char buffer5 [sizeof(header.headerLength)];
			in.read(buffer5, sizeof(header.headerLength));
			header.headerLength =  be_int2(*(unsigned short *)(&buffer5));
					
			//read key length
			char buffer6 [sizeof(header.keyLength)];
			in.read(buffer6, sizeof(header.keyLength));
			header.keyLength = be_int2(*(unsigned short *)(&buffer6));
			
			//read number of flow reads
			char buffer7 [sizeof(header.numFlowsPerRead)];
			in.read(buffer7, sizeof(header.numFlowsPerRead));
			header.numFlowsPerRead =  be_int2(*(unsigned short *)(&buffer7));
				
			//read format code
			char buffer8 [1];
			in.read(buffer8, 1);
			header.flogramFormatCode = (int)(buffer8[0]);
			
			//read flow chars
			char tempBuffer [header.numFlowsPerRead];
			in.read(tempBuffer, header.numFlowsPerRead); 
			header.flowChars = tempBuffer;
			if (header.flowChars.length() > header.numFlowsPerRead) { header.flowChars = header.flowChars.substr(0, header.numFlowsPerRead);  }
			
			//read key
			char tempBuffer2 [header.keyLength];
			in.read(tempBuffer2, header.keyLength);
			header.keySequence = tempBuffer2;
			if (header.keySequence.length() > header.keyLength) { header.keySequence = header.keySequence.substr(0, header.keyLength);  }
				
			/* Pad to 8 chars */
			int spotInFile = in.tellg();
			int spot = (spotInFile + 7)& ~7;  // ~ inverts
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
			char buffer [sizeof(header.headerLength)];
			in.read(buffer, sizeof(header.headerLength));
			header.headerLength = be_int2(*(unsigned short *)(&buffer));
						
			//read name length
			char buffer2 [sizeof(header.nameLength)];
			in.read(buffer2, sizeof(header.nameLength));
			header.nameLength = be_int2(*(unsigned short *)(&buffer2));

			//read num bases
			char buffer3 [sizeof(header.numBases)];
			in.read(buffer3, sizeof(header.numBases));
			header.numBases =  be_int4(*(unsigned int *)(&buffer3));
			
			//read clip qual left
			char buffer4 [sizeof(header.clipQualLeft)];
			in.read(buffer4, sizeof(header.clipQualLeft));
			header.clipQualLeft =  be_int2(*(unsigned short *)(&buffer4));
			
			//read clip qual right
			char buffer5 [sizeof(header.clipQualRight)];
			in.read(buffer5, sizeof(header.clipQualRight));
			header.clipQualRight =  be_int2(*(unsigned short *)(&buffer5));
			
			//read clipAdapterLeft
			char buffer6 [sizeof(header.clipAdapterLeft)];
			in.read(buffer6, sizeof(header.clipAdapterLeft));
			header.clipAdapterLeft = be_int2(*(unsigned short *)(&buffer6));

			//read clipAdapterRight
			char buffer7 [sizeof(header.clipAdapterRight)];
			in.read(buffer7, sizeof(header.clipAdapterRight));
			header.clipAdapterRight = be_int2(*(unsigned short *)(&buffer7));
		
			//read name
			char tempBuffer [header.nameLength];
			in.read(tempBuffer, header.nameLength);
			header.name = tempBuffer;
			if (header.name.length() > header.nameLength) { header.name = header.name.substr(0, header.nameLength);  }
			
			/* Pad to 8 chars */
			int spotInFile = in.tellg();
			int spot = (spotInFile + 7)& ~7;
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
				char buffer [sizeof(unsigned short)];
				in.read(buffer, (sizeof(unsigned short)));
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
			char tempBuffer[numBases];
			in.read(tempBuffer, numBases);
			read.bases = tempBuffer;
			if (read.bases.length() > numBases) { read.bases = read.bases.substr(0, numBases);  }

			//read flowgram
			read.qualScores.resize(numBases);
			for (int i = 0; i < numBases; i++) {  
				char temp[1];
				in.read(temp, 1);
				read.qualScores[i] = be_int1(*(unsigned char *)(&temp));
			}
		
			/* Pad to 8 chars */
			int spotInFile = in.tellg();
			int spot = (spotInFile + 7)& ~7;
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
int SffInfoCommand::printSeqData(ofstream& out, seqRead& read) {
	try {
		
		out << "FlowGram: ";
		for (int i = 0; i < read.flowgram.size(); i++) { out << setprecision(2) << (read.flowgram[i]/(float)100) << '\t';  }
		
		out << endl <<  "Flow Indexes: ";
		int sum = 0;
		for (int i = 0; i < read.flowIndex.size(); i++) {  sum +=  read.flowIndex[i];  out << sum << '\t'; }
		
		out << endl <<  "Bases: " << read.bases << endl << "Quality Scores: ";
		for (int i = 0; i < read.qualScores.size(); i++) {   out << read.qualScores[i] << '\t';  }
		out << endl << endl;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffInfoCommand", "printSeqData");
		exit(1);
	}
}
//**********************************************************************************************************************/
