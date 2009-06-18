/*
 *  screenseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "screenseqscommand.h"
#include "sequence.hpp"

//***************************************************************************************************************

ScreenSeqsCommand::ScreenSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta", "start", "end", "maxambig", "maxhomop", "minlength", "maxlength",
									"name", "group", "alignreport"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { cout << "fasta is a required parameter for the screen.seqs command." << endl; abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	

			alignreport = validParameter.validFile(parameters, "alignreport", true);
			if (alignreport == "not open") { abort = true; }
			else if (alignreport == "not found") { namefile = ""; }	
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "start", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, startPos); 
		
			temp = validParameter.validFile(parameters, "end", false);			if (temp == "not found") { temp = "-1"; }
			convert(temp, endPos);  

			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxAmbig);  

			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "-1"; }
			convert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "-1"; }
			convert(temp, maxLength); 
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ScreenSeqsCommand class Function ScreenSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ScreenSeqsCommand class function ScreenSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void ScreenSeqsCommand::help(){
	try {
		cout << "The screen.seqs command reads a fastafile and creates ....." << "\n";
		cout << "The screen.seqs command parameters are fasta, start, end, maxambig, maxhomop, minlength, maxlength, name, and group." << "\n";
		cout << "The fasta parameter is required." << "\n";
		cout << "The start parameter .... The default is -1." << "\n";
		cout << "The end parameter .... The default is -1." << "\n";
		cout << "The maxambig parameter .... The default is -1." << "\n";
		cout << "The maxhomop parameter .... The default is -1." << "\n";
		cout << "The minlength parameter .... The default is -1." << "\n";
		cout << "The maxlength parameter .... The default is -1." << "\n";
		cout << "The name parameter allows you to provide a namesfile, and the group parameter allows you to provide a groupfile." << "\n";
		cout << "The screen.seqs command should be in the following format: " << "\n";
		cout << "screen.seqs(fasta=yourFastaFile, name=youNameFile, group=yourGroupFIle, start=yourStart, end=yourEnd, maxambig=yourMaxambig,  " << "\n";
		cout << "maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  " << "\n";	
		cout << "Example screen.seqs(fasta=abrecovery.fasta, name=abrecovery.names, group=abrecovery.groups, start=..., end=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...)." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta)." << "\n" << "\n";

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ScreenSeqsCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ScreenSeqsCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//***************************************************************************************************************

ScreenSeqsCommand::~ScreenSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ScreenSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
				
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);
		
		set<string> badSeqNames;
		
		string goodSeqFile = getRootName(fastafile) + "good" + getExtension(fastafile);
		string badSeqFile = getRootName(fastafile) + "bad" + getExtension(fastafile);
		
		ofstream goodSeqOut;	openOutputFile(goodSeqFile, goodSeqOut);
		ofstream badSeqOut;		openOutputFile(badSeqFile, badSeqOut);		
		
		while(!inFASTA.eof()){
			Sequence currSeq(inFASTA);
			bool goodSeq = 1;		//	innocent until proven guilty
			if(goodSeq == 1 && startPos != -1 && startPos < currSeq.getStartPos())			{	goodSeq = 0;	}
			if(goodSeq == 1 && endPos != -1 && endPos > currSeq.getEndPos())				{	goodSeq = 0;	}
			if(goodSeq == 1 && maxAmbig != -1 && maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	}
			if(goodSeq == 1 && maxHomoP != -1 && maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	}
			if(goodSeq == 1 && minLength != -1 && minLength > currSeq.getNumBases())		{	goodSeq = 0;	}
			if(goodSeq == 1 && maxLength != -1 && maxLength < currSeq.getNumBases())		{	goodSeq = 0;	}
			
			if(goodSeq == 1){
				currSeq.printSequence(goodSeqOut);	
			}
			else{
				currSeq.printSequence(badSeqOut);	
				badSeqNames.insert(currSeq.getName());
			}
			gobble(inFASTA);
		}	
		if(namefile != "")		{	screenNameGroupFile(badSeqNames);	}
		if(groupfile != "")		{	screenGroupFile(badSeqNames);		}
		if(alignreport != "")	{	screenAlignReport(badSeqNames);		}
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ScreenSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ScreenSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}

//***************************************************************************************************************

void ScreenSeqsCommand::screenNameGroupFile(set<string> badSeqNames){

	ifstream inputNames;
	openInputFile(namefile, inputNames);
	set<string> badSeqGroups;
	string seqName, seqList, group;
	set<string>::iterator it;

	string goodNameFile = getRootName(namefile) + "good" + getExtension(namefile);
	string badNameFile = getRootName(namefile) + "bad" + getExtension(namefile);
	
	ofstream goodNameOut;	openOutputFile(goodNameFile, goodNameOut);
	ofstream badNameOut;	openOutputFile(badNameFile, badNameOut);		
	
	while(!inputNames.eof()){
		inputNames >> seqName >> seqList;
		it = badSeqNames.find(seqName);
		
		if(it != badSeqNames.end()){
			badSeqNames.erase(it);
			badNameOut << seqName << '\t' << seqList << endl;
			if(namefile != ""){
				int start = 0;
				for(int i=0;i<seqList.length();i++){
					if(seqList[i] == ','){
						badSeqGroups.insert(seqList.substr(start,i-start));
						start = i+1;
					}					
				}
				badSeqGroups.insert(seqList.substr(start,seqList.length()-start));
			}
		}
		else{
			goodNameOut << seqName << '\t' << seqList << endl;
		}
		gobble(inputNames);
	}
	inputNames.close();
	goodNameOut.close();
	badNameOut.close();
	
	if(groupfile != ""){
		
		ifstream inputGroups;
		openInputFile(groupfile, inputGroups);

		string goodGroupFile = getRootName(groupfile) + "good" + getExtension(groupfile);
		string badGroupFile = getRootName(groupfile) + "bad" + getExtension(groupfile);
		
		ofstream goodGroupOut;	openOutputFile(goodGroupFile, goodGroupOut);
		ofstream badGroupOut;	openOutputFile(badGroupFile, badGroupOut);		
		
		while(!inputGroups.eof()){
			inputGroups >> seqName >> group;

			it = badSeqGroups.find(seqName);
			
			if(it != badSeqGroups.end()){
				badSeqGroups.erase(it);
				badGroupOut << seqName << '\t' << group << endl;
			}
			else{
				goodGroupOut << seqName << '\t' << group << endl;
			}
			gobble(inputGroups);
		}
		inputGroups.close();
		goodGroupOut.close();
		badGroupOut.close();
	}
}

//***************************************************************************************************************

void ScreenSeqsCommand::screenGroupFile(set<string> badSeqNames){

	ifstream inputGroups;
	openInputFile(groupfile, inputGroups);
	string seqName, group;
	set<string>::iterator it;
	
	string goodGroupFile = getRootName(groupfile) + "good" + getExtension(groupfile);
	string badGroupFile = getRootName(groupfile) + "bad" + getExtension(groupfile);
	
	ofstream goodGroupOut;	openOutputFile(goodGroupFile, goodGroupOut);
	ofstream badGroupOut;	openOutputFile(badGroupFile, badGroupOut);		
	
	while(!inputGroups.eof()){
		inputGroups >> seqName >> group;
		it = badSeqNames.find(seqName);
		
		if(it != badSeqNames.end()){
			badSeqNames.erase(it);
			badGroupOut << seqName << '\t' << group << endl;
		}
		else{
			goodGroupOut << seqName << '\t' << group << endl;
		}
		gobble(inputGroups);
	}
	inputGroups.close();
	goodGroupOut.close();
	badGroupOut.close();
	
}

//***************************************************************************************************************

void ScreenSeqsCommand::screenAlignReport(set<string> badSeqNames){
	
	ifstream inputAlignReport;
	openInputFile(alignreport, inputAlignReport);
	string seqName, group;
	set<string>::iterator it;
	
	string goodAlignReportFile = getRootName(alignreport) + "good" + getExtension(alignreport);
	string badAlignReportFile = getRootName(alignreport) + "bad" + getExtension(alignreport);
	
	ofstream goodAlignReportOut;	openOutputFile(goodAlignReportFile, goodAlignReportOut);
	ofstream badAlignReportOut;		openOutputFile(badAlignReportFile, badAlignReportOut);		

	while (!inputAlignReport.eof())	{		//	need to copy header
		char c = inputAlignReport.get();
		goodAlignReportOut << c;
		badAlignReportOut << c;
		if (c == 10 || c == 13){	break;	}	
	}

	while(!inputAlignReport.eof()){
		inputAlignReport >> seqName;
		it = badSeqNames.find(seqName);
		string line;		
		while (!inputAlignReport.eof())	{		//	need to copy header
			char c = inputAlignReport.get();
			line += c;
			if (c == 10 || c == 13){	break;	}	
		}
		
		if(it != badSeqNames.end()){
			badSeqNames.erase(it);
			badAlignReportOut << seqName << '\t' << line;;
		}
		else{
			goodAlignReportOut << seqName << '\t' << line;
		}
		gobble(inputAlignReport);
	}
	inputAlignReport.close();
	goodAlignReportOut.close();
	badAlignReportOut.close();
	
}

//***************************************************************************************************************


