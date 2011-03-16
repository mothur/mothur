/*
 *  ryanscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/23/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clusterfragmentscommand.h"
#include "needlemanoverlap.hpp"

//**********************************************************************************************************************
//sort by unaligned
inline bool comparePriority(seqRNode first, seqRNode second) {  
	bool better = false;
	
	if (first.length > second.length) { 
		better = true;
	}else if (first.length == second.length) {
		if (first.numIdentical > second.numIdentical) {
			better = true;
		}
	}
	
	return better; 
}
//**********************************************************************************************************************
vector<string> ClusterFragmentsCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"fasta","name","diffs","percent","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterFragmentsCommand::ClusterFragmentsCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "ClusterFragmentsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClusterFragmentsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClusterFragmentsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterFragmentsCommand::ClusterFragmentsCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","name","diffs","percent","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the cluster.fragments command."); m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(fastafile); 	}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { abort = true; }	
			else {  readNameFile();  }
			
			string temp;
			temp = validParameter.validFile(parameters, "diffs", false);		if (temp == "not found"){	temp = "0";				}
			convert(temp, diffs); 
			
			temp = validParameter.validFile(parameters, "percent", false);		if (temp == "not found"){	temp = "0";				}
			convert(temp, percent);
			
		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "ClusterFragmentsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
ClusterFragmentsCommand::~ClusterFragmentsCommand(){}	
//**********************************************************************************************************************
void ClusterFragmentsCommand::help(){
	try {
		m->mothurOut("The cluster.fragments command groups sequences that are part of a larger sequence.\n");
		m->mothurOut("The cluster.fragments command outputs a new fasta and name file.\n");
		m->mothurOut("The cluster.fragments command parameters are fasta, name, diffs and percent. The fasta parameter is required. \n");
		m->mothurOut("The names parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n");
		m->mothurOut("The diffs parameter allows you to set the number of differences allowed, default=0. \n");
		m->mothurOut("The percent parameter allows you to set percentage of differences allowed, default=0. percent=2 means if the number of difference is less than or equal to two percent of the length of the fragment, then cluster.\n");
		m->mothurOut("You may use diffs and percent at the same time to say something like: If the number or differences is greater than 1 or more than 2% of the fragment length, don't merge. \n");
		m->mothurOut("The cluster.fragments command should be in the following format: \n");
		m->mothurOut("cluster.fragments(fasta=yourFastaFile, names=yourNamesFile) \n");
		m->mothurOut("Example cluster.fragments(fasta=amazon.fasta).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterFragmentsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int start = time(NULL);
		
		//reads fasta file and return number of seqs
		int numSeqs = readFASTA(); //fills alignSeqs and makes all seqs active
		
		if (m->control_pressed) { return 0; }
	
		if (numSeqs == 0) { m->mothurOut("Error reading fasta file...please correct."); m->mothurOutEndLine(); return 0;  }
		
		//sort seqs by length of unaligned sequence
		sort(alignSeqs.begin(), alignSeqs.end(), comparePriority);
	
		int count = 0;

		//think about running through twice...
		for (int i = 0; i < numSeqs; i++) {
			
			if (alignSeqs[i].active) {  //this sequence has not been merged yet
				
				string iBases = alignSeqs[i].seq.getUnaligned();
				
				//try to merge it with all smaller seqs
				for (int j = i+1; j < numSeqs; j++) {
					
					if (m->control_pressed) { return 0; }
					
					if (alignSeqs[j].active) {  //this sequence has not been merged yet
						
						string jBases = alignSeqs[j].seq.getUnaligned();
													
						if (isFragment(iBases, jBases)) {
							//merge
							alignSeqs[i].names += ',' + alignSeqs[j].names;
							alignSeqs[i].numIdentical += alignSeqs[j].numIdentical;

							alignSeqs[j].active = 0;
							alignSeqs[j].numIdentical = 0;
							count++;
						}
					}//end if j active
				}//end if i != j
			
				//remove from active list 
				alignSeqs[i].active = 0;
				
			}//end if active i
			if(i % 100 == 0)	{ m->mothurOut(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); m->mothurOutEndLine();	}
		}
		
		if(numSeqs % 100 != 0)	{ m->mothurOut(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); m->mothurOutEndLine();	}
	
		
		string fileroot = outputDir + m->getRootName(m->getSimpleName(fastafile));
		
		string newFastaFile = fileroot + "fragclust.fasta";
		string newNamesFile = fileroot + "names";
		
		if (m->control_pressed) { return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Total number of sequences before cluster.fragments was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
		m->mothurOut("cluster.fragments removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
		
		printData(newFastaFile, newNamesFile);
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
		
		if (m->control_pressed) { remove(newFastaFile.c_str()); remove(newNamesFile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(newFastaFile); m->mothurOutEndLine();	
		m->mothurOut(newNamesFile); m->mothurOutEndLine();	
		outputNames.push_back(newFastaFile);  outputNames.push_back(newNamesFile); outputTypes["fasta"].push_back(newFastaFile); outputTypes["name"].push_back(newNamesFile);
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************
bool ClusterFragmentsCommand::isFragment(string seq1, string seq2){
	try {
		bool fragment = false;
		
		//exact match
		int pos = seq1.find(seq2);
		if (pos != string::npos) { return true; }
		//no match, no diffs wanted
		else if ((diffs == 0) && (percent == 0)) { return false; }
		else { //try aligning and see if you can find it
			
			//find number of acceptable differences for this sequence fragment
			int totalDiffs;
			if (diffs == 0) { //you didnt set diffs you want a percentage
				totalDiffs = floor((seq2.length() * (percent / 100.0)));
			}else if (percent == 0) { //you didn't set percent you want diffs
				totalDiffs = diffs;
			}else if ((percent != 0) && (diffs != 0)) { //you want both, set total diffs to smaller of 2
				totalDiffs = diffs;
				int percentDiff = floor((seq2.length() * (percent / 100.0)));
				if (percentDiff < totalDiffs) { totalDiffs = percentDiff; }
			}
			
			Alignment* alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (seq1.length()+totalDiffs+1));
							
			//use needleman to align 
			alignment->align(seq2, seq1);
			string tempSeq2 = alignment->getSeqAAln();
			string temp = alignment->getSeqBAln();
			
			delete alignment;
			
			//chop gap ends
			int startPos = 0;
			int endPos = tempSeq2.length()-1;
			for (int i = 0; i < tempSeq2.length(); i++) {  if (isalpha(tempSeq2[i])) { startPos = i; break; } }
			for (int i = tempSeq2.length()-1; i >= 0; i--) {  if (isalpha(tempSeq2[i])) { endPos = i; break; } }
			
			//count number of diffs
			int numDiffs = 0;
			for (int i = startPos; i <= endPos; i++) {
				if (tempSeq2[i] != temp[i]) { numDiffs++; }
			}
			
			if (numDiffs <= totalDiffs) { fragment = true; }
			
		}
		
		return fragment;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "isFragment");
		exit(1);
	}
}
/**************************************************************************************************/
int ClusterFragmentsCommand::readFASTA(){
	try {
	
		ifstream inFasta;
		m->openInputFile(fastafile, inFasta);
		
		while (!inFasta.eof()) {
			
			if (m->control_pressed) { inFasta.close(); return 0; }
			
			Sequence seq(inFasta);  m->gobble(inFasta);
			
			if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file
				if (namefile != "") {
					itSize = sizes.find(seq.getName());
					
					if (itSize == sizes.end()) { m->mothurOut(seq.getName() + " is not in your names file, please correct."); m->mothurOutEndLine(); exit(1); }
					else{
						seqRNode tempNode(itSize->second, seq, names[seq.getName()], seq.getUnaligned().length());
						alignSeqs.push_back(tempNode);
					}	
				}else { //no names file, you are identical to yourself 
					seqRNode tempNode(1, seq, seq.getName(), seq.getUnaligned().length());
					alignSeqs.push_back(tempNode);
				}
			}
		}
		
		inFasta.close();
		return alignSeqs.size();
	}
	
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "readFASTA");
		exit(1);
	}
}
/**************************************************************************************************/
void ClusterFragmentsCommand::printData(string newfasta, string newname){
	try {
		ofstream outFasta;
		ofstream outNames;
		
		m->openOutputFile(newfasta, outFasta);
		m->openOutputFile(newname, outNames);
		
		for (int i = 0; i < alignSeqs.size(); i++) {
			if (alignSeqs[i].numIdentical != 0) {
				alignSeqs[i].seq.printSequence(outFasta); 
				outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;
			}
		}
		
		outFasta.close();
		outNames.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "printData");
		exit(1);
	}
}
/**************************************************************************************************/

void ClusterFragmentsCommand::readNameFile(){
	try {
		ifstream in;
		m->openInputFile(namefile, in);
		string firstCol, secondCol;
				
		while (!in.eof()) {
			in >> firstCol >> secondCol; m->gobble(in);
			names[firstCol] = secondCol;
			int size = 1;

			for(int i=0;i<secondCol.size();i++){
				if(secondCol[i] == ','){	size++;	}
			}
			sizes[firstCol] = size;
		}
		in.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "readNameFile");
		exit(1);
	}
}

/**************************************************************************************************/

