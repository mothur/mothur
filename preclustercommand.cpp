/*
 *  preclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "preclustercommand.h"

//**********************************************************************************************************************
inline bool comparePriority(seqPNode first, seqPNode second) {  return (first.numIdentical > second.numIdentical); }
//**********************************************************************************************************************

PreClusterCommand::PreClusterCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "name", "diffs"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { mothurOut("fasta is a required parameter for the pre.cluster command."); mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { abort = true; }	
			else {  readNameFile();  }
			
			string temp	= validParameter.validFile(parameters, "diffs", false);				if(temp == "not found"){	temp = "1"; }
			convert(temp, diffs); 
		}
				
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
PreClusterCommand::~PreClusterCommand(){}	
//**********************************************************************************************************************

void PreClusterCommand::help(){
	try {
		mothurOut("The pre.cluster command groups sequences that are within a given number of base mismatches.\n");
		mothurOut("The pre.cluster command outputs a new fasta and name file.\n");
		mothurOut("The pre.cluster command parameters are fasta, names and diffs. The fasta parameter is required. \n");
		mothurOut("The names parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n");
		mothurOut("The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n");
		mothurOut("The pre.cluster command should be in the following format: \n");
		mothurOut("pre.cluster(fasta=yourFastaFile, names=yourNamesFile, diffs=yourMaxDiffs) \n");
		mothurOut("Example pre.cluster(fasta=amazon.fasta, diffs=2).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************

int PreClusterCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//reads fasta file and return number of seqs
		int numSeqs = readSeqs(); //fills alignSeqs and makes all seqs active
	
		if (numSeqs == 0) { mothurOut("Error reading fasta file...please correct."); mothurOutEndLine(); return 0;  }
		if (diffs > length) { mothurOut("Error: diffs is greater than your sequence length."); mothurOutEndLine(); return 0;  }
		
		//clear sizes since you only needed this info to build the alignSeqs seqPNode structs
		sizes.clear();
	
		//sort seqs by number of identical seqs
		sort(alignSeqs.begin(), alignSeqs.end(), comparePriority);
	
		//go through active list and cluster everthing you can, until all nodes are inactive
		//taking advantage of the fact that maps are already sorted
		map<string, bool>::iterator itActive;
		map<string, bool>::iterator it2Active;
		int count = 0;
		
		for (int i = 0; i < alignSeqs.size(); i++) {
		
			//are you active
			itActive = active.find(alignSeqs[i].seq.getName());
			
			if (itActive != active.end()) {  //this sequence has not been merged yet
			
				//try to merge it with all smaller seqs
				for (int j = i; j < alignSeqs.size(); j++) {
					
					if (i != j) {
						//are you active
						it2Active = active.find(alignSeqs[j].seq.getName());
						if (it2Active != active.end()) {  //this sequence has not been merged yet
							//are you within "diff" bases
							int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
							
							if (mismatch <= diffs) {
								//merge
								names[alignSeqs[i].seq.getName()] += "," + names[alignSeqs[j].seq.getName()];
							
								//remove from active list
								active.erase(it2Active);
								
								//set numIdentical to 0, so you only print the representative seqs in the fasta file
								alignSeqs[j].numIdentical = 0;
								count++;
							}
						}//end if j active
					}//end if i != j
				}//end for loop
				
				//remove from active list 
				active.erase(itActive);
			}//end if active i
		}
	
		string newFastaFile = getRootName(fastafile) + "precluster" + getExtension(fastafile);
		string newNamesFile = getRootName(fastafile) + "precluster.names";
		
		printData(newFastaFile, newNamesFile);
		
		mothurOut("Total number of sequences before precluster was " + toString(alignSeqs.size()) + "."); mothurOutEndLine();
		mothurOut("pre.cluster removed " + toString(count) + " sequences."); mothurOutEndLine(); 
		
		return 0;
		
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::readSeqs(){
	try {
		ifstream inFasta;
		openInputFile(fastafile, inFasta);
		length = 0;
		map<string, string>::iterator it;
				
		while (!inFasta.eof()) {
			Sequence temp(inFasta);  //read seq
			
			if (temp.getName() != "") {  
				if (namefile != "") {
					//make sure fasta and name files match
					it = names.find(temp.getName());
					if (it == names.end()) { mothurOut(temp.getName() + " is not in your names file, please correct."); mothurOutEndLine(); exit(1); }
				}else { sizes[temp.getName()] = 1; }
				
				seqPNode tempNode(sizes[temp.getName()], temp);
				alignSeqs.push_back(tempNode); 
				active[temp.getName()] = true;
			}
			gobble(inFasta);
		}
		inFasta.close();
		
		if (alignSeqs.size() != 0) {  length = alignSeqs[0].seq.getAligned().length();  }
		
		return alignSeqs.size();
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "readSeqs");
		exit(1);
	}
}
/**************************************************************************************************/
void PreClusterCommand::readNameFile(){
	try {
		ifstream in;
		openInputFile(namefile, in);
		string firstCol, secondCol;
				
		while (!in.eof()) {
			in >> firstCol >> secondCol; gobble(in);
			names[firstCol] = secondCol;
			int size = 1;
			while (secondCol.find_first_of(',') != -1) { 
				size++;
				secondCol = secondCol.substr(secondCol.find_first_of(',')+1, secondCol.length());
			}
			sizes[firstCol] = size;
		}
		in.close();
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "readNameFile");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::calcMisMatches(string seq1, string seq2){
	try {
		int numBad = 0;
		
		for (int i = 0; i < seq1.length(); i++) {
			//do they match
			if (seq1[i] != seq2[i]) { numBad++; }
			if (numBad > diffs) { return length;  } //to far to cluster
		}
		
		return numBad;
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "calcMisMatches");
		exit(1);
	}
}
/**************************************************************************************************/
void PreClusterCommand::printData(string newfasta, string newname){
	try {
		ofstream outFasta;
		ofstream outNames;
		openOutputFile(newfasta, outFasta);
		openOutputFile(newname, outNames);
				
		map<string, string>::iterator itNames;
		
		for (int i = 0; i < alignSeqs.size(); i++) {
			if (alignSeqs[i].numIdentical != 0) {
				alignSeqs[i].seq.printSequence(outFasta); 
				
				itNames = names.find(alignSeqs[i].seq.getName());
				
				outNames << itNames->first << '\t' << itNames->second << endl;
			}
		}
		
		outFasta.close();
		outNames.close();
		
	}
	catch(exception& e) {
		errorOut(e, "PreClusterCommand", "printData");
		exit(1);
	}
}

/**************************************************************************************************/


