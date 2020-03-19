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
vector<string> ClusterFragmentsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pdiffs("diffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pdiffs);
		CommandParameter ppercent("percent", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppercent);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterFragmentsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The cluster.fragments command groups sequences that are part of a larger sequence.\n";
		helpString += "The cluster.fragments command outputs a new fasta and name or count file.\n";
		helpString += "The cluster.fragments command parameters are fasta, name, count, diffs and percent. The fasta parameter is required, unless you have a valid current file. \n";
		helpString += "The names parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The diffs parameter allows you to set the number of differences allowed, default=0. \n";
		helpString += "The percent parameter allows you to set percentage of differences allowed, default=0. percent=2 means if the number of difference is less than or equal to two percent of the length of the fragment, then cluster.\n";
		helpString += "You may use diffs and percent at the same time to say something like: If the number or differences is greater than 1 or more than 2% of the fragment length, don't merge. \n";
		helpString += "The cluster.fragments command should be in the following format: \n";
		helpString += "cluster.fragments(fasta=yourFastaFile, names=yourNamesFile) \n";
		helpString += "Example cluster.fragments(fasta=amazon.fasta).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterFragmentsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],fragclust.fasta"; } 
        else if (type == "name") {  pattern = "[filename],fragclust.names"; } 
        else if (type == "count") {  pattern = "[filename],fragclust.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFragmentsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ClusterFragmentsCommand::ClusterFragmentsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "ClusterFragmentsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterFragmentsCommand::ClusterFragmentsCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (!validParameter.isValidParameter(it2->first, myArray, it2->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
            outputTypes["count"] = tempOutNames;

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else { current->setFastaFile(fastafile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(fastafile); 	}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { namefile = ""; abort = true; }	
			else {  readNameFile(); current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { ct.readTable(countfile, true, false); current->setCountFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster.fragments command you must enter ONLY ONE of the following: count or name.\n");  abort = true; }
			
			string temp;
			temp = validParameter.valid(parameters, "diffs");		if (temp == "not found"){	temp = "0";				}
			util.mothurConvert(temp, diffs); 
			
			temp = validParameter.valid(parameters, "percent");		if (temp == "not found"){	temp = "0";				}
			util.mothurConvert(temp, percent);
			
			if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
			
		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterFragmentsCommand", "ClusterFragmentsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterFragmentsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		long start = time(NULL);
		
		//reads fasta file and return number of seqs
		int numSeqs = readFASTA(); //fills alignSeqs and makes all seqs active
		
		if (m->getControl_pressed()) { return 0; }
	
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
					
					if (m->getControl_pressed()) { return 0; }
					
					if (alignSeqs[j].active) {  //this sequence has not been merged yet
						
						string jBases = alignSeqs[j].seq.getUnaligned();
													
						if (isFragment(iBases, jBases)) {
                            if (countfile != "") {
                                ct.mergeCounts(alignSeqs[i].names, alignSeqs[j].names);
                            }else {
                                //merge
                                alignSeqs[i].names += ',' + alignSeqs[j].names;
                                alignSeqs[i].numIdentical += alignSeqs[j].numIdentical;
                            }
							alignSeqs[j].active = 0;
							alignSeqs[j].numIdentical = 0;
							count++;
						}
					}//end if j active
				}//end if i != j
			
				//remove from active list 
				alignSeqs[i].active = 0;
				
			}//end if active i
			if(i % 100 == 0)	{ m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
		}
		
		if(numSeqs % 100 != 0)	{ m->mothurOutJustToScreen(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n");	}
	
		
		string fileroot = outputDir + util.getRootName(util.getSimpleName(fastafile));
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
		string newFastaFile = getOutputFileName("fasta", variables);
		string newNamesFile = getOutputFileName("name", variables);
        if (countfile != "") { newNamesFile = getOutputFileName("count", variables); }
		
		if (m->getControl_pressed()) { return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Total number of sequences before cluster.fragments was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
		m->mothurOut("cluster.fragments removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
		
		printData(newFastaFile, newNamesFile);
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
		
		if (m->getControl_pressed()) { util.mothurRemove(newFastaFile); util.mothurRemove(newNamesFile); return 0; }
		
		m->mothurOut("\nOutput File Names: \n"); 
		m->mothurOut(newFastaFile); m->mothurOutEndLine();	
		m->mothurOut(newNamesFile); m->mothurOutEndLine();	
		outputNames.push_back(newFastaFile);  outputNames.push_back(newNamesFile); outputTypes["fasta"].push_back(newFastaFile); outputTypes["name"].push_back(newNamesFile);
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
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
			int totalDiffs = 0;
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
		util.openInputFile(fastafile, inFasta);
		
		while (!inFasta.eof()) {
			
			if (m->getControl_pressed()) { inFasta.close(); return 0; }
			
			Sequence seq(inFasta);  util.gobble(inFasta);
			
			if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file
				if (namefile != "") {
					itSize = sizes.find(seq.getName());
					
					if (itSize == sizes.end()) { m->mothurOut(seq.getName() + " is not in your names file, please correct."); m->mothurOutEndLine(); exit(1); }
					else{
						seqRNode tempNode(itSize->second, seq, names[seq.getName()], seq.getUnaligned().length());
						alignSeqs.push_back(tempNode);
					}
                }else if(countfile != "") {
                    seqRNode tempNode(ct.getNumSeqs(seq.getName()), seq, seq.getName(), seq.getUnaligned().length());
                    alignSeqs.push_back(tempNode);
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
		
		util.openOutputFile(newfasta, outFasta);
		if (countfile == "") {  util.openOutputFile(newname, outNames); }
		
		for (int i = 0; i < alignSeqs.size(); i++) {
			if (alignSeqs[i].numIdentical != 0) {
				alignSeqs[i].seq.printSequence(outFasta); 
				if (countfile == "") {  outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;  }
			}
		}
		
		outFasta.close();
		if (countfile == "") {  outNames.close(); }
        else { ct.printTable(newname); }
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
		util.openInputFile(namefile, in);
		string firstCol, secondCol;
				
		while (!in.eof()) {
			in >> firstCol >> secondCol; util.gobble(in);
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

