/*
 *  preclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "preclustercommand.h"
#include "deconvolutecommand.h"


//**********************************************************************************************************************
vector<string> PreClusterCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pdiffs("diffs", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter palign("align", "Multiple", "needleman-gotoh-blast-noalign", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
        CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
        CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
        CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);

        CommandParameter ptopdown("topdown", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(ptopdown);

		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PreClusterCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The pre.cluster command groups sequences that are within a given number of base mismatches.\n";
		helpString += "The pre.cluster command outputs a new fasta and name file.\n";
		helpString += "The pre.cluster command parameters are fasta, name, group, count, topdown, processors and diffs. The fasta parameter is required. \n";
		helpString += "The name parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The group parameter allows you to provide a group file so you can cluster by group. \n";
        helpString += "The count parameter allows you to provide a count file so you can cluster by group. \n";
		helpString += "The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n";
        helpString += "The topdown parameter allows you to specify whether to cluster from largest abundance to smallest or smallest to largest.  Default=T, meaning largest to smallest.\n";
        helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n";
        helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
        helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
        helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
        helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
		helpString += "The pre.cluster command should be in the following format: \n";
		helpString += "pre.cluster(fasta=yourFastaFile, names=yourNamesFile, diffs=yourMaxDiffs) \n";
		helpString += "Example pre.cluster(fasta=amazon.fasta, diffs=2).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PreClusterCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],precluster,[extension]"; } 
        else if (type == "name") {  pattern = "[filename],precluster.names"; } 
        else if (type == "count") {  pattern = "[filename],precluster.count_table"; }
        else if (type == "map") {  pattern =  "[filename],precluster.map"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
PreClusterCommand::PreClusterCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		outputTypes["map"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

PreClusterCommand::PreClusterCommand(string option) {
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
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["map"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { current->setFastaFile(fastafile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += util.hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { namefile = ""; abort = true; }	
			else {  current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") { groupfile =  "";  bygroup = false; }
			else if (groupfile == "not open") { abort = true; groupfile =  ""; }	
			else {   current->setGroupFile(groupfile); bygroup = true;  }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not found") { countfile =  "";   }
			else if (countfile == "not open") { abort = true; countfile =  ""; }	
			else {   
                current->setCountFile(countfile); 
                ct.readTable(countfile, true, false);
                if (ct.hasGroupInfo()) { bygroup = true; }
                else { bygroup = false;  }
            }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			string temp	= validParameter.valid(parameters, "diffs");		if(temp == "not found"){	temp = "1"; }
			util.mothurConvert(temp, diffs); 
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
            temp = validParameter.valid(parameters, "topdown");		if(temp == "not found"){  temp = "T"; }
			topdown = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "match");		if (temp == "not found"){	temp = "1.0";			}
            util.mothurConvert(temp, match);
            
            temp = validParameter.valid(parameters, "mismatch");		if (temp == "not found"){	temp = "-1.0";			}
            util.mothurConvert(temp, misMatch);
            if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }
            
            temp = validParameter.valid(parameters, "gapopen");		if (temp == "not found"){	temp = "-2.0";			}
            util.mothurConvert(temp, gapOpen);
            if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }
            
            temp = validParameter.valid(parameters, "gapextend");	if (temp == "not found"){	temp = "-1.0";			}
            util.mothurConvert(temp, gapExtend);
            if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }

            align = validParameter.valid(parameters, "align");		if (align == "not found"){	align = "needleman";	}
            
            method = "unaligned";
            
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
/**************************************************************************************************/
int calcMisMatches(string seq1, string seq2, preClusterData* params){
    try {
        int numBad = 0;
        
        if (params->method == "unaligned") {
            //align to eachother
            Sequence seqI("seq1", seq1);
            Sequence seqJ("seq2", seq2);
            
            //align seq2 to seq1 - less abundant to more abundant
            params->alignment->align(seqJ.getUnaligned(), seqI.getUnaligned());
            seq2 = params->alignment->getSeqAAln();
            seq1 = params->alignment->getSeqBAln();
            
            //chop gap ends
            int startPos = 0;
            int endPos = seq2.length()-1;
            for (int i = 0; i < seq2.length(); i++) {  if (isalpha(seq2[i])) { startPos = i; break; } }
            for (int i = seq2.length()-1; i >= 0; i--) {  if (isalpha(seq2[i])) { endPos = i; break; } }
            
            //count number of diffs
            for (int i = startPos; i <= endPos; i++) {
                if (seq2[i] != seq1[i]) { numBad++; }
                if (numBad > params->diffs) { return params->length;  } //to far to cluster
            }
            
        }else {
            //count diffs
            for (int i = 0; i < seq1.length(); i++) {
                //do they match
                if (seq1[i] != seq2[i]) { numBad++; }
                if (numBad > params->diffs) { return params->length;  } //to far to cluster
            }
        }
        return numBad;
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "calcMisMatches");
        exit(1);
    }
}
/**************************************************************************************************/
int process(string newMapFile, preClusterData* params){
    try {
        ofstream out;
        params->util.openOutputFile(newMapFile, out);
        
        //sort seqs by number of identical seqs
        if (params->topdown) { sort(params->alignSeqs.begin(), params->alignSeqs.end(), comparePriorityTopDown);  }
        else {  sort(params->alignSeqs.begin(), params->alignSeqs.end(), comparePriorityDownTop);  }
        
        int count = 0;
        long long numSeqs = params->alignSeqs.size();
        
        if (params->topdown) {
            //think about running through twice...
            for (int i = 0; i < numSeqs; i++) {
                
                if (params->alignSeqs[i].active) {  //this sequence has not been merged yet
                    
                    string chunk = params->alignSeqs[i].seq.getName() + "\t" + toString(params->alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + params->alignSeqs[i].seq.getAligned() + "\n";
                    
                    //try to merge it with all smaller seqs
                    for (int j = i+1; j < numSeqs; j++) {
                        
                        if (params->m->getControl_pressed()) { out.close(); return 0; }
                        
                        if (params->alignSeqs[j].active) {  //this sequence has not been merged yet
                            //are you within "diff" bases
                            int mismatch = params->length;
                            if (params->method == "unaligned") { mismatch = calcMisMatches(params->alignSeqs[i].seq.getAligned(), params->alignSeqs[j].seq.getAligned(), params); }
                            else { mismatch = calcMisMatches(params->alignSeqs[i].filteredSeq, params->alignSeqs[j].filteredSeq, params); }
                            
                            if (mismatch <= params->diffs) {
                                //merge
                                params->alignSeqs[i].names += ',' + params->alignSeqs[j].names;
                                params->alignSeqs[i].numIdentical += params->alignSeqs[j].numIdentical;
                                
                                chunk += params->alignSeqs[j].seq.getName() + "\t" + toString(params->alignSeqs[j].numIdentical) + "\t" + toString(mismatch) + "\t" + params->alignSeqs[j].seq.getAligned() + "\n";
                                
                                params->alignSeqs[j].active = 0;
                                params->alignSeqs[j].numIdentical = 0;
                                params->alignSeqs[j].diffs = mismatch;
                                count++;
                            }
                        }//end if j active
                    }//end for loop j
                    
                    //remove from active list
                    params->alignSeqs[i].active = 0;
                    
                    out << "ideal_seq_" << (i+1) << '\t' << params->alignSeqs[i].numIdentical << endl << chunk << endl;
                    
                }//end if active i
                if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
        }else {
            map<int, string> mapFile;
            map<int, int> originalCount;
            map<int, int>::iterator itCount;
            for (int i = 0; i < numSeqs; i++) { mapFile[i] = ""; originalCount[i] = params->alignSeqs[i].numIdentical; }
            
            //think about running through twice...
            for (int i = 0; i < numSeqs; i++) {
                
                //try to merge it into larger seqs
                for (int j = i+1; j < numSeqs; j++) {
                    
                    if (params->m->getControl_pressed()) { out.close(); return 0; }
                    
                    if (originalCount[j] > originalCount[i]) {  //this sequence is more abundant than I am
                        //are you within "diff" bases
                        int mismatch = params->length;
                        if (params->method == "unaligned") { mismatch = calcMisMatches(params->alignSeqs[i].seq.getAligned(), params->alignSeqs[j].seq.getAligned(), params); }
                        else { mismatch = calcMisMatches(params->alignSeqs[i].filteredSeq, params->alignSeqs[j].filteredSeq, params); }
                        
                        if (mismatch <= params->diffs) {
                            //merge
                            params->alignSeqs[j].names += ',' + params->alignSeqs[i].names;
                            params->alignSeqs[j].numIdentical += params->alignSeqs[i].numIdentical;
                            
                            mapFile[j] = params->alignSeqs[i].seq.getName() + "\t" + toString(params->alignSeqs[i].numIdentical) + "\t" + toString(mismatch) + "\t" + params->alignSeqs[i].seq.getAligned() + "\n" + mapFile[i];
                            params->alignSeqs[i].numIdentical = 0;
                            originalCount.erase(i);
                            mapFile[i] = "";
                            count++;
                            j+=numSeqs; //exit search, we merged this one in.
                        }
                    }//end abundance check
                }//end for loop j
                
                if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
            
            for (int i = 0; i < numSeqs; i++) {
                if (params->alignSeqs[i].numIdentical != 0) {
                    out << "ideal_seq_" << (i+1) << '\t' << params->alignSeqs[i].numIdentical << endl  << params->alignSeqs[i].seq.getName() + "\t" + toString(params->alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + params->alignSeqs[i].seq.getAligned() + "\n" << mapFile[i] << endl;
                }
            }
            
        }
        out.close();
        
        if(numSeqs % 100 != 0)	{ params->m->mothurOut(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count) + "\n"); 	}
        
        return count;
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "process");
        exit(1);
    }
}
/**************************************************************************************************/
void filterSeqs(preClusterData* params){
    try {
        string filterString = "";
        Filters F;
        
        F.setLength(params->length);
        F.initialize();
        F.setFilter(string(params->length, '1'));
        
        for (int i = 0; i < params->alignSeqs.size(); i++) { F.getFreqs(params->alignSeqs[i].seq); }
        
        F.setNumSeqs(params->alignSeqs.size());
        F.doVerticalAllBases();
        filterString = F.getFilter();
        
        //run filter
        for (int i = 0; i < params->alignSeqs.size(); i++) {
            if (params->m->getControl_pressed()) { break; }
            params->alignSeqs[i].filteredSeq = "";
            string align = params->alignSeqs[i].seq.getAligned();
            for(int j=0;j<params->length;j++){ if(filterString[j] == '1'){ params->alignSeqs[i].filteredSeq += align[j]; } }
        }
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "filterSeqs");
        exit(1);
    }
}
/**************************************************************************************************/
long long readFASTA(CountTable ct, preClusterData* params){
    try {
        map<string, string> nameMap;
        map<string, string>::iterator it;
        if (params->hasName) { params->util.readNames(params->namefile, nameMap); }
        
        ifstream inFasta;
        params->util.openInputFile(params->fastafile, inFasta);
        set<int> lengths;
        
        while (!inFasta.eof()) {
            
            if (params->m->getControl_pressed()) { inFasta.close(); return 0; }
            
            Sequence seq(inFasta);  params->util.gobble(inFasta);
            
            if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file
                if (params->hasName) {
                    it = nameMap.find(seq.getName());
                    
                    if (it == nameMap.end()) { params->m->mothurOut("[ERROR]: " + seq.getName() + " is not in your names file, please correct.\n"); exit(1); }
                    else{
                        string second = it->second;
                        int numReps = params->util.getNumNames(second);
                        seqPNode tempNode(numReps, seq, second);
                        params->alignSeqs.push_back(tempNode);
                        lengths.insert(seq.getAligned().length());
                    }
                }else { //no names file, you are identical to yourself
                    int numRep = 1;
                    if (params->hasCount) { numRep = ct.getNumSeqs(seq.getName()); }
                    seqPNode tempNode(numRep, seq, seq.getName());
                    params->alignSeqs.push_back(tempNode);
                    lengths.insert(seq.getAligned().length());
                }
            }
        }
        inFasta.close();
        
        params->length = *(lengths.begin());
        
        if (lengths.size() > 1) { params->method = "unaligned"; }
        else if (lengths.size() == 1) {  params->method = "aligned"; filterSeqs(params); }
        
        return params->alignSeqs.size();
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "readFASTA");
        exit(1);
    }
}
/**************************************************************************************************/
int PreClusterCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		long start = time(NULL);
        
		string fileroot = outputDir + util.getRootName(util.getSimpleName(fastafile));
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
		string newNamesFile = getOutputFileName("name",variables);
        string newCountFile = getOutputFileName("count",variables);
		string newMapFile = getOutputFileName("map",variables); //add group name if by group
        variables["[extension]"] = util.getExtension(fastafile);
		string newFastaFile = getOutputFileName("fasta", variables);
		outputNames.push_back(newFastaFile); outputTypes["fasta"].push_back(newFastaFile);
		if (countfile == "") { outputNames.push_back(newNamesFile); outputTypes["name"].push_back(newNamesFile); }
		else { outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile); }
		
		if (bygroup) {
			//clear out old files
			ofstream outFasta; util.openOutputFile(newFastaFile, outFasta); outFasta.close();
			ofstream outNames; util.openOutputFile(newNamesFile, outNames);  outNames.close();
			newMapFile = fileroot + "precluster.";
			
			createProcessesGroups(newFastaFile, newNamesFile, newMapFile);
			
			if (countfile != "") { 
                mergeGroupCounts(newCountFile, newNamesFile, newFastaFile);
            }else {
                //run unique.seqs for deconvolute results
                string inputString = "fasta=" + newFastaFile;
                if (namefile != "") { inputString += ", name=" + newNamesFile; }
                m->mothurOutEndLine(); 
                m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
                m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
                current->setMothurCalling(true);
                
                Command* uniqueCommand = new DeconvoluteCommand(inputString);
                uniqueCommand->execute();
                
                map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
                
                delete uniqueCommand;
                current->setMothurCalling(false);
                m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
                
                util.renameFile(filenames["fasta"][0], newFastaFile);
                util.renameFile(filenames["name"][0], newNamesFile); 
			}
            if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}	 return 0; }
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run pre.cluster."); m->mothurOutEndLine(); 
				
		}else {
            if (processors != 1) { m->mothurOut("When using running without group information mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }

            preClusterData* params = new preClusterData(fastafile, namefile, groupfile, countfile, NULL, NULL, newMapFile, nullVector);
            params->setVariables(0,0, diffs, topdown, method, align, match, misMatch, gapOpen, gapExtend);
            
            //reads fasta file and return number of seqs
			int numSeqs = readFASTA(ct, params); //fills alignSeqs and makes all seqs active
		
			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }
	
			if (numSeqs == 0) { m->mothurOut("Error reading fasta file...please correct."); m->mothurOutEndLine();  return 0;  }
			if (diffs > length) { m->mothurOut("Error: diffs is greater than your sequence length."); m->mothurOutEndLine();  return 0;  }
			
			int count = process(newMapFile, params);
			outputNames.push_back(newMapFile); outputTypes["map"].push_back(newMapFile);
			
			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }
			
			m->mothurOut("Total number of sequences before precluster was " + toString(params->alignSeqs.size()) + ".\n");
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences.\n\n");
			if (countfile != "") { newNamesFile = newCountFile; }
            print(newFastaFile, newNamesFile, params);
            			
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
		}
				
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}		
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
		m->errorOut(e, "PreClusterCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
void PreClusterCommand::print(string newfasta, string newname, preClusterData* params){
    try {
        ofstream outFasta;
        ofstream outNames;

        util.openOutputFile(newfasta, outFasta);
        util.openOutputFile(newname, outNames);

        if (countfile != "")  { outNames << "Representative_Sequence\ttotal\n";  }
        
        if (countfile != "") {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i].numIdentical != 0) {
                    params->alignSeqs[i].seq.printSequence(outFasta);
                    outNames << params->alignSeqs[i].seq.getName() << '\t' << params->alignSeqs[i].numIdentical << endl;
                }
            }
        }else {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i].numIdentical != 0) {
                    params->alignSeqs[i].seq.printSequence(outFasta);
                    outNames << params->alignSeqs[i].seq.getName() << '\t' << params->alignSeqs[i].names << endl;
                }
            }
        }
        outFasta.close();
        outNames.close();
        
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "printData");
        exit(1);
    }
}
/**************************************************************************************************/
long long loadSeqs(map<string, string>& thisName, vector<Sequence>& thisSeqs, map<string, int>& thisCount, string group, preClusterData* params){
    try {
        set<int> lengths;
        params->alignSeqs.clear();
        bool error = false;
        
        for (int i = 0; i < thisSeqs.size(); i++) {
            
            if (params->m->getControl_pressed()) { return 0; }
            
            if (params->hasName) {
                map<string, string>::iterator it = thisName.find(thisSeqs[i].getName());
                
                //should never be true since parser checks for this
                if (it == thisName.end()) { params->m->mothurOut("[ERROR]: " + thisSeqs[i].getName() + " is not in your names file, please correct.\n"); error = true; }
                else{
                    //get number of reps
                    int numReps = params->util.getNumNames(it->second);
                    seqPNode tempNode(numReps, thisSeqs[i], it->second);
                    params->alignSeqs.push_back(tempNode);
                    lengths.insert(thisSeqs[i].getAligned().length());
                }
            }else { //no names file, you are identical to yourself
                int numRep = 1;
                if (params->hasCount) {
                    map<string, int>::iterator it2 = thisCount.find(thisSeqs[i].getName());
                    
                    //should never be true since parser checks for this
                    if (it2 == thisCount.end()) { params->m->mothurOut("[ERROR]: " + thisSeqs[i].getName() + " is not in your count file, please correct.\n"); error = true; }
                    else { numRep = it2->second;  }
                }
                seqPNode tempNode(numRep, thisSeqs[i], thisSeqs[i].getName());
                params->alignSeqs.push_back(tempNode);
                lengths.insert(thisSeqs[i].getAligned().length());
            }
        }
        
        params->length = *(lengths.begin());
        
        if (lengths.size() > 1) { params->method = "unaligned"; }
        else if (lengths.size() == 1) {  params->method = "aligned"; filterSeqs(params); }
        
        //sanity check
        if (error) { params->m->setControl_pressed(true); }
        
        thisSeqs.clear();
        
        return params->alignSeqs.size();
    }
    
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "loadSeqs");
        exit(1);
    }
}
/**************************************************************************************************/
void printData(string group, preClusterData* params){
    try {
        if ((params->hasCount) && (group == ""))  { params->newNName->write("Representative_Sequence\ttotal\n");  }
        
        if (params->hasCount) {
            if (group != "") {
                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i].numIdentical != 0) {
                        params->alignSeqs[i].seq.printSequence(params->newFName);
                        params->newNName->write(group + '\t' + params->alignSeqs[i].seq.getName() + '\t' + params->alignSeqs[i].names + '\n');
                    }
                }
            }
            else {
                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i].numIdentical != 0) {
                        params->alignSeqs[i].seq.printSequence(params->newFName);
                        params->newNName->write(params->alignSeqs[i].seq.getName()  + '\t' + toString(params->alignSeqs[i].numIdentical) + '\n');
                    }
                }
            }
        }else {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i].numIdentical != 0) {
                    params->alignSeqs[i].seq.printSequence(params->newFName);
                    params->newNName->write(group + '\t' + params->alignSeqs[i].seq.getName() + '\t' + params->alignSeqs[i].names + '\n');
                }
            }
        }
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "printData");
        exit(1);
    }
}
/**************************************************************************************************/
long long driverGroups(preClusterData* params){
	try {
        vector<string> subsetGroups;
        for (int i = params->start; i < params->end; i++) {  subsetGroups.push_back(params->groups[i]);  }
        
        //parse fasta and name file by group
        SequenceCountParser* cparser = NULL;
        SequenceParser* parser = NULL;
        if (params->hasCount) {
            cparser = new SequenceCountParser(params->countfile, params->fastafile, subsetGroups);
        }else {
            if (params->hasName) { parser = new SequenceParser(params->groupfile, params->fastafile, params->namefile, subsetGroups);	}
            else				{ parser = new SequenceParser(params->groupfile, params->fastafile, subsetGroups);                      }
        }
        
		long long numSeqs = 0;
		
		//precluster each group
		for (int i = params->start; i < params->end; i++) {
			if (params->m->getControl_pressed()) { if (params->hasCount) { delete cparser; }else { delete parser; } return numSeqs; }
            
            params->m->mothurOut("\nProcessing group " + params->groups[i] + ":\n");
            
			time_t start = time(NULL);
			map<string, string> thisNameMap;
            vector<Sequence> thisSeqs;
			if (params->groupfile != "")        {  thisSeqs = parser->getSeqs(params->groups[i]);       }
            else if (params->hasCount)          { thisSeqs = cparser->getSeqs(params->groups[i]);       }
			
            if (params->hasName)                {  thisNameMap = parser->getNameMap(params->groups[i]); }
            
            map<string, int> thisCount;
            if (params->hasCount) { thisCount = cparser->getCountTable(params->groups[i]);  }
			numSeqs = loadSeqs(thisNameMap, thisSeqs, thisCount, params->groups[i], params);
			
			if (params->m->getControl_pressed()) {   return 0; }
			
            if (params->method == "aligned") { if (params->diffs > params->length) { params->m->mothurOut("[ERROR]: diffs is greater than your sequence length.\n"); params->m->setControl_pressed(true); return 0;  } }
			
            string extension = params->groups[i]+".map";
			long long count = process(params->newMName+extension, params);
			params->outputNames.push_back(params->newMName+extension); params->outputTypes["map"].push_back(params->newMName+extension);
			
			if (params->m->getControl_pressed()) {  return 0; }
			
			params->m->mothurOut("Total number of sequences before pre.cluster was " + toString(params->alignSeqs.size()) + ".\n");
			params->m->mothurOut("pre.cluster removed " + toString(count) + " sequences.\n\n");
			printData(params->groups[i], params);
			
			params->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(count) + " sequences.\n");
		}
		
        if (params->hasCount) { delete cparser; }else { delete parser; }
        
		return numSeqs;
	}
	catch(exception& e) {
		params->m->errorOut(e, "PreClusterCommand", "driverGroups");
		exit(1);
	}
}
/**************************************************************************************************/

int PreClusterCommand::mergeGroupCounts(string newcount, string newname, string newfasta){
	try {
		ifstream inNames;
        util.openInputFile(newname, inNames);
        
        string group, first, second;
        set<string> uniqueNames;
        while (!inNames.eof()) {
            if (m->getControl_pressed()) { break; }
            inNames >> group; util.gobble(inNames);
            inNames >> first; util.gobble(inNames);
            inNames >> second; util.gobble(inNames);
            
            vector<string> names;
            util.splitAtComma(second, names);
            
            uniqueNames.insert(first);
            
            int total = ct.getGroupCount(first, group);
            for (int i = 1; i < names.size(); i++) {
                total += ct.getGroupCount(names[i], group);
                ct.setAbund(names[i], group, 0);
            }
            ct.setAbund(first, group, total);
        }
        inNames.close();
        
        vector<string> namesOfSeqs = ct.getNamesOfSeqs();
        for (int i = 0; i < namesOfSeqs.size(); i++) {
            if (ct.getNumSeqs(namesOfSeqs[i]) == 0) {
                ct.remove(namesOfSeqs[i]);
            }
        }
        
        ct.printTable(newcount); 
        util.mothurRemove(newname);
        
        if (bygroup) { //if by group, must remove the duplicate seqs that are named the same
            ifstream in;
            util.openInputFile(newfasta, in);
            
            ofstream out;
            util.openOutputFile(newfasta+"temp", out);
            
            int count = 0;
            set<string> already;
            while(!in.eof()) {
                if (m->getControl_pressed()) { break; }
                
                Sequence seq(in); util.gobble(in);
                
                if (seq.getName() != "") {
                    count++;
                    if (already.count(seq.getName()) == 0) {
                        seq.printSequence(out);
                        already.insert(seq.getName());
                    }
                }
            }
            in.close();
            out.close();
            util.mothurRemove(newfasta);
            util.renameFile(newfasta+"temp", newfasta);
        }
		        return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "mergeGroupCounts");
		exit(1);
	}
}

/**************************************************************************************************/
void PreClusterCommand::createProcessesGroups(string newFName, string newNName, string newMFile) {
    try {
        //parse fasta and name file by group
        vector<string> groups;
        if (countfile != "") { CountTable ct; ct.testGroups(countfile, groups); }
        else { GroupMap gp; gp.readMap(groupfile); groups = gp.getNamesOfGroups(); }
        
        //sanity check
        if (groups.size() < processors) { processors = groups.size(); m->mothurOut("Reducing processors to " + toString(groups.size()) + ".\n"); }
        
        //divide the groups between the processors
        vector<linePair> lines;
        int remainingPairs = groups.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<preClusterData*> data;
        
        auto synchronizedFastaFile = std::make_shared<SynchronizedOutputFile>(newFName);
        auto synchronizedNameFile = std::make_shared<SynchronizedOutputFile>(newNName);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadFastaWriter = new OutputWriter(synchronizedFastaFile);
            OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);
            
            preClusterData* dataBundle = new preClusterData(fastafile, namefile, groupfile, countfile, threadFastaWriter, threadNameWriter, newMFile, groups);
            dataBundle->setVariables(lines[i+1].start, lines[i+1].end, diffs, topdown, method, align, match, misMatch, gapOpen, gapExtend);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverGroups, dataBundle));
        }
        
        OutputWriter* threadFastaWriter = new OutputWriter(synchronizedFastaFile);
        OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);
        
        preClusterData* dataBundle = new preClusterData(fastafile, namefile, groupfile, countfile, threadFastaWriter, threadNameWriter, newMFile, groups);
        dataBundle->setVariables(lines[0].start, lines[0].end, diffs, topdown, method, align, match, misMatch, gapOpen, gapExtend);

        driverGroups(dataBundle);
        
        outputNames.insert(outputNames.end(), dataBundle->outputNames.begin(), dataBundle->outputNames.end());
        outputTypes.insert(dataBundle->outputTypes.begin(), dataBundle->outputTypes.end());
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            delete data[i]->newFName;
            delete data[i]->newNName;
            
            outputNames.insert(outputNames.end(), data[i]->outputNames.begin(), data[i]->outputNames.end());
            outputTypes.insert(data[i]->outputTypes.begin(), data[i]->outputTypes.end());
            
            delete data[i];
            delete workerThreads[i];
        }
        delete threadFastaWriter;
        delete threadNameWriter;
        delete dataBundle;

    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "createProcessesGroups");
        exit(1);
    }
}
/**************************************************************************************************/


