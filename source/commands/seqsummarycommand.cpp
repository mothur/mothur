/*
 *  seqcoordcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqsummarycommand.h"
#include "counttable.h"
#include "summary.hpp"

//**********************************************************************************************************************
vector<string> SeqSummaryCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pfasta);
        CommandParameter psummary("summary", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(psummary);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SeqSummaryCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.seqs command reads a fastafile or summary file and summarizes the sequences.\n";
		helpString += "The summary.seqs command parameters are fasta, name, count, summary and processors, fasta or summary is required, unless you have a valid current fasta file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your fasta file. \n";
        helpString += "The count parameter allows you to enter a count file associated with your fasta file. \n";
		helpString += "The summary.seqs command should be in the following format: \n";
		helpString += "summary.seqs(fasta=yourFastaFile, processors=2) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SeqSummaryCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SeqSummaryCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
SeqSummaryCommand::SeqSummaryCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("summary.seqs");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
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
                
                it = parameters.find("summary");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["summary"] = inputDir + it->second;		}
                }
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") {  fastafile = "";  }
            else { m->setFastaFile(fastafile); }
			
            summaryfile = validParameter.validFile(parameters, "summary", true);
            if (summaryfile == "not open") { abort = true; }
            else if (summaryfile == "not found") {  summaryfile = "";  }
            else { m->setSummaryFile(summaryfile); }

            if ((summaryfile == "") && (fastafile == "")) {
                fastafile = m->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                else {
                    summaryfile = m->getSummaryFile();
                    if (summaryfile != "") { m->mothurOut("Using " + summaryfile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                    else { 	m->mothurOut("You have no current fastafile or summaryfile, one is required."); m->mothurOutEndLine(); abort = true; }
                }
            }
            
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    parser.getNameFile(files);
                }
            }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SeqSummaryCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        int start = time(NULL);
        
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string outputFile = getOutputFileName("summary",variables);
        
        string nameOrCount = countfile;
        if (namefile != "") { nameOrCount = namefile; }
        
        Summary sum;
        if (fastafile != "") {  sum.summarizeFasta(fastafile, nameOrCount, outputFile);  }
        else if (summaryfile != "") {  sum.summarizeFastaSummary(summaryfile, nameOrCount);  }
        
        vector<long long> starts = sum.getStart();
        vector<long long> ends = sum.getEnd();
        vector<long long> ambigs = sum.getAmbig();
        vector<long long> lengths = sum.getLength();
        vector<long long> homops = sum.getHomop();
        vector <long long> ptiles = sum.getDefaults();
        
		long long numSeqs = 0;
        long long size = sum.getTotalSeqs();
        long long numUniques = sum.getUniqueSeqs();
        
        if (m->control_pressed) {  m->mothurRemove(outputFile); return 0; }
        
        m->mothurOutEndLine();
        m->mothurOut("\t\tStart\tEnd\tNBases\tAmbigs\tPolymer\tNumSeqs"); m->mothurOutEndLine();
        m->mothurOut("Minimum:\t" + toString(starts[0]) + "\t" + toString(ends[0]) + "\t" + toString(lengths[0]) + "\t" + toString(ambigs[0]) + "\t" + toString(homops[0]) + "\t" + toString(ptiles[0])); m->mothurOutEndLine();
        m->mothurOut("2.5%-tile:\t" + toString(starts[1]) + "\t" + toString(ends[1]) + "\t" + toString(lengths[1]) + "\t" + toString(ambigs[1]) + "\t" + toString(homops[1]) + "\t" + toString(ptiles[1])); m->mothurOutEndLine();
        m->mothurOut("25%-tile:\t" + toString(starts[2]) + "\t" + toString(ends[2]) + "\t" + toString(lengths[2]) + "\t" + toString(ambigs[2]) + "\t" + toString(homops[2]) + "\t" + toString(ptiles[2])); m->mothurOutEndLine();
        m->mothurOut("Median: \t" + toString(starts[3]) + "\t" + toString(ends[3]) + "\t" + toString(lengths[3]) + "\t" + toString(ambigs[3]) + "\t" + toString(homops[3]) + "\t" + toString(ptiles[3])); m->mothurOutEndLine();
        m->mothurOut("75%-tile:\t" + toString(starts[4]) + "\t" + toString(ends[4]) + "\t" + toString(lengths[4]) + "\t" + toString(ambigs[4]) + "\t" + toString(homops[4]) + "\t" + toString(ptiles[4])); m->mothurOutEndLine();
        m->mothurOut("97.5%-tile:\t" + toString(starts[5]) + "\t" + toString(ends[5]) + "\t" + toString(lengths[5]) + "\t" + toString(ambigs[5]) + "\t" + toString(homops[5]) + "\t" + toString(ptiles[5])); m->mothurOutEndLine();
        m->mothurOut("Maximum:\t" + toString(starts[6]) + "\t" + toString(ends[6]) + "\t" + toString(lengths[6]) + "\t" + toString(ambigs[6]) + "\t" + toString(homops[6]) + "\t" + toString(ptiles[6])); m->mothurOutEndLine();
        m->mothurOut("Mean:\t" + toString(starts[7]) + "\t" + toString(ends[7]) + "\t" + toString(lengths[7]) + "\t" + toString(ambigs[7]) + "\t" + toString(homops[7])); m->mothurOutEndLine();
        if ((namefile == "") && (countfile == "") && (summaryfile == "")) {  m->mothurOut("# of Seqs:\t" + toString(numUniques)); m->mothurOutEndLine(); }
        else { m->mothurOut("# of unique seqs:\t" + toString(numUniques)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(size)); m->mothurOutEndLine(); }
        
        if (m->control_pressed) {  m->mothurRemove(outputFile); return 0; }
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        if (summaryfile == "") { m->mothurOut(outputFile); m->mothurOutEndLine();	outputNames.push_back(outputFile); outputTypes["summary"].push_back(outputFile);  }
        m->mothurOutEndLine();
        
        if ((namefile == "") && (countfile == "") && (summaryfile == "")) {  m->mothurOut("It took " + toString(time(NULL) - start) + " secs to summarize " + toString(numSeqs) + " sequences.\n");  }
        else{  m->mothurOut("It took " + toString(time(NULL) - start) + " secs to summarize " + toString(size) + " sequences.\n");   }
        
        //set fasta file as new current fastafile
        string current = "";
        itTypes = outputTypes.find("summary");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSummaryFile(current); }
        }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
