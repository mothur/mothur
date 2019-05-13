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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FastaReport", "none", "none","summary",false,true,true); parameters.push_back(pfasta);
        CommandParameter psummary("summary", "InputTypes", "", "", "FastaReport", "none", "none","",false,false,true); parameters.push_back(psummary);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "FastaReport", "none", "none","",false,false,true); parameters.push_back(pcontigsreport);
        CommandParameter palignreport("alignreport", "InputTypes", "", "", "FastaReport", "none", "none","",false,false,true); parameters.push_back(palignreport);
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
		helpString += "The summary.seqs command reads a fastafile, summary, contigsreport or alignreport file and summarizes it.\n";
		helpString += "The summary.seqs command parameters are fasta, name, count, summary, contigsreport, alignreport and processors, fasta, contigsreport, alignreport or summary is required, unless you have a valid current files.\n";
		helpString += "The name parameter allows you to enter a name file associated with your fasta file. \n";
        helpString += "The count parameter allows you to enter a count file associated with your fasta file. \n";
		helpString += "The summary.seqs command should be in the following format: \n";
		helpString += "summary.seqs(fasta=yourFastaFile, processors=2) \n";
			
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
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
                
                it = parameters.find("summary");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["summary"] = inputDir + it->second;		}
                }
                
                it = parameters.find("contigsreport");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["contigsreport"] = inputDir + it->second;		}
                }
				
                it = parameters.find("alignreport");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
                }
                
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") {  fastafile = "";  }
            else { current->setFastaFile(fastafile); }
			
            summaryfile = validParameter.validFile(parameters, "summary");
            if (summaryfile == "not open") { abort = true; }
            else if (summaryfile == "not found") {  summaryfile = "";  }
            else { current->setSummaryFile(summaryfile); }
            
            contigsfile = validParameter.validFile(parameters, "contigsreport");
            if (contigsfile == "not open") { abort = true; }
            else if (contigsfile == "not found") {  contigsfile = "";  }
            else { current->setContigsReportFile(contigsfile); }
            
            alignfile = validParameter.validFile(parameters, "alignreport");
            if (alignfile == "not open") { abort = true; }
            else if (alignfile == "not found") {  alignfile = "";  }

            if ((summaryfile == "") && (fastafile == "") && (contigsfile == "") && (alignfile == "")) {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                else {
                    summaryfile = current->getSummaryFile();
                    if (summaryfile != "") { m->mothurOut("Using " + summaryfile + " as input file for the summary parameter."); m->mothurOutEndLine(); }
                    else {
                        contigsfile = current->getContigsReportFile();
                        if (contigsfile != "") { m->mothurOut("Using " + contigsfile + " as input file for the contigsreport parameter."); m->mothurOutEndLine(); }
                        else { 	m->mothurOut("You have no current fasta, summary, contigsreport or alignreport file, one is required."); m->mothurOutEndLine(); abort = true; }
                    }
                }
            }
            
            if (((fastafile != "") && ((summaryfile != "") || (contigsfile != "") || (alignfile != ""))) ||
                ((summaryfile != "") && ((fastafile != "") || (contigsfile != "") || (alignfile != ""))) ||
                ((contigsfile != "") && ((summaryfile != "") || (fastafile != "") || (alignfile != ""))) ||
                ((alignfile != "") && ((summaryfile != "") || (contigsfile != "") || (fastafile != "")))) {
                m->mothurOut("[ERROR]: you may only use one of the following: fasta, summary, contigsreport or alignreport."); m->mothurOutEndLine(); abort = true;
            }
            
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += util.hasPath(fastafile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        long start = time(NULL);
        
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
		string outputFile = getOutputFileName("summary",variables);
        
        string nameOrCount = countfile;
        if (namefile != "") { nameOrCount = namefile; }
        
        Summary sum(processors);
        if (fastafile != "") {  sum.summarizeFasta(fastafile, nameOrCount, outputFile);  }
        else if (summaryfile != "") {  sum.summarizeFastaSummary(summaryfile, nameOrCount);  }
        else if (contigsfile != "") {  sum.summarizeContigsSummary(contigsfile, nameOrCount);  }
        else if (alignfile != "") {  sum.summarizeAlignSummary(alignfile, nameOrCount);  }
        else { m->mothurOut("[ERROR]: Unknown type: you may only use one of the following: fasta, summary, contigsreport or alignreport."); m->mothurOutEndLine(); m->setControl_pressed(true); }

        if (m->getControl_pressed()) {  util.mothurRemove(outputFile); return 0; }
        
        long long size = sum.getTotalSeqs();
        long long numUniques = sum.getUniqueSeqs();

        vector <long long> ptiles = sum.getDefaults();

        if ((fastafile != "") || (summaryfile != "")) {
            vector<long long> starts = sum.getStart();
            vector<long long> ends = sum.getEnd();
            vector<long long> ambigs = sum.getAmbig();
            vector<long long> lengths = sum.getLength();
            vector<long long> homops = sum.getHomop();
            
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
        }else if (contigsfile != "") {
            vector<long long> ostarts = sum.getOStart();
            vector<long long> oends = sum.getOEnd();
            vector<long long> length = sum.getLength();
            vector<long long> olengths = sum.getOLength();
            vector<long long> numns = sum.getNumNs();
            vector<long long> mismatches = sum.getMisMatches();
            
            m->mothurOutEndLine();
            m->mothurOut("\t\tLength\tOverlap_Length\tOverlap_Start\tOverlap_End\tMisMatches\tNum_Ns\tNumSeqs"); m->mothurOutEndLine();
            m->mothurOut("Minimum:\t" + toString(length[0]) + "\t" + toString(olengths[0]) + "\t" + toString(ostarts[0]) + "\t" + toString(oends[0]) + "\t" + toString(mismatches[0]) + "\t" + toString(numns[0]) + "\t" + toString(ptiles[0])); m->mothurOutEndLine();
            m->mothurOut("2.5%-tile:\t" + toString(length[1]) + "\t" + toString(olengths[1]) + "\t" + toString(ostarts[1]) + "\t" + toString(oends[1]) + "\t" + toString(mismatches[1]) + "\t" + toString(numns[1]) + "\t" + toString(ptiles[1])); m->mothurOutEndLine();
            m->mothurOut("25%-tile:\t" + toString(length[2]) + "\t" + toString(olengths[2]) + "\t" + toString(ostarts[2]) + "\t" + toString(oends[2]) + "\t" + toString(mismatches[2]) + "\t" + toString(numns[2]) + "\t" + toString(ptiles[2])); m->mothurOutEndLine();
            m->mothurOut("Median: \t" + toString(length[3]) + "\t" + toString(olengths[3]) + "\t" + toString(ostarts[3]) + "\t" + toString(oends[3]) + "\t" + toString(mismatches[3]) + "\t" + toString(numns[3]) + "\t" + toString(ptiles[3])); m->mothurOutEndLine();
            m->mothurOut("75%-tile:\t" + toString(length[4]) + "\t" + toString(olengths[4]) + "\t" + toString(ostarts[4]) + "\t" + toString(oends[4]) + "\t" + toString(mismatches[4]) + "\t" + toString(numns[4])  + "\t" + toString(ptiles[4])); m->mothurOutEndLine();
            m->mothurOut("97.5%-tile:\t" + toString(length[5]) + "\t" + toString(olengths[5]) + "\t" + toString(ostarts[5]) + "\t" + toString(oends[5]) + "\t" + toString(mismatches[5]) + "\t" + toString(numns[5])  + "\t" + toString(ptiles[5])); m->mothurOutEndLine();
            m->mothurOut("Maximum:\t" + toString(length[6]) + "\t" + toString(olengths[6]) + "\t" + toString(ostarts[6]) + "\t" + toString(oends[6]) + "\t" + toString(mismatches[6]) + "\t" + toString(numns[6]) + "\t" + toString(ptiles[6])); m->mothurOutEndLine();
            m->mothurOut("Mean:\t" + toString(length[7]) + "\t" + toString(olengths[7]) + "\t" + toString(ostarts[7]) + "\t" + toString(oends[7]) + "\t" + toString(mismatches[7]) + "\t" + toString(numns[7]) ); m->mothurOutEndLine();
        }else if (alignfile != "") {
            vector<long long> sims = sum.getSims();
            vector<long long> scores = sum.getScores();
            vector<long long> inserts = sum.getNumInserts();
            vector<long long> length = sum.getLength();
            
            m->mothurOutEndLine();
            m->mothurOut("\t\tLength\tSimBtwnQueryTemplate\tLongestInsert\tSearchScore\tNumSeqs"); m->mothurOutEndLine();
            m->mothurOut("Minimum:\t" + toString(length[0]) + "\t" + toString(sims[0]) + "\t" + toString(inserts[0]) + "\t" + toString(scores[0]) + "\t" +  toString(ptiles[0])); m->mothurOutEndLine();
            m->mothurOut("2.5%-tile:\t" + toString(length[1]) + "\t" + toString(sims[1]) + "\t" + toString(inserts[1]) + "\t" + toString(scores[1]) + "\t" + toString(ptiles[1])); m->mothurOutEndLine();
            m->mothurOut("25%-tile:\t" + toString(length[2]) + "\t" + toString(sims[2]) + "\t" + toString(inserts[2]) + "\t" + toString(scores[2]) + "\t" + toString(ptiles[2])); m->mothurOutEndLine();
            m->mothurOut("Median: \t" + toString(length[3]) + "\t" + toString(sims[3]) + "\t" + toString(inserts[3]) + "\t" + toString(scores[3]) + "\t" + toString(ptiles[3])); m->mothurOutEndLine();
            m->mothurOut("75%-tile:\t" + toString(length[4]) + "\t" + toString(sims[4]) + "\t" + toString(inserts[4]) + "\t" + toString(scores[4]) + "\t" + toString(ptiles[4])); m->mothurOutEndLine();
            m->mothurOut("97.5%-tile:\t" + toString(length[5]) + "\t" + toString(sims[5]) + "\t" + toString(inserts[5]) + "\t" + toString(scores[5]) + "\t" +  toString(ptiles[5])); m->mothurOutEndLine();
            m->mothurOut("Maximum:\t" + toString(length[6]) + "\t" + toString(sims[6]) + "\t" + toString(inserts[6]) + "\t" + toString(scores[6]) + "\t" +  toString(ptiles[6])); m->mothurOutEndLine();
            m->mothurOut("Mean:\t" + toString(length[7]) + "\t" + toString(sims[7]) + "\t" + toString(inserts[7]) + "\t" + toString(scores[7])); m->mothurOutEndLine();
        }
        
        if (m->getControl_pressed()) {  util.mothurRemove(outputFile); return 0; }
        
        if ((namefile == "") && (countfile == "") && (summaryfile == "")) {  m->mothurOut("# of Seqs:\t" + toString(numUniques)); m->mothurOutEndLine(); }
        else { m->mothurOut("# of unique seqs:\t" + toString(numUniques)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(size)); m->mothurOutEndLine(); }
        
        if (((namefile == "") && (countfile == "")) && (summaryfile == "")) {  m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to summarize " + toString(numUniques) + " sequences.\n");  }
        else{  m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to summarize " + toString(size) + " sequences.\n");   }
        
        m->mothurOut("\nOutput File Names:\n");
        if ((summaryfile == "") && (contigsfile == "") && (alignfile == "")) {
            m->mothurOut(outputFile); m->mothurOutEndLine();
            outputNames.push_back(outputFile); outputTypes["summary"].push_back(outputFile);
            //set fasta file as new current fastafile
            string currentName = "";
            itTypes = outputTypes.find("summary");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSummaryFile(currentName); }
            }
        }
        m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
