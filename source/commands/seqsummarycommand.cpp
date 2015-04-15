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

//**********************************************************************************************************************
vector<string> SeqSummaryCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pfasta);
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
		helpString += "The summary.seqs command reads a fastafile and summarizes the sequences.\n";
		helpString += "The summary.seqs command parameters are fasta, name, count and processors, fasta is required, unless you have a valid current fasta file.\n";
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
			else if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }	
			
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
        
		//set current fasta to fastafile
		m->setFastaFile(fastafile);
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string summaryFile = getOutputFileName("summary",variables);
				
		long long numSeqs = 0;
        long long size = 0;
		map<int, long long> startPosition;
		map<int, long long> endPosition;
		map<int, long long> seqLength;
		map<int, long long> ambigBases;
		map<int, long long> longHomoPolymer;
		
		if (namefile != "") { nameMap = m->readNames(namefile); }
        else if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, false, false);
            nameMap = ct.getNameMap();
            size = ct.getNumSeqs();
        }
		
		if (m->control_pressed) { return 0; }
			

			vector<unsigned long long> positions; 
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				positions = m->divideFile(fastafile, processors);
				for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
			#else
				positions = m->setFilePosFasta(fastafile, numSeqs); 
                if (positions.size() < processors) { processors = positions.size(); }
		
				//figure out how many sequences you have to process
				int numSeqsPerProcessor = numSeqs / processors;
				for (int i = 0; i < processors; i++) {
					int startIndex =  i * numSeqsPerProcessor;
					if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
					lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
				}
			#endif
			

			if(processors == 1){
				numSeqs = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
			}else{
				numSeqs = createProcessesCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile); 
			}
			
			if (m->control_pressed) {  return 0; }
			
		
        
        //set size
        if (countfile != "") {}//already set
        else if (namefile == "") { size = numSeqs;  }
        else { for (map<int, long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++) { size += it->second; } }
        
        long long ptile0_25	= 1+(long long)(size * 0.025); //number of sequences at 2.5%
        long long ptile25		= 1+(long long)(size * 0.250); //number of sequences at 25%
        long long ptile50		= 1+(long long)(size * 0.500);
        long long ptile75		= 1+(long long)(size * 0.750);
        long long ptile97_5	= 1+(long long)(size * 0.975);
        long long ptile100	= (long long)(size);
        vector<int> starts; starts.resize(7,0); vector<int> ends; ends.resize(7,0); vector<int> ambigs; ambigs.resize(7,0); vector<int> lengths; lengths.resize(7,0); vector<int> homops; homops.resize(7,0);
        
		//find means
		long long meanStartPosition, meanEndPosition, meanSeqLength, meanAmbigBases, meanLongHomoPolymer;
        meanStartPosition = 0; meanEndPosition = 0; meanSeqLength = 0; meanAmbigBases = 0; meanLongHomoPolymer = 0;
        //minimum
        if ((startPosition.begin())->first == -1) { starts[0] = 0; }
        else {starts[0] = (startPosition.begin())->first; }
        long long totalSoFar = 0;
        starts[0] = (startPosition.begin())->first;
        int lastValue = 0;
        for (map<int, long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanStartPosition += (value*it->second);
            totalSoFar += it->second;
            if (((totalSoFar <= ptile0_25) && (totalSoFar > 1)) || ((lastValue < ptile0_25) && (totalSoFar > ptile0_25))){  starts[1] = value;   } //save value
            if (((totalSoFar <= ptile25) && (totalSoFar > ptile0_25)) ||  ((lastValue < ptile25) && (totalSoFar > ptile25))) { starts[2] = value;  } //save value
            if (((totalSoFar <= ptile50) && (totalSoFar > ptile25)) ||  ((lastValue < ptile50) && (totalSoFar > ptile50))) {  starts[3] = value; } //save value
            if (((totalSoFar <= ptile75) && (totalSoFar > ptile50)) ||  ((lastValue < ptile75) && (totalSoFar > ptile75))) {  starts[4] = value; } //save value
            if (((totalSoFar <= ptile97_5) && (totalSoFar > ptile75)) ||  ((lastValue < ptile97_5) && (totalSoFar > ptile97_5))) {  starts[5] = value;  } //save value
            if ((totalSoFar <= ptile100) && (totalSoFar > ptile97_5)) {  starts[6] = value; } //save value
            lastValue = totalSoFar;
        }
        starts[6] = (startPosition.rbegin())->first;
        
        if ((endPosition.begin())->first == -1) { ends[0] = 0; }
        else {ends[0] = (endPosition.begin())->first; }
        totalSoFar = 0;
        ends[0] = (endPosition.begin())->first;
        lastValue = 0;
        for (map<int, long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanEndPosition += (value*it->second);
            totalSoFar += it->second;
            
            if (((totalSoFar <= ptile0_25) && (totalSoFar > 1)) || ((lastValue < ptile0_25) && (totalSoFar > ptile0_25))){  ends[1] = value;  } //save value
            if (((totalSoFar <= ptile25) && (totalSoFar > ptile0_25)) ||  ((lastValue < ptile25) && (totalSoFar > ptile25))) { ends[2] = value;  } //save value
            if (((totalSoFar <= ptile50) && (totalSoFar > ptile25)) ||  ((lastValue < ptile50) && (totalSoFar > ptile50))) {  ends[3] = value; } //save value
            if (((totalSoFar <= ptile75) && (totalSoFar > ptile50)) ||  ((lastValue < ptile75) && (totalSoFar > ptile75))) {  ends[4] = value; } //save value
            if (((totalSoFar <= ptile97_5) && (totalSoFar > ptile75)) ||  ((lastValue < ptile97_5) && (totalSoFar > ptile97_5))) {  ends[5] = value;  } //save value
            if ((totalSoFar <= ptile100) && (totalSoFar > ptile97_5)) {   ends[6] = value; } //save value
            lastValue = totalSoFar;
        }
        ends[6] = (endPosition.rbegin())->first;
        
        lengths[0] = (seqLength.begin())->first;
        totalSoFar = 0;
        lastValue = 0;
        for (map<int, long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++) {
            int value = it->first;
            meanSeqLength += (value*it->second);
            totalSoFar += it->second;
            
            if (((totalSoFar <= ptile0_25) && (totalSoFar > 1)) || ((lastValue < ptile0_25) && (totalSoFar > ptile0_25))){  lengths[1] = value;  } //save value
            if (((totalSoFar <= ptile25) && (totalSoFar > ptile0_25)) ||  ((lastValue < ptile25) && (totalSoFar > ptile25))) {   lengths[2] = value;  } //save value
            if (((totalSoFar <= ptile50) && (totalSoFar > ptile25)) ||  ((lastValue < ptile50) && (totalSoFar > ptile50))) {  lengths[3] = value; } //save value
            if (((totalSoFar <= ptile75) && (totalSoFar > ptile50)) ||  ((lastValue < ptile75) && (totalSoFar > ptile75))) {  lengths[4] = value; } //save value
            if (((totalSoFar <= ptile97_5) && (totalSoFar > ptile75)) ||  ((lastValue < ptile97_5) && (totalSoFar > ptile97_5))) {  lengths[5] = value;  } //save value
            if ((totalSoFar <= ptile100) && (totalSoFar > ptile97_5)) {  lengths[6] = value; } //save value
            lastValue = totalSoFar;
        }
        lengths[6] = (seqLength.rbegin())->first;
                
        ambigs[0] = (ambigBases.begin())->first;
        totalSoFar = 0;
        lastValue = 0;
        for (map<int, long long>::iterator it = ambigBases.begin(); it != ambigBases.end(); it++) {
            int value = it->first;
            meanAmbigBases += (value*it->second);
            totalSoFar += it->second;
            
            if (((totalSoFar <= ptile0_25) && (totalSoFar > 1)) || ((lastValue < ptile0_25) && (totalSoFar > ptile0_25))){  ambigs[1] = value;  } //save value
            if (((totalSoFar <= ptile25) && (totalSoFar > ptile0_25)) ||  ((lastValue < ptile25) && (totalSoFar > ptile25))) {   ambigs[2] = value;  } //save value
            if (((totalSoFar <= ptile50) && (totalSoFar > ptile25)) ||  ((lastValue < ptile50) && (totalSoFar > ptile50))) {  ambigs[3] = value; } //save value
            if (((totalSoFar <= ptile75) && (totalSoFar > ptile50)) ||  ((lastValue < ptile75) && (totalSoFar > ptile75))) {  ambigs[4] = value; } //save value
            if (((totalSoFar <= ptile97_5) && (totalSoFar > ptile75)) ||  ((lastValue < ptile97_5) && (totalSoFar > ptile97_5))) {  ambigs[5] = value;  } //save value
            if ((totalSoFar <= ptile100) && (totalSoFar > ptile97_5)) {  ambigs[6] = value; } //save value
            lastValue = totalSoFar;
        }
        ambigs[6] = (ambigBases.rbegin())->first;
        
        homops[0] = (longHomoPolymer.begin())->first;
        totalSoFar = 0;
        lastValue = 0;
        for (map<int, long long>::iterator it = longHomoPolymer.begin(); it != longHomoPolymer.end(); it++) {
            int value = it->first;
            meanLongHomoPolymer += (it->first*it->second);
            totalSoFar += it->second;
            
            if (((totalSoFar <= ptile0_25) && (totalSoFar > 1)) || ((lastValue < ptile0_25) && (totalSoFar > ptile0_25))){  homops[1] = value;  } //save value
            if (((totalSoFar <= ptile25) && (totalSoFar > ptile0_25)) ||  ((lastValue < ptile25) && (totalSoFar > ptile25))) {   homops[2] = value;  } //save value
            if (((totalSoFar <= ptile50) && (totalSoFar > ptile25)) ||  ((lastValue < ptile50) && (totalSoFar > ptile50))) {  homops[3] = value; } //save value
            if (((totalSoFar <= ptile75) && (totalSoFar > ptile50)) ||  ((lastValue < ptile75) && (totalSoFar > ptile75))) {  homops[4] = value; } //save value
            if (((totalSoFar <= ptile97_5) && (totalSoFar > ptile75)) ||  ((lastValue < ptile97_5) && (totalSoFar > ptile97_5))) {  homops[5] = value;  } //save value
            if ((totalSoFar <= ptile100) && (totalSoFar > ptile97_5)) {  homops[6] = value; } //save value
            lastValue = totalSoFar;
        }
        homops[6] = (longHomoPolymer.rbegin())->first;
        		      
        double meanstartPosition, meanendPosition, meanseqLength, meanambigBases, meanlongHomoPolymer;
                
		meanstartPosition = meanStartPosition / (double) size; meanendPosition = meanEndPosition /(double) size; meanlongHomoPolymer = meanLongHomoPolymer / (double) size; meanseqLength = meanSeqLength / (double) size; meanambigBases = meanAmbigBases /(double) size;
		
		//to compensate for blank sequences that would result in startPosition and endPostion equalling -1
		if (startPosition[0] == -1) {  startPosition[0] = 0;	}
		if (endPosition[0] == -1)	{  endPosition[0] = 0;		}
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("\t\tStart\tEnd\tNBases\tAmbigs\tPolymer\tNumSeqs"); m->mothurOutEndLine();
		m->mothurOut("Minimum:\t" + toString(starts[0]) + "\t" + toString(ends[0]) + "\t" + toString(lengths[0]) + "\t" + toString(ambigs[0]) + "\t" + toString(homops[0]) + "\t" + toString(1)); m->mothurOutEndLine();
		m->mothurOut("2.5%-tile:\t" + toString(starts[1]) + "\t" + toString(ends[1]) + "\t" + toString(lengths[1]) + "\t" + toString(ambigs[1]) + "\t" + toString(homops[1]) + "\t" + toString(ptile0_25)); m->mothurOutEndLine();
		m->mothurOut("25%-tile:\t" + toString(starts[2]) + "\t" + toString(ends[2]) + "\t" + toString(lengths[2]) + "\t" + toString(ambigs[2]) + "\t" + toString(homops[2]) + "\t" + toString(ptile25)); m->mothurOutEndLine();
		m->mothurOut("Median: \t" + toString(starts[3]) + "\t" + toString(ends[3]) + "\t" + toString(lengths[3]) + "\t" + toString(ambigs[3]) + "\t" + toString(homops[3]) + "\t" + toString(ptile50)); m->mothurOutEndLine();
		m->mothurOut("75%-tile:\t" + toString(starts[4]) + "\t" + toString(ends[4]) + "\t" + toString(lengths[4]) + "\t" + toString(ambigs[4]) + "\t" + toString(homops[4]) + "\t" + toString(ptile75)); m->mothurOutEndLine();
		m->mothurOut("97.5%-tile:\t" + toString(starts[5]) + "\t" + toString(ends[5]) + "\t" + toString(lengths[5]) + "\t" + toString(ambigs[5]) + "\t" + toString(homops[5]) + "\t" + toString(ptile97_5)); m->mothurOutEndLine();
		m->mothurOut("Maximum:\t" + toString(starts[6]) + "\t" + toString(ends[6]) + "\t" + toString(lengths[6]) + "\t" + toString(ambigs[6]) + "\t" + toString(homops[6]) + "\t" + toString(ptile100)); m->mothurOutEndLine();
		m->mothurOut("Mean:\t" + toString(meanstartPosition) + "\t" + toString(meanendPosition) + "\t" + toString(meanseqLength) + "\t" + toString(meanambigBases) + "\t" + toString(meanlongHomoPolymer)); m->mothurOutEndLine();
		if ((namefile == "") && (countfile == "")) {  m->mothurOut("# of Seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); }
		else { m->mothurOut("# of unique seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(size)); m->mothurOutEndLine(); }
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
		m->mothurOutEndLine();

        if ((namefile == "") && (countfile == "")) {  m->mothurOut("It took " + toString(time(NULL) - start) + " secs to summarize " + toString(numSeqs) + " sequences.\n");  }
        else{ m->mothurOut("It took " + toString(time(NULL) - start) + " secs to summarize " + toString(size) + " sequences.\n");   }
        
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
/**************************************************************************************/
 long long SeqSummaryCommand::driverCreateSummary(map<int, long long>& startPosition, map<int,  long long>& endPosition, map<int,  long long>& seqLength, map<int,  long long>& ambigBases, map<int,  long long>& longHomoPolymer, string filename, string sumFile, linePair* filePos) {
	try {
		
		ofstream outSummary;
		m->openOutputFile(sumFile, outSummary);
        
        ifstream in;
        m->openInputFile(filename, in);
        
        in.seekg(filePos->start);
		
		//print header if you are process 0
		if (filePos->start == 0) {
			outSummary << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl;
            m->zapGremlins(in); m->gobble(in);
		}
				
		bool done = false;
		int count = 0;
        
        
		while (!done) {
				
			if (m->control_pressed) { in.close(); outSummary.close(); return 1; }
            
            if (m->debug) { m->mothurOut("[DEBUG]: count = " + toString(count) + "\n");  }
            
			Sequence current(in); m->gobble(in);
           
			if (current.getName() != "") {
				
                if (m->debug) { m->mothurOut("[DEBUG]: " + current.getName() + '\t' + toString(current.getNumBases()) + "\n");  }
                
				int num = 1;
				if ((namefile != "") || (countfile != "")) {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(current.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: '" + current.getName() + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				int thisStartPosition = current.getStartPos();
                map<int, long long>::iterator it = startPosition.find(thisStartPosition);
                if (it == startPosition.end()) { startPosition[thisStartPosition] = num; } //first finding of this start position, set count.
                else { it->second += num; } //add counts
                
                int thisEndPosition = current.getEndPos();
                it = endPosition.find(thisEndPosition);
                if (it == endPosition.end()) { endPosition[thisEndPosition] = num; } //first finding of this end position, set count.
                else { it->second += num; } //add counts
                
                int thisSeqLength = current.getNumBases();
                it = seqLength.find(thisSeqLength);
                if (it == seqLength.end()) { seqLength[thisSeqLength] = num; } //first finding of this length, set count.
                else { it->second += num; } //add counts
                
                int thisAmbig = current.getAmbigBases();
                it = ambigBases.find(thisAmbig);
                if (it == ambigBases.end()) { ambigBases[thisAmbig] = num; } //first finding of this ambig, set count.
                else { it->second += num; } //add counts
                
                int thisHomoP = current.getLongHomoPolymer();
                it = longHomoPolymer.find(thisHomoP);
                if (it == longHomoPolymer.end()) { longHomoPolymer[thisHomoP] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts

				count++;
				outSummary << current.getName() << '\t';
				outSummary << thisStartPosition << '\t' << thisEndPosition << '\t';
				outSummary << thisSeqLength << '\t' << thisAmbig << '\t';
				outSummary << thisHomoP << '\t' << num << endl;
                
                if (m->debug) { m->mothurOut("[DEBUG]: " + current.getName() + '\t' + toString(num) + "\n");  }
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (in.eof()) { break; }
			#endif
		}
				
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "driverCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
 long long SeqSummaryCommand::createProcessesCreateSummary(map<int, long long>& startPosition, map<int, long long>& endPosition, map<int,  long long>& seqLength, map<int,  long long>& ambigBases, map<int,  long long>& longHomoPolymer, string filename, string sumFile) {
	try {
		int process = 1;
		int num = 0;
		processIDS.clear();
        bool recalc = false;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile + m->mothurGetpid(process) + ".temp", lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = fastafile + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << num << endl;
                out << startPosition.size() << endl;
				for (map<int,  long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << endPosition.size() << endl;
				for (map<int,  long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << seqLength.size() << endl;
				for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << ambigBases.size() << endl;
				for (map<int,  long long>::iterator it = ambigBases.begin(); it != ambigBases.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << longHomoPolymer.size() << endl;
				for (map<int,  long long>::iterator it = longHomoPolymer.begin(); it != longHomoPolymer.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(fastafile + (toString(processIDS[i]) + ".num.temp"));
                    m->mothurRemove(sumFile + (toString(processIDS[i]) + ".temp"));
                }
                recalc = true;
                break;
			}
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(fastafile + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(new linePair(positions[i], positions[(i+1)]));  }
            
            startPosition.clear();
            endPosition.clear();
            seqLength.clear();
            ambigBases.clear();
            longHomoPolymer.clear();
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile + m->mothurGetpid(process) + ".temp", lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = fastafile + m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << startPosition.size() << endl;
                    for (map<int,  long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << endPosition.size() << endl;
                    for (map<int,  long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << seqLength.size() << endl;
                    for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << ambigBases.size() << endl;
                    for (map<int,  long long>::iterator it = ambigBases.begin(); it != ambigBases.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << longHomoPolymer.size() << endl;
                    for (map<int,  long long>::iterator it = longHomoPolymer.begin(); it != longHomoPolymer.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

        
		//do your part
		num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile, lines[0]);

		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = fastafile + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
            long long  tempNum;
			in >> tempNum; m->gobble(in); num += tempNum;
			in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{
                 long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = startPosition.find(first);
                if (it == startPosition.end()) { startPosition[first] = second; } //first finding of this start position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{
                 long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = endPosition.find(first);
                if (it == endPosition.end()) { endPosition[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{
                 long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = seqLength.find(first);
                if (it == seqLength.end()) { seqLength[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{
                 long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = ambigBases.find(first);
                if (it == ambigBases.end()) { ambigBases[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{
                 long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = longHomoPolymer.find(first);
                if (it == longHomoPolymer.end()) { longHomoPolymer[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
							
			in.close();
			m->mothurRemove(tempFilename);
			
			m->appendFiles((sumFile + toString(processIDS[i]) + ".temp"), sumFile);
			m->mothurRemove((sumFile + toString(processIDS[i]) + ".temp"));
		}
		
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<seqSumData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
        
		bool hasNameMap = false;
        if ((namefile !="") || (countfile != "")) { hasNameMap = true; }
        
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
            string extension = "";
            if (i != 0) { extension = toString(i) + ".temp"; processIDS.push_back(i); }
			// Allocate memory for thread data.
			seqSumData* tempSum = new seqSumData(filename, (sumFile+extension), m, lines[i]->start, lines[i]->end, hasNameMap, nameMap);
			pDataArray.push_back(tempSum);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySeqSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
		num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, (sumFile+toString(processors-1)+".temp"), lines[processors-1]);
        processIDS.push_back(processors-1);

		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            for (map<int, long long>::iterator it = pDataArray[i]->startPosition.begin(); it != pDataArray[i]->startPosition.end(); it++)		{
                map<int, long long>::iterator itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->endPosition.begin(); it != pDataArray[i]->endPosition.end(); it++)		{
                map<int, long long>::iterator itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->seqLength.begin(); it != pDataArray[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->ambigBases.begin(); it != pDataArray[i]->ambigBases.end(); it++)		{
                map<int, long long>::iterator itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->longHomoPolymer.begin(); it != pDataArray[i]->longHomoPolymer.end(); it++)		{
                map<int, long long>::iterator itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
    
		//append files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((sumFile + toString(processIDS[i]) + ".temp"), sumFile);
			m->mothurRemove((sumFile + toString(processIDS[i]) + ".temp"));
		}
#endif		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**********************************************************************************************************************/



