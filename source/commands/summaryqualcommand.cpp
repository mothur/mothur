/*
 *  summaryqualcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/28/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "summaryqualcommand.h"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> SummaryQualCommand::setParameters(){	
	try {
		CommandParameter pqual("qfile", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pqual);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
        
        abort = false; calledHelp = false;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryQualCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.qual command reads a quality file and an optional name or count file, and summarizes the quality information.\n";
		helpString += "The summary.qual command parameters are qfile, name, count and processors. qfile is required, unless you have a valid current quality file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your quality file. \n";
        helpString += "The count parameter allows you to enter a count file associated with your quality file. \n";
		helpString += "The summary.qual command should be in the following format: \n";
		helpString += "summary.qual(qfile=yourQualityFile) \n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryQualCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],qual.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SummaryQualCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
SummaryQualCommand::SummaryQualCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { qualfile = ""; abort = true; }
			else if (qualfile == "not found") { 				
				qualfile = current->getQualFile(); 
				if (qualfile != "") { m->mothurOut("Using " + qualfile + " as input file for the qfile parameter.\n");  }
				else { 	m->mothurOut("You have no current quality file and the qfile parameter is required.\n");  abort = true; }
			}else { current->setQualFile(qualfile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name.\n");  abort = true; }
			
			if (outputdir == ""){	 outputdir += util.hasPath(qualfile);  }
			
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
            
			if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(qualfile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "SummaryQualCommand");
		exit(1);
	}
}
//***************************************************************************************************************
int SummaryQualCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		long start = time(NULL);
		long long numSeqs = 0;
        hasNameMap = false;
		
		vector<int> position;
		vector<int> averageQ;
		vector< vector<int> > scores;
				
		if (m->getControl_pressed()) { return 0; }
		
        if (namefile != "") { hasNameMap = true; nameMap = util.readNames(namefile); }
		else if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, false, false);
            nameMap = ct.getNameMap();
            hasNameMap = true;
        }
        
        numSeqs = createProcessesCreateSummary(position, averageQ, scores, qualfile);
		
		if (m->getControl_pressed()) {  return 0; }
		
		//print summary file
        map<string, string> variables; 
		variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(qualfile));
		string summaryFile = getOutputFileName("summary",variables);
		printQual(summaryFile, position, averageQ, scores);
		
		if (m->getControl_pressed()) {  util.mothurRemove(summaryFile); return 0; }
		
		//output results to screen
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		m->mothurOut("\nPosition\tNumSeqs\tAverageQ\n");
		for (int i = 0; i < position.size(); i+=100) {
			float average = averageQ[i] / (float) position[i];
			cout << i << '\t' << position[i] << '\t' << average << '\n';
			m->mothurOutJustToLog(toString(i) + "\t" + toString(position[i]) + "\t" + toString(average)+"\n");
		}
		
        outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
		m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numSeqs) + " sequences.\n\n");
		m->mothurOut("Output File Names: \n");
		m->mothurOut(summaryFile+"\n\n");
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct seqSumQualData {
    vector<int> position;
    vector<int> averageQ;
    vector< vector<int> > scores;
    string filename;
    unsigned long long start;
    unsigned long long end;
    int count, numSeqs;
    MothurOut* m;
    bool hasNameMap;
    map<string, int> nameMap;
    Utils util;
    
    ~seqSumQualData(){}
    seqSumQualData(string f, unsigned long long st, unsigned long long en, bool n, map<string, int> nam) {
        filename = f;
        m = MothurOut::getInstance();
        start = st;
        end = en;
        hasNameMap = n;
        nameMap = nam;
        count = 0;
    }
};

/**************************************************************************************/
void driverCreateSummary(seqSumQualData* params) {
	try {
		ifstream in;
		params->util.openInputFile(params->filename, in);
		
		in.seekg(params->start);
        
        //adjust start if null strings
        if (params->start == 0) {  params->util.zapGremlins(in); params->util.gobble(in);  }
		
		bool done = false;
		params->count = 0;
        int count = 0;
		
		while (!done) {
			
			if (params->m->getControl_pressed()) { in.close(); break; }
			
			QualityScores current(in); params->util.gobble(in);
			
			if (current.getName() != "") {
				
				int num = 1;
				if (params->hasNameMap) {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = params->nameMap.find(current.getName());
					
					if (it == params->nameMap.end()) { params->m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct.\n"); params->m->setControl_pressed(true); }
					else { num = it->second; }
				}
				
				vector<int> thisScores = current.getScores();
				
				//resize to num of positions setting number of seqs with that size to 1
				if (params->position.size() < thisScores.size()) { params->position.resize(thisScores.size(), 0); }
				if (params->averageQ.size() < thisScores.size()) { params->averageQ.resize(thisScores.size(), 0); }
				if (params->scores.size() < thisScores.size()) {
					params->scores.resize(thisScores.size());
					for (int i = 0; i < params->scores.size(); i++) { params->scores[i].resize(41, 0); }
				}
				
				//increase counts of number of seqs with this position
				//average is really the total, we will average in execute
				for (int i = 0; i < thisScores.size(); i++) { 
					params->position[i] += num;
					params->averageQ[i] += (thisScores[i] * num); //weighting for namesfile
                    if (thisScores[i] > 41) { params->m->mothurOut("[WARNING]: " + current.getName() + " has a quality scores of " + toString(thisScores[i]) + ", expecting values to be less than 40. Setting to 40.\n");  thisScores[i] = 40;
                    }
					else { params->scores[i][thisScores[i]] += num; }
				}
				
				params->count += num;   //totalSeqs
                count++;                //uniqueSeqs
			}
			
#if defined NON_WINDOWS
			unsigned long long pos = in.tellg();
			if ((pos == -1) || (pos >= params->end)) { break; }
#else
			if ((count == params->end) || (in.eof())) { break; }
#endif
		}
		
		in.close();
    }
	catch(exception& e) {
		params->m->errorOut(e, "SummaryQualCommand", "driverCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
long long SummaryQualCommand::createProcessesCreateSummary(vector<int>& position, vector<int>& averageQ, vector< vector<int> >& scores, string filename) {
	try {
        long long numSeqs = 0;
        vector<double> positions;
        vector<linePair> lines;
#if defined NON_WINDOWS
        positions = util.divideFile(filename, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        
            positions = util.setFilePosFasta(qualfile, numSeqs);
            if (numSeqs < processors) { processors = numSeqs; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        
#endif

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<seqSumQualData*> data;
        //string f, unsigned long long st, unsigned long long en, bool n, map<string, int> nam
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            seqSumQualData* dataBundle = new seqSumQualData(filename, lines[i+1].start, lines[i+1].end, hasNameMap, nameMap);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverCreateSummary, dataBundle));
        }
        
        seqSumQualData* dataBundle = new seqSumQualData(filename, lines[0].start, lines[0].end, hasNameMap, nameMap);
        
        driverCreateSummary(dataBundle);
        numSeqs = dataBundle->count;
        position = dataBundle->position;
        averageQ = dataBundle->averageQ;
        scores = dataBundle->scores;
        
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            numSeqs += data[i]->count;
            
            int tempNum = data[i]->position.size();
            if (position.size() < tempNum) { position.resize(tempNum, 0); }
            if (averageQ.size() < tempNum) { averageQ.resize(tempNum, 0); }
            if (scores.size() < tempNum) {
                scores.resize(tempNum);
                for (int i = 0; i < scores.size(); i++) { scores[i].resize(41, 0); }
            }
            
            for (int k = 0; k < tempNum; k++)			{		 position[k]    +=  data[i]->position[k];         }
            for (int k = 0; k < tempNum; k++)			{		 averageQ[k]    +=  data[i]->averageQ[k];         }
            for (int k = 0; k < tempNum; k++)			{	for (int j = 0; j < 41; j++) {  scores[k][j] += data[i]->scores[k][j];   }	}
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
int SummaryQualCommand::printQual(string sumFile, vector<int>& position, vector<int>& averageQ, vector< vector<int> >& scores) {
	try {
		ofstream out;
		util.openOutputFile(sumFile, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		outputNames.push_back(sumFile); outputTypes["summary"].push_back(sumFile);
		
		//print headings
		out << "Position\tnumSeqs\tAverageQ";
		for (int i = 0; i < 41; i++) { out << '\t' << "q" << i; }
		out << endl;
		
		for (int i = 0; i < position.size(); i++) {
			
			if (m->getControl_pressed()) { out.close(); return 0; }
			
			double average = averageQ[i] / (float) position[i];
			out << i << '\t' << position[i] << '\t' << average;
			
			for (int j = 0; j < 41; j++) {
				out  << '\t' << scores[i][j];
			}
			out << endl;
		}
		
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "printQual");
		exit(1);
	}
}

/**************************************************************************************/


