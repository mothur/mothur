/*
 *  filterseqscommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "filterseqscommand.h"
#include "sequence.hpp"


//**********************************************************************************************************************
vector<string> FilterSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-filter",false,true, true); parameters.push_back(pfasta);
		CommandParameter phard("hard", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(phard);
		CommandParameter ptrump("trump", "String", "", "*", "", "", "","",false,false, true); parameters.push_back(ptrump);
		CommandParameter psoft("soft", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psoft);
		CommandParameter pvertical("vertical", "Boolean", "", "T", "", "", "","",false,false, true); parameters.push_back(pvertical);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false, true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The filter.seqs command reads a file containing sequences and creates a .filter and .filter.fasta file.\n";
		helpString += "The filter.seqs command parameters are fasta, trump, soft, hard, processors and vertical. \n";
		helpString += "The fasta parameter is required, unless you have a valid current fasta file. You may enter several fasta files to build the filter from and filter, by separating their names with -'s.\n";
		helpString += "For example: fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The trump option will remove a column if the trump character is found at that position in any sequence of the alignment. Default=*, meaning no trump. \n";
		helpString += "A soft mask removes any column where the dominant base (i.e. A, T, G, C, or U) does not occur in at least a designated percentage of sequences. Default=0.\n";
		helpString += "The hard parameter allows you to enter a file containing the filter you want to use.\n";
		helpString += "The vertical parameter removes columns where all sequences contain a gap character. The default is T.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The filter.seqs command should be in the following format: \n";
		helpString += "filter.seqs(fasta=yourFastaFile, trump=yourTrump) \n";
		helpString += "Example filter.seqs(fasta=abrecovery.fasta, trump=.).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],filter.fasta"; } 
        else if (type == "filter") {  pattern =  "[filename],filter"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "FilterSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
FilterSeqsCommand::FilterSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["filter"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "FilterSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
FilterSeqsCommand::FilterSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;  recalced = false;
		filterFileName = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("filter.seqs");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["filter"] = tempOutNames;
		
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
				
				it = parameters.find("hard");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["hard"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fasta = validParameter.valid(parameters, "fasta");
			if (fasta == "not found") { 				
				fasta = current->getFastaFile(); 
				if (fasta != "") { 
                    fastafileNames.push_back(fasta);  
                    m->mothurOut("Using " + fasta + " as input file for the fasta parameter."); m->mothurOutEndLine();
                    string simpleName = util.getSimpleName(fasta);
                    filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
                }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else { 
				util.splitAtDash(fasta, fastafileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastafileNames.size(); i++) {
					
					bool ignore = false;
					if (fastafileNames[i] == "current") { 
						fastafileNames[i] = current->getFastaFile(); 
						if (fastafileNames[i] != "") {  m->mothurOut("Using " + fastafileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastafileNames.erase(fastafileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastafileNames[i], current->getLocations())) {
                            string simpleName = util.getSimpleName(fastafileNames[i]);
                            filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
                            current->setFastaFile(fastafileNames[i]);
                        }
                        else { fastafileNames.erase(fastafileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (fastafileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (!abort) {
				//if the user changes the output directory command factory will send this info to us in the output parameter 
				outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
					outputDir = "";	
					outputDir += util.hasPath(fastafileNames[0]); //if user entered a file with a path then preserve it	
				}
			}
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			hard = validParameter.validFile(parameters, "hard");				if (hard == "not found") { hard = ""; }
			else if (hard == "not open") { hard = ""; abort = true; }	

			temp = validParameter.valid(parameters, "trump");			if (temp == "not found") { temp = "*"; }
			trump = temp[0];
			
			temp = validParameter.valid(parameters, "soft");				if (temp == "not found") { soft = 0; }
			else {  soft = (float)atoi(temp.c_str()) / 100.0;  }
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp); 
			
			vertical = validParameter.valid(parameters, "vertical");		
			if (vertical == "not found") { 
				if ((hard == "") && (trump == '*') && (soft == 0)) { vertical = "T"; } //you have not given a hard file or set the trump char.
				else { vertical = "F";  }
			}
			
			numSeqs = 0;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "FilterSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		ifstream inFASTA;
		util.openInputFile(fastafileNames[0], inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.close();
		
		////////////create filter/////////////////
		m->mothurOut("Creating Filter... "); m->mothurOutEndLine();
		
		filter = createFilter();
		
		m->mothurOutEndLine();  m->mothurOutEndLine();
		
		if (m->getControl_pressed()) { outputTypes.clear(); return 0; }
		
		ofstream outFilter;
		
		//prevent giantic file name
        map<string, string> variables; 
        variables["[filename]"] = outputDir + filterFileName + ".";
		if (fastafileNames.size() > 3) { variables["[filename]"] = outputDir + "merge."; }
		string filterFile = getOutputFileName("filter", variables);  
		
		util.openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		outputNames.push_back(filterFile); outputTypes["filter"].push_back(filterFile);
				
		////////////run filter/////////////////
		
		m->mothurOut("Running Filter... "); m->mothurOutEndLine();
		
		filterSequences();
		
		m->mothurOutEndLine();	m->mothurOutEndLine();
					
		int filteredLength = 0;
		for(int i=0;i<alignmentLength;i++){
			if(filter[i] == '1'){	filteredLength++;	}
		}
		
		if (m->getControl_pressed()) {  outputTypes.clear(); for(int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }

		
		m->mothurOutEndLine();
		m->mothurOut("Length of filtered alignment: " + toString(filteredLength)); m->mothurOutEndLine();
		m->mothurOut("Number of columns removed: " + toString((alignmentLength-filteredLength))); m->mothurOutEndLine();
		m->mothurOut("Length of the original alignment: " + toString(alignmentLength)); m->mothurOutEndLine();
		m->mothurOut("Number of sequences used to construct filter: " + toString(numSeqs)); m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
int FilterSeqsCommand::filterSequences() {	
	try {
		
		numSeqs = 0;
		
		for (int s = 0; s < fastafileNames.size(); s++) {
			
            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafileNames[s]));
            string filteredFasta = getOutputFileName("fasta", variables);
            
            vector<unsigned long long> positions;
            if (savedPositions.size() != 0) { positions = savedPositions[s]; }
            else {
#if defined NON_WINDOWS
            positions = util.divideFile(fastafileNames[s], processors);
#else
            positions = util.setFilePosFasta(fastafileNames[s], numSeqs);
            if (numSeqs < processors) { processors = numSeqs; }
#endif
            }
            
            vector<linePair> lines;
		#if defined NON_WINDOWS
			for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
		#else
            
            long long numFSeqs = positions.size()-1;
            if (numFSeqs < processors) { processors = numFSeqs; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFSeqs / processors;
            for (int i = 0; i < processors; i++) {
                long long startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }

		#endif
            
            long long numFastaSeqs = createProcessesRunFilter(filter, fastafileNames[s], filteredFasta, lines);
            numSeqs += numFastaSeqs;
			outputNames.push_back(filteredFasta); outputTypes["fasta"].push_back(filteredFasta);
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "filterSequences");
		exit(1);
	}
}
/**************************************************************************************/
void driverRunFilter(filterRunData* params) {
	try {
		ofstream out;
		params->util.openOutputFile(params->outputFilename, out);
		
		ifstream in;
		params->util.openInputFile(params->filename, in);
				
		in.seekg(params->start);
        
        //adjust start if null strings
        if (params->start == 0) {  params->util.zapGremlins(in); params->util.gobble(in);  }

		bool done = false;
		params->count = 0;
	
		while (!done) {
				
				if (params->m->getControl_pressed()) { break; }
				
				Sequence seq(in); params->util.gobble(in);
				if (seq.getName() != "") {
					string align = seq.getAligned();
					string filterSeq = "";
                
					for(int j=0;j<params->alignmentLength;j++){ if(params->filter[j] == '1'){ filterSeq += align[j]; } }
					
					out << '>' << seq.getName() << endl << filterSeq << endl;
                }
				params->count++;
        
        
			#if defined NON_WINDOWS
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= params->end)) { break; }
			#else
				if (params->count == params->end) { break; }
			#endif
			
			//report progress
			if((params->count) % 100 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 	}
        }
		//report progress
		if((params->count) % 100 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 		}
		
		
		out.close();
		in.close();
		
	}
	catch(exception& e) {
		params->m->errorOut(e, "FilterSeqsCommand", "driverRunFilter");
		exit(1);
	}
}
/**************************************************************************************************/

long long FilterSeqsCommand::createProcessesRunFilter(string F, string filename, string filteredFastaName, vector<linePair> lines) {
	try {
        util.mothurRemove(filteredFastaName);
        long long num = 0;
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<filterRunData*> data;
    
        time_t start, end;
        time(&start);
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            filterRunData* dataBundle = new filterRunData(filter, filename, (filteredFastaName + extension), m, lines[i+1].start, lines[i+1].end, alignmentLength);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverRunFilter, dataBundle));
        }
        
        filterRunData* dataBundle = new filterRunData(filter, filename, filteredFastaName, m, lines[0].start, lines[0].end, alignmentLength);
        data.push_back(dataBundle);
        driverRunFilter(dataBundle);
        num = dataBundle->count;
        
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            num += data[i]->count;
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to filter " + toString(num) + " sequences.\n");
        
        //append and remove temp files
        for (int i=0;i<processors-1;i++) {
            util.appendFiles((filteredFastaName + toString(i+1) + ".temp"), filteredFastaName);
            util.mothurRemove((filteredFastaName + toString(i+1) + ".temp"));
        }

        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesRunFilter");
		exit(1);
	}
}
/**************************************************************************************/
string FilterSeqsCommand::createFilter() {	
	try {
		string filterString = "";			
		Filters F;
		
		if (soft != 0)			{  F.setSoft(soft);		}
		if (trump != '*')		{  F.setTrump(trump);	}
		
		F.setLength(alignmentLength);
		
		if(trump != '*' || util.isTrue(vertical) || soft != 0){ F.initialize(); }
		
		if(hard.compare("") != 0)	{	F.doHard(hard);		}
		else						{	F.setFilter(string(alignmentLength, '1'));	}
		
		numSeqs = 0;
		if(trump != '*' || util.isTrue(vertical) || soft != 0){
			for (int s = 0; s < fastafileNames.size(); s++) {
			
                numSeqs += createProcessesCreateFilter(F, fastafileNames[s]);
                
				if (m->getControl_pressed()) {  return filterString; }
			}
		}

		F.setNumSeqs(numSeqs);
		if(util.isTrue(vertical) == 1)	{	F.doVertical();	}
		if(soft != 0)                   {	F.doSoft();		}
		filterString = F.getFilter();
        
		return filterString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createFilter");
		exit(1);
	}
}
/**************************************************************************************/
void driverCreateFilter(filterData* params) {
	try {
        if (params->soft != 0)			{  params->F.setSoft(params->soft);		}
        if (params->trump != '*')		{  params->F.setTrump(params->trump);	}
        
        params->F.setLength(params->alignmentLength);
        
        if(params->trump != '*' || params->vertical || params->soft != 0){ params->F.initialize(); }
        
        if(params->hard.compare("") != 0)	{	params->F.doHard(params->hard);                             }
        else                                {	params->F.setFilter(string(params->alignmentLength, '1'));	}
        
		ifstream in;
		params->util.openInputFile(params->filename, in);
				
		in.seekg(params->start);
        
        //adjust start if null strings
        if (params->start == 0) {  params->util.zapGremlins(in); params->util.gobble(in);  }

		bool done = false;
		params->count = 0;
        bool error = false;
        
		while (!done) {
				
            if (params->m->getControl_pressed()) { break; }
					
			Sequence seq(in); params->util.gobble(in);
			if (seq.getName() != "") {
                    if (params->m->getDebug()) { params->m->mothurOutJustToScreen("[DEBUG]: " + seq.getName() + " length = " + toString(seq.getAligned().length()) + '\n'); }
                if (seq.getAligned().length() != params->alignmentLength) { params->m->mothurOut("[ERROR]: Sequences are not all the same length, please correct.\n"); error = true; if (!params->m->getDebug()) { params->m->setControl_pressed(true); }else{ params->m->mothurOutJustToLog("[DEBUG]: " + seq.getName() + " length = " + toString(seq.getAligned().length()) + '\n'); } }
					
					if(params->trump != '*')                    {	params->F.doTrump(seq);		}
					if(params->vertical || params->soft != 0)	{	params->F.getFreqs(seq);	}
					cout.flush();
					params->count++;
            }
			
			#if defined NON_WINDOWS
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= params->end)) { break; }
			#else
				if (params->count == params->end) { break; }
			#endif
			
			//report progress
			if((params->count) % 100 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 		}
		}
		//report progress
		if((params->count) % 100 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 	}
		in.close();
        
        if (error) { params->m->setControl_pressed(true); }
        
	}
	catch(exception& e) {
		params->m->errorOut(e, "FilterSeqsCommand", "driverCreateFilter");
		exit(1);
	}
}
/**************************************************************************************************/

long long FilterSeqsCommand::createProcessesCreateFilter(Filters& F, string filename) {
	try {
        vector<linePair> lines;
        vector<unsigned long long> positions;
        
#if defined NON_WINDOWS
        positions = util.divideFile(filename, processors);
        for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
#else
        long long numFastaSeqs = 0;
        positions = util.setFilePosFasta(filename, numFastaSeqs);
        if (numFastaSeqs < processors) { processors = numFastaSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numFastaSeqs / processors;
        for (int i = 0; i < processors; i++) {
            long long startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        //save the file positions so we can reuse them in the runFilter function
        if (!recalced) {  savedPositions.push_back(positions); }
        
		long long num = 0;
        bool doVertical = util.isTrue(vertical);

        //create array of worker threads
        vector<thread*> workerThreads;
        vector<filterData*> data;
        
        time_t start, end;
        time(&start);
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            filterData* dataBundle = new filterData(filename, m, lines[i+1].start, lines[i+1].end, alignmentLength, trump, doVertical, soft, hard, i+1);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverCreateFilter, dataBundle));
        }
        
        filterData* dataBundle = new filterData(filename, m, lines[0].start, lines[0].end, alignmentLength, trump, doVertical, soft, hard, 0);
        driverCreateFilter(dataBundle);
        
        num = dataBundle->count;
        F.mergeFilter(dataBundle->F.getFilter());
        
        for (int k = 0; k < alignmentLength; k++) {	 F.a[k] += dataBundle->F.a[k];       }
        for (int k = 0; k < alignmentLength; k++) {	 F.t[k] += dataBundle->F.t[k];       }
        for (int k = 0; k < alignmentLength; k++) {	 F.g[k] += dataBundle->F.g[k];       }
        for (int k = 0; k < alignmentLength; k++) {	 F.c[k] += dataBundle->F.c[k];       }
        for (int k = 0; k < alignmentLength; k++) {	 F.gap[k] += dataBundle->F.gap[k];   }

        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            num += data[i]->count;
            F.mergeFilter(data[i]->F.getFilter());
            
            for (int k = 0; k < alignmentLength; k++) {	 F.a[k] += data[i]->F.a[k];       }
            for (int k = 0; k < alignmentLength; k++) {	 F.t[k] += data[i]->F.t[k];       }
            for (int k = 0; k < alignmentLength; k++) {	 F.g[k] += data[i]->F.g[k];       }
            for (int k = 0; k < alignmentLength; k++) {	 F.c[k] += data[i]->F.c[k];       }
            for (int k = 0; k < alignmentLength; k++) {	 F.gap[k] += data[i]->F.gap[k];   }
            
            
            delete data[i];
            delete workerThreads[i];
        }
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to create filter for " + toString(num) + " sequences.\n");
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: filter = " + F.getFilter() + "\n\n");  }
        
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesCreateFilter");
		exit(1);
	}
}
/**************************************************************************************/
