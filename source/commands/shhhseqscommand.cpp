/*
 *  shhhseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "shhhseqscommand.h"



//**********************************************************************************************************************
vector<string> ShhhSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-map",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","name",false,true,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter psigma("sigma", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(psigma);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhhSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The shhh.seqs command reads a fasta and name file and ....\n";
		helpString += "The shhh.seqs command parameters are fasta, name, group, sigma and processors.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file associated with your fasta file. It is required. \n";
		helpString += "The group parameter allows you to provide a group file.  When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The sigma parameter ....  The default is 0.01. \n";
		helpString += "The shhh.seqs command should be in the following format: \n";
		helpString += "shhh.seqs(fasta=yourFastaFile, name=yourNameFile) \n";
		helpString += "Example: shhh.seqs(fasta=AD.align, name=AD.names) \n";
			
		return helpString;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhhSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],shhh_seqs.fasta";   }
        else if (type == "name")    {   pattern = "[filename],shhh_seqs.names";   }
        else if (type == "map")        {   pattern = "[filename],shhh_seqs.map";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ShhhSeqsCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************

ShhhSeqsCommand::ShhhSeqsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["map"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "ShhhSeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
ShhhSeqsCommand::ShhhSeqsCommand(string option) {
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
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { 			
				namefile = current->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current namefile and the name parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (namefile == "not open") { namefile =  ""; abort = true; }	
			else {  current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") { groupfile =  "";   }
			else if (groupfile == "not open") { abort = true; groupfile =  ""; }	
			else {   current->setGroupFile(groupfile);  }
			
			string temp	= validParameter.valid(parameters, "sigma");		if(temp == "not found"){	temp = "0.01"; }
			util.mothurConvert(temp, sigma); 
			sigma = 1/sigma;
            
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
                    
			if (namefile == "") {
				vector<string> files; files.push_back(fastafile);
				if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "ShhhSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int driver(seqNoise& noise,
           vector<string>& sequences,
           vector<string>& uniqueNames,
           vector<string>& redundantNames,
           vector<int>& seqFreq,
           string distFileName, string outputFileName, string nameFileName, string mapFileName, MothurOut* m, int sigma) {
    try {
        Utils util;
        double cutOff = 0.08;
        int minIter = 10;
        int maxIter = 1000;
        double minDelta = 1e-6;
        int numIters = 0;
        double maxDelta = MOTHURMAX;
        int numSeqs = sequences.size();
        
        //run cluster command
        string inputString = "phylip=" + distFileName + ", method=furthest, cutoff=0.08";
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: cluster(" + inputString + ")\n");
        
        Command* clusterCommand = new ClusterCommand(inputString);
        clusterCommand->execute();
        
        map<string, vector<string> > filenames = clusterCommand->getOutputFiles();
        string listFileName = filenames["list"][0];
        string rabundFileName = filenames["rabund"][0]; util.mothurRemove(rabundFileName);
        string sabundFileName = filenames["sabund"][0]; util.mothurRemove(sabundFileName);
        
        delete clusterCommand;
        m->mothurOut("/******************************************/\n");
        
        if (m->getControl_pressed()) { util.mothurRemove(distFileName); util.mothurRemove(listFileName); return 0; }
        
        vector<double> distances(numSeqs * numSeqs);
        noise.getDistanceData(distFileName, distances);
        util.mothurRemove(distFileName);
        if (m->getControl_pressed()) { util.mothurRemove(listFileName); return 0; }
        
        vector<int> otuData(numSeqs);
        vector<int> otuFreq;
        vector<vector<int> > otuBySeqLookUp;
        noise.getListData(listFileName, cutOff, otuData, otuFreq, otuBySeqLookUp);
        util.mothurRemove(listFileName);
        if (m->getControl_pressed()) { return 0; }
        
        int numOTUs = otuFreq.size();
        vector<double> weights(numOTUs, 0);
        vector<int> change(numOTUs, 1);
        vector<int> centroids(numOTUs, -1);
        vector<int> cumCount(numOTUs, 0);
        
        vector<double> tau(numSeqs, 1);
        vector<int> anP(numSeqs, 0);
        vector<int> anI(numSeqs, 0);
        vector<int> anN(numSeqs, 0);
        vector<vector<int> > aanI = otuBySeqLookUp;
        
        while(numIters < minIter || ((maxDelta > minDelta) && (numIters < maxIter))){
            
            if (m->getControl_pressed()) { return 0; }
            
            noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount); if (m->getControl_pressed()) { return 0; }
            maxDelta = noise.calcNewWeights(weights, seqFreq, anI, cumCount, anP, otuFreq, tau);  if (m->getControl_pressed()) { return 0; }
            
            noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (m->getControl_pressed()) { return 0; }
            noise.checkCentroids(weights, centroids); if (m->getControl_pressed()) { return 0; }
            
            otuFreq.assign(numOTUs, 0);
            
            int total = 0;
            
            for(int i=0;i<numSeqs;i++){
                if (m->getControl_pressed()) { return 0; }
                
                double offset = MOTHURMAX;
                double norm = 0.0000;
                double minWeight = 0.1;
                vector<double> currentTau(numOTUs);
                
                for(int j=0;j<numOTUs;j++){
                    if (m->getControl_pressed()) { return 0; }
                    if(weights[j] > minWeight && distances[i * numSeqs+centroids[j]] < offset){
                        offset = distances[i * numSeqs+centroids[j]];
                    }
                }
                
                for(int j=0;j<numOTUs;j++){
                    if (m->getControl_pressed()) { return 0; }
                    if(weights[j] > minWeight){
                        currentTau[j] = exp(sigma * (-distances[(i * numSeqs + centroids[j])] + offset)) * weights[j];
                        norm += currentTau[j];
                    }
                    else{ currentTau[j] = 0.0000; }
                }
                
                for(int j=0;j<numOTUs;j++){ currentTau[j] /= norm; }
                
                for(int j=0;j<numOTUs;j++){
                    if (m->getControl_pressed()) { return 0; }
                    
                    if(currentTau[j] > 1.0e-4){
                        int oldTotal = total;
                        total++;
                        
                        tau.resize(oldTotal+1);
                        tau[oldTotal] = currentTau[j];
                        otuBySeqLookUp[j][otuFreq[j]] = oldTotal;
                        aanI[j][otuFreq[j]] = i;
                        otuFreq[j]++;
                        
                    }
                }
                
                anP.resize(total);
                anI.resize(total);
            }
            
            numIters++;
        }
        
        noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);  if (m->getControl_pressed()) { return 0; }
        
        vector<double> percentage(numSeqs);
        noise.setUpOTUData(otuData, percentage, cumCount, tau, otuFreq, anP, anI);  if (m->getControl_pressed()) { return 0; }
        noise.finishOTUData(otuData, otuFreq, anP, anI, cumCount, otuBySeqLookUp, aanI, tau);  if (m->getControl_pressed()) { return 0; }
        
        change.assign(numOTUs, 1);
        noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (m->getControl_pressed()) { return 0; }
        
        
        vector<int> finalTau(numOTUs, 0);
        for(int i=0;i<numSeqs;i++){
            if (m->getControl_pressed()) { return 0; }
            finalTau[otuData[i]] += int(seqFreq[i]);
        }
        
        noise.writeOutput(outputFileName, nameFileName, mapFileName, finalTau, centroids, otuData, sequences, uniqueNames, redundantNames, seqFreq, distances);
        
        return 0;
        
    }catch(exception& e) {
        m->errorOut(e, "ShhhSeqsCommand", "driver");
        exit(1);
    }
}
//**********************************************************************************************************************
int ShhhSeqsCommand::execute() {
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (outputDir == "") { outputDir = util.hasPath(fastafile);  }//if user entered a file with a path then preserve it		
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
		string outputFileName = getOutputFileName("fasta",variables);
		string nameFileName = getOutputFileName("name",variables);
		string mapFileName = getOutputFileName("map",variables);
		
		if (groupfile != "") {
            mapFileName = outputDir + util.getRootName(util.getSimpleName(fastafile))  + "shhh.";
			vector<string> mapFileNames = createProcessesGroups(outputFileName, nameFileName, mapFileName);

			if (m->getControl_pressed()) {    return 0;	}	
			
			for (int j = 0; j < mapFileNames.size(); j++) { outputNames.push_back(mapFileNames[j]); outputTypes["map"].push_back(mapFileNames[j]); }
			
			//deconvolute results by running unique.seqs
			deconvoluteResults(outputFileName, nameFileName);
			
			if (m->getControl_pressed()) {   return 0;	}				
			
		}else{	
			vector<string> sequences;
			vector<string> uniqueNames;
			vector<string> redundantNames;
			vector<int> seqFreq;
			
			seqNoise noise;
			correctDist* correct = new correctDist(processors);
			
			//reads fasta and name file and loads them in order
			readData(correct, noise, sequences, uniqueNames, redundantNames, seqFreq);
			if (m->getControl_pressed()) { return 0; }
			
			//calc distances for cluster
			string distFileName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + "shhh.dist";
			correct->execute(distFileName);
			delete correct;
			
			if (m->getControl_pressed()) { util.mothurRemove(distFileName); return 0; }
			
			driver(noise, sequences, uniqueNames, redundantNames, seqFreq, distFileName, outputFileName, nameFileName, mapFileName, m, sigma);
			outputNames.push_back(mapFileName); outputTypes["map"].push_back(mapFileName);
		}
		
		if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
		
		outputNames.push_back(outputFileName); outputTypes["fasta"].push_back(outputFileName);
		outputNames.push_back(nameFileName); outputTypes["name"].push_back(nameFileName);
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::readData(correctDist* correct, seqNoise& noise, vector<string>& seqs, vector<string>& uNames, vector<string>& rNames, vector<int>& freq) {
	try {
		map<string, string> nameMap; 
		map<string, string>::iterator it;
		util.readNames(namefile, nameMap);
		bool error = false;
		
		ifstream in;
		util.openInputFile(fastafile, in);
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			Sequence seq(in); util.gobble(in);
			
			if (seq.getName() != "") {
				correct->addSeq(seq.getName(), seq.getAligned());
				
				it = nameMap.find(seq.getName());
				if (it != nameMap.end()) {
					noise.addSeq(seq.getAligned(), seqs);
					noise.addRedundantName(it->first, it->second, uNames, rNames, freq);
				}else {
					m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your namefile, please correct.");
					error = true;
				}
			}
		}
		in.close();
		
		if (error) { m->setControl_pressed(true); }
		
		return seqs.size();
		
	}catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "readData");
		exit(1);
	}
}
//**********************************************************************************************************************
int loadData(MothurOut* m, correctDist* correct, seqNoise& noise, vector<string>& seqs, vector<string>& uNames, vector<string>& rNames, vector<int>& freq, map<string, string>& nameMap, vector<Sequence>& sequences) {
	try {
		bool error = false;
		map<string, string>::iterator it;
		
		for (int i = 0; i < sequences.size(); i++) {
			
			if (m->getControl_pressed()) { return 0; }
			
			if (sequences[i].getName() != "") {
				correct->addSeq(sequences[i].getName(), sequences[i].getAligned());
				
				it = nameMap.find(sequences[i].getName());
				if (it != nameMap.end()) {
					noise.addSeq(sequences[i].getAligned(), seqs);
					noise.addRedundantName(it->first, it->second, uNames, rNames, freq);
				}else {
					m->mothurOut("[ERROR]: " + sequences[i].getName() + " is in your fasta file and not in your namefile, please correct.");
					error = true;
				}
			}
		}
				
		if (error) { m->setControl_pressed(true); }
		
		return seqs.size();
		
	}catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "loadData");
		exit(1);
	}
}
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct shhhseqsData {
    string fastafile;
    string namefile;
    string groupfile;
    string newFFile, newNFile, newMFile, outputDir, extension;
    MothurOut* m;
    int sigma, count;
    vector<string> groups;
    vector<string> mapfileNames;
    Utils util;
    
    shhhseqsData(){}
    shhhseqsData(string opd, string f, string n, string g, string nff,  string nnf, string nmf, vector<string> gr, int s, string ex) {
        outputDir = opd;
        fastafile = f;
        namefile = n;
        groupfile = g;
        newFFile = nff;
        newNFile = nnf;
        newMFile = nmf;
        m = MothurOut::getInstance();
        sigma = s;
        groups = gr;
        extension = ex;
        count=0;
    }
};
/**************************************************************************************************/
void driverShhSeqsGroups(shhhseqsData* params){
    try {
        //Parse sequences by group
        //Parse sequences by group
        vector<string> groups;
        map<string, vector<string> > group2Files;
        
        SequenceParser sparser(params->groupfile, params->fastafile, params->namefile, params->groups);
        groups = sparser.getNamesOfGroups();
        group2Files = sparser.getFiles();

        string fileroot = params->outputDir + params->util.getRootName(params->util.getSimpleName(params->fastafile));
        
        for (map<string, vector<string> >::iterator it = group2Files.begin(); it != group2Files.end(); it++) {
            long start = time(NULL);	 if (params->m->getControl_pressed()) {  break; }
            
            string thisGroup = it->first;
            
            string lowerCaseName = thisGroup;
            for (int j = 0; j < lowerCaseName.length(); j++) { lowerCaseName[j] = tolower(lowerCaseName[j]);    }
            
            if (lowerCaseName == "ignore") {   }
            else {
                params->m->mothurOut("\nProcessing group " + thisGroup + ":\n");
                
                map<string, string> thisNameMap;
                params->util.readNames(it->second[1], thisNameMap);

                vector<Sequence> thisSeqs;
                ifstream in; params->util.openInputFile(it->second[0], in);
                while (!in.eof()) {
                    if (params->m->getControl_pressed()) { break; }
                    
                    Sequence seq(in); params->util.gobble(in);
                    
                    if (seq.getName() != "") { thisSeqs.push_back(seq); }
                }
                in.close();

                vector<string> sequences;
                vector<string> uniqueNames;
                vector<string> redundantNames;
                vector<int> seqFreq;
                
                seqNoise noise;
                correctDist* correct = new correctDist(1); //we use one processor since we already split up the work load.
                
                //load this groups info in order
                loadData(params->m, correct, noise, sequences, uniqueNames, redundantNames, seqFreq, thisNameMap, thisSeqs);
                if (params->m->getControl_pressed()) { break; }
                
                //calc distances for cluster
                string distFileName = fileroot + thisGroup + ".shhh.dist";
                correct->execute(distFileName);
                delete correct;
                
                if (params->m->getControl_pressed()) { params->util.mothurRemove(distFileName); break; }
                
                driver(noise, sequences, uniqueNames, redundantNames, seqFreq, distFileName, params->newFFile+thisGroup, params->newNFile+thisGroup, params->newMFile+thisGroup+".map", params->m, params->sigma);
                
                if (params->m->getControl_pressed()) { break; }
                
                params->util.appendFiles(params->newFFile+thisGroup, params->newFFile+params->extension); params->util.mothurRemove(params->newFFile+thisGroup);
                params->util.appendFiles(params->newNFile+thisGroup, params->newNFile+params->extension); params->util.mothurRemove(params->newNFile+thisGroup);
                params->mapfileNames.push_back(params->newMFile+thisGroup+".map");
                
                params->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process group " + thisGroup + ".\n");
            }
        }
    }
    catch(exception& e) {
        params->m->errorOut(e, "ShhhSeqsCommand", "driverShhSeqsGroups");
        exit(1);
    }
}
/**************************************************************************************************/
vector<string> ShhhSeqsCommand::createProcessesGroups(string newFName, string newNName, string newMName) {
	try {
        GroupMap groupMap(groupfile); groupMap.readMap();vector<string> groups = groupMap.getNamesOfGroups();
		if (groups.size() < processors) { processors = groups.size(); m->mothurOut("Reducing processors to " + toString(groups.size()) + ".\n"); }
		
		//divide the groups between the processors
		vector<vector<string> > dividedGroupNames;
		int remainingPairs = groups.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            vector<string> thisProcessorsGroups;
            for (int i = startIndex; i < (startIndex+numPairs); i++) { thisProcessorsGroups.push_back(groups[i]);  }
            dividedGroupNames.push_back(thisProcessorsGroups); 
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        
        vector<std::thread*> workerThreads;
        vector<shhhseqsData*> data;
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            util.mothurRemove(newFName+extension);
            util.mothurRemove(newNName+extension);
            shhhseqsData* dataBundle = new shhhseqsData(outputDir, fastafile, namefile, groupfile, newFName, newNName, newMName, dividedGroupNames[i+1], sigma, extension);
            data.push_back(dataBundle);
            
            std::thread* thisThread = new thread(driverShhSeqsGroups, dataBundle);
            workerThreads.push_back(thisThread);
        }
        
        util.mothurRemove(newFName);
        util.mothurRemove(newNName);
        shhhseqsData* dataBundle = new shhhseqsData(outputDir, fastafile, namefile, groupfile, newFName, newNName, newMName, dividedGroupNames[0], sigma, "");
        driverShhSeqsGroups(dataBundle);
        vector<string> mapFileNames = dataBundle->mapfileNames;
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            mapFileNames.insert(mapFileNames.end(), data[i]->mapfileNames.begin(), data[i]->mapfileNames.end());
            util.appendFiles(data[i]->newFFile+data[i]->extension, newFName); util.mothurRemove(data[i]->newFFile+data[i]->extension);
            util.appendFiles(data[i]->newNFile+data[i]->extension, newNName); util.mothurRemove(data[i]->newNFile+data[i]->extension);
            
            delete data[i];
            delete workerThreads[i];
        }
		
		return mapFileNames;	
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "createProcessesGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::deconvoluteResults(string fastaFile, string nameFile){
	try {
		m->mothurOutEndLine(); m->mothurOut("Deconvoluting results:"); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		//use unique.seqs to create new name and fastafile
		string inputString = "fasta=" + fastaFile + ", name=" + nameFile;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
		current->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		current->setMothurCalling(false);
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		string newnameFile = filenames["name"][0];
		string newfastaFile = filenames["fasta"][0];
		
		util.mothurRemove(fastaFile); rename(newfastaFile.c_str(), fastaFile.c_str()); 
		if (nameFile != newnameFile) { util.mothurRemove(nameFile); rename(newnameFile.c_str(), nameFile.c_str()); }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "deconvoluteResults");
		exit(1);
	}
}
//**********************************************************************************************************************



