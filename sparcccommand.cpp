//
//  sparcccommand.cpp
//  Mothur
//
//  Created by SarahsWork on 5/10/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "sparcccommand.h"


//**********************************************************************************************************************
vector<string> SparccCommand::setParameters(){
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","outputType",false,true); parameters.push_back(pshared);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter psamplings("samplings", "Number", "", "20", "", "", "","",false,false,false); parameters.push_back(psamplings);
        CommandParameter piterations("iterations", "Number", "", "10", "", "", "","",false,false,false); parameters.push_back(piterations);
        CommandParameter ppermutations("permutations", "Number", "", "1000", "", "", "","",false,false,false); parameters.push_back(ppermutations);
        CommandParameter pmethod("method", "Multiple", "relabund-dirichlet", "dirichlet", "", "", "","",false,false); parameters.push_back(pmethod);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors); 
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SparccCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The sparcc command allows you to ....\n";
		helpString += "The sparcc command parameters are: shared, groups, label, samplings, iterations, permutations, processors and method.\n";
		helpString += "The samplings parameter is used to .... Default=20.\n";
        helpString += "The iterations parameter is used to ....Default=10.\n";
        helpString += "The permutations parameter is used to ....Default=1000.\n";
        helpString += "The method parameter is used to ....Options are relabund and dirichlet. Default=dirichlet.\n";
        helpString += "The default value for groups is all the groups in your sharedfile.\n";
		helpString += "The label parameter is used to analyze specific labels in your shared file.\n";
		helpString += "The sparcc command should be in the following format: sparcc(shared=yourSharedFile)\n";
		helpString += "sparcc(shared=final.an.shared)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SparccCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "corr") {  pattern = "[filename],[distance],sparcc_correlation"; }
        else if (type == "pvalue") {  pattern = "[filename],[distance],sparcc_pvalue"; }
        else if (type == "sparccrelabund") {  pattern = "[filename],[distance],sparcc_relabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SparccCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SparccCommand::SparccCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["corr"] = tempOutNames; 
		outputTypes["pvalue"] = tempOutNames;
		outputTypes["sparccrelabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "SparccCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SparccCommand::SparccCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["corr"] = tempOutNames; //filetypes should be things like: shared, fasta, accnos...
            outputTypes["pvalue"] = tempOutNames;
            outputTypes["sparccrelabund"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                 
				string path;
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
                        
			//check for parameters
            //get shared file, it is required
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile();
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
				outputDir = m->hasPath(sharedfile); //if user entered a file with a path then preserve it
			}
            
			normalizeMethod = validParameter.validFile(parameters, "method", false);
			if (normalizeMethod == "not found") { normalizeMethod = "dirichlet"; }
			if ((normalizeMethod == "dirichlet") || (normalizeMethod == "relabund")) { }
			else { m->mothurOut(normalizeMethod + " is not a valid method.  Valid methods are dirichlet and relabund."); m->mothurOutEndLine(); abort = true; }
            
            
            string temp = validParameter.validFile(parameters, "samplings", false);	if (temp == "not found"){	temp = "20";	}
			m->mothurConvert(temp, numSamplings);
            
            if(normalizeMethod == "relabund"){ numSamplings = 1; }
            
            temp = validParameter.validFile(parameters, "iterations", false);	if (temp == "not found"){	temp = "10";	}
			m->mothurConvert(temp, maxIterations);
            
            temp = validParameter.validFile(parameters, "permutations", false);	if (temp == "not found"){	temp = "1000";	}
			m->mothurConvert(temp, numPermutations);
            
            temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
    
            string groups = validParameter.validFile(parameters, "groups", false);
			if (groups == "not found") { groups = ""; }
			else { m->splitAtDash(groups, Groups); }
			m->setGroups(Groups);
            
            string label = validParameter.validFile(parameters, "label", false);
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "SparccCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SparccCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        int start = time(NULL);
        
        InputData input(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
        string lastLabel = lookup[0]->getLabel();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }return 0; }
            
            if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
                
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
            }
            
            if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup[0]->getLabel();
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                lookup = input.getSharedRAbundVectors(lastLabel);
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup);                
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
                
                //restore real lastlabel to save below
                lookup[0]->setLabel(saveLabel);
            }
            
            lastLabel = lookup[0]->getLabel();
            //prevent memory leak
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
            
            if (m->control_pressed) { return 0; }
            
            //get next line to process
            lookup = input.getSharedRAbundVectors();
        }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        //output error messages about any remaining user labels
        set<string>::iterator it;
        bool needToRun = false;
        for (it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it);
            if (processedLabels.count(lastLabel) != 1) {
                m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
                needToRun = true;
            }else {
                m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
            }
        }
        
        //run last label if you need to
        if (needToRun == true)  {
            for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            
            process(lookup);           
            
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
        }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " seconds to process.");
        m->mothurOutEndLine();
        m->mothurOutEndLine();
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
vector<vector<float> > SparccCommand::shuffleSharedVector(vector<vector<float> >& sharedVector){
    try {
        int numGroups = (int)sharedVector.size();
        int numOTUs = (int)sharedVector[0].size();
        
        vector<vector<float> > shuffledVector = sharedVector;
        
        for(int i=0;i<numGroups;i++){
            for(int j=0;j<numOTUs;j++){
                shuffledVector[i][j] = sharedVector[rand()%numGroups][j];
            }
        }
        
        return shuffledVector;
    }
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SparccCommand::process(vector<SharedRAbundVector*>& lookup){
	try {
        cout.setf(ios::fixed, ios::floatfield);
        cout.setf(ios::showpoint);
        
        vector<vector<float> > sharedVector;
        vector<string> otuNames = m->currentSharedBinLabels;
        
        //fill sharedVector to pass to CalcSparcc
        for (int i = 0; i < lookup.size(); i++) {
            vector<int> abunds = lookup[i]->getAbundances();
            vector<float> temp;
            for (int j = 0; j < abunds.size(); j++) { temp.push_back((float) abunds[j]); }
            sharedVector.push_back(temp);
        }
        int numOTUs = (int)sharedVector[0].size();
        int numGroups = lookup.size();
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        

        string relAbundFileName = getOutputFileName("sparccrelabund", variables);
        ofstream relAbundFile;
        m->openOutputFile(relAbundFileName, relAbundFile);
        outputNames.push_back(relAbundFileName);  outputTypes["sparccrelabund"].push_back(relAbundFileName);
        
        relAbundFile << "OTU\taveRelAbund\n";
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) { relAbundFile.close(); return 0; }
            
            double relAbund = 0.0000;
            for(int j=0;j<numGroups;j++){
                relAbund += sharedVector[j][i]/(double)lookup[j]->getNumSeqs();
            }
            relAbundFile << otuNames[i] <<'\t' << relAbund / (double) numGroups << endl;
        }
        relAbundFile.close();
        
        CalcSparcc originalData(sharedVector, maxIterations, numSamplings, normalizeMethod);
        vector<vector<float> > origCorrMatrix = originalData.getRho();
        
        string correlationFileName = getOutputFileName("corr", variables);
        ofstream correlationFile;
        m->openOutputFile(correlationFileName, correlationFile);
        outputNames.push_back(correlationFileName);  outputTypes["corr"].push_back(correlationFileName);
        correlationFile.setf(ios::fixed, ios::floatfield);
        correlationFile.setf(ios::showpoint);
        
        for(int i=0;i<numOTUs;i++){ correlationFile << '\t' << otuNames[i];    }   correlationFile << endl;
        for(int i=0;i<numOTUs;i++){
            correlationFile << otuNames[i];
            for(int j=0;j<numOTUs;j++){
                correlationFile << '\t' << origCorrMatrix[i][j];
            }
            correlationFile << endl;
        }
        
        
        if(numPermutations != 0){
            vector<vector<float> > pValues = createProcesses(sharedVector, origCorrMatrix);
            
            if (m->control_pressed) { return 0; }
            
            string pValueFileName = getOutputFileName("pvalue", variables);
            ofstream pValueFile;
            m->openOutputFile(pValueFileName, pValueFile);
            outputNames.push_back(pValueFileName);  outputTypes["pvalue"].push_back(pValueFileName);
            pValueFile.setf(ios::fixed, ios::floatfield);
            pValueFile.setf(ios::showpoint);
            
            for(int i=0;i<numOTUs;i++){ pValueFile << '\t' << otuNames[i];    }   pValueFile << endl;
            for(int i=0;i<numOTUs;i++){
                pValueFile << otuNames[i];
                for(int j=0;j<numOTUs;j++){
                    pValueFile << '\t' << pValues[i][j];
                }
                pValueFile << endl;
            }
        }


        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<vector<float> > SparccCommand::createProcesses(vector<vector<float> >& sharedVector, vector<vector<float> >& origCorrMatrix){
	try {
        int numOTUs = sharedVector[0].size();
        vector<vector<float> > pValues;
        
        if(processors == 1){
			pValues = driver(sharedVector, origCorrMatrix, numPermutations);
		}else{
            //divide iters between processors
			vector<int> procIters;
            int numItersPerProcessor = numPermutations / processors;
            
            //divide iters between processes
            for (int h = 0; h < processors; h++) {
                if(h == processors - 1){ numItersPerProcessor = numPermutations - h * numItersPerProcessor; }
                procIters.push_back(numItersPerProcessor);
            }
            
            vector<int> processIDS;
            int process = 1;
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            
			//loop through and create all the processes you want
			while (process != processors) {
				pid_t pid = fork();
				
				if (pid > 0) {
					processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					pValues = driver(sharedVector, origCorrMatrix, procIters[process]);
					
					//pass pvalues to parent
					ofstream out;
					string tempFile = m->mothurGetpid(process) + ".pvalues.temp";
					m->openOutputFile(tempFile, out);
					
					//pass values
					for (int i = 0; i < pValues.size(); i++) {
                        for (int j = 0; j < pValues[i].size(); j++) {
                            out << pValues[i][j] << '\t';
                        }
                        out << endl;
					}
					out << endl;
					
					out.close();
					
					exit(0);
				}else {
					m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
					for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
					exit(0);
				}
			}
			
			//do my part
			pValues = driver(sharedVector, origCorrMatrix, procIters[0]);
			
			//force parent to wait until all the processes are done
			for (int i=0;i<processIDS.size();i++) {
				int temp = processIDS[i];
				wait(&temp);
			}
			
			//combine results
			for (int i = 0; i < processIDS.size(); i++) {
				ifstream in;
				string tempFile =  toString(processIDS[i]) + ".pvalues.temp";
				m->openInputFile(tempFile, in);
				
				////// to do ///////////
				int numTemp; numTemp = 0;
				
                for (int j = 0; j < pValues.size(); j++) {
                    for (int k = 0; k < pValues.size(); k++) {
                        in >> numTemp; m->gobble(in);
                        pValues[j][k] += numTemp;
                    }
                    m->gobble(in);
				}
				in.close(); m->mothurRemove(tempFile);
			}
#else
            
            //fill in functions
            vector<sparccData*> pDataArray;
            DWORD   dwThreadIdArray[processors-1];
            HANDLE  hThreadArray[processors-1];
            
            //Create processor worker threads.
            for( int i=1; i<processors; i++ ){
                
                //make copy so we don't get access violations
                vector< vector<float> > copySharedVector = sharedVector;
                vector< vector<float> > copyOrig = origCorrMatrix;
                
                sparccData* temp = new sparccData(m, procIters[i], copySharedVector, copyOrig, numSamplings, maxIterations, numPermutations, normalizeMethod);
                pDataArray.push_back(temp);
                processIDS.push_back(i);
                
                hThreadArray[i-1] = CreateThread(NULL, 0, MySparccThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
            }
            
            //do my part
			pValues = driver(sharedVector, origCorrMatrix, procIters[0]);
            
            //Wait until all threads have terminated.
            WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
            
            //Close all thread handles and free memory allocations.
            for(int i=0; i < pDataArray.size(); i++){
                for (int j = 0; j < pDataArray[i]->pValues.size(); j++) {
                    for (int k = 0; k < pDataArray[i]->pValues[j].size(); k++) {
                        pValues[j][k] += pDataArray[i]->pValues[j][k];
                    }
                }
                
                CloseHandle(hThreadArray[i]);
                delete pDataArray[i];
            }
#endif
		}

        for(int i=0;i<numOTUs;i++){
            pValues[i][i] = 1;
            for(int j=0;j<i;j++){
                pValues[i][j] /= (double)numPermutations;
                pValues[j][i] = pValues[i][j];
            }
        }
        
        return pValues;
    }
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "createProcesses");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<vector<float> > SparccCommand::driver(vector<vector<float> >& sharedVector, vector<vector<float> >& origCorrMatrix, int numPerms){
	try {
        int numOTUs = sharedVector[0].size();
        vector<vector<float> > sharedShuffled = sharedVector;
        vector<vector<float> > pValues(numOTUs);
        for(int i=0;i<numOTUs;i++){ pValues[i].assign(numOTUs, 0);  }

        for(int i=0;i<numPerms;i++){
            if (m->control_pressed) { return pValues; }
            sharedShuffled = shuffleSharedVector(sharedVector);
            CalcSparcc permutedData(sharedShuffled, maxIterations, numSamplings, normalizeMethod);
            vector<vector<float> > permuteCorrMatrix = permutedData.getRho();
            
            for(int j=0;j<numOTUs;j++){
                for(int k=0;k<j;k++){
                    double randValue = permuteCorrMatrix[j][k];
                    double observedValue = origCorrMatrix[j][k];
                    if(observedValue >= 0 &&  randValue > observedValue)   { pValues[j][k]++; }//this method seems to deflate the
                    else if(observedValue < 0 && randValue < observedValue){ pValues[j][k]++; }//pvalues of small rho values
                }
            }
            if((i+1) % (int)(numPermutations * 0.05) == 0){ cout << i+1 << endl;  }
        }
        
        return pValues;
    }
	catch(exception& e) {
		m->errorOut(e, "SparccCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************

