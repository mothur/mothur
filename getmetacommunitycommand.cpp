//
//  getmetacommunitycommand.cpp
//  Mothur
//
//  Created by SarahsWork on 4/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "getmetacommunitycommand.h"
#include "qFinderDMM.h"

//**********************************************************************************************************************
vector<string> GetMetaCommunityCommand::setParameters(){
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","outputType",false,true); parameters.push_back(pshared);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pminpartitions("minpartitions", "Number", "", "5", "", "", "","",false,false,true); parameters.push_back(pminpartitions);
        CommandParameter pmaxpartitions("maxpartitions", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pmaxpartitions);
        CommandParameter poptimizegap("optimizegap", "Number", "", "3", "", "", "","",false,false,true); parameters.push_back(poptimizegap);
   		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMetaCommunityCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The get.metacommunity command parameters are shared, label, groups, minpartitions, maxpartitions and optimizegap. The shared file is required. \n";
        helpString += "The label parameter is used to analyze specific labels in your input. labels are separated by dashes.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your shared file you would like analyzed.  Group names are separated by dashes.\n";
		helpString += "The minpartitions parameter is used to .... Default=5.\n";
        helpString += "The maxpartitions parameter is used to .... Default=10.\n";
        helpString += "The optimizegap parameter is used to .... Default=3.\n";
		helpString += "The get.metacommunity command should be in the following format: get.metacommunity(shared=yourSharedFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMetaCommunityCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fit")              {  pattern = "[filename],[distance],mix.fit"; }
        else if (type == "relabund")    {  pattern = "[filename],[distance],[tag],mix.relabund"; }
        else if (type == "design")      {  pattern = "[filename],mix.design"; }
        else if (type == "matrix")      {  pattern = "[filename],[distance],[tag],mix.posterior"; }
        else if (type == "parameters")  {  pattern = "[filename],[distance],mix.parameters"; }
        else if (type == "summary")  {  pattern = "[filename],[distance],mix.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetMetaCommunityCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetMetaCommunityCommand::GetMetaCommunityCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["fit"] = tempOutNames;
        outputTypes["relabund"] = tempOutNames;
        outputTypes["matrix"] = tempOutNames;
        outputTypes["design"] = tempOutNames;
        outputTypes["parameters"] = tempOutNames;
        outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "GetMetaCommunityCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetMetaCommunityCommand::GetMetaCommunityCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        allLines=true;
		
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
            outputTypes["fit"] = tempOutNames;
            outputTypes["relabund"] = tempOutNames;
            outputTypes["matrix"] = tempOutNames;
            outputTypes["design"] = tempOutNames;
            outputTypes["parameters"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                string path;
                it = parameters.find("shared");
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
                       
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
            
            string temp = validParameter.validFile(parameters, "minpartitions", false);	if (temp == "not found"){	temp = "5";      }
			m->mothurConvert(temp, minpartitions);
            
            temp = validParameter.validFile(parameters, "maxpartitions", false);        if (temp == "not found"){	temp = "10";	 }
			m->mothurConvert(temp, maxpartitions);
            
            temp = validParameter.validFile(parameters, "optimizegap", false);          if (temp == "not found"){	temp = "3";	 }
			m->mothurConvert(temp, optimizegap);
            
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
		m->errorOut(e, "GetMetaCommunityCommand", "GetMetaCommunityCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetMetaCommunityCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        InputData input(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
        string lastLabel = lookup[0]->getLabel();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
            
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
        
        if (m->control_pressed) {  return 0; }
        
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
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetMetaCommunityCommand::process(vector<SharedRAbundVector*>& thislookup){
	try {
        
        double minLaplace = 1e10;
        int minPartition = 0;
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = thislookup[0]->getLabel();
		string outputFileName = getOutputFileName("fit", variables);
        outputNames.push_back(outputFileName); outputTypes["fit"].push_back(outputFileName);
        
		ofstream fitData;
		m->openOutputFile(outputFileName, fitData);
        fitData.setf(ios::fixed, ios::floatfield);
        fitData.setf(ios::showpoint);
        cout.setf(ios::fixed, ios::floatfield);
        cout.setf(ios::showpoint);

        vector< vector<int> > sharedMatrix;
        for (int i = 0; i < thislookup.size(); i++) { sharedMatrix.push_back(thislookup[i]->getAbundances()); }
        
        m->mothurOut("K\tNLE\t\tlogDet\tBIC\t\tAIC\t\tLaplace\n");
        fitData << "K\tNLE\tlogDet\tBIC\tAIC\tLaplace" << endl;
        
        for(int numPartitions=1;numPartitions<=maxpartitions;numPartitions++){
            
            if (m->control_pressed) { break; }
            
            qFinderDMM findQ(sharedMatrix, numPartitions);
            
            double laplace = findQ.getLaplace();
            m->mothurOut(toString(numPartitions) + '\t');
            cout << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
            m->mothurOutJustToLog(toString(findQ.getNLL()) + '\t' + toString(findQ.getLogDet()) + '\t');
            cout << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace;
            m->mothurOutJustToLog(toString(findQ.getBIC()) + '\t' + toString(findQ.getAIC()) + '\t' + toString(laplace));
            
            fitData << numPartitions << '\t';
            fitData << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
            fitData << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace << endl;
            
            if(laplace < minLaplace){
                minPartition = numPartitions;
                minLaplace = laplace;
                m->mothurOut("***");
            }
            m->mothurOutEndLine();
            
            variables["[tag]"] = toString(numPartitions);
            string relabund = getOutputFileName("relabund", variables);
            outputNames.push_back(relabund); outputTypes["relabund"].push_back(relabund);
            string matrix = getOutputFileName("matrix", variables);
            outputNames.push_back(matrix); outputTypes["matrix"].push_back(matrix);
            
            findQ.printZMatrix(matrix, m->getGroups());
            findQ.printRelAbund(relabund, m->currentBinLabels);
            
            if(optimizegap != -1 && (numPartitions - minPartition) >= optimizegap && numPartitions >= minpartitions){ break;  }
        }
        fitData.close();
        
        //minPartition = 4;
        
        if (m->control_pressed) { return 0; }
        
        generateSummaryFile(minpartitions, outputTypes["relabund"][0], outputTypes["relabund"][outputTypes["relabund"].size()-1], outputTypes["matrix"][outputTypes["matrix"].size()-1], thislookup[0]->getLabel());

        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/

vector<double> GetMetaCommunityCommand::generateDesignFile(int numPartitions, string input){
    try {
        vector<double> piValues(numPartitions, 0);
        
        ifstream postFile;
        m->openInputFile(input, postFile);//((fileRoot + toString(numPartitions) + "mix.posterior").c_str()); //matrix file
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(input));
		string outputFileName = getOutputFileName("design", variables);
        ofstream designFile;
        m->openOutputFile(outputFileName, designFile);
        outputNames.push_back(outputFileName); outputTypes["design"].push_back(outputFileName);
        
        
        vector<string> titles(numPartitions);
        
        for(int i=0;i<numPartitions;i++){   postFile >> titles[i];  }
        
        double posterior;
        string sampleName;
        int numSamples = 0;
        
        while(postFile){
            
            if (m->control_pressed) { break; }
            
            double maxPosterior = 0.0000;
            int maxPartition = -1;
            
            postFile >> sampleName;
            
            for(int i=0;i<numPartitions;i++){
                
                postFile >> posterior;
                if(posterior > maxPosterior){
                    maxPosterior = posterior;
                    maxPartition = i;
                }
                piValues[i] += posterior;
                
            }
            
            designFile << sampleName << '\t' << titles[maxPartition] << endl;
            
            numSamples++;
            m->gobble(postFile);
        }
        for(int i=0;i<numPartitions;i++){
            piValues[i] /= (double)numSamples;
        }
        
        
        postFile.close();
        designFile.close();
        
        return piValues;
    }
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "generateDesignFile");
		exit(1);
	}
}

/**************************************************************************************************/

inline bool summaryFunction(summaryData i, summaryData j){ return i.difference > j.difference;   }

/**************************************************************************************************/
int GetMetaCommunityCommand::generateSummaryFile(int numPartitions, string reference, string partFile, string designInput, string label){
    try {
        vector<summaryData> summary;
        
        vector<double> pMean(numPartitions, 0);
        vector<double> pLCI(numPartitions, 0);
        vector<double> pUCI(numPartitions, 0);
        
        string name, header;
        double mean, lci, uci;
        
        
        vector<double> piValues = generateDesignFile(numPartitions, designInput);
        
        ifstream referenceFile;
        m->openInputFile(reference, referenceFile); //((fileRoot + label + ".1mix.relabund").c_str());
        ifstream partitionFile;
        m->openInputFile(partFile, partitionFile); //((fileRoot + toString(numPartitions) + "mix.relabund").c_str());
        
        header = m->getline(referenceFile);
        header = m->getline(partitionFile);
        stringstream head(header);
        string dummy, label;
        head >> dummy;
        vector<string> thetaValues(numPartitions, "");
        for(int i=0;i<numPartitions;i++){
            head >> label >> dummy >> dummy;
            thetaValues[i] = label.substr(label.find_last_of('_')+1);
        }
        
        
        vector<double> partitionDiff(numPartitions, 0.0000);
        
        while(referenceFile){
            
            if (m->control_pressed) { break; }
            
            referenceFile >> name >> mean >> lci >> uci;
            
            summaryData tempData;
            tempData.name = name;
            tempData.refMean = mean;
            
            double difference = 0.0000;
            
            partitionFile >> name;
            for(int j=0;j<numPartitions;j++){
                partitionFile >> pMean[j] >> pLCI[j] >> pUCI[j];
                difference += abs(mean - pMean[j]);
                partitionDiff[j] += abs(mean - pMean[j]);;
            }
            
            tempData.partMean = pMean;
            tempData.partLCI = pLCI;
            tempData.partUCI = pUCI;
            tempData.difference = difference;
            summary.push_back(tempData);
            
            m->gobble(referenceFile);
            m->gobble(partitionFile);
        }
        referenceFile.close();
        partitionFile.close();
        
        if (m->control_pressed) { return 0; }
        
        int numOTUs = (int)summary.size();
        
        sort(summary.begin(), summary.end(), summaryFunction);
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = label;
		string outputFileName = getOutputFileName("parameters", variables);
        outputNames.push_back(outputFileName); outputTypes["parameters"].push_back(outputFileName);
        
        ofstream parameterFile;
        m->openOutputFile(outputFileName, parameterFile); //((fileRoot + "mix.parameters").c_str());
        parameterFile.setf(ios::fixed, ios::floatfield);
        parameterFile.setf(ios::showpoint);
        
        double totalDifference =  0.0000;
        parameterFile << "Part\tDif2Ref_i\ttheta_i\tpi_i\n";
        for(int i=0;i<numPartitions;i++){
            if (m->control_pressed) { break; }
            parameterFile << i+1 << '\t' << setprecision(2) << partitionDiff[i] << '\t' << thetaValues[i] << '\t' << piValues[i] << endl;
            totalDifference += partitionDiff[i];
        }
        parameterFile.close();
        
        if (m->control_pressed) { return 0; }
        
        string summaryFileName = getOutputFileName("summary", variables);
        outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
        
        ofstream summaryFile;
        m->openOutputFile(summaryFileName, summaryFile); //((fileRoot + "mix.summary").c_str());
        summaryFile.setf(ios::fixed, ios::floatfield);
        summaryFile.setf(ios::showpoint);
        
        
        summaryFile << "OTU\tP0.mean";
        for(int i=0;i<numPartitions;i++){
            summaryFile << "\tP" << i+1 << ".mean\tP" << i+1 << ".lci\tP" << i+1 << ".uci";
        }
        summaryFile << "\tDifference\tCumFraction" << endl;
        
        double cumDiff = 0.0000;
        
        for(int i=0;i<numOTUs;i++){
            if (m->control_pressed) { break; }
            summaryFile << summary[i].name << setprecision(2) << '\t' << summary[i].refMean;
            for(int j=0;j<numPartitions;j++){
                summaryFile  << '\t' << summary[i].partMean[j] << '\t' << summary[i].partLCI[j] << '\t' << summary[i].partUCI[j];
            }
            
            cumDiff += summary[i].difference/totalDifference;
            summaryFile << '\t' << summary[i].difference << '\t' << cumDiff << endl;
        }
        summaryFile.close();
        
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "GetMetaCommunityCommand", "generateSummaryFile");
		exit(1);
	}
    
}
//**********************************************************************************************************************

