//
//  sffmultiplecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sffmultiplecommand.h"
#include "sffinfocommand.h"
#include "seqsummarycommand.h"
#include "trimflowscommand.h"
#include "shhhercommand.h"
#include "trimseqscommand.h"


//**********************************************************************************************************************
vector<string> SffMultipleCommand::setParameters(){	
	try {		
		CommandParameter pfile("file", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfile);
        
        //sffinfo
		CommandParameter ptrim("trim", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(ptrim);
        
        //trim.flows
 		CommandParameter pmaxhomop("maxhomop", "Number", "", "9", "", "", "",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pmaxflows("maxflows", "Number", "", "450", "", "", "",false,false); parameters.push_back(pmaxflows);
		CommandParameter pminflows("minflows", "Number", "", "450", "", "", "",false,false); parameters.push_back(pminflows);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ptdiffs);
		CommandParameter psignal("signal", "Number", "", "0.50", "", "", "",false,false); parameters.push_back(psignal);
		CommandParameter pnoise("noise", "Number", "", "0.70", "", "", "",false,false); parameters.push_back(pnoise);
		CommandParameter porder("order", "String", "", "", "", "", "",false,false); parameters.push_back(porder);

        //shhh.flows
        CommandParameter plookup("lookup", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(plookup);
		CommandParameter pcutoff("cutoff", "Number", "", "0.01", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pmaxiter("maxiter", "Number", "", "1000", "", "", "",false,false); parameters.push_back(pmaxiter);
        CommandParameter plarge("large", "Number", "", "-1", "", "", "",false,false); parameters.push_back(plarge);
		CommandParameter psigma("sigma", "Number", "", "60", "", "", "",false,false); parameters.push_back(psigma);
		CommandParameter pmindelta("mindelta", "Number", "", "0.000001", "", "", "",false,false); parameters.push_back(pmindelta);
        
        //trim.seqs parameters
        CommandParameter pallfiles("allfiles", "Boolean", "", "t", "", "", "",false,false); parameters.push_back(pallfiles);
        CommandParameter pflip("flip", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pflip);
		CommandParameter pmaxambig("maxambig", "Number", "", "-1", "", "", "",false,false); parameters.push_back(pmaxambig);
		CommandParameter pminlength("minlength", "Number", "", "0", "", "", "",false,false); parameters.push_back(pminlength);
		CommandParameter pmaxlength("maxlength", "Number", "", "0", "", "", "",false,false); parameters.push_back(pmaxlength);
		CommandParameter pkeepforward("keepforward", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pkeepforward);
		CommandParameter pqtrim("qtrim", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pqtrim);
		CommandParameter pqthreshold("qthreshold", "Number", "", "0", "", "", "",false,false); parameters.push_back(pqthreshold);
		CommandParameter pqaverage("qaverage", "Number", "", "0", "", "", "",false,false); parameters.push_back(pqaverage);
		CommandParameter prollaverage("rollaverage", "Number", "", "0", "", "", "",false,false); parameters.push_back(prollaverage);
		CommandParameter pqwindowaverage("qwindowaverage", "Number", "", "0", "", "", "",false,false); parameters.push_back(pqwindowaverage);
		CommandParameter pqstepsize("qstepsize", "Number", "", "1", "", "", "",false,false); parameters.push_back(pqstepsize);
		CommandParameter pqwindowsize("qwindowsize", "Number", "", "50", "", "", "",false,false); parameters.push_back(pqwindowsize);
		CommandParameter pkeepfirst("keepfirst", "Number", "", "0", "", "", "",false,false); parameters.push_back(pkeepfirst);
		CommandParameter premovelast("removelast", "Number", "", "0", "", "", "",false,false); parameters.push_back(premovelast);

        
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffMultipleCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sff.multiple command reads a file containing sff filenames and optional oligos filenames. It runs the files through sffinfo, trim.flows, shhh.flows and trim.seqs combining the results.\n";
		helpString += "The sff.multiple command parameters are: ";
        vector<string> parameters = setParameters();
        for (int i = 0; i < parameters.size()-1; i++) {
            helpString += parameters[i] + ", ";
        }
        helpString += parameters[parameters.size()-1] + ".\n";
		helpString += "The file parameter allows you to enter the a file containing the list of sff files and optional oligos files.\n";
        helpString += "The trim parameter allows you to indicate if you would like a sequences and quality scores generated by sffinfo trimmed to the clipQualLeft and clipQualRight values.  Default=True. \n";
        helpString += "The maxambig parameter allows you to set the maximum number of ambigious bases allowed. The default is -1.\n";
		helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. \n";
		helpString += "The minlength parameter allows you to set and minimum sequence length. \n";
		helpString += "The maxlength parameter allows you to set and maximum sequence length. \n";
		helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The qfile parameter allows you to provide a quality file.\n";
		helpString += "The qthreshold parameter allows you to set a minimum quality score allowed. \n";
		helpString += "The qaverage parameter allows you to set a minimum average quality score allowed. \n";
		helpString += "The qwindowsize parameter allows you to set a number of bases in a window. Default=50.\n";
		helpString += "The qwindowaverage parameter allows you to set a minimum average quality score allowed over a window. \n";
		helpString += "The rollaverage parameter allows you to set a minimum rolling average quality score allowed over a window. \n";
		helpString += "The qstepsize parameter allows you to set a number of bases to move the window over. Default=1.\n";
		helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";
		helpString += "The keepforward parameter allows you to indicate whether you want the forward primer removed or not. The default is F, meaning remove the forward primer.\n";
		helpString += "The qtrim parameter will trim sequence from the point that they fall below the qthreshold and put it in the .trim file if set to true. The default is T.\n";
		helpString += "The keepfirst parameter trims the sequence to the first keepfirst number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements. \n";
		helpString += "The removelast removes the last removelast number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements.\n";

		helpString += "Example sff.multiple(file=mySffOligosFile.txt, trim=F).\n";
		helpString += "Note: No spaces between parameter labels (i.e. file), '=' and parameters (i.e.mySffOligosFile.txt).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffMultipleCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "getOutputFileNameTag");
		exit(1);
	}
}


//**********************************************************************************************************************
SffMultipleCommand::SffMultipleCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["flow"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "SffMultipleCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

SffMultipleCommand::SffMultipleCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
            map<string,string>::iterator it;
            
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["flow"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
			}
            
			filename = validParameter.validFile(parameters, "file", true);
            if (filename == "not open") { filename = ""; abort = true; }
            else if (filename == "not found") { filename = "";  }
			
			string temp;
			temp = validParameter.validFile(parameters, "trim", false);					if (temp == "not found"){	temp = "T";				}
			trim = m->isTrue(temp); 
            
            temp = validParameter.validFile(parameters, "minflows", false);	if (temp == "not found") { temp = "450"; }
			m->mothurConvert(temp, minFlows);  
            
			temp = validParameter.validFile(parameters, "maxflows", false);	if (temp == "not found") { temp = "450"; }
			m->mothurConvert(temp, maxFlows);  
            
            temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found"){	temp = "9";		}
			m->mothurConvert(temp, maxHomoP);  
            
			temp = validParameter.validFile(parameters, "signal", false);		if (temp == "not found"){	temp = "0.50";	}
			m->mothurConvert(temp, signal);  
            
			temp = validParameter.validFile(parameters, "noise", false);		if (temp == "not found"){	temp = "0.70";	}
			m->mothurConvert(temp, noise);  
            
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, pdiffs);
			
            temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, ldiffs);
            
            temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, sdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
            
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
			flowOrder = validParameter.validFile(parameters, "order", false);
			if (flowOrder == "not found"){ flowOrder = "TACG";		}
			else if(flowOrder.length() != 4){
				m->mothurOut("The value of the order option must be four bases long\n");
			}
            
            temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found"){	temp = "0.01";		}
			m->mothurConvert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "mindelta", false);	if (temp == "not found"){	temp = "0.000001";	}
			m->mothurConvert(temp, minDelta); 
            
			temp = validParameter.validFile(parameters, "maxiter", false);	if (temp == "not found"){	temp = "1000";		}
			m->mothurConvert(temp, maxIters); 
            
            temp = validParameter.validFile(parameters, "large", false);	if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, largeSize); 
            if (largeSize != 0) { large = true; }
            else { large = false;  }
            if (largeSize < 0) {  m->mothurOut("The value of the large cannot be negative.\n"); }
            
			temp = validParameter.validFile(parameters, "sigma", false);if (temp == "not found")	{	temp = "60";		}
			m->mothurConvert(temp, sigma); 
            
            temp = validParameter.validFile(parameters, "flip", false);
			if (temp == "not found")    {	flip = 0;	}
			else {  flip = m->isTrue(temp);		}
			
			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxAmbig);  
                       
			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, maxLength);
						
			temp = validParameter.validFile(parameters, "qthreshold", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, qThreshold);
			
			temp = validParameter.validFile(parameters, "qtrim", false);		if (temp == "not found") { temp = "t"; }
			qtrim = m->isTrue(temp);
            
			temp = validParameter.validFile(parameters, "rollaverage", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qRollAverage);
            
			temp = validParameter.validFile(parameters, "qwindowaverage", false);if (temp == "not found") { temp = "0"; }
			convert(temp, qWindowAverage);
            
			temp = validParameter.validFile(parameters, "qwindowsize", false);	if (temp == "not found") { temp = "50"; }
			convert(temp, qWindowSize);
            
			temp = validParameter.validFile(parameters, "qstepsize", false);	if (temp == "not found") { temp = "1"; }
			convert(temp, qWindowStep);
            
			temp = validParameter.validFile(parameters, "qaverage", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);
            
			temp = validParameter.validFile(parameters, "keepfirst", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);
            
			temp = validParameter.validFile(parameters, "removelast", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "keepforward", false);		if (temp == "not found") { temp = "F"; }
			keepforward = m->isTrue(temp);
	          
            numFPrimers = 0;
			numRPrimers = 0;
            numLinkers = 0;
            numSpacers = 0;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "SffMultipleCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		vector<string> sffFiles, oligosFiles;
        readFile(sffFiles, oligosFiles);
        
        if (m->control_pressed) { return 0; }
        
        if (sffFiles.size() < processors) { processors = sffFiles.size(); }
    
        if (processors == 1) { driver(sffFiles, oligosFiles, 0, sffFiles.size()); }
        else { createProcesses(sffFiles, oligosFiles); } 
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
		
		//report output filenames
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::readFile(vector<string>& sffFiles, vector<string>& oligosFiles){
	try {
        
        ifstream in;
        m->openInputFile(filename, in);
        
        string oligos, sff;
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            in >> sff;
            
            //ignore file pairing
            if(sff[0] == '#'){ while (!in.eof())	{	char c = in.get();  if (c == 10 || c == 13){	break;	}	} m->gobble(in); }
            else { //check for oligos file
                oligos = "";
            
                // get rest of line in case there is a oligos filename
                while (!in.eof())	{	
                    char c = in.get(); 
                    if (c == 10 || c == 13){	break;	}
                    else if (c == 32 || c == 9){;} //space or tab
                    else { 	oligos += c;  }
                } 
            }
            m->gobble(in);
            
            sffFiles.push_back(sff);
            oligosFiles.push_back(oligos); //will push a blank if there is not an oligos for this sff file
        }
        in.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::driver(vector<string> sffFiles, vector<string> oligosFiles, int start, int end){
    try {
        int count = 0;
        for (int i = start; i < end; i++) {
            string sff = sffFiles[i];
            string oligos = oligosFiles[i];
            
            //run sff.info
            string inputString = "sff=" + sff + ", flow=T";
            if (trim) { inputString += ", trim=T"; }
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            m->mothurOut("Running command: sffinfo(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* sffCommand = new SffInfoCommand(inputString);
            sffCommand->execute();
            
            map<string, vector<string> > filenames = sffCommand->getOutputFiles();
            
            delete sffCommand;
            m->mothurCalling = false;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            
            //run summary.seqs on the fasta file
            string fastaFile = "";
            map<string, vector<string> >::iterator it = filenames.find("fasta");
            if (it != filenames.end()) {  if ((it->second).size() != 0) { fastaFile = (it->second)[0];  } }
            else {  m->mothurOut("[ERROR]: sffinfo did not create a fasta file, quitting.\n"); m->control_pressed = true; break;  }
            
            inputString = "fasta=" + fastaFile;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            m->mothurOut("Running command: summary.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* summarySeqsCommand = new SeqSummaryCommand(inputString);
            summarySeqsCommand->execute();
            
            delete summarySeqsCommand;
            m->mothurCalling = false;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            
            count++;
        }
        
        return count;
    }
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::createProcesses(vector<string> sffFiles, vector<string> oligosFiles){
    try {
        vector<int> processIDS;
		int process = 1;
		int num = 0;
				
		//divide the groups between the processors
		vector<linePair> lines;
        vector<int> numFilesToComplete;
		int numFilesPerProcessor = sffFiles.size() / processors;
		for (int i = 0; i < processors; i++) {
			int startIndex =  i * numFilesPerProcessor;
			int endIndex = (i+1) * numFilesPerProcessor;
			if(i == (processors - 1)){	endIndex = sffFiles.size(); 	}
			lines.push_back(linePair(startIndex, endIndex));
            numFilesToComplete.push_back((endIndex-startIndex));
		}
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(sffFiles, oligosFiles, lines[process].start, lines[process].end);
                
                //pass numSeqs to parent
				ofstream out;
				string tempFile = toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
                
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		num = driver(sffFiles, oligosFiles, lines[0].start, lines[0].end);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
#else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the sffMultiplesData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		/*
         vector<shhhFlowsData*> pDataArray; 
         DWORD   dwThreadIdArray[processors-1];
         HANDLE  hThreadArray[processors-1]; 
         
         //Create processor worker threads.
         for( int i=0; i<processors-1; i++ ){
         // Allocate memory for thread data.
         string extension = "";
         if (i != 0) { extension = toString(i) + ".temp"; }
         
         shhhFlowsData* tempFlow = new shhhFlowsData(filenames, (compositeFASTAFileName + extension), (compositeNamesFileName + extension), outputDir, flowOrder, jointLookUp, singleLookUp, m, lines[i].start, lines[i].end, cutoff, sigma, minDelta, maxIters, i);
         pDataArray.push_back(tempFlow);
         processIDS.push_back(i);
         
         hThreadArray[i] = CreateThread(NULL, 0, ShhhFlowsThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
         }
         
         //using the main process as a worker saves time and memory
         //do my part
         driver(filenames, compositeFASTAFileName, compositeNamesFileName, lines[processors-1].start, lines[processors-1].end);
         
         //Wait until all threads have terminated.
         WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
         
         //Close all thread handles and free memory allocations.
         for(int i=0; i < pDataArray.size(); i++){
         for(int j=0; j < pDataArray[i]->outputNames.size(); j++){ outputNames.push_back(pDataArray[i]->outputNames[j]); }
         CloseHandle(hThreadArray[i]);
         delete pDataArray[i];
         }
         */
#endif
        
        for (int i=0;i<processIDS.size();i++) { 
            ifstream in;
			string tempFile = toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { 
                int tempNum = 0; 
                in >> tempNum; 
                if (tempNum != numFilesToComplete[i+1]) {
                    m->mothurOut("[ERROR]: main process expected " + toString(processIDS[i]) + " to complete " + toString(numFilesToComplete[i+1]) + " files, and it only reported completing " + toString(tempNum) + ". This will cause file mismatches.  The flow files may be too large to process with multiple processors. \n");
                }
            }
			in.close(); m->mothurRemove(tempFile);
        }
        
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************




