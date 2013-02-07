//
//  newcommandtemplate.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "newcommandtemplate.h"

//**********************************************************************************************************************
vector<string> NewCommand::setParameters(){	
	try {
		//eaxamples of each type of parameter. more info on the types of parameters can be found in commandparameter.h
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pprocessors);
        
        //files that have dependancies
        CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none","outputType",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName","outputType",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "ColumnName","outputType",false,false); parameters.push_back(pcolumn);		
        //files that do not have dependancies - fasta is set to not be required whereas shared is set to be required
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","outputType",false,false); parameters.push_back(pfasta);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","outputType",false,true); parameters.push_back(pshared);		

        
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		
        //choose more than one multiple options
        CommandParameter pcalc("calc", "Multiple", "jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-morisitahorn-braycurtis", "jest-thetayc", "", "", "","",true,false); parameters.push_back(pcalc);
        //choose only one multiple options
        CommandParameter pdistance("distance", "Multiple", "column-lt-square", "column", "", "", "","",false,false); parameters.push_back(pdistance);
        
        CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptiming);
        
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
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
string NewCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The new command allows you to ....\n";
		helpString += "The new command parameters are: ....\n";
		helpString += "The whatever parameter is used to ....\n";
		helpString += "The new command should be in the following format: \n";
		helpString += "new(...)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string NewCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fileType1") {  pattern = "[filename],tag1"; }
        else if (type == "fileType2") {  pattern = "[filename],tag2"; }
        else if (type == "fileType3") {  pattern = "[filename],tag3"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "NewCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
NewCommand::NewCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["fileType1"] = tempOutNames; //filetypes should be things like: shared, fasta, accnos...
		outputTypes["fileType2"] = tempOutNames;
		outputTypes["FileType3"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "NewCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
NewCommand::NewCommand(string option)  {
	try {
////////////////////////////////////////////////////////
/////////////////// start leave alone block ////////////
////////////////////////////////////////////////////////
		abort = false; calledHelp = false;   
		
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
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
///////////////////////////////////////////////////////////////
//////////////// stop leave alone block ///////////////////////
///////////////////////////////////////////////////////////////
                
//edit file types below to include only the types you added as parameters
                
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
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
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
///////////////////////////////////////////////////////////////////////////////
/////////// example of getting filenames and checking dependancies ////////////
// the validParameter class will make sure file exists, fill with correct    //
// and name is current is given ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
            
            
            ///variables for examples below that you will most likely want to put in the header for 
            //use by the other class functions.
            string phylipfile, columnfile, namefile, fastafile, sharedfile, method, countfile;
            int processors;
            bool useTiming, allLines;
            vector<string> Estimators, Groups;
            set<string> labels;
            //if allLines is used it should be initialized to 1 above.
            
            
			//check for parameters
            phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else { 	m->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {   m->setColumnFile(columnfile);	}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
            
            //get fastafile - it is not required
            fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort=true;  }
			else if (fastafile == "not found") {  fastafile = "";  }
			if (fastafile != "") { m->setFastaFile(fastafile); }

			
			if ((phylipfile == "") && (columnfile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {   m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file before you can use the cluster command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
			
			if (columnfile != "") {
				if (namefile == "") { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format."); m->mothurOutEndLine(); 
						abort = true; 
					}	
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

			
//////////////////////////////////////////////////////////////////////
////////// example of getting other types of parameters //////////////
//////////////////////////////////////////////////////////////////////
            
			//use only one Mutliple type
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, and weighted."); m->mothurOutEndLine(); abort = true; }
            
            //use more than one multiple type. do not check to make sure the entry is valid.
            string calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
                if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			m->splitAtDash(calc, Estimators);
            
            //Boolean type - m->isTrue looks for t, true, f or false and is case insensitive
			string timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
            useTiming = m->isTrue(timing);
			
            //Number type - mothurConvert makes sure the convert can happen to avoid a crash.
            string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
            //Groups must be checked later to make sure they are valid. SharedUtilities has functions of check the validity, just make to so m->setGroups() after the checks.  If you are using these with a shared file no need to check the SharedRAbundVector class will call SharedUtilites for you, kinda nice, huh?
            string groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { m->splitAtDash(groups, Groups); }
			m->setGroups(Groups);
            
            //Commonly used to process list, rabund, sabund, shared and relabund files.  Look at "smart distancing" examples below in the execute function.
            string label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
            
            //if your command has a namefile as an option, you may want ot check to see if there is a current namefile
            //saved by mothur that is associated with the other files you are using as inputs.  
            //You can do so by adding the files associated with the namefile to the files vector and then asking parser to check.  
            //This saves our users headaches over file mismatches because they forgot to include the namefile, :)
            if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    parser.getNameFile(files);
                }
            }

			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "NewCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int NewCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        // reading and processing a shared file code example
        // Note: As long as you set groups and labels as shown in the constructor, you can use this code without modification other than adding your function call which is passed the lookup vector.
        // The classes used below will handle the checking of groups to make sure they are valid and returning only the groups you selected.  The while loop implements mothur "smart distancing" so as long as you filled label as shown above in the constructor the code below will handle bad labels or labels not included in the sharedfile.
         
         //Reads sharefile, binLabels are stored in m->currentBinLabels, lookup will be filled with groups in m->getGroups() or all groups in file if m->getGroups is empty. If groups are selected, some bins maybe eliminated if they only contained seqs from groups not included. No need to worry about the details of this, SharedRAbundVector takes care of it.  Just make sure to use m->currentBinLabels if you are outputting OTU labels so that if otus are eliminated you still have the correct names.
        
        /*
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
         
         ///////////////////////////////////////////////////////////////////////////////////
         //// Call your function to process specific distance in sharedfile, ie lookup /////
         ///////////////////////////////////////////////////////////////////////////////////
         
         processedLabels.insert(lookup[0]->getLabel());
         userLabels.erase(lookup[0]->getLabel());
         }
         
         if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
         string saveLabel = lookup[0]->getLabel();
         
         for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
         lookup = input.getSharedRAbundVectors(lastLabel);
         m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
         
         ///////////////////////////////////////////////////////////////////////////////////
         //// Call your function to process specific distance in sharedfile, ie lookup /////
         ///////////////////////////////////////////////////////////////////////////////////

         
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
         
         ///////////////////////////////////////////////////////////////////////////////////
         //// Call your function to process specific distance in sharedfile, ie lookup /////
         ///////////////////////////////////////////////////////////////////////////////////

         
         for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
         }
         */
        
        
        
        //if you make a new file or a type that mothur keeps track of the current version, you can update it with something like the following.
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
		}
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "NewCommand");
		exit(1);
	}
}
//**********************************************************************************************************************


