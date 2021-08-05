//
//  newcommandtemplate.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

//
#include "newcommandtemplate.h"

// Test Change.
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
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        //set output file types
        vector<string> tempOutNames;
        outputTypes["fileType1"] = tempOutNames; //filetypes should be things like: shared, fasta, accnos...
        outputTypes["fileType2"] = tempOutNames;
        outputTypes["FileType3"] = tempOutNames;
        
        //set abort and called Help
        abort = false; calledHelp = false;
		
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "NewCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
NewCommand::NewCommand(string option) : Command()  {
	try {
////////////////////////////////////////////////////////
/////////////////// start leave alone block ////////////
////////////////////////////////////////////////////////
		  
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			//valid paramters for this command
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
            ///variables for examples below that you will most likely want to put in the header for 
            //use by the other class functions.
            string phylipfile, columnfile, namefile, fastafile, sharedfile, method, countfile;
            int processors;
            bool useTiming, allLines;
            vector<string> Estimators, Groups;
            set<string> labels;
            //if allLines is used it should be initialized to 1 above.
            
            
			//check for parameters
            phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else { 	current->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {   current->setColumnFile(columnfile);	}
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
            
            //get fastafile - it is not required
            fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort=true;  }
			else if (fastafile == "not found") {  fastafile = "";  }
			if (fastafile != "") { current->setFastaFile(fastafile); }

			
			if ((phylipfile == "") && (columnfile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = current->getColumnFile(); 
				if (columnfile != "") {   m->mothurOut("Using " + columnfile + " as input file for the column parameter.\n");  }
				else { 
					phylipfile = current->getPhylipFile(); 
					if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter.\n");  }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file before you can use the cluster command.\n");  
						abort = true;
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip or column.\n");  abort = true; }
			
			if (columnfile != "") {
				if (namefile == "") { 
					namefile = current->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format.\n");  
						abort = true; 
					}	
				}
			}
            
            //get shared file, it is required
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required.\n");  abort = true; }
			}else { current->setSharedFile(sharedfile); }
            
			
//////////////////////////////////////////////////////////////////////
////////// example of getting other types of parameters //////////////
//////////////////////////////////////////////////////////////////////
            
			//use only one Mutliple type
			method = validParameter.valid(parameters, "method");
			if (method == "not found") { method = "average"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, and weighted.\n");  abort = true; }
            
            //use more than one multiple type. do not check to make sure the entry is valid.
            string calc = validParameter.valid(parameters, "calc");			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
                if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			util.splitAtDash(calc, Estimators);
            
            //Boolean type - m->isTrue looks for t, true, f or false and is case insensitive
			string timing = validParameter.valid(parameters, "timing");
			if (timing == "not found") { timing = "F"; }
            useTiming = util.isTrue(timing);
			
            //Number type - mothurConvert makes sure the convert can happen to avoid a crash.
            string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
            
            //Groups must be checked later to make sure they are valid. SharedUtilities has functions of check the validity, just make to so m->setGroups() after the checks.  If you are using these with a shared file no need to check the SharedRAbundVector class will call SharedUtilites for you, kinda nice, huh?
            string groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else {
                util.splitAtDash(groups, Groups);
                    if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
            }
            
            //Commonly used to process list, rabund, sabund, shared and relabund files.  Look at "smart distancing" examples below in the execute function.
            string label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
            
            //if your command has a namefile as an option, you may want ot check to see if there is a current namefile
            //saved by mothur that is associated with the other files you are using as inputs.  
            //You can do so by adding the files associated with the namefile to the files vector and then asking parser to check.  
            //This saves our users headaches over file mismatches because they forgot to include the namefile, :)
            if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        /*
         InputData input(inputFileName, format, Groups);
         set<string> processedLabels;
         set<string> userLabels = labels;
         string lastLabel = "";
         
         if (format == "relabund") {
             SharedRAbundFloatVectors* lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
         Groups = lookup->getNamesGroups();
             
             while (lookup != NULL) {
                 
                 if (m->getControl_pressed()) { delete lookup; break; }
                 
         //////// myfunction(lookup); - call your function to process relabund data ////////////////////
                 
                delete lookup;
                 
                 lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
             }

         }else if (format == "sharedfile") {
         
             SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
         Groups = lookup->getNamesGroups();
             
             while (lookup != NULL) {
                 
                 if (m->getControl_pressed()) { delete lookup; break; }
                 
                 //////// myfunction(lookup); - call your function to process shared data ////////////////////
                delete lookup;
                 
                 lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
             }
             
         }else if (format == "list") {
             ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
                    
             while (list != NULL) {
                        
                 if (m->getControl_pressed()) { delete list; break; }
                        
                 //////// myfunction(list); - call your function to process list data //////////////////// delete list;
                       
                 list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
             }
             
         }else if (format == "rabund") {
             RAbundVector* rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
                    
             while (rabund != NULL) {
                        
                 if (m->getControl_pressed()) { delete rabund; break; }
                        
                 //////// myfunction(rabund); - call your function to process list data //////////////////// delete rabund;
                       
                 rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
             }
             
         }
        */
        
        //if you make a new file or a type that mothur keeps track of the current version, you can update it with something like the following.
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; current->setFastaFile(currentFasta); }
		}
		
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "NewCommand");
		exit(1);
	}
}
//**********************************************************************************************************************


