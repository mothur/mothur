//
//  sortseqscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 2/3/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sortseqscommand.h"
#include "sequence.hpp"
#include "qualityscores.h"

//**********************************************************************************************************************
vector<string> SortSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none","fasta",false,false); parameters.push_back(pfasta);
        CommandParameter pflow("flow", "InputTypes", "", "", "none", "FNGLT", "none","flow",false,false); parameters.push_back(pflow);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none","name",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none","count",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none","group",false,false); parameters.push_back(pgroup);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none","taxonomy",false,false); parameters.push_back(ptaxonomy);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none","qfile",false,false); parameters.push_back(pqfile);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(plarge);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SortSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sort.seqs command puts the sequences in the same order for the following file types: accnos fasta, name, group, count, taxonomy, flow or quality file.\n";
        helpString += "The sort.seqs command parameters are accnos, fasta, name, group, count, taxonomy, flow, qfile and large.\n";
        helpString += "The accnos file allows you to specify the order you want the files in.  If none is provided, mothur will use the order of the first file it reads.\n";
        helpString += "The large parameters is used to indicate your files are too large to fit in RAM.\n";
		helpString += "The sort.seqs command should be in the following format: sort.seqs(fasta=yourFasta).\n";
		helpString += "Example sort.seqs(fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SortSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "flow")        {   pattern = "[filename],sorted,[extension]";    }
        else if (type == "qfile")      {   pattern = "[filename],sorted,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SortSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SortSeqsCommand::SortSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["flow"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "SortSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SortSeqsCommand::SortSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["flow"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
                
                it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
            
			
			//check for parameters
            accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { accnosfile = ""; abort = true; }
			else if (accnosfile == "not found") {  accnosfile = "";  }	
			else { current->setAccnosFile(accnosfile); }
            
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			else { current->setFastaFile(fastafile); }
            
            flowfile = validParameter.validFile(parameters, "flow");
			if (flowfile == "not open") { flowfile = ""; abort = true; }
			else if (flowfile == "not found") {  flowfile = "";  }	
			else { current->setFlowFile(flowfile); }
            
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { current->setNameFile(namefile); } 
            
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }
			else { current->setGroupFile(groupfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { current->setTaxonomyFile(taxfile); }
			
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") {  qualfile = "";  }			
			else { current->setQualFile(qualfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
            string temp = validParameter.valid(parameters, "large");		if (temp == "not found") { temp = "f"; }
			large = util.isTrue(temp);
            
			if ((fastafile == "") && (namefile == "") && (countfile == "") && (groupfile == "") && (taxfile == "") && (flowfile == "") && (qualfile == ""))  { m->mothurOut("You must provide at least one of the following: fasta, name, group, count, taxonomy, flow or quality."); m->mothurOutEndLine(); abort = true; }
			
            if (countfile == "") {
                if ((fastafile != "") && (namefile == "")) {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "SortSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SortSeqsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//read through the correct file and output lines you want to keep
        if (accnosfile != "")		{		
            vector<string> temp;
            util.readAccnos(accnosfile, temp);
            for (int i = 0; i < temp.size(); i++) {  names[temp[i]] = i;  }
            m->mothurOut("\nUsing " + accnosfile + " to determine the order. It contains " + toString(temp.size()) + " representative sequences.\n");	
        }
        
		if (fastafile != "")		{		readFasta();	}
        if (flowfile != "")         {		readFlow();     }
        if (qualfile != "")			{		readQual();		}
        if (namefile != "")			{		readName();		}
		if (groupfile != "")		{		readGroup();	}
        if (countfile != "")		{		readCount();	}
        if (taxfile != "")			{		readTax();		}
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string currentName = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
			}
			
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
			}
			
			itTypes = outputTypes.find("qfile");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
			}	
            
            itTypes = outputTypes.find("flow");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFlowFile(currentName); }
			}
            
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
			}
		}
		
		return 0;		
	}
    
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int SortSeqsCommand::readFasta(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(fastafile, in);
		string name;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    Sequence currSeq(in);
                    name = currSeq.getName();
                    
                    if (name != "") {
                        numNamesInFile++;
                        map<string, int>::iterator it = names.find(name);
                        if (it == names.end()) { 
                            names[name] = numNames; numNames++;
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();
                out.close();
                
                int numLeft = names.size();
                if (numNamesInFile < numLeft) { numLeft = numNamesInFile; }
                
                int size = 1000; //assume that user can hold 1000 seqs in memory
                if (numLeft < size) { size = numLeft; }
                int times = 0;
                
                vector<Sequence> seqs; seqs.resize(size);
                for (int i = 0; i < seqs.size(); i++) { seqs[i].setName(""); } //this is so if some of the seqs are missing we dont print out garbage
                
                while (numLeft > 0) {
                    
                    ifstream in2;
                    util.openInputFile(fastafile, in2);
                    
                    if (m->getControl_pressed()) { in2.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->getControl_pressed()) { in2.close();   util.mothurRemove(outputFileName);  return 0; }
                        
                        //stop reading if we already found the seqs we are looking for
                        if (found >= needToFind) { break; }
                        
                        Sequence currSeq(in2);
                        name = currSeq.getName();
                        
                        if (name != "") {
                            map<string, int>::iterator it = names.find(name);
                            if (it != names.end()) { //we found it, so put it in the vector in the right place.
                                //is it in the set of seqs we are looking for this time around
                                int thisSeqsPlace = it->second;
                                thisSeqsPlace -= (times * size);
                                if ((thisSeqsPlace < size) && (thisSeqsPlace >= 0)) {
                                    seqs[thisSeqsPlace] = currSeq; 
                                    found++;
                                }
                            }else { m->mothurOut("[ERROR]: in logic of readFasta function.\n"); m->setControl_pressed(true); }
                        }
                        util.gobble(in2);
                    }
                    in2.close();	

                    ofstream out2;
                    util.openOutputFileAppend(outputFileName, out2);
                    
                    int output = seqs.size();
                    if (numLeft < seqs.size()) { output = numLeft; }
                        
                    for (int i = 0; i < output; i++) {
                        if (seqs[i].getName() != "") { seqs[i].printSequence(out2); }
                    }
                    out2.close();
                    
                    times++;
                    numLeft -= output;
                }
                
                m->mothurOut("Ordered " + toString(numNamesInFile) + " sequences from " + fastafile + ".\n");
            }else {
                
                vector<Sequence> seqs; seqs.resize(names.size());
                for (int i = 0; i < seqs.size(); i++) { seqs[i].setName(""); } //this is so if some of the seqs are missing we dont print out garbage
                
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    Sequence currSeq(in);
                    name = currSeq.getName();
                    
                    if (name != "") {
                        map<string, int>::iterator it = names.find(name);
                        if (it != names.end()) { //we found it, so put it in the vector in the right place.
                            seqs[it->second] = currSeq;  
                        }else { //if we cant find it then add it to the end
                            names[name] = seqs.size();
                            seqs.push_back(currSeq);
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();	
                
                int count = 0;
                for (int i = 0; i < seqs.size(); i++) {
                    if (seqs[i].getName() != "") {
                        seqs[i].printSequence(out); count++;
                    }
                }
                out.close();
                
                m->mothurOut("Ordered " + toString(count) + " sequences from " + fastafile + ".\n");
            }
                        
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                Sequence currSeq(in);
                name = currSeq.getName();
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    currSeq.printSequence(out);
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + fastafile + " to determine the order. It contains " + toString(count) + " sequences.\n");
        }
				
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int SortSeqsCommand::readFlow(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(flowfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(flowfile));
        variables["[extension]"] = util.getExtension(flowfile);
		string outputFileName = getOutputFileName("flow", variables);
		outputTypes["flow"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(flowfile, in);
        int numFlows;
		string name;
        
        in >> numFlows; util.gobble(in);
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    in >> name;	
                    string rest = util.getline(in);
                    
                    if (name != "") {
                        numNamesInFile++;
                        map<string, int>::iterator it = names.find(name);
                        if (it == names.end()) { 
                            names[name] = numNames; numNames++;
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();
                out.close();
                
                int numLeft = names.size();
                if (numNamesInFile < numLeft) { numLeft = numNamesInFile; }
                
                int size = 1000; //assume that user can hold 1000 seqs in memory
                if (numLeft < size) { size = numLeft; }
                int times = 0;
                
                vector<string> seqs; seqs.resize(size, "");
                
                while (numLeft > 0) {
                    
                    ifstream in2;
                    util.openInputFile(flowfile, in2); in2 >> numFlows; util.gobble(in2);
                    
                    if (m->getControl_pressed()) { in2.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->getControl_pressed()) { in2.close();   util.mothurRemove(outputFileName);  return 0; }
                        
                        //stop reading if we already found the seqs we are looking for
                        if (found >= needToFind) { break; }
                        
                        in2 >> name;	
                        string rest = util.getline(in2);
                        
                        if (name != "") {
                            map<string, int>::iterator it = names.find(name);
                            if (it != names.end()) { //we found it, so put it in the vector in the right place.
                                //is it in the set of seqs we are looking for this time around
                                int thisSeqsPlace = it->second;
                                thisSeqsPlace -= (times * size);
                                if ((thisSeqsPlace < size) && (thisSeqsPlace >= 0)) {
                                    seqs[thisSeqsPlace] = (name +'\t' + rest); 
                                    found++;
                                }
                            }else { m->mothurOut("[ERROR]: in logic of readFlow function.\n"); m->setControl_pressed(true); }
                        }
                        util.gobble(in2);
                    }
                    in2.close();	
                    
                    ofstream out2;
                    util.openOutputFileAppend(outputFileName, out2);
                    
                    int output = seqs.size();
                    if (numLeft < seqs.size()) { output = numLeft; }
                    
                    for (int i = 0; i < output; i++) {
                        if (seqs[i] != "") {
                            out2 << seqs[i] << endl;
                        }
                    }
                    out2.close();
                    
                    times++;
                    numLeft -= output;
                }
                
                m->mothurOut("Ordered " + toString(numNamesInFile) + " flows from " + flowfile + ".\n");
            }else {
                
                vector<string> seqs; seqs.resize(names.size(), "");
                
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    in >> name;	
                    string rest = util.getline(in);
                    
                    if (name != "") {
                        map<string, int>::iterator it = names.find(name);
                        if (it != names.end()) { //we found it, so put it in the vector in the right place.
                            seqs[it->second] = (name + '\t' + rest);  
                        }else { //if we cant find it then add it to the end
                            names[name] = seqs.size();
                            seqs.push_back((name + '\t' + rest));
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();	
                
                int count = 0;
                for (int i = 0; i < seqs.size(); i++) {
                    if (seqs[i] != "") {
                        out << seqs[i] << endl;
                        count++;
                    }
                }
                out.close();
                
                m->mothurOut("Ordered " + toString(count) + " flows from " + flowfile + ".\n");
            }
            
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> name;	
                string rest = util.getline(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << rest << endl;
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + flowfile + " to determine the order. It contains " + toString(count) + " flows.\n");
        }
        
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readFlow");
		exit(1);
	}
}

//**********************************************************************************************************************
int SortSeqsCommand::readQual(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(qualfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(qualfile));
        variables["[extension]"] = util.getExtension(qualfile);
		string outputFileName = getOutputFileName("qfile", variables);
        outputTypes["qfile"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(qualfile, in);
		string name;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    QualityScores currQual;
                    currQual = QualityScores(in); 
                    name = currQual.getName();
                    
                    if (name != "") {
                        numNamesInFile++;
                        map<string, int>::iterator it = names.find(name);
                        if (it == names.end()) { 
                            names[name] = numNames; numNames++;
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();
                out.close();
                
                int numLeft = names.size();
                if (numNamesInFile < numLeft) { numLeft = numNamesInFile; }
                
                int size = 1000; //assume that user can hold 1000 seqs in memory
                if (numLeft < size) { size = numLeft; }
                int times = 0;

                
                vector<QualityScores> seqs; seqs.resize(size);
                for (int i = 0; i < seqs.size(); i++) { seqs[i].setName(""); } //this is so if some of the seqs are missing we dont print out garbage
                
                while (numLeft > 0) {
                    
                    ifstream in2;
                    util.openInputFile(qualfile, in2);
                    
                    if (m->getControl_pressed()) { in2.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->getControl_pressed()) { in2.close();   util.mothurRemove(outputFileName);  return 0; }
                        
                        //stop reading if we already found the seqs we are looking for
                        if (found >= needToFind) { break; }
                        
                        QualityScores currQual;
                        currQual = QualityScores(in2); 
                        name = currQual.getName();
                        
                        if (name != "") {
                            map<string, int>::iterator it = names.find(name);
                            if (it != names.end()) { //we found it, so put it in the vector in the right place.
                                //is it in the set of seqs we are looking for this time around
                                int thisSeqsPlace = it->second;
                                thisSeqsPlace -= (times * size);
                                if ((thisSeqsPlace < size) && (thisSeqsPlace >= 0)) {
                                    seqs[thisSeqsPlace] = currQual; 
                                    found++;
                                }
                            }else { m->mothurOut("[ERROR]: in logic of readQual function.\n"); m->setControl_pressed(true); }
                        }
                        util.gobble(in2);
                    }
                    in2.close();	
                    
                    ofstream out2;
                    util.openOutputFileAppend(outputFileName, out2);
                    
                    int output = seqs.size();
                    if (numLeft < seqs.size()) { output = numLeft; }
                    
                    for (int i = 0; i < output; i++) {
                        if (seqs[i].getName() != "") {
                            seqs[i].printQScores(out2);
                        }
                    }
                    out2.close();
                    
                    times++;
                    numLeft -= output;
                }
                
                 m->mothurOut("Ordered " + toString(numNamesInFile) + " sequences from " + qualfile + ".\n");
                
            }else {
                
                vector<QualityScores> seqs; seqs.resize(names.size());
                for (int i = 0; i < seqs.size(); i++) { seqs[i].setName(""); } //this is so if some of the seqs are missing we dont print out garbage
                
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    QualityScores currQual;
                    currQual = QualityScores(in); 
                    name = currQual.getName();
                    
                    if (name != "") {
                        map<string, int>::iterator it = names.find(name);
                        if (it != names.end()) { //we found it, so put it in the vector in the right place.
                            seqs[it->second] = currQual;  
                        }else { //if we cant find it then add it to the end
                            names[name] = seqs.size();
                            seqs.push_back(currQual);
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    util.gobble(in);
                }
                in.close();	
                
                int count = 0;
                for (int i = 0; i < seqs.size(); i++) {
                    if (seqs[i].getName() != "") { seqs[i].printQScores(out); count++; }
                }
                out.close();
                
                m->mothurOut("Ordered " + toString(count) + " sequences from " + qualfile + ".\n");
            }
            
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                QualityScores currQual;
                currQual = QualityScores(in);  
                               
                util.gobble(in);
                
                if (currQual.getName() != "") {
                    //if this name is in the accnos file
                    names[currQual.getName()] = count;
                    count++;
                    currQual.printQScores(out);
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + qualfile + " to determine the order. It contains " + toString(count) + " sequences.\n");
        }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readQual");
		exit(1);
	}
}
//**********************************************************************************************************************
int SortSeqsCommand::readName(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[extension]"] = util.getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);
        outputTypes["name"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
        
                vector<string> seqs; seqs.resize(names.size(), "");
                
                while(!in.eof()){
                    if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                    
                    in >> firstCol;		util.gobble(in);		
                    in >> secondCol;    util.gobble(in);
                    
                    if (firstCol != "") {
                        map<string, int>::iterator it = names.find(firstCol);
                        if (it != names.end()) { //we found it, so put it in the vector in the right place.
                            seqs[it->second] = firstCol + '\t' + secondCol;  
                        }else { //if we cant find it then add it to the end
                            names[firstCol] = seqs.size();
                            seqs.push_back((firstCol + '\t' + secondCol));
                            m->mothurOut(firstCol + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                }
                in.close();	
                
                int count = 0;
                for (int i = 0; i < seqs.size(); i++) {
                    if (seqs[i] != "") { out << seqs[i] << endl; count++; }
                }
                out.close();
                
                m->mothurOut("Ordered " + toString(count) + " sequences from " + namefile + ".\n");
            
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> firstCol;		util.gobble(in);		
                in >> secondCol;    util.gobble(in);
                
                if (firstCol != "") {
                    //if this name is in the accnos file
                    names[firstCol] = count;
                    count++;
                    out << firstCol << '\t' << secondCol << endl;
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + namefile + " to determine the order. It contains " + toString(count) + " representative sequences.\n");
        }
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readName");
		exit(1);
	}
}
//**********************************************************************************************************************
int SortSeqsCommand::readCount(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
		string outputFileName = getOutputFileName("count", variables);
        outputTypes["count"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(countfile, in);
		string firstCol, rest;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            vector<string> seqs; seqs.resize(names.size(), "");
            
            string headers = util.getline(in); util.gobble(in);
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> firstCol;		util.gobble(in);		
                rest = util.getline(in);    util.gobble(in);
                
                if (firstCol != "") {
                    map<string, int>::iterator it = names.find(firstCol);
                    if (it != names.end()) { //we found it, so put it in the vector in the right place.
                        seqs[it->second] = firstCol + '\t' + rest;  
                    }else { //if we cant find it then add it to the end
                        names[firstCol] = seqs.size();
                        seqs.push_back((firstCol + '\t' + rest));
                        m->mothurOut(firstCol + " was not in the contained the file which determined the order, adding it to the end.\n");
                    }
                }
            }
            in.close();	
            
            int count = 0;
            out << headers << endl;
            for (int i = 0; i < seqs.size(); i++) {
                if (seqs[i] != "") { out << seqs[i] << endl; count++; }
            }
            out.close();
            
            m->mothurOut("Ordered " + toString(count) + " sequences from " + countfile + ".\n");
            
        }else { //read in file to fill names
            int count = 0;
            
            string headers = util.getline(in); util.gobble(in);
            out << headers << endl;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> firstCol;		util.gobble(in);		
                rest = util.getline(in);  util.gobble(in);
                
                if (firstCol != "") {
                    //if this name is in the accnos file
                    names[firstCol] = count;
                    count++;
                    out << firstCol << '\t' << rest << endl;
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + countfile + " to determine the order. It contains " + toString(count) + " representative sequences.\n");
        }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
int SortSeqsCommand::readGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(groupfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);
        outputTypes["group"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(groupfile, in);
		string name, group;
		
		if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            vector<string> seqs; seqs.resize(names.size(), "");
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> name;		util.gobble(in);		
                in >> group;    util.gobble(in);
                
                if (name != "") {
                    map<string, int>::iterator it = names.find(name);
                    if (it != names.end()) { //we found it, so put it in the vector in the right place.
                        seqs[it->second] = name + '\t' + group;  
                    }else { //if we cant find it then add it to the end
                        names[name] = seqs.size();
                        seqs.push_back((name + '\t' + group));
                        m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                    }
                }
            }
            in.close();	
            
            int count = 0;
            for (int i = 0; i < seqs.size(); i++) {
                if (seqs[i] != "") { out << seqs[i] << endl; count++; }
            }
            out.close();
            
            m->mothurOut("Ordered " + toString(count) + " sequences from " + groupfile + ".\n");
            
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> name;		util.gobble(in);		
                in >> group;    util.gobble(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << group << endl;
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + groupfile + " to determine the order. It contains " + toString(count) + " sequences.\n");
        }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int SortSeqsCommand::readTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);

        outputTypes["taxonomy"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(taxfile, in);
		string name, tax;
		
		if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            vector<string> seqs; seqs.resize(names.size(), "");
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> name; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                if (name != "") {
                    map<string, int>::iterator it = names.find(name);
                    if (it != names.end()) { //we found it, so put it in the vector in the right place.
                        seqs[it->second] = name + '\t' + tax;  
                    }else { //if we cant find it then add it to the end
                        names[name] = seqs.size();
                        seqs.push_back((name + '\t' + tax));
                        m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                    }
                }
            }
            in.close();	
            
            int count = 0;
            for (int i = 0; i < seqs.size(); i++) {
                if (seqs[i] != "") { out << seqs[i] << endl; count++; }
            }
            out.close();
            
            m->mothurOut("Ordered " + toString(count) + " sequences from " + taxfile + ".\n");
            
        }else { //read in file to fill names
            int count = 0;
            
            while(!in.eof()){
                if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
                
                in >> name; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << tax << endl;
                }
                util.gobble(in);
            }
            in.close();	
            out.close();
            
            m->mothurOut("\nUsing " + taxfile + " to determine the order. It contains " + toString(count) + " sequences.\n");
        }
        
		return 0;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************





