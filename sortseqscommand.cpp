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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pfasta);
        CommandParameter pflow("flow", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pflow);
		CommandParameter pname("name", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pgroup);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(ptaxonomy);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pqfile);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(plarge);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(paccnos);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The sort.seqs command puts the sequences in the same order for the following file types: accnos fasta, name, group, taxonomy, flow or quality file.\n";
        helpString += "The sort.seqs command parameters are accnos, fasta, name, group, taxonomy, flow, qfile and large.\n";
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
string SortSeqsCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "fasta")            {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else if (type == "taxonomy")    {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else if (type == "name")        {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else if (type == "group")       {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else if (type == "flow")        {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else if (type == "qfile")       {   outputFileName =  "sorted" + m->getExtension(inputName);   }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "SortSeqsCommand", "getOutputFileNameTag");
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
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
                
                it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
			}
            
			
			//check for parameters
            accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { accnosfile = ""; abort = true; }
			else if (accnosfile == "not found") {  accnosfile = "";  }	
			else { m->setAccnosFile(accnosfile); }
            
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			else { m->setFastaFile(fastafile); }
            
            flowfile = validParameter.validFile(parameters, "flow", true);
			if (flowfile == "not open") { flowfile = ""; abort = true; }
			else if (flowfile == "not found") {  flowfile = "";  }	
			else { m->setFlowFile(flowfile); }
            
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { m->setNameFile(namefile); } 
            
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }
			else { m->setGroupFile(groupfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { m->setTaxonomyFile(taxfile); }
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") {  qualfile = "";  }			
			else { m->setQualFile(qualfile); }
			
            string temp = validParameter.validFile(parameters, "large", false);		if (temp == "not found") { temp = "f"; }
			large = m->isTrue(temp);
            
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (taxfile == "") && (flowfile == "") && (qualfile == ""))  { m->mothurOut("You must provide at least one of the following: fasta, name, group, taxonomy, flow or quality."); m->mothurOutEndLine(); abort = true; }
			
			if ((fastafile != "") && (namefile == "")) {
				vector<string> files; files.push_back(fastafile);
				parser.getNameFile(files);
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
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read through the correct file and output lines you want to keep
        if (accnosfile != "")		{		
            vector<string> temp;
            m->readAccnos(accnosfile, temp);
            for (int i = 0; i < temp.size(); i++) {  names[temp[i]] = i;  }
            m->mothurOut("\nUsing " + accnosfile + " to determine the order. It contains " + toString(temp.size()) + " representative sequences.\n");	
        }
        
		if (fastafile != "")		{		readFasta();	}
        if (flowfile != "")         {		readFlow();     }
        if (qualfile != "")			{		readQual();		}
        if (namefile != "")			{		readName();		}
		if (groupfile != "")		{		readGroup();	}
        if (taxfile != "")			{		readTax();		}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string current = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
			}
			
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
			}
			
			itTypes = outputTypes.find("qfile");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
			}	
            
            itTypes = outputTypes.find("flow");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFlowFile(current); }
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + getOutputFileNameTag("fasta", fastafile);
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
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
                    m->gobble(in);
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
                    m->openInputFile(fastafile, in2);
                    
                    if (m->control_pressed) { in2.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->control_pressed) { in2.close();   m->mothurRemove(outputFileName);  return 0; }
                        
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
                            }else { m->mothurOut("[ERROR]: in logic of readFasta function.\n"); m->control_pressed = true; }
                        }
                        m->gobble(in2);
                    }
                    in2.close();	

                    ofstream out2;
                    m->openOutputFileAppend(outputFileName, out2);
                    
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
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
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
                    m->gobble(in);
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                Sequence currSeq(in);
                name = currSeq.getName();
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    currSeq.printSequence(out);
                }
                m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(flowfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(flowfile)) + getOutputFileNameTag("flow", flowfile);
		outputTypes["flow"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(flowfile, in);
        int numFlows;
		string name;
        
        in >> numFlows; m->gobble(in);
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    in >> name;	
                    string rest = m->getline(in);
                    
                    if (name != "") {
                        numNamesInFile++;
                        map<string, int>::iterator it = names.find(name);
                        if (it == names.end()) { 
                            names[name] = numNames; numNames++;
                            m->mothurOut(name + " was not in the contained the file which determined the order, adding it to the end.\n");
                        }
                    }
                    m->gobble(in);
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
                    m->openInputFile(flowfile, in2); in2 >> numFlows; m->gobble(in2);
                    
                    if (m->control_pressed) { in2.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->control_pressed) { in2.close();   m->mothurRemove(outputFileName);  return 0; }
                        
                        //stop reading if we already found the seqs we are looking for
                        if (found >= needToFind) { break; }
                        
                        in2 >> name;	
                        string rest = m->getline(in2);
                        
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
                            }else { m->mothurOut("[ERROR]: in logic of readFlow function.\n"); m->control_pressed = true; }
                        }
                        m->gobble(in2);
                    }
                    in2.close();	
                    
                    ofstream out2;
                    m->openOutputFileAppend(outputFileName, out2);
                    
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
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    in >> name;	
                    string rest = m->getline(in);
                    
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
                    m->gobble(in);
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> name;	
                string rest = m->getline(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << rest << endl;
                }
                m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(qualfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(qualfile)) + getOutputFileNameTag("qfile", qualfile);
        outputTypes["qfile"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(qualfile, in);
		string name;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            if (large) { //if the file is too large to fit in memory we can still process it, but the io will be very time consuming.
                //read through the file looking for 1000 seqs at a time. Once we find them output them and start looking for the next 1000.
                //this way we only store 1000 seqs in memory at a time.
                
                int numNames = names.size();
                int numNamesInFile = 0;
                
                //to make sure we dont miss any seqs, add any seqs that are not in names but in the file to the end of names
                while(!in.eof()){
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
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
                    m->gobble(in);
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
                    m->openInputFile(qualfile, in2);
                    
                    if (m->control_pressed) { in2.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    int found = 0;
                    int needToFind = size;
                    if (numLeft < size) { needToFind = numLeft; }
                    
                    while(!in2.eof()){
                        if (m->control_pressed) { in2.close();   m->mothurRemove(outputFileName);  return 0; }
                        
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
                            }else { m->mothurOut("[ERROR]: in logic of readQual function.\n"); m->control_pressed = true; }
                        }
                        m->gobble(in2);
                    }
                    in2.close();	
                    
                    ofstream out2;
                    m->openOutputFileAppend(outputFileName, out2);
                    
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
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
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
                    m->gobble(in);
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                QualityScores currQual;
                currQual = QualityScores(in);  
                               
                m->gobble(in);
                
                if (currQual.getName() != "") {
                    //if this name is in the accnos file
                    names[currQual.getName()] = count;
                    count++;
                    currQual.printQScores(out);
                }
                m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(namefile)) + getOutputFileNameTag("name", namefile); 
        outputTypes["name"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
        if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
        
                vector<string> seqs; seqs.resize(names.size(), "");
                
                while(!in.eof()){
                    if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                    
                    in >> firstCol;		m->gobble(in);		
                    in >> secondCol;    m->gobble(in);
                    
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> firstCol;		m->gobble(in);		
                in >> secondCol;    m->gobble(in);
                
                if (firstCol != "") {
                    //if this name is in the accnos file
                    names[firstCol] = count;
                    count++;
                    out << firstCol << '\t' << secondCol << endl;
                }
                m->gobble(in);
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
int SortSeqsCommand::readGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + getOutputFileNameTag("group", groupfile);	
        outputTypes["group"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            vector<string> seqs; seqs.resize(names.size(), "");
            
            while(!in.eof()){
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> name;		m->gobble(in);		
                in >> group;    m->gobble(in);
                
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> name;		m->gobble(in);		
                in >> group;    m->gobble(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << group << endl;
                }
                m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + getOutputFileNameTag("taxonomy", taxfile); 
        outputTypes["taxonomy"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		if (names.size() != 0) {//this is not the first file we are reading so we need to use the order we already have
            
            vector<string> seqs; seqs.resize(names.size(), "");
            
            while(!in.eof()){
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> name;		m->gobble(in);		
                in >> tax;    m->gobble(in);
                
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
                if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
                
                in >> name;		m->gobble(in);		
                in >> tax;    m->gobble(in);
                
                if (name != "") {
                    //if this name is in the accnos file
                    names[name] = count;
                    count++;
                    out << name << '\t' << tax << endl;
                }
                m->gobble(in);
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





