//
//  removedistscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/29/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "removedistscommand.h"

//**********************************************************************************************************************
vector<string> RemoveDistsCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "PhylipColumn", "none","phylip",false,false,true); parameters.push_back(pphylip);
        CommandParameter pcolumn("column", "InputTypes", "", "", "none", "PhylipColumn", "none","column",false,false,true); parameters.push_back(pcolumn);	
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(paccnos);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
		
        vector<string> tempOutNames;
        outputTypes["column"] = tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveDistsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.dists command removes distances from a phylip or column file related to groups or sequences listed in an accnos file.\n";
		helpString += "The remove.dists command parameters are accnos, phylip and column.\n";
		helpString += "The remove.dists command should be in the following format: get.dists(accnos=yourAccnos, phylip=yourPhylip).\n";
		helpString += "Example remove.dists(accnos=final.accnos, phylip=final.an.thetayc.0.03.lt.ave.dist).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveDistsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "phylip")           {   pattern = "[filename],pick,[extension]";    }
        else if (type == "column")      {   pattern = "[filename],pick,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveDistsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RemoveDistsCommand::RemoveDistsCommand(string option)  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = current->getAccnosFile(); 
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You have no valid accnos file and accnos is required."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else { current->setAccnosFile(accnosfile); }	
			
			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else { 	current->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  current->setColumnFile(columnfile);	}
			
			if ((phylipfile == "") && (columnfile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = current->getColumnFile(); 
				if (columnfile != "") {  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = current->getPhylipFile(); 
					if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "RemoveDistsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveDistsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//get names you want to keep
		names = util.readAccnos(accnosfile);
		
		if (m->getControl_pressed()) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (phylipfile != "")		{		readPhylip();		}
		if (columnfile != "")		{		readColumn();       }
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
		
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string currentName = "";
			itTypes = outputTypes.find("phylip");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setPhylipFile(currentName); }
			}
			
			itTypes = outputTypes.find("column");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); }
			}
        }
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveDistsCommand::readPhylip(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(phylipfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(phylipfile));
        variables["[extension]"] = util.getExtension(phylipfile);
		string outputFileName = getOutputFileName("phylip", variables);
		
        ifstream in;
        util.openInputFile(phylipfile, in);
        
        float distance;
        int square, nseqs; 
        string name;
        unsigned int row;
        set<unsigned int> rows; //converts names in names to a index
        row = 0;
        
        string numTest;
        in >> numTest >> name;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        //not one we want to remove
        if (names.count(name) == 0) { rows.insert(row); }
        row++;
        
        //is the matrix square?
        char d;
        while((d=in.get()) != EOF){
            
            if(isalnum(d)){
                square = 1;
                in.putback(d);
                for(int i=0;i<nseqs;i++){
                    in >> distance;
                }
                break;
            }
            if(d == '\n'){
                square = 0;
                break;
            }
        }
        
        //map name to row/column        
        if(square == 0){
            for(int i=1;i<nseqs;i++){
                in >> name;  
                if (names.count(name) == 0) { rows.insert(row); }
                row++;
                
                for(int j=0;j<i;j++){
                    if (m->getControl_pressed()) {  in.close(); return 0;  }
                    in >> distance;
                }
            }
        }
        else{
            for(int i=1;i<nseqs;i++){
                in >> name;  
                if (names.count(name) == 0) { rows.insert(row);  }
                row++;
                for(int j=0;j<nseqs;j++){
                    if (m->getControl_pressed()) {  in.close(); return 0;  }
                    in >> distance;
                }
            }
        }
        in.close();
        
        if (m->getControl_pressed()) {  return 0; }
        
        //read through file only printing rows and columns of seqs in names
        ifstream inPhylip;
        util.openInputFile(phylipfile, inPhylip);
        
        inPhylip >> numTest;
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        outputTypes["phylip"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        out << (nseqs-names.size()) << endl;
            
        unsigned int count = 0;
        unsigned int keptCount = 0;
		if(square == 0){
            for(int i=0;i<nseqs;i++){
                inPhylip >> name;  
                bool ignoreRow = false;
                
                if (names.count(name) != 0) { ignoreRow = true; count++; }
                else{ out << name; keptCount++; }
                
                for(int j=0;j<i;j++){
                    if (m->getControl_pressed()) {  inPhylip.close(); out.close();  return 0;  }
                    inPhylip >> distance;
                    if (!ignoreRow) {
                        //is this a column we want
                        if(rows.count(j) != 0) {  out << '\t' << distance;  }
                    }
                }
                if (!ignoreRow) { out << endl; }
            }
        }
        else{
            for(int i=0;i<nseqs;i++){
                inPhylip >> name; 
                
                bool ignoreRow = false;
                
                if (names.count(name) != 0) { ignoreRow = true; count++; }
                else{ out << name; keptCount++; }
                
                for(int j=0;j<nseqs;j++){
                    if (m->getControl_pressed()) {  inPhylip.close(); out.close(); return 0;  }
                    inPhylip >> distance;
                    if (!ignoreRow) {
                        //is this a column we want
                        if(rows.count(j) != 0) {  out << '\t' << distance;  }
                    }
                }
                if (!ignoreRow) { out << endl; }
            }
        }
        inPhylip.close();
		out.close();
		
		if (keptCount == 0) {  m->mothurOut("Your file contains ONLY distances related to groups or sequences listed in the accnos file."); m->mothurOutEndLine();  }
        else if (count != names.size()) {
            m->mothurOut("[WARNING]: Your accnos file contains " + toString(names.size()) + " groups or sequences, but I only found " + toString(count) + " of them in the phylip file."); m->mothurOutEndLine();
            //rewrite with new number
            util.renameFile(outputFileName, outputFileName+".temp");
            ofstream out2;
            util.openOutputFile(outputFileName, out2);
            out2 << keptCount << endl;
            
            ifstream in3;
            util.openInputFile(outputFileName+".temp", in3);
            in3 >> nseqs; util.gobble(in3);
            char buffer[4096];        
            while (!in3.eof()) {
                in3.read(buffer, 4096);
                out2.write(buffer, in3.gcount());
            }
            in3.close();
            out2.close();
            util.mothurRemove(outputFileName+".temp");
        }
		
		m->mothurOut("Removed " + toString(count) + " groups or sequences from your phylip file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "readPhylip");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveDistsCommand::readColumn(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(columnfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(columnfile));
        variables["[extension]"] = util.getExtension(columnfile);
		string outputFileName = getOutputFileName("column", variables);
        outputTypes["column"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
        ifstream in;
        util.openInputFile(columnfile, in);
        
        set<string> removeNames;
        string firstName, secondName;
        float distance;
        bool wrote = false;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { out.close(); in.close(); return 0; }
            
            in >> firstName >> secondName >> distance; util.gobble(in);
            
            //is either names in the accnos file
            if (names.count(firstName) != 0)       { 
                removeNames.insert(firstName);  
                if (names.count(secondName) != 0)  { removeNames.insert(secondName);      }   }
            else if (names.count(secondName) != 0) { 
                removeNames.insert(secondName); 
                if (names.count(firstName) != 0)   { removeNames.insert(firstName);     }   }
            else {
                wrote = true;
                out << firstName << '\t' << secondName << '\t' << distance << endl;
            }
        }
		in.close();
		out.close();
        
        if (!wrote) {  m->mothurOut("Your file contains ONLY distances related to groups or sequences listed in the accnos file."); m->mothurOutEndLine();  }
        else if (removeNames.size() != names.size()) {
            m->mothurOut("[WARNING]: Your accnos file contains " + toString(names.size()) + " groups or sequences, but I only found " + toString(removeNames.size()) + " of them in the column file."); m->mothurOutEndLine();
        }
		
		m->mothurOut("Removed " + toString(removeNames.size()) + " groups or sequences from your column file."); m->mothurOutEndLine();
        
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveDistsCommand", "readColumn");
		exit(1);
	}
}
//**********************************************************************************************************************


