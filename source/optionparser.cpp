/*
 *  optionparser.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "optionparser.h"

/***********************************************************************/

OptionParser::OptionParser(string option, vector<string> parametersAllowedByThisCommand) {
	try {
		m = MothurOut::getInstance();
        current = CurrentFile::getInstance();
        
        ValidParameters validParameter;
        
		if (option != "") {
			
			string key, value;		
			//reads in parameters and values
			while((option.find_first_of(',') != -1)) {  //while there are parameters
				util.splitAtComma(value, option);
				util.splitAtEquals(key, value);
				if ((key == "candidate") || (key == "query")) { key = "fasta"; }
				if (key == "template") { key = "reference"; }
				key = util.splitWhiteSpace(key).front();
                //if value is wrapped in '' preserve spaces
                if ((value[0] == '\'') && (value[(value.length()-1)] == '\'')) {  value = value.substr(1); value = value.substr(0, (value.length()-1)); }
                else {
                    value = util.splitWhiteSpace(value).front();
                }
				
                if (!validParameter.isValidParameter(key, parametersAllowedByThisCommand, value)) {} //ignore invalid parameters
                else { parameters[key] = value; }
			}
			
			//in case there is no comma and to get last parameter after comma
			util.splitAtEquals(key, option);
			if ((key == "candidate") || (key == "query")) { key = "fasta"; }
			if (key == "template") { key = "reference"; }
            key = util.splitWhiteSpace(key).front();
            //if value is wrapped in '' preserve spaces
            if ((option[0] == '\'') && (option[(option.length()-1)] == '\'')) {  option = option.substr(1); option = option.substr(0, (option.length()-1)); }
            else {
                option = util.splitWhiteSpace(option).front();
            }
            
			if (!validParameter.isValidParameter(key, parametersAllowedByThisCommand, option)) {} //ignore invalid parameters
            else { parameters[key] = option; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "OptionParser", "OptionParser");
		exit(1);
	}
}
/***********************************************************************/

OptionParser::OptionParser(string option, map<string, string>& copy) {
	try {
		m = MothurOut::getInstance();
        current = CurrentFile::getInstance();
		if (option != "") {
			
			string key, value;		
			//reads in parameters and values
			while((option.find_first_of(',') != -1)) {  //while there are parameters
				util.splitAtComma(value, option);
				util.splitAtEquals(key, value);
				if ((key == "candidate") || (key == "query")) { key = "fasta"; }
				if (key == "template") { key = "reference"; }
				key = util.splitWhiteSpace(key).front();
                //if value is wrapped in '' preserve spaces
                if ((value[0] == '\'') && (value[(value.length()-1)] == '\'')) {  value = value.substr(1); value = value.substr(0, (value.length()-1)); }
                else {
                    value = util.splitWhiteSpace(value).front();
                }
				parameters[key] = value;
			}
			
			//in case there is no comma and to get last parameter after comma
			util.splitAtEquals(key, option);
			if ((key == "candidate") || (key == "query")) { key = "fasta"; }
			if (key == "template") { key = "reference"; }
			key = util.splitWhiteSpace(key).front();
            //if value is wrapped in '' preserve spaces
            if ((option[0] == '\'') && (option[(option.length()-1)] == '\'')) {  option = option.substr(1); option = option.substr(0, (option.length()-1)); }
            else {
                option = util.splitWhiteSpace(option).front();
            }
			parameters[key] = option;
		}
        
        copy = parameters;
	}
	catch(exception& e) {
		m->errorOut(e, "OptionParser", "OptionParser");
		exit(1);
	}
}
/***********************************************************************/

map<string, string> OptionParser::getParameters() {	
    try {
        //loop through parameters and look for "current" so you can return the appropriate file
        //doing it here to avoid code duplication in each of the commands
        
        map<string, string>::iterator it;
        
        for (it = parameters.begin(); it != parameters.end();) {
            
            if (it->second == "current") {
                
                //look for file types
                if (it->first == "fasta") {
                    it->second = current->getFastaFile();
                }else if (it->first == "qfile") {
                    it->second = current->getQualFile();
                }else if (it->first == "phylip") {
                    it->second = current->getPhylipFile();
                }else if (it->first == "column") {
                    it->second = current->getColumnFile();
                }else if (it->first == "list") {
                    it->second = current->getListFile();
                }else if (it->first == "rabund") {
                    it->second = current->getRabundFile();
                }else if (it->first == "clr") {
                    it->second = current->getCLRFile();
                }else if (it->first == "sabund") {
                    it->second = current->getSabundFile();
                }else if (it->first == "name") {
                    it->second = current->getNameFile();
                }else if (it->first == "group") {
                    it->second = current->getGroupFile();
                }else if (it->first == "order") {
                    it->second = current->getOrderFile();
                }else if (it->first == "ordergroup") {
                    it->second = current->getOrderGroupFile();
                }else if (it->first == "tree") {
                    it->second = current->getTreeFile();
                }else if (it->first == "shared") {
                    it->second = current->getSharedFile();
                }else if (it->first == "relabund") {
                    it->second = current->getRelAbundFile();
                }else if (it->first == "design") {
                    it->second = current->getDesignFile();
                }else if (it->first == "sff") {
                    it->second = current->getSFFFile();
                }else if (it->first == "flow") {
                    it->second = current->getFlowFile();
                }else if (it->first == "oligos") {
                    it->second = current->getOligosFile();
                }else if (it->first == "accnos") {
                    it->second = current->getAccnosFile();
                }else if (it->first == "taxonomy") {
                    it->second = current->getTaxonomyFile();
                }else if (it->first == "constaxonomy") {
                    it->second = current->getConsTaxonomyFile();
                }else if (it->first == "contigsreport") {
                    it->second = current->getContigsReportFile();
                }else if (it->first == "biom") {
                    it->second = current->getBiomFile();
                }else if (it->first == "count") {
                    it->second = current->getCountFile();
                }else if (it->first == "summary") {
                    it->second = current->getSummaryFile();
                }else if (it->first == "file") {
                    it->second = current->getFileFile();
                }else if (it->first == "sample") {
                    it->second = current->getSampleFile();
                }else {
                    m->mothurOut("[ERROR]: mothur does not save a current file for " + it->first); m->mothurOutEndLine();
                }
                
                if (it->second == "") { //no file was saved for that type, warn and remove from parameters
                    m->mothurOut("[WARNING]: no file was saved for " + it->first + " parameter.\n"); 
                    parameters.erase(it++);
                }else {
                    m->mothurOut("Using " + it->second + " as input file for the " + it->first + " parameter.\n");
                    it++;
                }
            }else{ it++; }
        }
        
        string inputDir = "";
        it = parameters.find("inputdir");
        if (it != parameters.end()){
            inputDir = it->second;
            
            set<string> fileTypes;
            fillFileTypes(fileTypes);
            
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (fileTypes.count(it->first) != 0) {
                    string path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {
                        it->second = inputDir + it->second;
                    }
                }
            }
        }
        return parameters;
    }
    catch(exception& e) {
        m->errorOut(e, "OptionParser", "getParameters");
        exit(1);
    }
}
/***********************************************************************/

void OptionParser::fillFileTypes(set<string>& fileTypes) {
    try {
        fileTypes.insert("fasta");
        fileTypes.insert("qfile");
        fileTypes.insert("phylip");
        fileTypes.insert("column");
        fileTypes.insert("list");
        fileTypes.insert("rabund");
        fileTypes.insert("clr");
        fileTypes.insert("sabund");
        
        fileTypes.insert("name");
        fileTypes.insert("group");
        fileTypes.insert("order");
        fileTypes.insert("ordergroup");
        fileTypes.insert("tree");
        fileTypes.insert("shared");
        fileTypes.insert("relabund");
        fileTypes.insert("design");
        
        fileTypes.insert("sff");
        fileTypes.insert("sfftxt");
        fileTypes.insert("flow");
        fileTypes.insert("oligos");
        fileTypes.insert("accnos");
        fileTypes.insert("taxonomy");
        fileTypes.insert("constaxonomy");
        fileTypes.insert("contigsreport");
        fileTypes.insert("biom");
        
        fileTypes.insert("count");
        fileTypes.insert("summary");
        fileTypes.insert("file");
        fileTypes.insert("sample");
        
        fileTypes.insert("list");
        fileTypes.insert("rabund");
        fileTypes.insert("clr");
        fileTypes.insert("sabund");
        fileTypes.insert("reference");
        fileTypes.insert("conservation");
        fileTypes.insert("quantile");
        fileTypes.insert("reffasta");
        fileTypes.insert("refcolumn");
        fileTypes.insert("reflist");
        fileTypes.insert("refname");
        fileTypes.insert("refcount");
        fileTypes.insert("reftaxonomy");
        
        fileTypes.insert("axes");
        fileTypes.insert("metadata");
        fileTypes.insert("refname");
        fileTypes.insert("repfasta");
        fileTypes.insert("oldfasta");
        fileTypes.insert("hard");
        fileTypes.insert("alignreport");
        fileTypes.insert("report");
        
        fileTypes.insert("corraxes");
        fileTypes.insert("otucorr");
        fileTypes.insert("accnos");
        fileTypes.insert("phylip1");
        fileTypes.insert("phylip2");
        fileTypes.insert("picrust");
        
        fileTypes.insert("ffastq");
        fileTypes.insert("rfastq");
        fileTypes.insert("ffasta");
        fileTypes.insert("rfasta");
        fileTypes.insert("fqfile");
        fileTypes.insert("rqfile");
        fileTypes.insert("findex");
        fileTypes.insert("rindex");
        fileTypes.insert("error");
        fileTypes.insert("xml");
        fileTypes.insert("ecoli");
        fileTypes.insert("map");
        fileTypes.insert("lookup");
        fileTypes.insert("project");
        fileTypes.insert("mimark");
        
        fileTypes.insert("vsearch");
        fileTypes.insert("blast");
        fileTypes.insert("uchime");
        fileTypes.insert("prefetch");
        fileTypes.insert("fasterq-dump");
        fileTypes.insert("input");
        
        //uchime, vsearch, prefetch  fasterq-dump, and blast are not included
       
    }
    catch(exception& e) {
        m->errorOut(e, "OptionParser", "fillFileTypes");
        exit(1);
    }
}
/***********************************************************************/
//pass a vector of filenames that may match the current namefile.  
//this function will look at each one, if the rootnames match, mothur will warn 
//the user that they may have neglected to provide a namefile.
//stops when it finds a match.
bool OptionParser::getNameFile(vector<string> files) {
    try {
        string namefile = current->getNameFile();
        bool match = false;
        
        if (namefile != "") {
            string temp = util.getRootName(util.getSimpleName(namefile));
            vector<string> rootName;
            util.splitAtChar(temp, rootName, '.');
            
            for (int i = 0; i < files.size(); i++) {
                temp = util.getRootName(util.getSimpleName(files[i]));
                vector<string> root;
                util.splitAtChar(temp, root, '.');
                
                int smallest = rootName.size();
                if (root.size() < smallest) { smallest = root.size(); }
                
                int numMatches = 0;
                for(int j = 0; j < smallest; j++) {
                    if (root[j] == rootName[j]) { numMatches++; }
                }
                
                if (smallest > 0) {
                    if ((numMatches >= (smallest-2)) && (root[0] == rootName[0])) {
                        m->mothurOut("[WARNING]: This command can take a namefile and you did not provide one. The current namefile is " + namefile + " which seems to match " + files[i] + ".");
                        m->mothurOutEndLine();
                        match = true;
                        break;
                    }
                }
            }
        }
        
        return match;
    }
    catch(exception& e) {
        m->errorOut(e, "OptionParser", "getNameFile");
        exit(1);
    }
}


/***********************************************************************/
