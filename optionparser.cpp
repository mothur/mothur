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

OptionParser::OptionParser(string option) {
	try {
		m = MothurOut::getInstance();
		if (option != "") {
			
			string key, value;		
			//reads in parameters and values
			while((option.find_first_of(',') != -1)) {  //while there are parameters
				m->splitAtComma(value, option);
				m->splitAtEquals(key, value);
				if ((key == "candidate") || (key == "query")) { key = "fasta"; }
				if (key == "template") { key = "reference"; }
				parameters[key] = value;
			}
			
			//in case there is no comma and to get last parameter after comma
			m->splitAtEquals(key, option);
			if ((key == "candidate") || (key == "query")) { key = "fasta"; }
			if (key == "template") { key = "reference"; }
			parameters[key] = option;
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
		if (option != "") {
			
			string key, value;		
			//reads in parameters and values
			while((option.find_first_of(',') != -1)) {  //while there are parameters
				m->splitAtComma(value, option);
				m->splitAtEquals(key, value);
				if ((key == "candidate") || (key == "query")) { key = "fasta"; }
				if (key == "template") { key = "reference"; }
				parameters[key] = value;
			}
			
			//in case there is no comma and to get last parameter after comma
			m->splitAtEquals(key, option);
			if ((key == "candidate") || (key == "query")) { key = "fasta"; }
			if (key == "template") { key = "reference"; }
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
                        it->second = m->getFastaFile();
                    }else if (it->first == "qfile") {
                        it->second = m->getQualFile();
                    }else if (it->first == "phylip") {
                        it->second = m->getPhylipFile();
                    }else if (it->first == "column") {
                        it->second = m->getColumnFile();
                    }else if (it->first == "list") {
                        it->second = m->getListFile();
                    }else if (it->first == "rabund") {
                        it->second = m->getRabundFile();
                    }else if (it->first == "sabund") {
                        it->second = m->getSabundFile();
                    }else if (it->first == "name") {
                        it->second = m->getNameFile();
                    }else if (it->first == "group") {
                        it->second = m->getGroupFile();
                    }else if (it->first == "order") {
                        it->second = m->getOrderFile();
                    }else if (it->first == "ordergroup") {
                        it->second = m->getOrderGroupFile();
                    }else if (it->first == "tree") {
                        it->second = m->getTreeFile();
                    }else if (it->first == "shared") {
                        it->second = m->getSharedFile();
                    }else if (it->first == "relabund") {
                        it->second = m->getRelAbundFile();
                    }else if (it->first == "design") {
                        it->second = m->getDesignFile();
                    }else if (it->first == "sff") {
                        it->second = m->getSFFFile();
                    }else if (it->first == "flow") {
                            it->second = m->getFlowFile();
                    }else if (it->first == "oligos") {
                        it->second = m->getOligosFile();
                    }else if (it->first == "accnos") {
                        it->second = m->getAccnosFile();
                    }else if (it->first == "taxonomy") {
                        it->second = m->getTaxonomyFile();
                    }else if (it->first == "biom") {
                        it->second = m->getBiomFile();
                    }else if (it->first == "count") {
                            it->second = m->getCountTableFile();
                    }else if (it->first == "summary") {
                            it->second = m->getSummaryFile();
                    }else {
                        m->mothurOut("[ERROR]: mothur does not save a current file for " + it->first); m->mothurOutEndLine();
                    }
                    
                    if (it->second == "") { //no file was saved for that type, warn and remove from parameters
                        m->mothurOut("[WARNING]: no file was saved for " + it->first + " parameter."); m->mothurOutEndLine();
                        parameters.erase(it++);
                    }else {
                        m->mothurOut("Using " + it->second + " as input file for the " + it->first + " parameter."); m->mothurOutEndLine();
                        it++;
                    }
                }else{ it++; }
            }
		
		return parameters;	
	}
	catch(exception& e) {
		m->errorOut(e, "OptionParser", "getParameters");
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
		string namefile = m->getNameFile();
		bool match = false;
		
		if ((namefile != "")&&(!m->mothurCalling)) {
			string temp = m->getRootName(m->getSimpleName(namefile));
			vector<string> rootName;
			m->splitAtChar(temp, rootName, '.');
			
			for (int i = 0; i < files.size(); i++) {
				temp = m->getRootName(m->getSimpleName(files[i]));
				vector<string> root;
				m->splitAtChar(temp, root, '.');
				
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
