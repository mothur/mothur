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
				}else if (it->first == "oligos") {
					it->second = m->getOligosFile();
				}else if (it->first == "accnos") {
					it->second = m->getAccnosFile();
				}else if (it->first == "taxonomy") {
					it->second = m->getTaxonomyFile();
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
