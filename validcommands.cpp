/*
 *  validcommands.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validcommands.h"

/***********************************************************************/

ValidCommands::ValidCommands() {
	try {
	
		commands["read.dist"]			= "read.dist"; 
		commands["read.otu"]			= "read.otu";
		commands["read.tree"]			= "read.tree"; 
		commands["bin.seqs"]			= "bin.seqs"; 
		commands["get.oturep"]			= "get.oturep";
		commands["cluster"]				= "cluster"; 
		commands["deconvolute"]			= "deconvolute"; 
		commands["dist.seqs"]			= "dist.seqs";
		commands["dist.shared"]			= "dist.shared";
		commands["collect.single"]		= "collect.single"; 
		commands["collect.shared"]		= "collect.shared"; 
		commands["rarefaction.single"]	= "rarefaction.single"; 
		commands["rarefaction.shared"]	= "rarefaction.shared"; 
		commands["summary.single"]		= "summary.single"; 
		commands["summary.shared"]		= "summary.shared"; 
		commands["parsimony"]			= "parsimony";
		commands["unifrac.weighted"]	= "unifrac.weighted"; 
		commands["unifrac.unweighted"]	= "unifrac.unweighted"; 
		commands["libshuff"]			= "libshuff";
		commands["tree.shared"]			= "tree.shared";
		commands["heatmap"]				= "heatmap";
		commands["venn"]				= "venn";
		commands["get.group"]           = "get.group";
		commands["get.label"]           = "get.label";
		commands["get.line"]            = "get.line";
		commands["get.sabund"]          = "get.sabund";
		commands["get.rabund"]          = "get.rabund";
		commands["bootstrap.shared"]	= "bootstrap.shared";
		commands["concensus"]			= "concensus";
		commands["help"]				= "help"; 
		commands["filter.seqs"]			= "filter.seqs";
		commands["align.seqs"]			= "align.seqs";
		commands["summary.seqs"]		= "summary.seqs";
		commands["screen.seqs"]			= "screen.seqs";
		commands["quit"]				= "quit"; 

				
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidCommands class Function ValidCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidCommands class function ValidCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ValidCommands::~ValidCommands() {}

/***********************************************************************/
bool ValidCommands::isValidCommand(string command) {
	try {	
	
		//is the command in the map
		if ((commands.find(command)) != (commands.end())) {
			return true;
		}else{
			cout << command << " is not a valid command in Mothur.  Valid commands are ";
			for (it = commands.begin(); it != commands.end(); it++) {
				cout << it->first << ", ";
			}
			cout << endl;
			return false;
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidCommands class Function isValidCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidCommands class function isValidCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/
void ValidCommands::printCommands(ostream& out) {
	try {	
		out << "Valid commands are ";
		for (it = commands.begin(); it != commands.end(); it++) {
			out << it->first << ", ";
		}
		out << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidCommands class Function printCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidCommands class function printCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}




