/*
 *  validparameter.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validparameter.h"

/***********************************************************************/

ValidParameters::ValidParameters() {
	try {
		initialReaddist();
		initialReadotu();
		initialReadtree();
		initialCluster();
		initialDeconvolute();
		initialParsimony();
		initialCollectsingle();
		initialCollectshared();
		initialRarefactsingle();
		initialRarefactshared();
		initialSummarysingle();
		initialSummaryshared();
		initialUnifracweighted();
		initialUnifracunweighted();
		initialLibshuff();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ValidParameters::~ValidParameters() {}

/***********************************************************************/
bool ValidParameters::isValidParameter(string parameter, string command) {
	try {	
	
		if (command == "read.dist") {
			//is it valid
			if ((readdist.find(parameter)) != (readdist.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = readdist.begin(); it != readdist.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "read.otu") {
			//is it valid
			if ((readotu.find(parameter)) != (readotu.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = readotu.begin(); it != readotu.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		//are you looking for a calculator for a rarefaction parameter
		}else if (command == "read.tree") {
			//is it valid
			if ((readtree.find(parameter)) != (readtree.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = readtree.begin(); it != readtree.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		//are you looking for a calculator for a summary parameter
		}else if (command == "cluster") {
			//is it valid
			if ((cluster.find(parameter)) != (cluster.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = cluster.begin(); it != cluster.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }		//are you looking for a calculator for a sharedsummary parameter
		}else if (command == "deconvolute") {
			//is it valid
			if ((deconvolute.find(parameter)) != (deconvolute.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = deconvolute.begin(); it != deconvolute.end(); it++) {
					cout << it->first;
				}
				cout << endl;
				return false; }
		}else if (command == "parsimony") {
			//is it valid
			if ((parsimony.find(parameter)) != (parsimony.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = parsimony.begin(); it != parsimony.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "collect.single") {
			//is it valid
			if ((collectsingle.find(parameter)) != (collectsingle.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = collectsingle.begin(); it != collectsingle.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "collect.shared") {
			//is it valid
			if ((collectshared.find(parameter)) != (collectshared.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = collectshared.begin(); it != collectshared.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "rarefaction.single") {
			//is it valid
			if ((rarefactsingle.find(parameter)) != (rarefactsingle.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = rarefactsingle.begin(); it != rarefactsingle.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "rarefaction.shared") {
			//is it valid
			if ((rarefactshared.find(parameter)) != (rarefactshared.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = rarefactshared.begin(); it != rarefactshared.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "summary.single") {
			//is it valid
			if ((summarysingle.find(parameter)) != (summarysingle.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = summarysingle.begin(); it != summarysingle.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "summary.shared") {
			//is it valid
			if ((summaryshared.find(parameter)) != (summaryshared.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = summaryshared.begin(); it != summaryshared.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "unifrac.weighted") {
			//is it valid
			if ((unifracweighted.find(parameter)) != (unifracweighted.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = unifracweighted.begin(); it != unifracweighted.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "unifrac.unweighted") {
			//is it valid
			if ((unifracunweighted.find(parameter)) != (unifracunweighted.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = unifracunweighted.begin(); it != unifracunweighted.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		}else if (command == "libshuff") {
			//is it valid
			if ((libshuff.find(parameter)) != (libshuff.end())) {
				return true;
			}else { 
				cout << parameter << " is not a valid parameter for the " + command + " command. Valid parameters are ";
				for (it = libshuff.begin(); it != libshuff.end(); it++) {
					cout << it->first << ", ";
				}
				cout << endl;
				return false; }
		//not a valid paramter
		}else if (command == "help") { cout << parameter << " is not a valid parameter for the " + command + " command. There are no vaild parameters." << endl;  
		}else if (command == "quit") { cout << parameter << " is not a valid parameter for the " + command + " command. There are no vaild parameters." << endl; 
		}else if (command == "get.group") { cout << parameter << " is not a valid parameter for the " + command + " command. There are no vaild parameters." << endl; 
		}else if (command == "get.label") { cout << parameter << " is not a valid parameter for the " + command + " command. There are no vaild parameters." << endl; 
		}else if (command == "get.line") { cout << parameter << " is not a valid parameter for the " + command + " command. There are no vaild parameters." << endl; }
		
		return false; 
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/
void ValidParameters::initialReaddist() {
	readdist["phylip"]		= "phylip";
	readdist["column"]		= "column";
	readdist["name"]		= "name"; 
	readdist["group"]		= "group"; 
	readdist["cutoff"]		= "cutoff"; 
	readdist["precision"]	= "precision"; 
}
/***********************************************************************/
void ValidParameters::initialReadotu() {
	readotu["list"]			= "list"; 
	readotu["rabund"]		= "rabund"; 
	readotu["sabund"]		= "sabund";
	readotu["shared"]		= "shared";
	readotu["group"]		= "group"; 
	readotu["order"]		= "order"; 
	readotu["label"]		= "label"; 
	readotu["line"]			= "line";
}
/***********************************************************************/
void ValidParameters::initialReadtree() {
	readtree["group"]	= "group"; 
	readtree["tree"]	= "tree";
}
/***********************************************************************/
void ValidParameters::initialCluster() {
	cluster["cutoff"]		= "cutoff"; 
	cluster["method"]		= "method";
	cluster["precision"]	= "precision"; 
}
/***********************************************************************/
void ValidParameters::initialDeconvolute() {
	deconvolute["fasta"]	= "fasta"; 
}
/***********************************************************************/
void ValidParameters::initialParsimony() {
	parsimony["iters"]		= "iters"; 
	parsimony["random"]		= "random";
	parsimony["groups"]		= "groups";
}
/***********************************************************************/
void ValidParameters::initialCollectsingle() {
	collectsingle["label"]		= "label"; 
	collectsingle["line"]		= "line";
	collectsingle["freq"]		= "freq"; 
	collectsingle["calc"]		= "calc";
}
/***********************************************************************/
void ValidParameters::initialCollectshared() {
	collectshared["label"]		= "label"; 
	collectshared["line"]		= "line";
	collectshared["freq"]		= "freq"; 
	collectshared["calc"]		= "calc";
	collectshared["jumble"]		= "jumble";
}
/***********************************************************************/
void ValidParameters::initialRarefactsingle() {
	rarefactsingle["label"]		= "label"; 
	rarefactsingle["line"]		= "line";
	rarefactsingle["freq"]		= "freq"; 
	rarefactsingle["calc"]		= "calc";
	rarefactsingle["iters"]		= "iters"; 
}
/***********************************************************************/
void ValidParameters::initialRarefactshared() {
	rarefactshared["label"]		= "label"; 
	rarefactshared["line"]		= "line";
	rarefactshared["jumble"]	= "jumble";
	rarefactshared["calc"]		= "calc";
	rarefactshared["iters"]		= "iters"; 
}
/***********************************************************************/
void ValidParameters::initialSummarysingle() {
	summarysingle["label"]		= "label"; 
	summarysingle["line"]		= "line";
	summarysingle["calc"]		= "calc";
}
/***********************************************************************/
void ValidParameters::initialSummaryshared() {
	summaryshared["label"]		= "label"; 
	summaryshared["line"]		= "line";
	summaryshared["calc"]		= "calc";
	summaryshared["jumble"]		= "jumble";

}
/***********************************************************************/
void ValidParameters::initialUnifracweighted() {
	unifracweighted["iters"]		= "iters"; 
	unifracweighted["groups"]		= "groups";
}
/***********************************************************************/
void ValidParameters::initialUnifracunweighted() {
	unifracunweighted["iters"]		= "iters"; 
	unifracunweighted["groups"]		= "groups";
}
/***********************************************************************/
void ValidParameters::initialLibshuff() {
	libshuff["cutoff"]		= "cutoff"; 
	libshuff["iters"]		= "iters"; 
	libshuff["groups"]		= "groups";
	libshuff["step"]		= "step";
	libshuff["form"]		= "form";
}
/***********************************************************************/
