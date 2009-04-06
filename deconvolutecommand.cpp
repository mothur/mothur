/*
 *  deconvolute.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "deconvolutecommand.h"

/**************************************************************************************/
int DeconvoluteCommand::execute() {	
	try {
		globaldata = GlobalData::getInstance();
	
		//prepare filenames and open files
		filename = globaldata->getFastaFile();
		outputFileName = (getRootName(filename) + "names");
		outFastafile = (getRootName(filename) + "unique.fasta");
		
		openInputFile(filename, in);
		openOutputFile(outputFileName, out);
		openOutputFile(outFastafile, outFasta);
	
		//constructor reads in file and store internally
		fastamap = new FastaMap();
	
		//two columns separated by tabs sequence name and then sequence
		fastamap->readFastaFile(in);
		
		//print out new names file 
		//file contains 2 columns separated by tabs.  the first column is the groupname(name of first sequence found.
		//the second column is the list of names of identical sequences separated by ','.
		fastamap->printNamesFile(out);
		fastamap->printCondensedFasta(outFasta);
	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DeconvoluteCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DeconvoluteCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************/
