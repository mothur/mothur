/*
 *  database.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "sequence.hpp"
#include "database.hpp"

/**************************************************************************************************/

Database::Database(string fastaFileName){		//	This assumes that the template database is in fasta format, may 
												//	need to alter this in the future?

	ifstream fastaFile;
	openInputFile(fastaFileName, fastaFile);

	cout << endl << "Reading in the " << fastaFileName << " template sequences...\t";	cout.flush();

	//all of this is elsewhere already!
	numSeqs=count(istreambuf_iterator<char>(fastaFile),istreambuf_iterator<char>(), '>');	//	count the number of
	fastaFile.seekg(0);																		//	sequences
	
	templateSequences.resize(numSeqs);
	
	string seqName, sequence;
	for(int i=0;i<numSeqs;i++){
		fastaFile >> seqName;
		seqName = seqName.substr(1);
		char letter;
		string aligned;
		
		while(fastaFile && (letter=fastaFile.get()) != '>'){
			if(isprint(letter)){
				letter = toupper(letter);
				if(letter == 'U'){letter = 'T';}
				aligned += letter;
			}
		}
		templateSequences[i] = Sequence(seqName, aligned);
		fastaFile.putback(letter);
	}
	
	fastaFile.close();
	//all of this is elsewhere already!
	
	cout << "DONE." << endl;	cout.flush();

}
/**************************************************************************************************/

Database::~Database(){														
	try {
		
		//for (int i = 0; i < templateSequences.size(); i++) {  delete templateSequences[i];    }
		templateSequences.clear();

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Database class Function ~Database. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Database class function ~Database. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**************************************************************************************************/

float Database::getSearchScore()	{	return searchScore;		}	//	we're assuming that the search is already done

/**************************************************************************************************/
