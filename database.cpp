/*
 *  database.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include "mothur.h"
#include "sequence.hpp"
#include "database.hpp"

/**************************************************************************************************/

Database::Database(string fastaFileName){
	
	ifstream fastaFile(fastaFileName.c_str());
	if(!fastaFile) {
		cerr << "Error: Could not open " << fastaFileName << endl;
		exit(1);
	}
	cout << endl << "Reading in the " << fastaFileName << " template sequences...\t";	cout.flush();

	numSeqs=count(istreambuf_iterator<char>(fastaFile),istreambuf_iterator<char>(), '>');
	fastaFile.seekg(0);
	
	templateSequences.resize(numSeqs);
	
	string seqName, sequence;
	for(int i=0;i<numSeqs;i++){
		templateSequences[i] = new Sequence();
		
		fastaFile >> seqName;
		templateSequences[i]->setName(seqName);
		
		char letter;
		string aligned;
		
		while(fastaFile && (letter=fastaFile.get()) != '>'){
			if(isprint(letter)){
				letter = toupper(letter);
				aligned += letter;
			}
		}
		templateSequences[i]->setAligned(aligned);
		templateSequences[i]->setUnaligned(aligned);
		fastaFile.putback(letter);
	}
	
	fastaFile.close();
	cout << "DONE." << endl;	cout.flush();

}

/**************************************************************************************************/
