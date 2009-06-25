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
	
	mothurOutEndLine();
	mothurOut("Reading in the " + fastaFileName + " template sequences...\t");	cout.flush();

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
	
	mothurOut("DONE.");
	mothurOutEndLine();	cout.flush();

}
/**************************************************************************************************/

Database::~Database(){														

		templateSequences.clear();
}

/**************************************************************************************************/

float Database::getSearchScore()	{	return searchScore;		}	//	we're assuming that the search is already done

/**************************************************************************************************/
