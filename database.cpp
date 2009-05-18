/*
 *  database.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;


#include "sequence.hpp"
#include "database.hpp"

/**************************************************************************************************/

Database::Database(string fastaFileName){		//	This assumes that the template database is in fasta format, may 
												//	need to alter this in the future?

	ifstream fastaFile;
	openInputFile(fastaFileName, fastaFile);

	cout << endl << "Reading in the " << fastaFileName << " template sequences...\t";	cout.flush();

	numSeqs=count(istreambuf_iterator<char>(fastaFile),istreambuf_iterator<char>(), '>');	//	count the number of
	fastaFile.seekg(0);																		//	sequences
	
	templateSequences.resize(numSeqs);
	
	string seqName, sequence;
	for(int i=0;i<numSeqs;i++){
		templateSequences[i] = new Sequence();		//	We're storing the aligned template sequences as a vector of
													//	pointers to Sequence objects
		fastaFile >> seqName;
		templateSequences[i]->setName(seqName);
		
		char letter;
		string aligned;
		
		while(fastaFile && (letter=fastaFile.get()) != '>'){
			if(isprint(letter)){
				letter = toupper(letter);
				if(letter == 'U'){letter = 'T';}
				aligned += letter;
			}
		}
		templateSequences[i]->setAligned(aligned);	//	Obviously, we need the fully aligned template sequence
		templateSequences[i]->setUnaligned(aligned);//	We will also need the unaligned sequence for pairwise alignments
		fastaFile.putback(letter);					//	and database searches
	}
	
	fastaFile.close();
	cout << "DONE." << endl;	cout.flush();

}

/**************************************************************************************************/

float Database::getSearchScore()	{	return searchScore;		}	//	we're assuming that the search is already done

/**************************************************************************************************/
