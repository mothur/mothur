/*
 *  blastdb.cpp
 *  
 *
 *  Created by Pat Schloss on 12/22/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "database.hpp"
#include "sequence.hpp"
#include "blastdb.hpp"

/**************************************************************************************************/

BlastDB::BlastDB(string fastaFileName, float gO, float gE, float m, float mM) : Database(fastaFileName), 
gapOpen(gO), gapExtend(gE), match(m), misMatch(mM) {
	
	globaldata = GlobalData::getInstance();
	
	cout << "Generating the temporary BLAST database...\t";	cout.flush();

	int randNumber = rand();
	dbFileName = toString(randNumber) + ".template.unaligned.fasta";
	queryFileName = toString(randNumber) + ".candidate.unaligned.fasta";
	blastFileName = toString(randNumber) + ".blast";


	ofstream unalignedFastaFile;
	openOutputFile(dbFileName, unalignedFastaFile);				
	
	for(int i=0;i<numSeqs;i++){									//	generating a fasta file with unaligned template
		unalignedFastaFile << '>' << i << endl;					//	sequences, which will be input to formatdb
		unalignedFastaFile << templateSequences[i]->getUnaligned() << endl;
	}
	unalignedFastaFile.close();
	
	path = globaldata->argv;
	path = path.substr(0, (path.find_last_of('m')));
	
	string formatdbCommand = path + "blast/bin/formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability
	system(formatdbCommand.c_str());								//	to get the right sequence names, i think. -p F
																	//	option tells formatdb that seqs are DNA, not prot
	cout << "DONE." << endl << endl;	cout.flush();
	emptySequence = new Sequence();
	emptySequence->setName("no_match");
	emptySequence->setUnaligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
	emptySequence->setAligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");


}

/**************************************************************************************************/

BlastDB::~BlastDB(){
	remove(queryFileName.c_str());				//	let's clean stuff up and remove the temp files
	remove(dbFileName.c_str());					//	let's clean stuff up and remove the temp files
	remove(blastFileName.c_str());				//	let's clean stuff up and remove the temp files
}

/**************************************************************************************************/

Sequence* BlastDB::findClosestSequence(Sequence* candidate){

	ofstream queryFile;
	openOutputFile(queryFileName, queryFile);
	queryFile << '>' << candidate->getName() << endl;
	queryFile << candidate->getUnaligned() << endl;
	queryFile.close();
	
	
//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
//	long.  With this setting, it seems comparable in speed to the suffix tree approach.

	string blastCommand = path + "blast/bin/blastall -p blastn -d " + dbFileName + " -b 1 -m 8 -W 28";
	blastCommand += (" -i " + queryFileName + " -o " + blastFileName);
	system(blastCommand.c_str());
	
	ifstream m8FileHandle;
	openInputFile(blastFileName, m8FileHandle);
	
	string dummy;
	int templateAccession;
	gobble(m8FileHandle);
	if(!m8FileHandle.eof()){
		m8FileHandle >> dummy >> templateAccession >> searchScore;
		m8FileHandle.close();
		return templateSequences[templateAccession];
	}
	else{
		searchScore = 0.00;
		return emptySequence;
	}
}

/**************************************************************************************************/
