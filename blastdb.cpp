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

BlastDB::BlastDB(float gO, float gE, float m, float mM) : Database(), 
gapOpen(gO), gapExtend(gE), match(m), misMatch(mM) {
	
	globaldata = GlobalData::getInstance();
	count = 0;

	int randNumber = rand();
	dbFileName = toString(randNumber) + ".template.unaligned.fasta";
	queryFileName = toString(randNumber) + ".candidate.unaligned.fasta";
	blastFileName = toString(randNumber) + ".blast";

}
/**************************************************************************************************/

BlastDB::BlastDB() : Database() {
	
	globaldata = GlobalData::getInstance();
	count = 0;

	int randNumber = rand();
	dbFileName = toString(randNumber) + ".template.unaligned.fasta";
	queryFileName = toString(randNumber) + ".candidate.unaligned.fasta";
	blastFileName = toString(randNumber) + ".blast";

}

/**************************************************************************************************/

BlastDB::~BlastDB(){
	remove(queryFileName.c_str());				//	let's clean stuff up and remove the temp files
	remove(dbFileName.c_str());					//	let's clean stuff up and remove the temp files
	remove(blastFileName.c_str());				//	let's clean stuff up and remove the temp files
}
/**************************************************************************************************/
//assumes you have added all the template sequences using the addSequence function and run generateDB.
vector<int> BlastDB::findClosestSequences(Sequence* seq, int n) {
	try{
		vector<int> topMatches;
		
		ofstream queryFile;
		openOutputFile(queryFileName, queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();
				
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
		
		string blastCommand = path + "blast/bin/blastall -p blastn -d " + dbFileName + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);;
		blastCommand += (" -i " + queryFileName + " -o " + blastFileName);
		system(blastCommand.c_str());
		
		ifstream m8FileHandle;
		openInputFile(blastFileName, m8FileHandle);
		
		string dummy;
		int templateAccession;
		gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore;
			
			//get rest of junk in line
			while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}	} 
			
			gobble(m8FileHandle);
			topMatches.push_back(templateAccession);
		}
		m8FileHandle.close();
		
		return topMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "findClosestSequences");
		exit(1);
	}

}
/**************************************************************************************************/
//assumes you have added all the template sequences using the addSequence function and run generateDB.
vector<int> BlastDB::findClosestMegaBlast(Sequence* seq, int n) {
	try{
		vector<int> topMatches;
		
		ofstream queryFile;
		openOutputFile(queryFileName, queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();
				
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
	
		string blastCommand = path + "blast/bin/megablast -e 1e-10 -d " + dbFileName + " -m 8 -b " + toString(n) + " -v " + toString(n); //-W 28 -p blastn
		blastCommand += (" -i " + queryFileName + " -o " + blastFileName);
		system(blastCommand.c_str());
		
		ifstream m8FileHandle;
		openInputFile(blastFileName, m8FileHandle, "no error");
	
		string dummy;
		int templateAccession;
		gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore;
			
			//get rest of junk in line
			while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}	} 
			
			gobble(m8FileHandle);
			topMatches.push_back(templateAccession);
//cout << templateAccession << endl;
		}
		m8FileHandle.close();
//cout << "\n\n" ;		
		return topMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "findClosest");
		exit(1);
	}
}
/**************************************************************************************************/
void BlastDB::addSequence(Sequence seq) {
	try {
	
		ofstream unalignedFastaFile;
		openOutputFileAppend(dbFileName, unalignedFastaFile);				
	
		//	generating a fasta file with unaligned template
		unalignedFastaFile << '>' << count << endl;					//	sequences, which will be input to formatdb
		unalignedFastaFile << seq.getUnaligned() << endl;
		unalignedFastaFile.close();
	
		count++;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "addSequence");
		exit(1);
	}
}
/**************************************************************************************************/
void BlastDB::generateDB() {
	try {
	
		//m->mothurOut("Generating the temporary BLAST database...\t");	cout.flush();
		
		path = globaldata->argv;
		path = path.substr(0, (path.find_last_of('m')));
	
		string formatdbCommand = path + "blast/bin/formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability
		system(formatdbCommand.c_str());								//	to get the right sequence names, i think. -p F
																	//	option tells formatdb that seqs are DNA, not prot
		//m->mothurOut("DONE."); m->mothurOutEndLine();	m->mothurOutEndLine(); cout.flush();
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "generateDB");
		exit(1);
	}
}
#ifdef USE_MPI	
/**************************************************************************************************/
int BlastDB::MPISend(int receiver) {
	try {
		
		//send gapOpen - float
		MPI_Send(&gapOpen, 1, MPI_FLOAT, receiver, 2001, MPI_COMM_WORLD); 

		//send gapExtend - float
		MPI_Send(&gapExtend, 1, MPI_FLOAT, receiver, 2001, MPI_COMM_WORLD); 
		
		//send match - float
		MPI_Send(&match, 1, MPI_FLOAT, receiver, 2001, MPI_COMM_WORLD); 
		
		//send mismatch - float
		MPI_Send(&misMatch, 1, MPI_FLOAT, receiver, 2001, MPI_COMM_WORLD); 
									
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "MPISend");
		exit(1);
	}
}
/**************************************************************************************************/
int BlastDB::MPIRecv(int sender) {
	try {
		MPI_Status status;
		
		//receive gapOpen - float
		MPI_Recv(&gapOpen, 1, MPI_FLOAT, sender, 2001, MPI_COMM_WORLD, &status);
		
		//receive gapExtend - float
		MPI_Recv(&gapExtend, 1, MPI_FLOAT, sender, 2001, MPI_COMM_WORLD, &status);
				
		//receive match - float
		MPI_Recv(&match, 1, MPI_FLOAT, sender, 2001, MPI_COMM_WORLD, &status);
		
		//receive mismatch - float
		MPI_Recv(&misMatch, 1, MPI_FLOAT, sender, 2001, MPI_COMM_WORLD, &status);		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "MPIRecv");
		exit(1);
	}
}
#endif	
/**************************************************************************************************/

/**************************************************************************************************/

