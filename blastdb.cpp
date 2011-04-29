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
	
	count = 0;

	int randNumber = rand();
	dbFileName = toString(randNumber) + ".template.unaligned.fasta";
	queryFileName = toString(randNumber) + ".candidate.unaligned.fasta";
	blastFileName = toString(randNumber) + ".blast";

}
/**************************************************************************************************/

BlastDB::BlastDB() : Database() {
	try {
		count = 0;

		int randNumber = rand();
		dbFileName = toString(randNumber) + ".template.unaligned.fasta";
		queryFileName = toString(randNumber) + ".candidate.unaligned.fasta";
		blastFileName = toString(randNumber) + ".blast";
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "BlastDB");
		exit(1);
	}
}

/**************************************************************************************************/

BlastDB::~BlastDB(){
	try{
		remove(queryFileName.c_str());				//	let's clean stuff up and remove the temp files
		remove(dbFileName.c_str());					//	let's clean stuff up and remove the temp files
		remove((dbFileName+".nsq").c_str());					//	let's clean stuff up and remove the temp files
		remove((dbFileName+".nsi").c_str());					//	let's clean stuff up and remove the temp files
		remove((dbFileName+".nsd").c_str());					//	let's clean stuff up and remove the temp files
		remove((dbFileName+".nin").c_str());					//	let's clean stuff up and remove the temp files
		remove((dbFileName+".nhr").c_str());					//	let's clean stuff up and remove the temp files
		remove(blastFileName.c_str());				//	let's clean stuff up and remove the temp files
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "~BlastDB");
		exit(1);
	}
}
/**************************************************************************************************/
//assumes you have added all the template sequences using the addSequence function and run generateDB.
vector<int> BlastDB::findClosestSequences(Sequence* seq, int n) {
	try{
		vector<int> topMatches;
		
		ofstream queryFile;
		m->openOutputFile((queryFileName+seq->getName()), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();

				
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
		
		string blastCommand = path + "blast/bin/blastall -p blastn -d " + dbFileName + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);;
		blastCommand += (" -i " + (queryFileName+seq->getName()) + " -o " + blastFileName+seq->getName());
		system(blastCommand.c_str());
		
		ifstream m8FileHandle;
		m->openInputFile(blastFileName+seq->getName(), m8FileHandle, "no error");
		
		string dummy;
		int templateAccession;
		m->gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore;
			
			//get rest of junk in line
			while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}	} 
			
			m->gobble(m8FileHandle);
			topMatches.push_back(templateAccession);
		}
		m8FileHandle.close();
		remove((queryFileName+seq->getName()).c_str());
		remove((blastFileName+seq->getName()).c_str());

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
		Scores.clear();
		
		ofstream queryFile;

		m->openOutputFile((queryFileName+seq->getName()), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();
				
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
	
		string blastCommand = path + "blast/bin/megablast -e 1e-10 -d " + dbFileName + " -m 8 -b " + toString(n) + " -v " + toString(n); //-W 28 -p blastn
		blastCommand += (" -i " + (queryFileName+seq->getName()) + " -o " + blastFileName+seq->getName());
		system(blastCommand.c_str());

		ifstream m8FileHandle;
		m->openInputFile(blastFileName+seq->getName(), m8FileHandle, "no error");
	
		string dummy;
		int templateAccession;
		m->gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore;
			//cout << templateAccession << '\t' << searchScore << endl;
			
			//get rest of junk in line
			while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}	} 
			
			m->gobble(m8FileHandle);
			topMatches.push_back(templateAccession);
			Scores.push_back(searchScore);
//cout << templateAccession << endl;
		}
		m8FileHandle.close();
		remove((queryFileName+seq->getName()).c_str());
		remove((blastFileName+seq->getName()).c_str());
//cout << "\n" ;		
		return topMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "findClosestMegaBlast");
		exit(1);
	}
}
/**************************************************************************************************/
void BlastDB::addSequence(Sequence seq) {
	try {
	
		ofstream unalignedFastaFile;
		m->openOutputFileAppend(dbFileName, unalignedFastaFile);				
	
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
		
		path = m->argv;
		string tempPath = path;
		for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
		path = path.substr(0, (tempPath.find_last_of('m')));
	
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
/**************************************************************************************************/

/**************************************************************************************************/

