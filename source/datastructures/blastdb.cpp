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

BlastDB::BlastDB(string tag, float gO, float gE, float mm, float mM, string b, int tid) : Database(), 
gapOpen(gO), gapExtend(gE), match(mm), misMatch(mM) {
	try {
		count = 0;
		path = b;
		threadID = tid;
        
        Utils util;
		int randNumber = util.getRandomNumber();
		string pid = toString(threadID);
		
        if (m->getDebug()) { m->mothurOut("[DEBUG]: tag = " + tag + "\t pid = " + pid + "\n"); }
        
		dbFileName = tag + pid + toString(randNumber) + ".template.unaligned.fasta";
		queryFileName = tag + pid + toString(randNumber) + ".candidate.unaligned.fasta";
		blastFileName = tag + pid + toString(randNumber) + ".blast";
		
		//make sure blast exists in the write place
		if (path == "") {  path = current->getBlastPath() + "blast" + PATH_SEPARATOR + "bin" + PATH_SEPARATOR; }

		//test to make sure formatdb exists
		ifstream in;
        string formatdbCommand = path + "formatdb" + EXECUTABLE_EXT;
		formatdbCommand = util.getFullPathName(formatdbCommand);
		bool ableToOpen = util.openInputFile(formatdbCommand, in, "no error"); in.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " + formatdbCommand + " file does not exist. mothur requires formatdb.exe.\n"); m->setControl_pressed(true); }
		
		//test to make sure formatdb exists
		ifstream in2;
        string blastCommand = path + "blastall" + EXECUTABLE_EXT;
		blastCommand = util.getFullPathName(blastCommand);
		ableToOpen = util.openInputFile(blastCommand, in2, "no error"); in2.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe.\n");  m->setControl_pressed(true);  }

		//test to make sure formatdb exists
		ifstream in3;
        string megablastCommand = path + "megablast" + EXECUTABLE_EXT;
		megablastCommand = util.getFullPathName(megablastCommand);
		ableToOpen = util.openInputFile(megablastCommand, in3, "no error"); in3.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " +  megablastCommand + " file does not exist. mothur requires megablast.exe.\n");  m->setControl_pressed(true);  }
        
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "BlastDB");
		exit(1);
	}
}
/**************************************************************************************************/

BlastDB::BlastDB(string b, int tid) : Database() {
	try {
		count = 0;
		path = b;
		threadID = tid;
        Utils util;
		
		//make sure blast exists in the write place
		if (path == "") {  path = current->getBlastPath() + "blast" + PATH_SEPARATOR + "bin" + PATH_SEPARATOR; }
		
		int randNumber = util.getRandomNumber();;
		string pid = toString(threadID);
		dbFileName = pid + toString(randNumber) + ".template.unaligned.fasta";
		queryFileName = pid + toString(randNumber) + ".candidate.unaligned.fasta";
		blastFileName = pid + toString(randNumber) + ".blast";
		
        //test to make sure formatdb exists
        ifstream in;
        string formatdbCommand = path + "formatdb" + EXECUTABLE_EXT;
        formatdbCommand = util.getFullPathName(formatdbCommand);
        bool ableToOpen = util.openInputFile(formatdbCommand, in, "no error"); in.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " +  formatdbCommand + " file does not exist. mothur requires formatdb.exe.\n");  m->setControl_pressed(true);  }
		
        //test to make sure formatdb exists
        ifstream in2;
        string blastCommand = path + "blastall" + EXECUTABLE_EXT;
        blastCommand = util.getFullPathName(blastCommand);
        ableToOpen = util.openInputFile(blastCommand, in2, "no error"); in2.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe.\n");  m->setControl_pressed(true); }
		
        //test to make sure formatdb exists
        ifstream in3;
        string megablastCommand = path + "megablast" + EXECUTABLE_EXT;
        megablastCommand = util.getFullPathName(megablastCommand);
        ableToOpen = util.openInputFile(megablastCommand, in3, "no error"); in3.close();
		if(!ableToOpen) {	m->mothurOut("[ERROR]: " + megablastCommand + " file does not exist. mothur requires megablast.exe.\n");  m->setControl_pressed(true);  }
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "BlastDB");
		exit(1);
	}
}

/**************************************************************************************************/

BlastDB::~BlastDB(){
	try{
        Utils util;
		util.mothurRemove(queryFileName);				//	let's clean stuff up and remove the temp files
		util.mothurRemove(dbFileName);					//	let's clean stuff up and remove the temp files
		util.mothurRemove((dbFileName+".nsq"));					//	let's clean stuff up and remove the temp files
		util.mothurRemove((dbFileName+".nsi"));					//	let's clean stuff up and remove the temp files
		util.mothurRemove((dbFileName+".nsd"));					//	let's clean stuff up and remove the temp files
		util.mothurRemove((dbFileName+".nin"));					//	let's clean stuff up and remove the temp files
		util.mothurRemove((dbFileName+".nhr"));					//	let's clean stuff up and remove the temp files
		util.mothurRemove(blastFileName.c_str());				//	let's clean stuff up and remove the temp files
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "~BlastDB");
		exit(1);
	}
}
/**************************************************************************************************/
//assumes you have added all the template sequences using the addSequence function and run generateDB.
vector<int> BlastDB::findClosestSequences(Sequence* seq, int n, vector<float>& scores) const {
	try{
		vector<int> topMatches;
		
		ofstream queryFile;
        Utils util;
		string pid = scrubName(seq->getName());
		
        util.openOutputFile((queryFileName+pid), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();

		//lock_guard<std::mutex> guard(mutex);
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
		
		string blastCommand;
		#if defined NON_WINDOWS
		
			blastCommand = path + "blastall -p blastn -d " + dbFileName + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);
			blastCommand += (" -i " + (queryFileName+pid) + " -o " + blastFileName+pid);
		#else
			blastCommand =  "\"" + path + "blastall\" -p blastn -d " + "\"" + dbFileName + "\"" + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);
			blastCommand += (" -i " + (queryFileName+pid+toString(randNumber)) + " -o " + blastFileName+pid+toString(randNumber));
			//wrap entire string in ""
			blastCommand = "\"" + blastCommand + "\"";
		#endif
		
		system(blastCommand.c_str());
		
		ifstream m8FileHandle;
		util.openInputFile(blastFileName+pid, m8FileHandle, "no error");
		
		string dummy;
		int templateAccession;
		util.gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
            float searchScore;
			m8FileHandle >> dummy >> templateAccession >> searchScore;
			
			//get rest of junk in line
			while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}	} 
			
			util.gobble(m8FileHandle);
			topMatches.push_back(templateAccession);
            scores.push_back(searchScore);
		}
		m8FileHandle.close();
		util.mothurRemove((queryFileName+pid));
		util.mothurRemove((blastFileName+pid));

		return topMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "findClosestSequences");
		exit(1);
	}

}
/**************************************************************************************************/
//assumes you have added all the template sequences using the addSequence function and run generateDB.
vector<int> BlastDB::findClosestMegaBlast(Sequence* seq, int n, int minPerID) {
	try{
		vector<int> topMatches;
		float numBases, mismatch, gap, startQuery, endQuery, startRef, endRef, score, searchScore;
		
		
		ofstream queryFile;
		int randNumber = util.getRandomNumber();
		string pid = scrubName(seq->getName());
		
        Utils util; util.openOutputFile((queryFileName+pid+toString(randNumber)), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();

		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
//7000004128189528left	0	100		66	0	0	1	66	61	126	1e-31	 131	
		string blastCommand;
		#if defined NON_WINDOWS
			blastCommand = path + "megablast -e 1e-10 -d " + dbFileName + " -m 8 -b " + toString(n) + " -v " + toString(n); //-W 28 -p blastn
			blastCommand += (" -i " + (queryFileName+pid+toString(randNumber)) + " -o " + blastFileName+pid+toString(randNumber));
		#else
		//blastCommand = path + "blast\\bin\\megablast -e 1e-10 -d " + dbFileName + " -m 8 -b " + toString(n) + " -v " + toString(n); //-W 28 -p blastn
		//blastCommand += (" -i " + (queryFileName+toString(randNumber)) + " -o " + blastFileName+toString(randNumber));

			blastCommand =  "\"" + path + "megablast\" -e 1e-10 -d " + "\"" + dbFileName + "\"" + " -m 8 -b " + toString(n) + " -v " + toString(n); //-W 28 -p blastn
			blastCommand += (" -i " + (queryFileName+pid+toString(randNumber)) + " -o " + blastFileName+pid+toString(randNumber));
			//wrap entire string in ""
			blastCommand = "\"" + blastCommand + "\"";

		#endif
		system(blastCommand.c_str());

		ifstream m8FileHandle;
		util.openInputFile(blastFileName+pid+toString(randNumber), m8FileHandle, "no error");
	
		string dummy, eScore;
		int templateAccession;
		util.gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;

			util.gobble(m8FileHandle);
			if (searchScore >= minPerID) { 
				topMatches.push_back(templateAccession);
				//Scores.push_back(searchScore);
			}

		}
		m8FileHandle.close();
		util.mothurRemove((queryFileName+pid+toString(randNumber)));
		util.mothurRemove((blastFileName+pid+toString(randNumber)));
		
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
        Utils util; util.openOutputFileAppend(dbFileName, unalignedFastaFile);
	
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
			
		string formatdbCommand;
		
		#if defined NON_WINDOWS
			formatdbCommand = path + "formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability
		#else
			//formatdbCommand = path + "blast\\bin\\formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability

			formatdbCommand = "\"" + path + "formatdb\" -p F -o T -i " + "\"" +  dbFileName + "\"";
			//wrap entire string in ""
			formatdbCommand = "\"" + formatdbCommand + "\"";
		#endif
		
		system(formatdbCommand.c_str());								//	to get the right sequence names, i think. -p F
																	//	option tells formatdb that seqs are DNA, not prot
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "generateDB");
		exit(1);
	}
}
/**************************************************************************************************/
string BlastDB::scrubName(string seqName) const {
	try {
		
		string cleanName = "";
		
		for (int i = 0; i < seqName.length(); i++) {
			if (isalnum(seqName[i])) { cleanName += seqName[i]; }
			else { cleanName += "_";  }
		}
		
		return cleanName;
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "scrubName");
		exit(1);
	}
}
/**************************************************************************************************/
