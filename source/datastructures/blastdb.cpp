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

		int randNumber = m->getRandomNumber();;
		//int randNumber = 12345;
		string pid = m->mothurGetpid(threadID);
		
        if (m->debug) { m->mothurOut("[DEBUG]: tag = " + tag + "\t pid = " + pid + "\n"); }
        
		dbFileName = tag + pid + toString(randNumber) + ".template.unaligned.fasta";
		queryFileName = tag + pid + toString(randNumber) + ".candidate.unaligned.fasta";
		blastFileName = tag + pid + toString(randNumber) + ".blast";
		
		//make sure blast exists in the write place
		if (path == "") {
			path = m->getBlastPath();
			//string tempPath = path;
			//for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
			//path = path.substr(0, (tempPath.find_last_of('m')));
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			path += "blast/bin/";	
#else
			path += "blast\\bin\\";
#endif			
		}
		
		
		string formatdbCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		formatdbCommand = path + "formatdb";	//	format the database, -o option gives us the ability
#else
		formatdbCommand = path + "formatdb.exe";
#endif
		
		//test to make sure formatdb exists
		ifstream in;
		formatdbCommand = m->getFullPathName(formatdbCommand);
		int ableToOpen = m->openInputFile(formatdbCommand, in, "no error"); in.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + formatdbCommand + " file does not exist. mothur requires formatdb.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		string blastCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		blastCommand = path + "blastall";	//	format the database, -o option gives us the ability
#else
		blastCommand = path + "blastall.exe";
		//wrap entire string in ""
		//blastCommand = "\"" + blastCommand + "\"";
#endif
		
		//test to make sure formatdb exists
		ifstream in2;
		blastCommand = m->getFullPathName(blastCommand);
		ableToOpen = m->openInputFile(blastCommand, in2, "no error"); in2.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		
		string megablastCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		megablastCommand = path + "megablast";	//	format the database, -o option gives us the ability
#else
		megablastCommand = path + "megablast.exe";
#endif
		
		//test to make sure formatdb exists
		ifstream in3;
		megablastCommand = m->getFullPathName(megablastCommand);
		ableToOpen = m->openInputFile(megablastCommand, in3, "no error"); in3.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " +  megablastCommand + " file does not exist. mothur requires megablast.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
        
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
		
		//make sure blast exists in the write place
		if (path == "") {
			path = m->getBlastPath();
			//string tempPath = path;
			//for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
			//path = path.substr(0, (tempPath.find_last_of('m')));
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			path += "blast/bin/";	
#else
			path += "blast\\bin\\";
#endif			
		}
		
		int randNumber = m->getRandomNumber();;
		string pid = m->mothurGetpid(threadID);
		dbFileName = pid + toString(randNumber) + ".template.unaligned.fasta";
		queryFileName = pid + toString(randNumber) + ".candidate.unaligned.fasta";
		blastFileName = pid + toString(randNumber) + ".blast";
		
		string formatdbCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		formatdbCommand = path + "formatdb";	//	format the database, -o option gives us the ability
#else
		formatdbCommand = path + "formatdb.exe";
		//wrap entire string in ""
		//formatdbCommand = "\"" + formatdbCommand + "\"";
#endif
		
		//test to make sure formatdb exists
		ifstream in;
		formatdbCommand = m->getFullPathName(formatdbCommand);
		int ableToOpen = m->openInputFile(formatdbCommand, in, "no error"); in.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " +  formatdbCommand + " file does not exist. mothur requires formatdb.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		string blastCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		blastCommand = path + "blastall";	//	format the database, -o option gives us the ability
#else
		blastCommand = path + "blastall.exe";
		//wrap entire string in ""
		//blastCommand = "\"" + blastCommand + "\"";
#endif
		
		//test to make sure formatdb exists
		ifstream in2;
		blastCommand = m->getFullPathName(blastCommand);
		ableToOpen = m->openInputFile(blastCommand, in2, "no error"); in2.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		
		string megablastCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		megablastCommand = path + "megablast";	//	format the database, -o option gives us the ability
#else
		megablastCommand = path + "megablast.exe";
		//wrap entire string in ""
		//megablastCommand = "\"" + megablastCommand + "\"";
#endif
		
		//test to make sure formatdb exists
		ifstream in3;
		megablastCommand = m->getFullPathName(megablastCommand);
		ableToOpen = m->openInputFile(megablastCommand, in3, "no error"); in3.close();
		if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + megablastCommand + " file does not exist. mothur requires megablast.exe."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "BlastDB", "BlastDB");
		exit(1);
	}
}

/**************************************************************************************************/

BlastDB::~BlastDB(){
	try{
		m->mothurRemove(queryFileName);				//	let's clean stuff up and remove the temp files
		m->mothurRemove(dbFileName);					//	let's clean stuff up and remove the temp files
		m->mothurRemove((dbFileName+".nsq"));					//	let's clean stuff up and remove the temp files
		m->mothurRemove((dbFileName+".nsi"));					//	let's clean stuff up and remove the temp files
		m->mothurRemove((dbFileName+".nsd"));					//	let's clean stuff up and remove the temp files
		m->mothurRemove((dbFileName+".nin"));					//	let's clean stuff up and remove the temp files
		m->mothurRemove((dbFileName+".nhr"));					//	let's clean stuff up and remove the temp files
		m->mothurRemove(blastFileName.c_str());				//	let's clean stuff up and remove the temp files
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
		int randNumber = m->getRandomNumber();;
		string pid = scrubName(seq->getName());
		
		m->openOutputFile((queryFileName+pid+toString(randNumber)), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();

				
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
		
		string blastCommand;
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
			blastCommand = path + "blastall -p blastn -d " + dbFileName + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);
			blastCommand += (" -i " + (queryFileName+pid+toString(randNumber)) + " -o " + blastFileName+pid+toString(randNumber));
		#else
			blastCommand =  "\"" + path + "blastall\" -p blastn -d " + "\"" + dbFileName + "\"" + " -m 8 -W 28 -v " + toString(n) + " -b " + toString(n);
			blastCommand += (" -i " + (queryFileName+pid+toString(randNumber)) + " -o " + blastFileName+pid+toString(randNumber));
			//wrap entire string in ""
			blastCommand = "\"" + blastCommand + "\"";
		#endif
		
		system(blastCommand.c_str());
		
		ifstream m8FileHandle;
		m->openInputFile(blastFileName+pid+toString(randNumber), m8FileHandle, "no error");
		
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
		m->mothurRemove((queryFileName+pid+toString(randNumber)));
		m->mothurRemove((blastFileName+pid+toString(randNumber)));

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
		float numBases, mismatch, gap, startQuery, endQuery, startRef, endRef, score;
		Scores.clear();
		
		ofstream queryFile;
		int randNumber = m->getRandomNumber();
		string pid = scrubName(seq->getName());
		
		m->openOutputFile((queryFileName+pid+toString(randNumber)), queryFile);
		queryFile << '>' << seq->getName() << endl;
		queryFile << seq->getUnaligned() << endl;
		queryFile.close();
//		cout << seq->getUnaligned() << endl;
		//	the goal here is to quickly survey the database to find the closest match.  To do this we are using the default
		//	wordsize used in megablast.  I'm sure we're sacrificing accuracy for speed, but anyother way would take way too
		//	long.  With this setting, it seems comparable in speed to the suffix tree approach.
//7000004128189528left	0	100		66	0	0	1	66	61	126	1e-31	 131	
		string blastCommand;
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
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
		m->openInputFile(blastFileName+pid+toString(randNumber), m8FileHandle, "no error");
	
		string dummy, eScore;
		int templateAccession;
		m->gobble(m8FileHandle);
		
		while(!m8FileHandle.eof()){
			m8FileHandle >> dummy >> templateAccession >> searchScore >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
//			cout << dummy << '\t' << templateAccession << '\t' << searchScore << '\t' << numBases << '\t' << mismatch << '\t' << gap << '\t' << startQuery << '\t' << endQuery << '\t' << startRef << '\t' << endRef << '\t' << eScore << '\t' << score << endl; 
			
			//get rest of junk in line
			//while (!m8FileHandle.eof())	{	char c = m8FileHandle.get(); if (c == 10 || c == 13){	break;	}else{ cout << c; }	} //
				//cout << endl;
			m->gobble(m8FileHandle);
			if (searchScore >= minPerID) { 
				topMatches.push_back(templateAccession);
				Scores.push_back(searchScore);
			}
//cout << templateAccession << endl;
		}
		m8FileHandle.close();
		m->mothurRemove((queryFileName+pid+toString(randNumber)));
		m->mothurRemove((blastFileName+pid+toString(randNumber)));
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
			
		string formatdbCommand;
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			formatdbCommand = path + "formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability
		#else
			//formatdbCommand = path + "blast\\bin\\formatdb -p F -o T -i " + dbFileName;	//	format the database, -o option gives us the ability

			formatdbCommand = "\"" + path + "formatdb\" -p F -o T -i " + "\"" +  dbFileName + "\"";
			//wrap entire string in ""
			formatdbCommand = "\"" + formatdbCommand + "\"";
		#endif
		//cout << formatdbCommand << endl;
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
string BlastDB::scrubName(string seqName) {
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

/**************************************************************************************************/

