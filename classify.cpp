/*
 *  classify.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "classify.h"
#include "sequence.hpp"
#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "blastdb.hpp"
#include "distancedb.hpp"

/**************************************************************************************************/
void Classify::generateDatabaseAndNames(string tfile, string tempFile, string method, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch)  {		
	try {	
		taxFile = tfile;
		readTaxonomy(taxFile);	
		
		templateFile = tempFile;	
		
		int start = time(NULL);
		int numSeqs = 0;
		
		m->mothurOut("Generating search database...    "); cout.flush();
#ifdef USE_MPI	
			int pid;
			vector<long> positions;
		
			MPI_Status status; 
			MPI_File inMPI;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

			//char* inFileName = new char[tempFile.length()];
			//memcpy(inFileName, tempFile.c_str(), tempFile.length());
			
			char inFileName[1024];
			strcpy(inFileName, tempFile.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			//delete inFileName;

			if (pid == 0) { //only one process needs to scan file
				positions = setFilePosFasta(tempFile, numSeqs); //fills MPIPos, returns numSeqs

				//send file positions to all processes
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
				MPI_Bcast(&positions[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
			}else{
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD); //get numSeqs
				positions.resize(numSeqs);
				MPI_Bcast(&positions[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //get file positions
			}
			
			//create database
			if(method == "kmer")			{	database = new KmerDB(tempFile, kmerSize);			}
			else if(method == "suffix")		{	database = new SuffixDB(numSeqs);								}
			else if(method == "blast")		{	database = new BlastDB(gapOpen, gapExtend, match, misMatch);	}
			else if(method == "distance")	{	database = new DistanceDB();	}
			else {
				m->mothurOut(method + " is not a valid search option. I will run the command using kmer, ksize=8."); m->mothurOutEndLine();
				database = new KmerDB(tempFile, 8);
			}

			//read file 
			for(int i=0;i<numSeqs;i++){
				//read next sequence
				int length = positions[i+1] - positions[i];
				char* buf4 = new char[length];
				MPI_File_read_at(inMPI, positions[i], buf4, length, MPI_CHAR, &status);
				
				string tempBuf = buf4;
				if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
				delete buf4;
				istringstream iss (tempBuf,istringstream::in);
				
				Sequence temp(iss);  
				if (temp.getName() != "") {
					names.push_back(temp.getName());
					database->addSequence(temp);	
				}
			}
			
			database->generateDB();
			MPI_File_close(&inMPI);
	#else
		
		//need to know number of template seqs for suffixdb
		if (method == "suffix") {
			ifstream inFASTA;
			openInputFile(tempFile, inFASTA);
			numSeqs = count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
		}

		bool needToGenerate = true;
		string kmerDBName;
		if(method == "kmer")			{	
			database = new KmerDB(tempFile, kmerSize);			
			
			kmerDBName = tempFile.substr(0,tempFile.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTest(kmerDBName.c_str());
			if(kmerFileTest){	needToGenerate = false;		}
		}
		else if(method == "suffix")		{	database = new SuffixDB(numSeqs);								}
		else if(method == "blast")		{	database = new BlastDB(gapOpen, gapExtend, match, misMatch);	}
		else if(method == "distance")	{	database = new DistanceDB();	}
		else {
			m->mothurOut(method + " is not a valid search option. I will run the command using kmer, ksize=8.");
			m->mothurOutEndLine();
			database = new KmerDB(tempFile, 8);
		}
		
		if (needToGenerate) {
			ifstream fastaFile;
			openInputFile(tempFile, fastaFile);
			
			while (!fastaFile.eof()) {
				Sequence temp(fastaFile);
				gobble(fastaFile);
			
				names.push_back(temp.getName());
								
				database->addSequence(temp);	
			}
			fastaFile.close();

			database->generateDB();
			
		}else if ((method == "kmer") && (!needToGenerate)) {	
			ifstream kmerFileTest(kmerDBName.c_str());
			database->readKmerDB(kmerFileTest);	
			
			ifstream fastaFile;
			openInputFile(tempFile, fastaFile);
			
			while (!fastaFile.eof()) {
				Sequence temp(fastaFile);
				gobble(fastaFile);

				names.push_back(temp.getName());
			}
			fastaFile.close();
		}
#endif		
		database->setNumSeqs(names.size());
		
		m->mothurOut("DONE."); m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " seconds generate search database. "); m->mothurOutEndLine();

	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "generateDatabaseAndNames");
		exit(1);
	}
}
/**************************************************************************************************/
Classify::Classify() {		m = MothurOut::getInstance();	database = NULL;	}
/**************************************************************************************************/

int Classify::readTaxonomy(string file) {
	try {
		
		phyloTree = new PhyloTree();
		string name, taxInfo;
		
		m->mothurOutEndLine();
		m->mothurOut("Reading in the " + file + " taxonomy...\t");	cout.flush();

#ifdef USE_MPI	
		int pid, num;
		vector<long> positions;
		
		MPI_Status status; 
		MPI_File inMPI;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

		//char* inFileName = new char[file.length()];
		//memcpy(inFileName, file.c_str(), file.length());
		
		char inFileName[1024];
		strcpy(inFileName, file.c_str());

		MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
		//delete inFileName;

		if (pid == 0) {
			positions = setFilePosEachLine(file, num);
			
			//send file positions to all processes
			MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
			MPI_Bcast(&positions[0], (num+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
		}else{
			MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD); //get numSeqs
			positions.resize(num);
			MPI_Bcast(&positions[0], (num+1), MPI_LONG, 0, MPI_COMM_WORLD); //get file positions
		}
	
		//read file 
		for(int i=0;i<num;i++){
			//read next sequence
			int length = positions[i+1] - positions[i];
			char* buf4 = new char[length];

			MPI_File_read_at(inMPI, positions[i], buf4, length, MPI_CHAR, &status);

			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
			delete buf4;

			istringstream iss (tempBuf,istringstream::in);
			iss >> name >> taxInfo;
			taxonomy[name] = taxInfo;
			phyloTree->addSeqToTree(name, taxInfo);
		}
		
		MPI_File_close(&inMPI);
#else				
		ifstream inTax;
		openInputFile(file, inTax);
	
		//read template seqs and save
		while (!inTax.eof()) {
			inTax >> name >> taxInfo;
			
			taxonomy[name] = taxInfo;
			
			phyloTree->addSeqToTree(name, taxInfo);
		
			gobble(inTax);
		}
		inTax.close();
#endif	
	
		phyloTree->assignHeirarchyIDs(0);
		
		phyloTree->setUp(file);
		
		m->mothurOut("DONE.");
		m->mothurOutEndLine();	cout.flush();
		
		return phyloTree->getNumSeqs();
	
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "readTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

vector<string> Classify::parseTax(string tax) {
	try {
		vector<string> taxons;
		
		tax = tax.substr(0, tax.length()-1);  //get rid of last ';'
	
		//parse taxonomy
		string individual;
		while (tax.find_first_of(';') != -1) {
			individual = tax.substr(0,tax.find_first_of(';'));
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			taxons.push_back(individual);
			
		}
		//get last one
		taxons.push_back(tax);
		
		return taxons;
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "parseTax");
		exit(1);
	}
}
/**************************************************************************************************/

