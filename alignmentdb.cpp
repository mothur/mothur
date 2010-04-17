/*
 *  alignmentdb.cpp
 *  Mothur
 *
 *  Created by westcott on 11/4/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "alignmentdb.h"
#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "blastdb.hpp"


/**************************************************************************************************/
AlignmentDB::AlignmentDB(string fastaFileName, string s, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch){		//	This assumes that the template database is in fasta format, may 
	try {											//	need to alter this in the future?
		m = MothurOut::getInstance();
		longest = 0;
		method = s;
		bool needToGenerate = true;
		
		m->mothurOutEndLine();
		m->mothurOut("Reading in the " + fastaFileName + " template sequences...\t");	cout.flush();
		
		#ifdef USE_MPI	
			int pid;
			vector<long> positions;
		
			MPI_Status status; 
			MPI_File inMPI;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
	
			char inFileName[1024];
			strcpy(inFileName, fastaFileName.c_str());
	
			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			
			if (pid == 0) {
				positions = setFilePosFasta(fastaFileName, numSeqs); //fills MPIPos, returns numSeqs

				//send file positions to all processes
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
				MPI_Bcast(&positions[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
			}else{
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD); //get numSeqs
				positions.resize(numSeqs+1);
				MPI_Bcast(&positions[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //get file positions
			}
		
			//read file 
			for(int i=0;i<numSeqs;i++){
				
				if (m->control_pressed) {  templateSequences.clear(); break;  }
				
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
					templateSequences.push_back(temp);
					//save longest base
					if (temp.getUnaligned().length() > longest)  { longest = temp.getUnaligned().length()+1; }
				}
			}
		
			MPI_File_close(&inMPI);
		
	#else
		ifstream fastaFile;
		openInputFile(fastaFileName, fastaFile);

		while (!fastaFile.eof()) {
			Sequence temp(fastaFile);  gobble(fastaFile);
			
			if (m->control_pressed) {  templateSequences.clear(); break;  }
			
			if (temp.getName() != "") {
				templateSequences.push_back(temp);
				//save longest base
				if (temp.getUnaligned().length() > longest)  { longest = temp.getUnaligned().length()+1; }
			}
		}
		fastaFile.close();
		
	#endif
	
		numSeqs = templateSequences.size();
		//all of this is elsewhere already!
		
		m->mothurOut("DONE.");
		m->mothurOutEndLine();	cout.flush();
		
		//in case you delete the seqs and then ask for them
		emptySequence = Sequence();
		emptySequence.setName("no_match");
		emptySequence.setUnaligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		emptySequence.setAligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		
		
		string kmerDBName;
		if(method == "kmer")			{	
			search = new KmerDB(fastaFileName, kmerSize);			
			
			#ifdef USE_MPI
			#else
				kmerDBName = fastaFileName.substr(0,fastaFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
				ifstream kmerFileTest(kmerDBName.c_str());
			
				if(kmerFileTest){	needToGenerate = false;		}
			#endif
		}
		else if(method == "suffix")		{	search = new SuffixDB(numSeqs);								}
		else if(method == "blast")		{	search = new BlastDB(gapOpen, gapExtend, match, misMatch);	}
		else {
			m->mothurOut(method + " is not a valid search option. I will run the command using kmer, ksize=8.");
			m->mothurOutEndLine();
			search = new KmerDB(fastaFileName, 8);
		}
		
		if (!(m->control_pressed)) {
			if (needToGenerate) {
				//add sequences to search 
				for (int i = 0; i < templateSequences.size(); i++) {
					search->addSequence(templateSequences[i]);
					
					if (m->control_pressed) {  templateSequences.clear(); break;  }
				}
				
				if (m->control_pressed) {  templateSequences.clear();  }
				
				search->generateDB();
				
			}else if ((method == "kmer") && (!needToGenerate)) {
				ifstream kmerFileTest(kmerDBName.c_str());
				search->readKmerDB(kmerFileTest);	
			}
			
			search->setNumSeqs(numSeqs);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "AlignmentDB", "AlignmentDB");
		exit(1);
	}
}
/**************************************************************************************************/
AlignmentDB::AlignmentDB(string s){		 
	try {											
		m = MothurOut::getInstance();
		method = s;
		
		if(method == "suffix")		{	search = new SuffixDB();	}
		else if(method == "blast")	{	search = new BlastDB();		}
		else						{	search = new KmerDB();		}

				
		//in case you delete the seqs and then ask for them
		emptySequence = Sequence();
		emptySequence.setName("no_match");
		emptySequence.setUnaligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		emptySequence.setAligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignmentDB", "AlignmentDB");
		exit(1);
	}
}
/**************************************************************************************************/
AlignmentDB::~AlignmentDB() {  delete search;	}
/**************************************************************************************************/
Sequence AlignmentDB::findClosestSequence(Sequence* seq) {
	try{
	
		vector<int> spot = search->findClosestSequences(seq, 1);

		if (spot.size() != 0)	{		return templateSequences[spot[0]];	}
		else					{		return emptySequence;				}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignmentDB", "findClosestSequence");
		exit(1);
	}
}
#ifdef USE_MPI	
/**************************************************************************************************/
int AlignmentDB::MPISend(int receiver) {
	try {
		
		//send numSeqs - int
		MPI_Send(&numSeqs, 1, MPI_INT, receiver, 2001, MPI_COMM_WORLD); 
									
		//send longest - int
		MPI_Send(&longest, 1, MPI_INT, receiver, 2001, MPI_COMM_WORLD); 
	
		//send templateSequences
		for (int i = 0; i < templateSequences.size(); i++) {
			templateSequences[i].MPISend(receiver);
		}
		
		//send Database
		search->MPISend(receiver);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignmentDB", "MPISend");
		exit(1);
	}
}
/**************************************************************************************************/
int AlignmentDB::MPIRecv(int sender) {
	try {
		MPI_Status status;
		//receive numSeqs - int
		MPI_Recv(&numSeqs, 1, MPI_INT, sender, 2001, MPI_COMM_WORLD, &status);
		
		//receive longest - int
		MPI_Recv(&longest, 1, MPI_INT, sender, 2001, MPI_COMM_WORLD, &status);

		//receive templateSequences
		templateSequences.resize(numSeqs);
		for (int i = 0; i < templateSequences.size(); i++) {
			templateSequences[i].MPIRecv(sender);
		}

		//receive Database
		search->MPIRecv(sender);
	
		for (int i = 0; i < templateSequences.size(); i++) {
			search->addSequence(templateSequences[i]);
		}
		search->generateDB();
		search->setNumSeqs(numSeqs);

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignmentDB", "MPIRecv");
		exit(1);
	}
}
#endif
/**************************************************************************************************/






