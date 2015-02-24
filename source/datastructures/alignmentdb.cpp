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
#include "referencedb.h"

/**************************************************************************************************/
AlignmentDB::AlignmentDB(string fastaFileName, string s, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch, int tid){		//	This assumes that the template database is in fasta format, may 
	try {											//	need to alter this in the future?
		m = MothurOut::getInstance();
		longest = 0;
		method = s;
		bool needToGenerate = true;
		ReferenceDB* rdb = ReferenceDB::getInstance();
		bool silent = false;
		threadID = tid;
		
		if (fastaFileName == "saved-silent") {
			fastaFileName = "saved"; silent = true;
		}
		
		if (fastaFileName == "saved") {
			int start = time(NULL);
			
			if (!silent) { m->mothurOutEndLine();  m->mothurOut("Using sequences from " + rdb->getSavedReference() + " that are saved in memory.");	m->mothurOutEndLine(); }

			for (int i = 0; i < rdb->referenceSeqs.size(); i++) {
				templateSequences.push_back(rdb->referenceSeqs[i]);
				//save longest base
				if (rdb->referenceSeqs[i].getUnaligned().length() >= longest)  { longest = (rdb->referenceSeqs[i].getUnaligned().length()+1); }
			}
			fastaFileName = rdb->getSavedReference();
			
			numSeqs = templateSequences.size();
			if (!silent) { m->mothurOut("It took " + toString(time(NULL) - start) + " to load " + toString(rdb->referenceSeqs.size()) + " sequences.");m->mothurOutEndLine();  }
            
		}else {
			int start = time(NULL);
			m->mothurOutEndLine();
			m->mothurOut("Reading in the " + fastaFileName + " template sequences...\t");	cout.flush();
			//bool aligned = false;
            int tempLength = 0;
            
			#ifdef USE_MPI	
				int pid, processors;
				vector<unsigned long long> positions;
			
				MPI_Status status; 
				MPI_File inMPI;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				int tag = 2001;
		
				char inFileName[1024];
				strcpy(inFileName, fastaFileName.c_str());
		
				MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				
				if (pid == 0) {
					positions = m->setFilePosFasta(fastaFileName, numSeqs); //fills MPIPos, returns numSeqs

					//send file positions to all processes
					for(int i = 1; i < processors; i++) { 
						MPI_Send(&numSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&positions[0], (numSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					}
				}else{
					MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					positions.resize(numSeqs+1);
					MPI_Recv(&positions[0], (numSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
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
						
						if (rdb->save) { rdb->referenceSeqs.push_back(temp); }
						
						//save longest base
						if (temp.getUnaligned().length() >= longest)  { longest = temp.getUnaligned().length()+1; }
                        if (tempLength != 0) {
                            if (tempLength != temp.getAligned().length()) { m->mothurOut("[ERROR]: template is not aligned, aborting.\n"); m->control_pressed=true; }
                        }else { tempLength = temp.getAligned().length(); }
					}
				}
				
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
				MPI_File_close(&inMPI);
			
		#else
			ifstream fastaFile;
			m->openInputFile(fastaFileName, fastaFile);

			while (!fastaFile.eof()) {
				Sequence temp(fastaFile);  m->gobble(fastaFile);
				
				if (m->control_pressed) {  templateSequences.clear(); break;  }
				
				if (temp.getName() != "") {
					templateSequences.push_back(temp);
					
					if (rdb->save) { rdb->referenceSeqs.push_back(temp); }
					
					//save longest base
					if (temp.getUnaligned().length() >= longest)  { longest = (temp.getUnaligned().length()+1); }
                    
                    if (tempLength != 0) {
                        if (tempLength != temp.getAligned().length()) { m->mothurOut("[ERROR]: template is not aligned, aborting.\n"); m->control_pressed=true; }
                    }else { tempLength = temp.getAligned().length(); }
				}
			}
			fastaFile.close();
		#endif
		
			numSeqs = templateSequences.size();
			//all of this is elsewhere already!
			
			m->mothurOut("DONE.");
			m->mothurOutEndLine();	cout.flush();
			m->mothurOut("It took " + toString(time(NULL) - start) + " to read  " + toString(templateSequences.size()) + " sequences."); m->mothurOutEndLine();  

		}
		
		
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
				
				if(kmerFileTest){	
					bool GoodFile = m->checkReleaseVersion(kmerFileTest, m->getVersion());
					if (GoodFile) {  needToGenerate = false;	}
				}
			#endif
		}
		else if(method == "suffix")		{	search = new SuffixDB(numSeqs);								}
		else if(method == "blast")		{	search = new BlastDB(fastaFileName.substr(0,fastaFileName.find_last_of(".")+1), gapOpen, gapExtend, match, misMatch, "", threadID);	}
		else {
			method = "kmer";
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
		else if(method == "blast")	{	search = new BlastDB("", 0);		}
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
/**************************************************************************************************/






