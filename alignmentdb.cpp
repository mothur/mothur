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
AlignmentDB::AlignmentDB(string fastaFileName, string method, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch){		//	This assumes that the template database is in fasta format, may 
	try {											//	need to alter this in the future?
		longest = 0;
		cout << longest;
		ifstream fastaFile;
		openInputFile(fastaFileName, fastaFile);
		
		mothurOutEndLine();
		mothurOut("Reading in the " + fastaFileName + " template sequences...\t");	cout.flush();
		
		while (!fastaFile.eof()) {
			Sequence temp(fastaFile);  gobble(fastaFile);
			
			if (temp.getName() != "") {
				templateSequences.push_back(temp);
				//save longest base
				if (temp.getUnaligned().length() > longest)  { longest = temp.getUnaligned().length()+1; }
			}
		}
		
		numSeqs = templateSequences.size();
		
		fastaFile.close();
		//all of this is elsewhere already!
		
		mothurOut("DONE.");
		mothurOutEndLine();	cout.flush();
		
		//in case you delete the seqs and then ask for them
		emptySequence = Sequence();
		emptySequence.setName("no_match");
		emptySequence.setUnaligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		emptySequence.setAligned("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		
		bool needToGenerate = true;
		string kmerDBName;
		if(method == "kmer")			{	
			search = new KmerDB(fastaFileName, kmerSize);			
			
			kmerDBName = fastaFileName.substr(0,fastaFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTest(kmerDBName.c_str());
			
			if(kmerFileTest){	needToGenerate = false;		}
		}
		else if(method == "suffix")		{	search = new SuffixDB(numSeqs);								}
		else if(method == "blast")		{	search = new BlastDB(gapOpen, gapExtend, match, misMatch);	}
		else {
			mothurOut(method + " is not a valid search option. I will run the command using kmer, ksize=8.");
			mothurOutEndLine();
			search = new KmerDB(fastaFileName, 8);
		}
		
		if (needToGenerate) {
		
			//add sequences to search 
			for (int i = 0; i < templateSequences.size(); i++) {
				search->addSequence(templateSequences[i]);
			}
			search->generateDB();
			
		}else if ((method == "kmer") && (!needToGenerate)) {
			ifstream kmerFileTest(kmerDBName.c_str());
			search->readKmerDB(kmerFileTest);	
		}
		
		search->setNumSeqs(numSeqs);
	}
	catch(exception& e) {
		errorOut(e, "AlignmentDB", "AlignmentDB");
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
		errorOut(e, "AlignmentDB", "findClosestSequence");
		exit(1);
	}
}
/**************************************************************************************************/




