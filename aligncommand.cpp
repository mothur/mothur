/*
 *  aligncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 *	This version of nast does everything I think that the greengenes nast server does and then some.  I have added the 
 *	feature of allowing users to define their database, kmer size for searching, alignment penalty values and alignment 
 *	method.  This latter feature is perhaps most significant.  nastPlus enables a user to use either a Needleman-Wunsch 
 *	(non-affine gap penalty) or Gotoh (affine gap penalty) pairwise alignment algorithm.  This is significant because it
 *	allows for a global alignment and not the local alignment provided by bLAst.  Furthermore, it has the potential to
 *	provide a better alignment because of the banding method employed by blast (I'm not sure about this).
 *
 *	to compile type:
 *		make
 *
 *	for basic instructions on how to run nastPlus type:
 *		./nastPlus
 */

#include "aligncommand.h"
#include "sequence.hpp"

#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "database.hpp"
#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "blastdb.hpp"

#include "nast.hpp"
#include "nastreport.hpp"


//**********************************************************************************************************************
AlignCommand::AlignCommand(){
	try {
		globaldata = GlobalData::getInstance();
		if(globaldata->getFastaFile() == "" && globaldata->getPhylipFile() == "" && globaldata->getNexusFile() == "" && globaldata->getClustalFile() == ""){
			cout << "you forgot a template file" << endl;
		}
		openInputFile(globaldata->getCandidateFile(), in);
		
		convert(globaldata->getKSize(), kmerSize);
		convert(globaldata->getMatch(), match);
		convert(globaldata->getMismatch(), misMatch);
		convert(globaldata->getGapopen(), gapOpen);
		convert(globaldata->getGapextend(), gapExtend);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function AlignCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function AlignCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

AlignCommand::~AlignCommand(){
	
}

//**********************************************************************************************************************

int AlignCommand::execute(){
	try {
		srand( (unsigned)time( NULL ) );  //needed to assign names to temporary files
		
		Database* templateDB;
		if(globaldata->getSearch() == "kmer")			{	templateDB = new KmerDB(globaldata->getFastaFile() , kmerSize);	}
		else if(globaldata->getSearch() == "suffix")	{	templateDB = new SuffixDB(globaldata->getFastaFile());			}
		else if(globaldata->getSearch() == "blast")		{	templateDB = new BlastDB(globaldata->getFastaFile(), gapOpen, gapExtend, match, misMatch);	}
		else {
			cout << globaldata->getSearch() << " is not a valid search option. I will run the command using kmer, ksize=8." << endl;
			templateDB = new KmerDB(globaldata->getFastaFile(), kmerSize);
		}
	
		Alignment* alignment;
		if(globaldata->getAlign() == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 3000);	}
		else if(globaldata->getAlign() == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);			}
		else if(globaldata->getAlign() == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(globaldata->getAlign() == "noalign")	{	alignment = new NoAlign();													}
		else {
			cout << globaldata->getAlign() << " is not a valid alignment option. I will run the command using needleman." << endl;
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);
		}
				
		int numFastaSeqs=count(istreambuf_iterator<char>(in),istreambuf_iterator<char>(), '>');
		in.seekg(0);
	
		candidateFileName = globaldata->getCandidateFile();
		string candidateAligngmentFName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align";
		ofstream candidateAlignmentFile;
		openOutputFile(candidateAligngmentFName, candidateAlignmentFile);

		string candidateReportFName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align.report";
		NastReport report(candidateReportFName);

		cout << "We are going to align the " << numFastaSeqs << " sequences in " << candidateFileName << "..." << endl;
		cout.flush();
	
		int start = time(NULL);

		for(int i=0;i<numFastaSeqs;i++){

			Sequence* candidateSeq = new Sequence(in);
			report.setCandidate(candidateSeq);

			Sequence* templateSeq = templateDB->findClosestSequence(candidateSeq);
			report.setTemplate(templateSeq);
			report.setSearchParameters(globaldata->getSearch(), templateDB->getSearchScore());
			
				
			Nast nast(alignment, candidateSeq, templateSeq);
			report.setAlignmentParameters(globaldata->getAlign(), alignment);
			report.setNastParameters(nast);

			candidateAlignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
			candidateAlignmentFile.flush();

			if((i+1) % 100 == 0){
				cout << "It has taken " << time(NULL) - start << " secs to align " << i+1 << " sequences" << endl;
			}
			report.print();
		
			delete candidateSeq;		
		}
		
		cout << "It took " << time(NULL) - start << " secs to align " << numFastaSeqs << " sequences" << endl;
		cout << endl;

		delete templateDB;
		delete alignment;

				
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************