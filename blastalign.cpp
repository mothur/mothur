/*
 *  blastalign.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is a basic alignment method that gets the blast program to do the heavy lifting.  In the future, we should
 *	probably incorporate NCBI's library so that we don't have to call on a user-supplied executable.  This is a child
 *	of the Alignment class, which requires a constructor and align method.
 *
 */


#include "alignment.hpp"
#include "blastalign.hpp"


//**************************************************************************************************/

BlastAlignment::BlastAlignment(float go, float ge, float m, float mm) : 
			match(m),				//	This is the score to award for two nucleotides matching (match >= 0)
			mismatch(mm)			//	This is the penalty to assess for a mismatch (mismatch <= 0)
{
	globaldata = GlobalData::getInstance();
	path = globaldata->argv;
	path = path.substr(0, (path.find_last_of('m')));
	
	gapOpen = abs(go);				//	This is the penalty to assess for opening a gap (gapOpen >= 0)
	gapExtend = abs(ge);				//	This is the penalty to assess for extending a gap (gapExtend >= 0)
		
	int randNumber = rand();
	candidateFileName = toString(randNumber) + ".candidate";
	templateFileName = toString(randNumber) + ".template";
	blastFileName = toString(randNumber) + ".pairwise";
}

//**************************************************************************************************/

BlastAlignment::~BlastAlignment(){		//	The desctructor should clean up by removing the temporary 
	remove(candidateFileName.c_str());	//	files used to run bl2seq
	remove(templateFileName.c_str());
	remove(blastFileName.c_str());
}

//**************************************************************************************************/

void BlastAlignment::align(string seqA, string seqB){	//Use blastn to align the two sequences

	ofstream candidateFile(candidateFileName.c_str());	//	Write the sequence to be aligned to a temporary candidate seq file
	candidateFile << ">candidate" << endl << seqA << endl;
	candidateFile.close();
	
	ofstream templateFile(templateFileName.c_str());	//	Write the unaligned template sequence to a temporary candidate seq file
	templateFile << ">template" << endl << seqB << endl;
	templateFile.close();
	
	//	The blastCommand assumes that we have DNA sequences (blastn) and that they are fairly similar (-e 0.001) and
	//	that we don't want to apply any kind of complexity filtering (-F F)
	string blastCommand = path + "blast/bin/bl2seq -p blastn -i " + candidateFileName + " -j " + templateFileName + " -e 0.0001 -F F -o " + blastFileName + " -W 11";
	blastCommand += " -r " + toString(match) + " -q " + toString(mismatch);
	blastCommand +=	" -G " + toString(gapOpen) + " -E " + toString(gapExtend);
	
	system(blastCommand.c_str());	//	Here we assume that "bl2seq" is in the users path or in the same folder as
									//	this executable
	setPairwiseSeqs();
}

/**************************************************************************************************/

void BlastAlignment::setPairwiseSeqs(){	//	This method call assigns the blast generated alignment
															//	to the pairwise entry in the Sequence class for the 
															//	candidate and template Sequence objects
	ifstream blastFile;
	openInputFile(blastFileName, blastFile);
	
	seqAaln = "";
	seqBaln = "";
	
	int candidateLength, templateLength;
	char d;
	
	string candidateName, templateName;
	
	while(d=blastFile.get() != '='){}
	blastFile >> candidateName;					//	Get the candidate sequence name from flatfile
	
	while(d=blastFile.get() != '('){}
	blastFile >> candidateLength;				//	Get the candidate sequence length from flatfile
	
	while(d=blastFile.get()){
		if(d == '>'){
			blastFile >> templateName;			//	Get the template sequence name from flatfile
			break;
		}
		else if(d == '*'){									//	We go here if there is no significant match
			
			seqAstart = 0;
			seqBstart = 0;
			seqAend = 0;
			seqBend = 0;
			pairwiseLength = 0;
			
//			string dummy;
//			while(dummy != "query:"){	mothurOut(dummy, ""); mothurOutEndLine(); blastFile >> dummy;	}
//			blastFile >> seqBend;
//			mothurOut(toString(seqBend), ""); mothurOutEndLine();
//			for(int i=0;i<seqBend;i++){
//				seqAaln += 'Z';
//				seqBaln += 'X';
//			}
//			pairwiseLength = 0;
			return;
		}
	}
	
	while(d=blastFile.get() != '='){}
	blastFile >> templateLength;				//	Get the template sequence length from flatfile
		
	while(d=blastFile.get() != 'Q'){}			//	Suck up everything else until we get to the start of the alignment
	int queryStart, sbjctStart, queryEnd, sbjctEnd;
	string queryLabel, sbjctLabel, query, sbjct;

	blastFile >> queryLabel;	queryLabel = 'Q' + queryLabel;

	
	while(queryLabel == "Query:"){
		blastFile >> queryStart >> query >> queryEnd;
		
		while(d=blastFile.get() != 'S'){};
		
		blastFile >> sbjctLabel >> sbjctStart >> sbjct >> sbjctEnd;
		
		if(seqAaln == ""){
			seqAstart = queryStart;
			seqBstart = sbjctStart;
		}

		seqAaln += query;					//	concatenate each line of the sequence to what we already have
		seqBaln += sbjct;					//	for the query and template (subject) sequence
		
		blastFile >> queryLabel;
	}
	seqAend = queryEnd;
	seqBend = sbjctEnd;
	pairwiseLength = seqAaln.length();

	for(int i=1;i<seqBstart;i++){				//	Since the alignments don't always start at (1, 1), we need to pad
		seqAaln = 'Z' + seqAaln;				//	the sequences so that they start at the same point
		seqBaln = 'X' + seqBaln;
	}
	
	for(int i=seqBend+1;i<=templateLength;i++){	//	since the sequences don't necessarily end at the same point, we
		seqAaln += 'Z';							//	again need ot pad the sequences so that they extend to the length
		seqBaln += 'X';							//	of the template sequence
	}
	blastFile.close();
}

//**************************************************************************************************/
