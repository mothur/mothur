#ifndef TRIMOLIGOS_H
#define TRIMOLIGOS_H

/*
 *  trimoligos.h
 *  Mothur
 *
 *  Created by westcott on 9/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"
#include "qualityscores.h"


class TrimOligos {
	
	public:
		TrimOligos(int,int, map<string, int>, map<string, int>, vector<string>); //pdiffs, bdiffs, primers, barcodes, revPrimers
		~TrimOligos();
	
		int stripBarcode(Sequence&, int&);	
		int stripBarcode(Sequence&, QualityScores&, int&);
	
		int stripForward(Sequence&, int&);
		int stripForward(Sequence&, QualityScores&, int&);
	
		bool stripReverse(Sequence&);
		bool stripReverse(Sequence&, QualityScores&);
				
	
	private:
		int pdiffs, bdiffs;
	
		map<string, int> barcodes;
		map<string, int> primers;
		vector<string> revPrimer;
	
		MothurOut* m;
	
		bool compareDNASeq(string, string);				
		int countDiffs(string, string);			
};

#endif

