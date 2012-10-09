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

struct oligosPair {
	string forward;
	string reverse;
	
	oligosPair() { forward = ""; reverse = "";  }
	oligosPair(string f, string r) : forward(f), reverse(r) {}
	~oligosPair() {}
};

class TrimOligos {
	
	public:
        TrimOligos(int,int, map<string, int>, map<string, int>, vector<string>); //pdiffs, bdiffs, primers, barcodes, revPrimers
        TrimOligos(int,int, int, int, map<string, int>, map<string, int>, vector<string>, vector<string>, vector<string>); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimers, linker, spacer
        TrimOligos(int,int, int, int, map<int, oligosPair>, map<int, oligosPair>, vector<string>, vector<string>); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, linker, spacer
		~TrimOligos();
	
		int stripBarcode(Sequence&, int&);	
		int stripBarcode(Sequence&, QualityScores&, int&);
        int stripBarcode(Sequence&, Sequence&, QualityScores&, QualityScores&, int&);
    	
		int stripForward(Sequence&, int&);
		int stripForward(Sequence&, QualityScores&, int&, bool);
        int stripForward(Sequence&, Sequence&, QualityScores&, QualityScores&, int&);
	
		bool stripReverse(Sequence&);
		bool stripReverse(Sequence&, QualityScores&);
    
        bool stripLinker(Sequence&);
        bool stripLinker(Sequence&, QualityScores&);
    
        bool stripSpacer(Sequence&);
        bool stripSpacer(Sequence&, QualityScores&);
				
	
	private:
		int pdiffs, bdiffs, ldiffs, sdiffs;
	
		map<string, int> barcodes;
		map<string, int> primers;
		vector<string> revPrimer;
        vector<string> linker;
        vector<string> spacer;
        map<int, oligosPair> ibarcodes;
        map<int, oligosPair> iprimers;
	
		MothurOut* m;
	
		bool compareDNASeq(string, string);				
		int countDiffs(string, string);			
};

#endif

