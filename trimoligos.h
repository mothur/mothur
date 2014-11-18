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
        TrimOligos(int,int, int, int, map<string, int>, map<string, int>, vector<string>, vector<string>, vector<string>); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimers, linker, spacer
        TrimOligos(int,int, int, int, map<int, oligosPair>, map<int, oligosPair>); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes
		~TrimOligos();
	
    
        //codes : 10 means sequence shorter than barcode, 100 means multiple matches, 1e6 no matches, 0 found match.
        //[0] = forward for paired, [1] = code for forward, [2] = reverse for paired, [3] = code for reverse
        //[0] = forward for single, [1] = code for forward
		vector<int> stripBarcode(Sequence&, int&);	
		vector<int> stripBarcode(Sequence&, QualityScores&, int&);
        vector<int> stripBarcode(Sequence&, Sequence&, QualityScores&, QualityScores&, int&);
        vector<int> stripBarcode(Sequence&, Sequence&, int&);
    	
		vector<int> stripForward(Sequence&, int&);
		vector<int> stripForward(Sequence&, QualityScores&, int&, bool);
        vector<int> stripForward(Sequence&, Sequence&, QualityScores&, QualityScores&, int&);
        vector<int> stripForward(Sequence&, Sequence&, int&);
	
		vector<int> stripReverse(Sequence&);
		vector<int> stripReverse(Sequence&, QualityScores&);
    
        bool stripLinker(Sequence&);
        bool stripLinker(Sequence&, QualityScores&);
    
        bool stripSpacer(Sequence&);
        bool stripSpacer(Sequence&, QualityScores&);
    
        //seq, primerStart, primerEnd
        bool findForward(Sequence&, int&, int&);
        bool findReverse(Sequence&, int&, int&);
    
        string reverseOligo(string);
        string getCodeValue(int, int);
	
	private:
		int pdiffs, bdiffs, ldiffs, sdiffs;
        bool paired;
	
		map<string, int> barcodes;
		map<string, int> primers;
		vector<string> revPrimer;
        vector<string> linker;
        vector<string> spacer;
        map<string, vector<int> > ifbarcodes;
        map<string, vector<int> > ifprimers;
        map<string, vector<int> > irbarcodes;
        map<string, vector<int> > irprimers;
        map<int, oligosPair> ipbarcodes;
        map<int, oligosPair> ipprimers;
    
        int maxFBarcodeLength, maxRBarcodeLength, maxFPrimerLength, maxRPrimerLength, maxLinkerLength, maxSpacerLength;
	
		MothurOut* m;
	
		bool compareDNASeq(string, string);				
		int countDiffs(string, string);
        
        vector<int> stripPairedBarcode(Sequence& seq, QualityScores& qual, int& group);
        vector<int> stripPairedPrimers(Sequence& seq, QualityScores& qual, int& group, bool);
        vector<int> stripPairedBarcode(Sequence& seq,  int& group);
        vector<int> stripPairedPrimers(Sequence& seq,  int& group);

};

#endif

