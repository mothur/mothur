/*
 *  trimoligos.cpp
 *  Mothur
 *
 *  Created by westcott on 9/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "trimoligos.h"
#include "alignment.hpp"
#include "needlemanoverlap.hpp"


/********************************************************************/
//strip, pdiffs, bdiffs, primers, barcodes, revPrimers
TrimOligos::TrimOligos(int p, int b, int l, int s, map<string, int> pr, map<string, int> br, vector<string> r, vector<string> lk, vector<string> sp){
	try {
		m = MothurOut::getInstance();
        paired = false;
        trashCode = "";
		
		pdiffs = p;
		bdiffs = b;
        ldiffs = l;
        sdiffs = s;
		
		barcodes = br;
		primers = pr;
		revPrimer = r;
        linker = lk;
        spacer = sp;
        maxFBarcodeLength = 0;
        for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
            if(it->first.length() > maxFBarcodeLength){
                maxFBarcodeLength = it->first.length();
            }
        }
        maxFPrimerLength = 0;
        map<string,int>::iterator it; 
        for(it=primers.begin();it!=primers.end();it++){
            if(it->first.length() > maxFPrimerLength){
                maxFPrimerLength = it->first.length();
            }
        }
        
        maxLinkerLength = 0;
        for(int i = 0; i < linker.size(); i++){
            if(linker[i].length() > maxLinkerLength){
                maxLinkerLength = linker[i].length();
            }
        }
        
        maxSpacerLength = 0;
        for(int i = 0; i < spacer.size(); i++){
            if(spacer[i].length() > maxSpacerLength){
                maxSpacerLength = spacer[i].length();
            }
        }
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "TrimOligos");
		exit(1);
	}
}
/********************************************************************/
//strip, pdiffs, bdiffs, primers, barcodes, revPrimers
TrimOligos::TrimOligos(int p, int b, int l, int s, map<int, oligosPair> pr, map<int, oligosPair> br){
	try {
		m = MothurOut::getInstance();
		
		pdiffs = p;
		bdiffs = b;
        ldiffs = l;
        sdiffs = s;
        paired = true;
        trashCode = "";
		
        maxFBarcodeLength = 0;
        for(map<int,oligosPair>::iterator it=br.begin();it!=br.end();it++){
            string forward = it->second.forward;
            map<string, vector<int> >::iterator itForward = ifbarcodes.find(forward);
            
            if(forward.length() > maxFBarcodeLength){ maxFBarcodeLength = forward.length();  }
            
            if (itForward == ifbarcodes.end()) {
                vector<int> temp; temp.push_back(it->first);
                ifbarcodes[forward] = temp;
            }else { itForward->second.push_back(it->first); }
        }

        maxFPrimerLength = 0;
        for(map<int,oligosPair>::iterator it=pr.begin();it!=pr.end();it++){
            string forward = it->second.forward;
            map<string, vector<int> >::iterator itForward = ifprimers.find(forward);
            
            if(forward.length() > maxFPrimerLength){ maxFPrimerLength = forward.length();  }
            
            if (itForward == ifprimers.end()) {
                vector<int> temp; temp.push_back(it->first);
                ifprimers[forward] = temp;
            }else { itForward->second.push_back(it->first); }
        }
        
        maxRBarcodeLength = 0;
        for(map<int,oligosPair>::iterator it=br.begin();it!=br.end();it++){
            string forward = it->second.reverse;
            map<string, vector<int> >::iterator itForward = irbarcodes.find(forward);
            
            if(forward.length() > maxRBarcodeLength){ maxRBarcodeLength = forward.length();  }
            
            if (itForward == irbarcodes.end()) {
                vector<int> temp; temp.push_back(it->first);
                irbarcodes[forward] = temp;
            }else { itForward->second.push_back(it->first); }
        }
        
        maxRPrimerLength = 0;
        for(map<int,oligosPair>::iterator it=pr.begin();it!=pr.end();it++){
            string forward = it->second.reverse;
            map<string, vector<int> >::iterator itForward = irprimers.find(forward);
            
            if(forward.length() > maxRPrimerLength){ maxRPrimerLength = forward.length();  }
            
            if (itForward == irprimers.end()) {
                vector<int> temp; temp.push_back(it->first);
                irprimers[forward] = temp;
            }else { itForward->second.push_back(it->first); }
        }

        ipbarcodes = br;
        ipprimers = pr;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "TrimOligos");
		exit(1);
	}
}
/********************************************************************/
//strip, pdiffs, bdiffs, primers, barcodes, revPrimers
TrimOligos::TrimOligos(int p, int b, map<string, int> pr, map<string, int> br, vector<string> r){
	try {
		m = MothurOut::getInstance();
		
		pdiffs = p;
		bdiffs = b;
		
		barcodes = br;
		primers = pr;
		revPrimer = r;
        paired = false;
        trashCode = "";
        
        maxFBarcodeLength = 0;
        for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
            string oligo = it->first;
            if(oligo.length() > maxFBarcodeLength){
                maxFBarcodeLength = oligo.length();
            }
        }
        maxRBarcodeLength = maxFBarcodeLength;
        
        maxFPrimerLength = 0;
        for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
            string oligo = it->first;
            if(oligo.length() > maxFPrimerLength){
                maxFPrimerLength = oligo.length();
            }
        }
        maxRPrimerLength = maxFPrimerLength;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "TrimOligos");
		exit(1);
	}
}
/********************************************************************/
TrimOligos::~TrimOligos() {}
//********************************************************************/
bool TrimOligos::findForward(Sequence& seq, int& primerStart, int& primerEnd){
	try {
		
		string rawSequence = seq.getUnaligned();
        trashCode = "";
		
		for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
			string oligo = it->first;
			
			if(rawSequence.length() < oligo.length()) {  break;  }
			
			//search for primer
            int olength = oligo.length();
            for (int j = 0; j < rawSequence.length()-olength; j++){
                if (m->control_pressed) {  primerStart = 0; primerEnd = 0; return false; }
                string rawChunk = rawSequence.substr(j, olength);
                if(compareDNASeq(oligo, rawChunk)) {
                    primerStart = j;
                    primerEnd = primerStart + olength;
                    return true;
                }
                
            }
        }	
		
        primerStart = 0; primerEnd = 0;
        //if you don't want to allow for diffs
		if (pdiffs == 0) { trashCode = "f"; return false;  }
		else { //try aligning and see if you can find it
						
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			
            Alignment* alignment;
            if (primers.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));  }
            else{ alignment = NULL; } 
            
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
                
                string prim = it->first;
                //search for primer
                int olength = prim.length();
                if (rawSequence.length() < olength) {} //ignore primers too long for this seq
                else{
                    for (int j = 0; j < rawSequence.length()-olength; j++){
                        
                        string oligo = it->first;

                        if (m->control_pressed) {  primerStart = 0; primerEnd = 0; return false; }
				
                        string rawChunk = rawSequence.substr(j, olength+pdiffs);
                        
                        //use needleman to align first primer.length()+numdiffs of sequence to each barcode
                        alignment->alignPrimer(oligo, rawChunk);
                        oligo = alignment->getSeqAAln();
                        string temp = alignment->getSeqBAln();
				
                        int alnLength = oligo.length();
				
                        for(int i=oligo.length()-1;i>=0;i--){
                            if(oligo[i] != '-'){	alnLength = i+1;	break;	}
                        }
                        oligo = oligo.substr(0,alnLength);
                        temp = temp.substr(0,alnLength);
				
                        int numDiff = countDiffs(oligo, temp);
				
                        if(numDiff < minDiff){
                            minDiff = numDiff;
                            minCount = 1;
                            primerStart = j;
                            primerEnd = primerStart + alnLength;
                        }else if(numDiff == minDiff){ minCount++; }
                    }
                }
			}
			
            if (alignment != NULL) {  delete alignment;  }
            
			if(minDiff > pdiffs)	{	primerStart = 0; primerEnd = 0; trashCode = "f"; return false;		}	//no good matches
			else if(minCount > 1)	{	primerStart = 0; primerEnd = 0; trashCode = "f"; return false;	}	//can't tell the difference between multiple primers
			else{  trashCode = ""; return true; }
		}
        
        trashCode = "f";
        primerStart = 0; primerEnd = 0;
		return false;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripForward");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::findReverse(Sequence& seq, int& primerStart, int& primerEnd){
	try {
		trashCode = "";
		string rawSequence = seq.getUnaligned();
		
		for(int i=0;i<revPrimer.size();i++){
			string oligo = revPrimer[i];
			if(rawSequence.length() < oligo.length()) {  break;  }
			
			//search for primer
            int olength = oligo.length();
            for (int j = rawSequence.length()-olength; j >= 0; j--){
                if (m->control_pressed) {  primerStart = 0; primerEnd = 0; return false; }
                string rawChunk = rawSequence.substr(j, olength);
                
                if(compareDNASeq(oligo, rawChunk)) {
                    primerStart = j;
                    primerEnd = primerStart + olength;
                    return true;
                }
                
            }
		}	
		
        trashCode = "r";
        primerStart = 0; primerEnd = 0;
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "findReverse");
		exit(1);
	}
}
//*******************************************************************/

int TrimOligos::stripBarcode(Sequence& seq, QualityScores& qual, int& group){
	try {
		trashCode = "";
        if (paired) {  int success = stripPairedBarcode(seq, qual, group);  return success; }
        
		string rawSequence = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the barcode
		for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				
				success = 0;
                trashCode = "";
				return success;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if (bdiffs == 0) { trashCode = "b"; return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (barcodes.size() > 0) {alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFBarcodeLength+bdiffs+1));   }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxFBarcodeLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10; trashCode = "b";
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > bdiffs)	{	trashCode = "b"; success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	trashCode = "b"; success = bdiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				group = minGroup;
				seq.setUnaligned(rawSequence.substr(minPos));
    
				if(qual.getName() != ""){
					qual.trimQScores(minPos, -1);
				}
				success = minDiff;
                trashCode = "";
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripBarcode");
		exit(1);
	}
}
//*******************************************************************/
int TrimOligos::stripBarcode(Sequence& forwardSeq, Sequence& reverseSeq, int& group){
	try {
        trashCode = "";
		//look for forward barcode
		string rawFSequence = forwardSeq.getUnaligned();
        string rawRSequence = reverseSeq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the forward barcode
		for(map<int,oligosPair>::iterator it=ipbarcodes.begin();it!=ipbarcodes.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            
			if(rawFSequence.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode += "b";
                break;
			}
            if(rawRSequence.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode += "d";
				break;	
			}
			
			if(compareDNASeq(foligo, rawFSequence.substr(0,foligo.length()))) {
                if (compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length()))) {
                    group = it->first;
                    forwardSeq.setUnaligned(rawFSequence.substr(foligo.length()));
                    reverseSeq.setUnaligned(rawRSequence.substr(0,(rawRSequence.length()-roligo.length())));
                    success = 0; trashCode = "";
                    return success;
                }
            }else {
                trashCode = "b";
                if (!compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length()))) { trashCode += "d";  }
            }
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if (bdiffs == 0) {  return success;  }
		else { //try aligning and see if you can find it
            trashCode = "";
			Alignment* alignment;
			if (ifbarcodes.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFBarcodeLength+bdiffs+1));   }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos; 
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
                1   Sarah   Westcott
                2   John    Westcott
                3   Anna    Westcott
                4   Sarah   Schloss
                5   Pat     Schloss
                6   Gail    Brown
                7   Pat     Moore
             
                only want to look for forward = Sarah, John, Anna, Pat, Gail
                                      reverse = Westcott, Schloss, Brown, Moore
             
                but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;         
			for(map<string, vector<int> >::iterator it=ifbarcodes.begin();it!=ifbarcodes.end();it++){
				string oligo = it->first;
				
				if(rawFSequence.length() < maxFBarcodeLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawFSequence.substr(0,oligo.length()+bdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawFSequence.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
               
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = bdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > bdiffs)	{	trashCode = "b"; minFGroup.clear(); success = minDiff;		}	//no good matches
			//else{
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irbarcodes.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRBarcodeLength+bdiffs+1));   }
                else{ alignment = NULL; } 
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                for(map<string, vector<int> >::iterator it=irbarcodes.begin();it!=irbarcodes.end();it++){
                    string oligo = it->first;
                    //cout << "before = " << oligo << '\t' << rawRSequence.substr(0,oligo.length()+bdiffs) << endl;
                    if(rawRSequence.length() < maxRBarcodeLength){	//let's just assume that the barcodes are the same length
                        success = bdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+bdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = bdiffs + 100; }
                    
                    //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);                    
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);  
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                for (int i = 0; i < minFGroup.size(); i++) { 
                    cout << i << '\t';
                    for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                    cout << endl;
                }
                cout << endl;
                for (int i = 0; i < minRGroup.size(); i++) { 
                    cout << i << '\t';
                    for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                    cout << endl;
                }
                cout << endl;*/
                if(minDiff > bdiffs)	{	trashCode += "d"; success = minDiff;		}	//no good matches
                else  {  
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) { 
                        foundMatch = true;
                        group = matches[0]; 
                        fStart = minFPos[0];
                        rStart = minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        forwardSeq.setUnaligned(rawFSequence.substr(fStart));
                        reverseSeq.setUnaligned(rawRSequence.substr(rStart));
                        success = minDiff;
                    }else { if (trashCode == "") { trashCode = "bd"; } success = bdiffs + 100;	}
                }
			//}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
        //catchall for case I didn't think of
        if ((trashCode == "") && (success > bdiffs)) { trashCode = "bd"; }
        
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripIBarcode");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripBarcode(Sequence& forwardSeq, Sequence& reverseSeq, QualityScores& forwardQual, QualityScores& reverseQual, int& group){
	try {
        trashCode = "";
        
		//look for forward barcode
		string rawFSequence = forwardSeq.getUnaligned();
        string rawRSequence = reverseSeq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the forward barcode
		for(map<int,oligosPair>::iterator it=ipbarcodes.begin();it!=ipbarcodes.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            
			if(rawFSequence.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode = "b";
				break;	
			}
            if(rawRSequence.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode += "d";
				break;	
			}
			
			if(compareDNASeq(foligo, rawFSequence.substr(0,foligo.length()))) {
                if (compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length()))) {
                    group = it->first;
                    forwardSeq.setUnaligned(rawFSequence.substr(foligo.length()));
                    reverseSeq.setUnaligned(rawRSequence.substr(0,(rawRSequence.length()-roligo.length())));
                    forwardQual.trimQScores(foligo.length(), -1);
                    reverseQual.trimQScores(-1, rawRSequence.length()-roligo.length());
                    success = 0; trashCode = "";
                    return success;
                }
            }else {
                trashCode = "b";
                if (!compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length()))) {  trashCode += "d";
                }
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if (bdiffs == 0) { return success;  }
		else { //try aligning and see if you can find it
            trashCode = "";
			Alignment* alignment;
			if (ifbarcodes.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFBarcodeLength+bdiffs+1));   }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos; 
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
             1   Sarah   Westcott
             2   John    Westcott
             3   Anna    Westcott
             4   Sarah   Schloss
             5   Pat     Schloss
             6   Gail    Brown
             7   Pat     Moore
             
             only want to look for forward = Sarah, John, Anna, Pat, Gail
             reverse = Westcott, Schloss, Brown, Moore
             
             but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;         
			for(map<string, vector<int> >::iterator it=ifbarcodes.begin();it!=ifbarcodes.end();it++){
				string oligo = it->first;
				
				if(rawFSequence.length() < maxFBarcodeLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawFSequence.substr(0,oligo.length()+bdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawFSequence.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = bdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > bdiffs)	{	trashCode = "b"; minFGroup.clear(); success = minDiff;		}	//no good matches
			//else{
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irbarcodes.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRBarcodeLength+bdiffs+1));   }
                else{ alignment = NULL; } 
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                for(map<string, vector<int> >::iterator it=irbarcodes.begin();it!=irbarcodes.end();it++){
                    string oligo = it->first;
                    //cout << "before = " << oligo << '\t' << rawRSequence.substr(0,oligo.length()+bdiffs) << endl;
                    if(rawRSequence.length() < maxRBarcodeLength){	//let's just assume that the barcodes are the same length
                        success = bdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+bdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = bdiffs + 100; }
                    
                    //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);                    
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);  
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                 for (int i = 0; i < minFGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;
                 for (int i = 0; i < minRGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;*/
                if(minDiff > bdiffs)	{	trashCode = "d"; success = minDiff;		}	//no good matches
                else  {  
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) { 
                        foundMatch = true;
                        group = matches[0]; 
                        fStart = minFPos[0];
                        rStart = minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        forwardSeq.setUnaligned(rawFSequence.substr(fStart));
                        reverseSeq.setUnaligned(rawRSequence.substr(rStart));
                        forwardQual.trimQScores(fStart, -1);
                        reverseQual.trimQScores(rStart, -1);
                        success = minDiff;
                    }else { if (trashCode == "") { trashCode = "bd"; } success = bdiffs + 100;	}
                }
			//}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
        //catchall for case I didn't think of
        if ((trashCode == "") && (success > bdiffs)) { trashCode = "bd"; }
        
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripIBarcode");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripPairedBarcode(Sequence& seq, QualityScores& qual, int& group){
	try {
        trashCode = "";
        
		//look for forward barcode
		string rawSeq = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		//cout << seq.getName() << endl;
		//can you find the forward barcode
		for(map<int,oligosPair>::iterator it=ipbarcodes.begin();it!=ipbarcodes.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            //cout << it->first << '\t' << foligo << '\t' << roligo << endl;
            //cout << it->first << '\t' << rawSeq.substr(0,foligo.length()) << '\t' << rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()) << endl;
			if(rawSeq.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode = "b";
				break;
			}
            if(rawSeq.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode += "d";
				break;
			}
			
			if(compareDNASeq(foligo, rawSeq.substr(0,foligo.length()))) { //find forward
                if (compareDNASeq(roligo, rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()))) { //find reverse
                    group = it->first;
                    string trimmedSeq = rawSeq.substr(foligo.length()); //trim forward barcode
                    seq.setUnaligned(trimmedSeq.substr(0,(trimmedSeq.length()-roligo.length()))); //trim reverse barcode
                    if(qual.getName() != ""){
                        qual.trimQScores(-1, rawSeq.length()-roligo.length());
                        qual.trimQScores(foligo.length(), -1);
                    }
                    success = 0;
                    trashCode = "";
                    return success;
                }
            }else {
                trashCode = "b";
				if (!compareDNASeq(roligo, rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()))) { trashCode += "d"; }
			}
		}
		//cout << "success=" << success << endl;
		//if you found the barcode or if you don't want to allow for diffs
		if (bdiffs == 0) { return success;  }
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (ifbarcodes.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFBarcodeLength+bdiffs+1));   }
			else{ alignment = NULL; }
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos;
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
             1   Sarah   Westcott
             2   John    Westcott
             3   Anna    Westcott
             4   Sarah   Schloss
             5   Pat     Schloss
             6   Gail    Brown
             7   Pat     Moore
             
             only want to look for forward = Sarah, John, Anna, Pat, Gail
             reverse = Westcott, Schloss, Brown, Moore
             
             but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;
			for(map<string, vector<int> >::iterator it=ifbarcodes.begin();it!=ifbarcodes.end();it++){
				string oligo = it->first;
				
				if(rawSeq.length() < maxFBarcodeLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawSeq.substr(0,oligo.length()+bdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSeq.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = bdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > bdiffs)	{	trashCode = "b"; minFGroup.clear(); success = minDiff;		}	//no good matches
			//else{
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irbarcodes.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRBarcodeLength+bdiffs+1));   }
                else{ alignment = NULL; }
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                string rawRSequence = reverseOligo(seq.getUnaligned());
                //cout << irbarcodes.size() << '\t' << maxRBarcodeLength <<  endl;
                for(map<string, vector<int> >::iterator it=irbarcodes.begin();it!=irbarcodes.end();it++){
                    string oligo = reverseOligo(it->first);
                    //cout << "r before = " << reverseOligo(oligo) << '\t' << reverseOligo(rawRSequence.substr(0,oligo.length()+bdiffs)) << endl;
                    if(rawRSequence.length() < maxRBarcodeLength){	//let's just assume that the barcodes are the same length
                        success = bdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+bdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = bdiffs + 100; }
                    
                    //cout << "r after = " << reverseOligo(oligo) << '\t' << reverseOligo(temp) << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                 for (int i = 0; i < minFGroup.size(); i++) {
                 cout << i << '\t';
                 for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;
                 for (int i = 0; i < minRGroup.size(); i++) {
                 cout << i << '\t';
                 for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;*/
                if(minDiff > bdiffs)	{	trashCode += "d"; success = minDiff;		}	//no good matches
                else  {
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) {
                        foundMatch = true;
                        group = matches[0];
                        fStart = minFPos[0];
                        rStart = rawSeq.length() - minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        string trimmedSeq = rawSeq.substr(0, rStart); //trim reverse barcode
                        seq.setUnaligned(trimmedSeq.substr(fStart)); //trim forward barcode
                        if(qual.getName() != ""){
                            qual.trimQScores(-1, rStart);
                            qual.trimQScores(fStart, -1);
                        }
                        success = minDiff;
                        //cout << "barcode = " << ipbarcodes[group].forward << '\t' << ipbarcodes[group].reverse << endl;
                    }else { if (trashCode == "") { trashCode = "bd"; } success = bdiffs + 100;	}
                }
			//}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
        //catchall for case I didn't think of
        if ((trashCode == "") && (success > bdiffs)) { trashCode = "bd"; }

		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripPairedBarcode");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripPairedPrimers(Sequence& seq, QualityScores& qual, int& group, bool keepForward){
	try {
        trashCode = "";
		//look for forward 
		string rawSeq = seq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
        //cout << seq.getName() << endl;
		//can you find the forward
		for(map<int,oligosPair>::iterator it=ipprimers.begin();it!=ipprimers.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            
            //cout << it->first << '\t' << foligo << '\t' << roligo << endl;
            //cout << it->first << '\t' << rawSeq.substr(0,foligo.length()) << '\t' << rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()) << endl;
			if(rawSeq.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode = "p";
				break;
			}
            if(rawSeq.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
                trashCode += "q";
				break;
			}
			
            if(compareDNASeq(foligo, rawSeq.substr(0,foligo.length()))) { //find forward
                if (compareDNASeq(roligo, rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()))) { //find reverse
                    group = it->first;
                    string trimmedSeq = rawSeq.substr(foligo.length()); //trim forward barcode
                    seq.setUnaligned(trimmedSeq.substr(0,(trimmedSeq.length()-roligo.length()))); //trim reverse barcode
                    if(qual.getName() != ""){
                        qual.trimQScores(-1, rawSeq.length()-roligo.length());
                        qual.trimQScores(foligo.length(), -1);
                    }
                    success = 0;
                    trashCode = "";
                    return success;
                }
            }else {
                trashCode = "b";
				if (!compareDNASeq(roligo, rawSeq.substr(rawSeq.length()-roligo.length(),roligo.length()))) { trashCode += "d"; }
			}

        }
		//cout << "success=" << success << endl;
		//if you found the barcode or if you don't want to allow for diffs
		if (pdiffs == 0) { return success;  }
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (ifprimers.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));   }
			else{ alignment = NULL; }
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos;
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
             1   Sarah   Westcott
             2   John    Westcott
             3   Anna    Westcott
             4   Sarah   Schloss
             5   Pat     Schloss
             6   Gail    Brown
             7   Pat     Moore
             
             only want to look for forward = Sarah, John, Anna, Pat, Gail
             reverse = Westcott, Schloss, Brown, Moore
             
             but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;
			for(map<string, vector<int> >::iterator it=ifprimers.begin();it!=ifprimers.end();it++){
				string oligo = it->first;
				
				if(rawSeq.length() < maxFPrimerLength){	//let's just assume that the barcodes are the same length
					success = pdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawSeq.substr(0,oligo.length()+bdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSeq.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = pdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > pdiffs)	{	trashCode = "p"; minFGroup.clear(); success = minDiff;		}	//no good matches
			//else{
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irprimers.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRPrimerLength+pdiffs+1));   }
                else{ alignment = NULL; }
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                string rawRSequence = reverseOligo(seq.getUnaligned());
                
                for(map<string, vector<int> >::iterator it=irprimers.begin();it!=irprimers.end();it++){
                    string oligo = reverseOligo(it->first);
                    //cout << "r before = " << reverseOligo(oligo) << '\t' << reverseOligo(rawRSequence.substr(0,oligo.length()+pdiffs)) << endl;
                    if(rawRSequence.length() < maxRPrimerLength){	//let's just assume that the barcodes are the same length
                        success = pdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+pdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = pdiffs + 100; }
                    
                    //cout << "r after = " << reverseOligo(oligo) << '\t' << reverseOligo(temp) << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                 for (int i = 0; i < minFGroup.size(); i++) {
                 cout << i << '\t';
                 for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;
                 for (int i = 0; i < minRGroup.size(); i++) {
                 cout << i << '\t';
                 for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;*/
                if(minDiff > pdiffs)	{ trashCode += "q"; 	success = minDiff;		}	//no good matches
                else  {
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) {
                        foundMatch = true;
                        group = matches[0];
                        fStart = minFPos[0];
                        rStart = rawSeq.length() - minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        if (!keepForward) {
                            string trimmedSeq = rawSeq.substr(0, rStart); //trim reverse barcode
                            seq.setUnaligned(trimmedSeq.substr(fStart)); //trim forward barcode
                            if(qual.getName() != ""){
                                qual.trimQScores(-1, rStart);
                                qual.trimQScores(fStart, -1);
                            }
                        }
                        success = minDiff;
                        //cout << "primer = " << ipprimers[group].forward << '\t' << ipprimers[group].reverse << endl;
                    }else { if (trashCode == "") { trashCode = "pq"; } success = pdiffs + 100;	}
                }
			//}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
        //catchall for case I didn't think of
        if ((trashCode == "") && (success > pdiffs)) { trashCode = "pq"; }
        
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripPairedPrimers");
		exit(1);
	}
	
}

//*******************************************************************/
int TrimOligos::stripForward(Sequence& forwardSeq, Sequence& reverseSeq, QualityScores& forwardQual, QualityScores& reverseQual, int& group){
	try {
		//look for forward barcode
		string rawFSequence = forwardSeq.getUnaligned();
        string rawRSequence = reverseSeq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
		
		//can you find the forward barcode
		for(map<int,oligosPair>::iterator it=ipprimers.begin();it!=ipprimers.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            
			if(rawFSequence.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
            if(rawRSequence.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if((compareDNASeq(foligo, rawFSequence.substr(0,foligo.length()))) && (compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length())))) {
				group = it->first;
				forwardSeq.setUnaligned(rawFSequence.substr(foligo.length()));
                reverseSeq.setUnaligned(rawRSequence.substr(0,(rawRSequence.length()-roligo.length())));
                forwardQual.trimQScores(foligo.length(), -1);
                reverseQual.trimQScores(-1, rawRSequence.length()-roligo.length());
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) { return success;  }
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (ifprimers.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));   }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos; 
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
             1   Sarah   Westcott
             2   John    Westcott
             3   Anna    Westcott
             4   Sarah   Schloss
             5   Pat     Schloss
             6   Gail    Brown
             7   Pat     Moore
             
             only want to look for forward = Sarah, John, Anna, Pat, Gail
             reverse = Westcott, Schloss, Brown, Moore
             
             but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;         
			for(map<string, vector<int> >::iterator it=ifprimers.begin();it!=ifprimers.end();it++){
				string oligo = it->first;
				
				if(rawFSequence.length() < maxFPrimerLength){	//let's just assume that the barcodes are the same length
					success = pdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawFSequence.substr(0,oligo.length()+pdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawFSequence.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = pdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
			else{	
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irbarcodes.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRPrimerLength+pdiffs+1));   }
                else{ alignment = NULL; } 
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                for(map<string, vector<int> >::iterator it=irprimers.begin();it!=irprimers.end();it++){
                    string oligo = it->first;
                    //cout << "before = " << oligo << '\t' << rawRSequence.substr(0,oligo.length()+pdiffs) << endl;
                    if(rawRSequence.length() < maxRPrimerLength){	//let's just assume that the barcodes are the same length
                        success = pdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+pdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = pdiffs + 100; }
                    
                    //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);                    
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);  
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                 for (int i = 0; i < minFGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;
                 for (int i = 0; i < minRGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;*/
                if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
                else  {  
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) { 
                        foundMatch = true;
                        group = matches[0]; 
                        fStart = minFPos[0];
                        rStart = minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        forwardSeq.setUnaligned(rawFSequence.substr(fStart));
                        reverseSeq.setUnaligned(rawRSequence.substr(rStart));
                        forwardQual.trimQScores(fStart, -1);
                        reverseQual.trimQScores(rStart, -1);
                        success = minDiff;
                    }else { success = pdiffs + 100;	}
                }
			}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripIForward");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripForward(Sequence& forwardSeq, Sequence& reverseSeq, int& group){
	try {
		//look for forward barcode
		string rawFSequence = forwardSeq.getUnaligned();
        string rawRSequence = reverseSeq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
		
		//can you find the forward barcode
		for(map<int,oligosPair>::iterator it=ipprimers.begin();it!=ipprimers.end();it++){
			string foligo = it->second.forward;
            string roligo = it->second.reverse;
            
			if(rawFSequence.length() < foligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
            if(rawRSequence.length() < roligo.length()){	//let's just assume that the barcodes are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if((compareDNASeq(foligo, rawFSequence.substr(0,foligo.length()))) && (compareDNASeq(roligo, rawRSequence.substr((rawRSequence.length()-roligo.length()),roligo.length())))) {
				group = it->first;
				forwardSeq.setUnaligned(rawFSequence.substr(foligo.length()));
                reverseSeq.setUnaligned(rawRSequence.substr(0,(rawRSequence.length()-roligo.length())));
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) { return success;  }
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (ifprimers.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));   }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			vector< vector<int> > minFGroup;
			vector<int> minFPos; 
            
            //the pair can have duplicates, but we only want to search for the unique forward and reverses, example
            /*
             1   Sarah   Westcott
             2   John    Westcott
             3   Anna    Westcott
             4   Sarah   Schloss
             5   Pat     Schloss
             6   Gail    Brown
             7   Pat     Moore
             
             only want to look for forward = Sarah, John, Anna, Pat, Gail
             reverse = Westcott, Schloss, Brown, Moore
             
             but if best match forward = 4, and reverse = 1, we want to count as a valid match because forward 1 and forward 4 are the same.  so both barcodes map to same group.
             */
            //cout << endl << forwardSeq.getName() << endl;         
			for(map<string, vector<int> >::iterator it=ifprimers.begin();it!=ifprimers.end();it++){
				string oligo = it->first;
				
				if(rawFSequence.length() < maxFPrimerLength){	//let's just assume that the barcodes are the same length
					success = pdiffs + 10;
					break;
				}
				//cout << "before = " << oligo << '\t' << rawFSequence.substr(0,oligo.length()+pdiffs) << endl;
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawFSequence.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
                int numDiff = countDiffs(oligo, temp);
				
                if (alnLength == 0) { numDiff = pdiffs + 100; }
                //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minFGroup.clear();
                    minFGroup.push_back(it->second);
					int tempminFPos = 0;
                    minFPos.clear();
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}else if(numDiff == minDiff){
					minFGroup.push_back(it->second);
                    int tempminFPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							tempminFPos++;
						}
					}
                    minFPos.push_back(tempminFPos);
				}
			}
            
			//cout << minDiff << '\t' << minCount << '\t' << endl;
			if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
			else{	
                //check for reverse match
                if (alignment != NULL) {  delete alignment;  }
                
                if (irbarcodes.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxRPrimerLength+pdiffs+1));   }
                else{ alignment = NULL; } 
                
                //can you find the barcode
                minDiff = 1e6;
                minCount = 1;
                vector< vector<int> > minRGroup;
                vector<int> minRPos;
                
                for(map<string, vector<int> >::iterator it=irprimers.begin();it!=irprimers.end();it++){
                    string oligo = it->first;
                    //cout << "before = " << oligo << '\t' << rawRSequence.substr(0,oligo.length()+pdiffs) << endl;
                    if(rawRSequence.length() < maxRPrimerLength){	//let's just assume that the barcodes are the same length
                        success = pdiffs + 10;
                        break;
                    }
                    
                    //use needleman to align first barcode.length()+numdiffs of sequence to each barcode
                    alignment->alignPrimer(oligo, rawRSequence.substr(0,oligo.length()+pdiffs));
                    oligo = alignment->getSeqAAln();
                    string temp = alignment->getSeqBAln();
                    
                    int alnLength = oligo.length();
                    for(int i=oligo.length()-1;i>=0;i--){ if(oligo[i] != '-'){	alnLength = i+1;	break;	} }
                    oligo = oligo.substr(0,alnLength);
                    temp = temp.substr(0,alnLength);
                    int numDiff = countDiffs(oligo, temp);
                    if (alnLength == 0) { numDiff = pdiffs + 100; }
                    
                    //cout << "after = " << oligo << '\t' << temp << '\t' << numDiff << endl;
                    if(numDiff < minDiff){
                        minDiff = numDiff;
                        minCount = 1;
                        minRGroup.clear();
                        minRGroup.push_back(it->second);
                        int tempminRPos = 0;
                        minRPos.clear();
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);                    
                    }else if(numDiff == minDiff){
                        int tempminRPos = 0;
                        for(int i=0;i<alnLength;i++){
                            if(temp[i] != '-'){
                                tempminRPos++;
                            }
                        }
                        minRPos.push_back(tempminRPos);  
                        minRGroup.push_back(it->second);
                    }
                    
                }
                
                /*cout << minDiff << '\t' << minCount << '\t' << endl;
                 for (int i = 0; i < minFGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minFGroup[i].size(); j++) {  cout << minFGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;
                 for (int i = 0; i < minRGroup.size(); i++) { 
                 cout << i << '\t';
                 for (int j = 0; j < minRGroup[i].size(); j++) {  cout << minRGroup[i][j] << " "; }
                 cout << endl;
                 }
                 cout << endl;*/
                if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
                else  {  
                    bool foundMatch = false;
                    vector<int> matches;
                    for (int i = 0; i < minFGroup.size(); i++) {
                        for (int j = 0; j < minFGroup[i].size(); j++) {
                            for (int k = 0; k < minRGroup.size(); k++) {
                                if (m->inUsersGroups(minFGroup[i][j], minRGroup[k])) { matches.push_back(minFGroup[i][j]); k+= minRGroup.size(); }
                            }
                        }
                    }
                    
                    int fStart = 0;
                    int rStart = 0;
                    if (matches.size() == 1) { 
                        foundMatch = true;
                        group = matches[0]; 
                        fStart = minFPos[0];
                        rStart = minRPos[0];
                    }
                    
                    //we have an acceptable match for the forward and reverse, but do they match?
                    if (foundMatch) {
                        forwardSeq.setUnaligned(rawFSequence.substr(fStart));
                        reverseSeq.setUnaligned(rawRSequence.substr(rStart));
                        success = minDiff;
                    }else { success = pdiffs + 100;	}
                }
			}
			
			if (alignment != NULL) {  delete alignment;  }
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripIForward");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripBarcode(Sequence& seq, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the barcode
		for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((bdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
            Alignment* alignment;
			if (barcodes.size() > 0) {  alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFBarcodeLength+bdiffs+1));  }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxFBarcodeLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				 
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				 
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > bdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = bdiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				group = minGroup;
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripBarcode");
		exit(1);
	}
	
}

/********************************************************************/
int TrimOligos::stripForward(Sequence& seq, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
		
		//can you find the primer
		for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the primers are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) {	return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (primers.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));   }
            else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxFPrimerLength){	
					success = pdiffs + 100;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = pdiffs + 10;	}	//can't tell the difference between multiple primers
			else{													//use the best match
				group = minGroup;
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripForward");
		exit(1);
	}
}
//*******************************************************************/
int TrimOligos::stripForward(Sequence& seq, QualityScores& qual, int& group, bool keepForward){
	try {
        
        if (paired) {  int success = stripPairedPrimers(seq, qual, group, keepForward);  return success; }
        
		string rawSequence = seq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
		
		//can you find the primer
		for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the primers are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				if (!keepForward) { seq.setUnaligned(rawSequence.substr(oligo.length())); }
				if(qual.getName() != ""){
					if (!keepForward) {  qual.trimQScores(oligo.length(), -1); }
				}
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (primers.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxFPrimerLength+pdiffs+1));  }
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxFPrimerLength){	
					success = pdiffs + 100;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = pdiffs + 10;	}	//can't tell the difference between multiple primers
			else{													//use the best match
				group = minGroup;
				if (!keepForward) { seq.setUnaligned(rawSequence.substr(minPos)); }
				if(qual.getName() != ""){
					if (!keepForward) { qual.trimQScores(minPos, -1); }
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripForward");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::stripReverse(Sequence& seq, QualityScores& qual){
	try {
		string rawSequence = seq.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		for(int i=0;i<revPrimer.size();i++){
			string oligo = revPrimer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = 0;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(rawSequence.length()-oligo.length(),oligo.length()))){
				seq.setUnaligned(rawSequence.substr(0,rawSequence.length()-oligo.length()));
				if(qual.getName() != ""){
					qual.trimQScores(-1, rawSequence.length()-oligo.length());
				}
				success = 1;
				break;
			}
		}	
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripReverse");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::stripReverse(Sequence& seq){
	try {
		
		string rawSequence = seq.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		for(int i=0;i<revPrimer.size();i++){
			string oligo = revPrimer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = 0;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(rawSequence.length()-oligo.length(),oligo.length()))){
				seq.setUnaligned(rawSequence.substr(0,rawSequence.length()-oligo.length()));
				success = 1;
				break;
			}
		}	
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripReverse");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::stripLinker(Sequence& seq, QualityScores& qual){
	try {
		string rawSequence = seq.getUnaligned();
		bool success = ldiffs + 1;	//guilty until proven innocent
		
		for(int i=0;i<linker.size();i++){
			string oligo = linker[i];
			
			if(rawSequence.length() < oligo.length()){
				success = ldiffs + 10;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				success = 0;
				break;
			}
		}
        
        //if you found the linker or if you don't want to allow for diffs
		if ((ldiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (linker.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLinkerLength+ldiffs+1));   }	
			else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < linker.size(); i++){
				string oligo = linker[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLinkerLength){	//let's just assume that the barcodes are the same length
					success = ldiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+ldiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > ldiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = ldiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				seq.setUnaligned(rawSequence.substr(minPos));
				
				if(qual.getName() != ""){
					qual.trimQScores(minPos, -1);
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}

       
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripLinker");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::stripLinker(Sequence& seq){
	try {
		
		string rawSequence = seq.getUnaligned();
		bool success = ldiffs +1;	//guilty until proven innocent
		
		for(int i=0;i<linker.size();i++){
			string oligo = linker[i];
			
			if(rawSequence.length() < oligo.length()){
				success = ldiffs +10;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 0;
				break;
			}
		}	
		
        //if you found the linker or if you don't want to allow for diffs
		if ((ldiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (linker.size() > 0) {alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLinkerLength+ldiffs+1));  }
            else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < linker.size(); i++){
				string oligo = linker[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLinkerLength){	//let's just assume that the barcodes are the same length
					success = ldiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+ldiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > ldiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = ldiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}

		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripLinker");
		exit(1);
	}
}

//******************************************************************/
bool TrimOligos::stripSpacer(Sequence& seq, QualityScores& qual){
	try {
		string rawSequence = seq.getUnaligned();
		bool success = sdiffs+1;	//guilty until proven innocent
		
		for(int i=0;i<spacer.size();i++){
			string oligo = spacer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = sdiffs+10;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				success = 0;
				break;
			}
		}
        
        //if you found the spacer or if you don't want to allow for diffs
		if ((sdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (spacer.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxSpacerLength+sdiffs+1));  }
            else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < spacer.size(); i++){
				string oligo = spacer[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxSpacerLength){	//let's just assume that the barcodes are the same length
					success = sdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+sdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > sdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = sdiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				seq.setUnaligned(rawSequence.substr(minPos));
				
				if(qual.getName() != ""){
					qual.trimQScores(minPos, -1);
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
        

		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripSpacer");
		exit(1);
	}
}
//******************************************************************/
bool TrimOligos::stripSpacer(Sequence& seq){
	try {
		
		string rawSequence = seq.getUnaligned();
		bool success = sdiffs+1;	//guilty until proven innocent
		
		for(int i=0;i<spacer.size();i++){
			string oligo = spacer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = sdiffs+10;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 0;
				break;
			}
		}	
		
        //if you found the spacer or if you don't want to allow for diffs
		if ((sdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			Alignment* alignment;
			if (spacer.size() > 0) { alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxSpacerLength+sdiffs+1));  }
            else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < spacer.size(); i++){
				string oligo = spacer[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxSpacerLength){	//let's just assume that the barcodes are the same length
					success = sdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->alignPrimer(oligo, rawSequence.substr(0,oligo.length()+sdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > sdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = sdiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}

		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripSpacer");
		exit(1);
	}
}

//******************************************************************/
bool TrimOligos::compareDNASeq(string oligo, string seq){
	try {
		bool success = 1;
		int length = oligo.length();
		
		for(int i=0;i<length;i++){
			
			if(oligo[i] != seq[i]){
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C')	{	success = 0; 	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))				{	success = 0;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))					{	success = 0;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	success = 0;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	success = 0;	}			
				
				if(success == 0)	{	break;	 }
			}
			else{
				success = 1;
			}
		}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "compareDNASeq");
		exit(1);
	}
	
}
//********************************************************************/
int TrimOligos::countDiffs(string oligo, string seq){
	try {
		
		int length = oligo.length();
		int countDiffs = 0;
		
		for(int i=0;i<length;i++){
			
			if(oligo[i] != seq[i]){
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C' || oligo[i] == '-' || oligo[i] == '.')	{	countDiffs++; 	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))				{	countDiffs++;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))					{	countDiffs++;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))					{	countDiffs++;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))					{	countDiffs++;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	countDiffs++;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	countDiffs++;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	countDiffs++;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	countDiffs++;	}	
			}
			
		}
		
		return countDiffs;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "countDiffs");
		exit(1);
	}
}
//********************************************************************/
string TrimOligos::reverseOligo(string oligo){
	try {
        string reverse = "";
        
        for(int i=oligo.length()-1;i>=0;i--){
            
            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}
            
            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}
            
            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}
            
            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}
            
            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}
            
            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}
            
            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}
            
            else if(oligo[i] == '-'){	reverse += '-';	} 
            else						{	reverse += 'N';	}
        }
        
        
        return reverse;
    }
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "reverseOligo");
		exit(1);
	}
}

/********************************************************************/



