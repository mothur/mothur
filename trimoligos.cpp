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
TrimOligos::TrimOligos(int p, int b, int l, int s, map<string, int> pr, map<string, int> br, map<string, int> rbr, vector<string> r, vector<string> lk, vector<string> sp){
	try {
		m = MothurOut::getInstance();
		
		pdiffs = p;
		bdiffs = b;
        ldiffs = l;
        sdiffs = s;
		
		barcodes = br;
        rbarcodes = rbr;
		primers = pr;
		revPrimer = r;
        linker = lk;
        spacer = sp;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "TrimOligos");
		exit(1);
	}
}
/********************************************************************/
//strip, pdiffs, bdiffs, primers, barcodes, revPrimers
TrimOligos::TrimOligos(int p, int b, int l, int s, map<string, int> pr, map<string, int> br, vector<string> r, vector<string> lk, vector<string> sp){
	try {
		m = MothurOut::getInstance();
		
		pdiffs = p;
		bdiffs = b;
        ldiffs = l;
        sdiffs = s;
		
		barcodes = br;
		primers = pr;
		revPrimer = r;
        linker = lk;
        spacer = sp;
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
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "TrimOligos");
		exit(1);
	}
}
/********************************************************************/
TrimOligos::~TrimOligos() {}
//*******************************************************************/
int TrimOligos::stripBarcode(Sequence& seq, QualityScores& qual, int& group){
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
				
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((bdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (barcodes.size() > 0) {
				map<string,int>::iterator it; 
				
				for(it=barcodes.begin();it!=barcodes.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+bdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+bdiffs));
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
		m->errorOut(e, "TrimOligos", "stripBarcode");
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (barcodes.size() > 0) {
				map<string,int>::iterator it=barcodes.begin();
				
				for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+bdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+bdiffs));
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
//*******************************************************************/
int TrimOligos::stripRBarcode(Sequence& seq, QualityScores& qual, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the barcode
		for(map<string,int>::iterator it=rbarcodes.begin();it!=rbarcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
            
			if(compareDNASeq(oligo, rawSequence.substr((rawSequence.length()-oligo.length()),oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(0,(rawSequence.length()-oligo.length())));
				
				if(qual.getName() != ""){
					qual.trimQScores(-1, rawSequence.length()-oligo.length());
				}
				
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((bdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (rbarcodes.size() > 0) {
				map<string,int>::iterator it; 
				
				for(it=rbarcodes.begin();it!=rbarcodes.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+bdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=rbarcodes.begin();it!=rbarcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr((rawSequence.length()-(oligo.length()+bdiffs)),oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
     
				int alnLength = oligo.length();
				
				for(int i=0;i<alnLength;i++){
					if(oligo[i] != '-'){	alnLength = i;	break;	}
				}
				oligo = oligo.substr(alnLength);
				temp = temp.substr(alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=alnLength-1;i>=0;i--){
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
				seq.setUnaligned(rawSequence.substr(0, (rawSequence.length()-minPos)));
                
				if(qual.getName() != ""){
					qual.trimQScores(-1, (rawSequence.length()-minPos));
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripRBarcode");
		exit(1);
	}
	
}
//*******************************************************************/
int TrimOligos::stripRBarcode(Sequence& seq, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the barcode
		for(map<string,int>::iterator it=rbarcodes.begin();it!=rbarcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
            
			if(compareDNASeq(oligo, rawSequence.substr((rawSequence.length()-oligo.length()),oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(0,(rawSequence.length()-oligo.length())));
				
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((bdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (rbarcodes.size() > 0) {
				map<string,int>::iterator it; 
				
				for(it=rbarcodes.begin();it!=rbarcodes.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+bdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=rbarcodes.begin();it!=rbarcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr((rawSequence.length()-(oligo.length()+bdiffs)),oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
                
				int alnLength = oligo.length();
				
				for(int i=0;i<alnLength;i++){
					if(oligo[i] != '-'){	alnLength = i;	break;	}
				}
				oligo = oligo.substr(alnLength);
				temp = temp.substr(alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=alnLength-1;i>=0;i--){
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
				seq.setUnaligned(rawSequence.substr(0, (rawSequence.length()-minPos)));
                
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripRBarcode");
		exit(1);
	}
	
}
//********************************************************************/
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (primers.size() > 0) {
				map<string,int>::iterator it; 
				
				for(it=primers.begin();it!=primers.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+pdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	
					success = pdiffs + 100;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+pdiffs));
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (primers.size() > 0) {
				map<string,int>::iterator it; 
				
				for(it=primers.begin();it!=primers.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+pdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	
					success = pdiffs + 100;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+pdiffs));
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (linker.size() > 0) {
				for(int i = 0; i < linker.size(); i++){
					if(linker[i].length() > maxLength){
						maxLength = linker[i].length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+ldiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < linker.size(); i++){
				string oligo = linker[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = ldiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+ldiffs));
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (linker.size() > 0) {
				for(int i = 0; i < linker.size(); i++){
					if(linker[i].length() > maxLength){
						maxLength = linker[i].length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+ldiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < linker.size(); i++){
				string oligo = linker[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = ldiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+ldiffs));
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (spacer.size() > 0) {
				for(int i = 0; i < spacer.size(); i++){
					if(spacer[i].length() > maxLength){
						maxLength = spacer[i].length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+sdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < spacer.size(); i++){
				string oligo = spacer[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = sdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+sdiffs));
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
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (spacer.size() > 0) {
				for(int i = 0; i < spacer.size(); i++){
					if(spacer[i].length() > maxLength){
						maxLength = spacer[i].length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+sdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minPos = 0;
			
			for(int i = 0; i < spacer.size(); i++){
				string oligo = spacer[i];
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = sdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+sdiffs));
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
/********************************************************************/



