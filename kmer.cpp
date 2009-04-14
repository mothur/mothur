/*
 *  kmer.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;


#include "kmer.hpp"

/**************************************************************************************************/

Kmer::Kmer(int size) : kmerSize(size) {
	
	int power4s[9] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536 };
	maxKmer = power4s[kmerSize]+1;// (int)pow(4.,k)+1;
	
}

/**************************************************************************************************/

string Kmer::getKmerString(string sequence){
	int length = sequence.length();
	int nKmers = length - kmerSize + 1;
	vector<int> counts(maxKmer, 0);
	
	for(int i=0;i<nKmers;i++){
		int kmerNumber = getKmerNumber(sequence, i);
		counts[kmerNumber]++;
	}
	
	string kmerString = "";
	for(int i=0;i<maxKmer;i++){
		kmerString += getASCII(counts[i]);
	}		
	
	return kmerString;	
}
	
/**************************************************************************************************/

int Kmer::getKmerNumber(string sequence, int index){
	int power4s[9] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536 };
	int kmer = 0;
	for(int i=0;i<kmerSize;i++){
		if(toupper(sequence[i+index]) == 'A')		{	kmer += (0 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'C')	{	kmer += (1 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'G')	{	kmer += (2 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'U')	{	kmer += (3 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'T')	{	kmer += (3 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'N')	{	return (int)power4s[kmerSize];			}
	}
	return kmer;	
}
	
/**************************************************************************************************/
	
string Kmer::getKmerBases(int kmerNumber){
	int power4s[9] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536 };
	
	string kmer = "";
	
	if(kmerNumber == power4s[kmerSize]){//pow(4.,7)){	
		for(int i=0;i<kmerSize;i++){
			kmer += 'N';
		}
	}
	else{
		for(int i=0;i<kmerSize;i++){
			int nt = (int)(kmerNumber / (float)power4s[i]) % 4;
			if(nt == 0)		{	kmer = 'A' + kmer;	}
			else if(nt == 1){	kmer = 'C' + kmer;	}
			else if(nt == 2){	kmer = 'G' + kmer;	}
			else if(nt == 3){	kmer = 'T' + kmer;	}
		}
	}
	return kmer;
}

/**************************************************************************************************/

char Kmer::getASCII(int number)		{	return (char)(33+number);			}

/**************************************************************************************************/

int Kmer::getNumber(char character)	{	return ((int)(character-'!'));		}

/**************************************************************************************************/
