/*
 *  kmer.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "kmer.hpp"

/**************************************************************************************************/

Kmer::Kmer(int size) : kmerSize(size) {	//	The constructor sets the size of the kmer
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
																		//	No reason to waste the time of recalculating
	maxKmer = power4s[kmerSize]+1;// (int)pow(4.,k)+1;					//	powers of 4 everytime through.  We need an
																		//	extra kmer if we get a non-ATGCU base
}

/**************************************************************************************************/

string Kmer::getKmerString(string sequence){	//	Calculate kmer for each position in the sequence, count the freq
	int length = sequence.length();				//	of each kmer, and convert it to an ascii character with base '!'.
	int nKmers = length - kmerSize + 1;			//	Export the string of characters as a string
	vector<int> counts(maxKmer, 0);
	
	for(int i=0;i<nKmers;i++){					//	Go though sequence and get the number between 0 and maxKmer for that
		int kmerNumber = getKmerNumber(sequence, i);//	kmer.  Increase the count for the kmer in the counts vector
		counts[kmerNumber]++;
	}
	
	string kmerString = "";						
	for(int i=0;i<maxKmer;i++){					//	Scan through the vector of counts and convert each element of the 
		kmerString += getASCII(counts[i]);		//	vector to an ascii character that is equal to or greater than '!'
	}		
	
	return kmerString;	
}
/**************************************************************************************************/

vector< map<int, int> > Kmer::getKmerCounts(string sequence){	//	Calculate kmer for each position in the sequence, save info in a map
	int length = sequence.length();				//	so you know at each spot in the sequence what kmers were found
	int nKmers = length - kmerSize + 1;			//	
	vector< map<int, int> > counts; counts.resize(nKmers);  // a map kmer counts for each spot
	 map<int, int>::iterator it;
	
	for(int i=0;i<nKmers;i++){					//	Go though sequence and get the number between 0 and maxKmer for that
		if (i == 0) {
			int kmerNumber = getKmerNumber(sequence, i);//	kmer. 
			counts[i][kmerNumber] = 1;			// add this kmer if not already there
		}else { 
			//your count is everything that came before and whatever you find now
			counts[i] = counts[i-1];
			int kmerNumber = getKmerNumber(sequence, i);//	kmer.  
			
			it = counts[i].find(kmerNumber);
			if (it!= counts[i].end()) {   counts[i][kmerNumber]++;  }  //increment number of times you have seen this kmer
			else {   counts[i][kmerNumber] = 1;   }  // add this kmer since not already there
		}
	}
	
	return counts;	
}
	
	
/**************************************************************************************************/

int Kmer::getKmerNumber(string sequence, int index){
	
//	Here we convert a kmer to a number between 0 and maxKmer.  For example, AAAA would equal 0 and TTTT would equal 255.
//	If there's an N in the kmer, it is set to 256 (if we are looking at 4mers).  The largest we can look at are 8mers,
//	this could be easily increased.
//	
//	Example:   ATGCAT (kSize = 6)
//		  i:   012345
//
//		Score	= 0*4^(6-0-1) + 3*4^(6-1-1) + 2*4^(6-2-1) + 1*4^(6-3-1) + 0*4^(6-4-1) + 3*4^(6-5-1)
//				= 0*4^5	+	3*4^4	+	2*4^3	+	1*4^2	+	0*4^1	+	3*4^0
//				= 0 + 3*256 + 2*64 + 1*16 + 0*4 + 3*1
//				= 0 + 768 + 128 + 16 + 0 + 3
//				= 915
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
	
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
	
//	Here we convert the kmer number into the kmer in terms of bases.
//
//	Example:	Score = 915 (for a 6-mer)
//				Base6 = (915 / 4^0) % 4 = 915 % 4 = 3 => T	[T]
//				Base5 = (915 / 4^1) % 4 = 228 % 4 = 0 => A	[AT]
//				Base4 = (915 / 4^2) % 4 = 57 % 4 = 1 => C	[CAT]
//				Base3 = (915 / 4^3) % 4 = 14 % 4 = 2 => G	[GCAT]
//				Base2 = (915 / 4^4) % 4 = 3 % 4 = 3 => T	[TGCAT]
//				Base1 = (915 / 4^5) % 4 = 0 % 4 = 0 => A	[ATGCAT] -> this checks out with the previous method
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
	
	string kmer = "";
	
	if(kmerNumber == power4s[kmerSize]){//pow(4.,7)){	//	if the kmer number is the same as the maxKmer then it must
		for(int i=0;i<kmerSize;i++){					//	have had an N in it and so we'll just call it N x kmerSize
			kmer += 'N';
		}
	}
	else{
		for(int i=0;i<kmerSize;i++){
			int nt = (int)(kmerNumber / (float)power4s[i]) % 4;		//	the '%' operator returns the remainder 
			if(nt == 0)		{	kmer = 'A' + kmer;	}				//	from int-based division ]
			else if(nt == 1){	kmer = 'C' + kmer;	}
			else if(nt == 2){	kmer = 'G' + kmer;	}
			else if(nt == 3){	kmer = 'T' + kmer;	}
		}
	}
	return kmer;
}
/**************************************************************************************************/

int Kmer::getReverseKmerNumber(int kmerNumber){
		
	string kmerString = getKmerBases(kmerNumber);
	
	//get Reverse
	string reverse = "";
	for(int i=kmerString.length()-1;i>=0;i--){
		if(kmerString[i] == 'A')		{	reverse += 'T';	}
		else if(kmerString[i] == 'T'){	reverse += 'A';	}
		else if(kmerString[i] == 'G'){	reverse += 'C';	}
		else if(kmerString[i] == 'C'){	reverse += 'G';	}
		else						{	reverse += 'N';	}
	}
	
	int reverseNumber = getKmerNumber(reverse, 0);
	
	return reverseNumber;
	
}

/**************************************************************************************************/

char Kmer::getASCII(int number)		{	return (char)(33+number);			}	// '!' is the first printable char and
																				// has the int value of 33
/**************************************************************************************************/

int Kmer::getNumber(char character)	{	return ((int)(character-'!'));		}	// '!' has the value of 33

/**************************************************************************************************/
