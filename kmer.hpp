#ifndef KMER_HPP
#define KMER_HPP

/*
 *  kmer.hpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

class Kmer {
	
public:
	Kmer(int);
	string getKmerString(string);
	int getKmerNumber(string, int);
	
	
private:
	string getKmerBases(int);
	
	char getASCII(int);
	int getNumber(char);
	int kmerSize;
	int maxKmer;
	int nKmers;
};

/**************************************************************************************************/

#endif
