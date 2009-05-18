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

#include "mothur.h"

/**************************************************************************************************/

class Kmer {
	
public:
	Kmer(int);
	string getKmerString(string);
	int getKmerNumber(string, int);
	string getKmerBases(int);
	
	
private:
	char getASCII(int);
	int getNumber(char);
	int kmerSize;
	int maxKmer;
	int nKmers;
};

/**************************************************************************************************/


#endif
