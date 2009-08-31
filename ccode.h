#ifndef CCODE_H
#define CCODE_H

/*
 *  ccode.h
 *  Mothur
 *
 *  Created by westcott on 8/24/09.
 *  Copyright 2009 Schloss LAB. All rights reserved.
 *
 */

#include "chimera.h"
#include "dist.h"
#include "decalc.h"

//This class was created using the algorythms described in the 
// "Evaluating putative chimeric sequences from PCR-amplified products" paper 
//by Juan M. Gonzalez, Johannes Zimmerman and Cesareo Saiz-Jimenez.

/***********************************************************/

class Ccode : public Chimera {
	
	public:
		Ccode(string, string);	
		~Ccode();
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string c)		{}
		void setQuantiles(string q) {}
		
		
	private:
	
		Dist* distCalc;
		DeCalculator* decalc;
		int iters;
		string fastafile, templateFile;
		
		
		vector<linePair*> lines;
		vector<linePair*> templateLines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		
		vector< vector<Sequence*> > closest;  //bestfit[0] is a vector of sequence at are closest to queryseqs[0]...
		
		vector< vector<Sequence*> > findClosest(int, int, int); 
		void removeSeqs(vector<Sequence*>);  //removes sequences from closest that are to different of too similar to eachother. 
		
		void createProcessesClosest();
		
};

/***********************************************************/

#endif


