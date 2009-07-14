#ifndef PINTAIL_H
#define PINTAIL_H

/*
 *  pintail.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimera.h"
#include "dist.h"

//This class was created using the algorythms described in the 
// "At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.

/***********************************************************/

class Pintail : public Chimera {
	
	public:
		Pintail(string);	
		~Pintail();
		
		void getChimeras();
		void print(ostream&);
		
		
	private:
	
		struct linePair {
			int start;
			int end;
			linePair(int i, int j) : start(i), end(j) {}
		};

	
		Dist* distCalculator;
		string fastafile;
		int iters;
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		
		map<Sequence*, Sequence*> bestfit;  //maps a query sequence to its most similiar sequence in the template
		map<Sequence*, Sequence*>::iterator itBest;
		
		map<Sequence*, vector<float> > obsDistance;  //maps a query sequence to its observed distance at each window
		map<Sequence*, vector<float> > expectedDistance;  //maps a query sequence to its expected distance at each window
		map<Sequence*, vector<float> >::iterator itObsDist;
		map<Sequence*, vector<float> >::iterator itExpDist;
		
		vector<float> averageProbability;			//Qav
		map<Sequence*, float> seqCoef;				//maps a sequence to its coefficient
		map<Sequence*, float> DE;					//maps a sequence to its deviation
		map<Sequence*, float>::iterator itCoef;	
		
		vector<Sequence*> readSeqs(string);
		vector<float> findQav(vector<float>);
		vector<float> calcFreq(vector<Sequence*>);
		map<Sequence*, float> getCoef(vector<float>);
		
		void findPairs(int, int);
		void calcObserved(int, int);
		void calcExpected(int, int);
		void calcDE(int, int);
	
		void createProcessesPairs();
		void createProcessesObserved();
		void createProcessesExpected();
		void createProcessesDE();
		
};

/***********************************************************/

#endif

