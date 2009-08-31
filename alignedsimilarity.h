#ifndef ALIGNSIM_H
#define ALIGNSIM_H

/*
 *  alignedsimilarity.h
 *  Mothur
 *
 *  Created by westcott on 8/18/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "chimera.h"
#include "dist.h"

//This class was created using the algorythms described in the 
//"Evaluation of Nearest Neighbor Methods for Detection of Chimeric Small-Subunit rRna Sequences" paper 
//by J.F. Robison-Cox 1, M.M. Bateson 2 and D.M. Ward 2


/***********************************************************/

class AlignSim : public Chimera {
	
	public:
		AlignSim(string, string);	
		~AlignSim();
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string){};
		void setQuantiles(string q) { quanfile = q; };
		
		
	private:
		
		Dist* distCalc;
		vector<linePair*> lines;
		vector<linePair*> templateLines;
		
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		
		vector< vector<sim> > IS;  //IS[0] is the vector os sim values for each window for querySeqs[0]
		vector< vector<sim> > templateIS;  //templateIS[0] is the vector os sim values for each window for templateSeqs[0]
		vector<int> windowBreak;
		
		string fastafile, templateFile;
		int iters;
		
		vector< vector<sim> > findIS(int, int, vector<Sequence*>);
		vector<Sequence*> findClosestSides(Sequence*, int, vector<int>&, int);
		int findNumMatchedBases(Sequence, Sequence);
		vector<int> findWindows();
		
		vector< vector<float> > quantile;
		
		vector< vector<sim> > createProcessesIS(vector<Sequence*>, vector<linePair*>);
		
};

/***********************************************************/

#endif

