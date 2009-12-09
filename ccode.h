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

/***********************************************************/
//This class was created using the algorythms described in the 
// "Evaluating putative chimeric sequences from PCR-amplified products" paper 
//by Juan M. Gonzalez, Johannes Zimmerman and Cesareo Saiz-Jimenez.

/***********************************************************/

class Ccode : public Chimera {
	
	public:
		Ccode(string, string);	
		~Ccode();
		
		int getChimeras();
		void print(ostream&);
		
		void setCons(string c)		{}
		void setQuantiles(string q) {}
		
		
	private:
	
		Dist* distCalc;
		DeCalculator* decalc;
		int iters;
		string fastafile, templateFile;
		
		
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		vector< map<int, int> > spotMap;
		map<int, int>::iterator it;
		
		vector< vector<int> > windows; //windows[0] is the vector of window breaks for querySeqs[0]
		vector<int> windowSizes;  //windowSizes[0] is the size of the windows for querySeqs[0]
		vector< map<int, int> > trim;  //trim[0] is the map containing the starting and ending positions for querySeqs[0] set of seqs
		vector< vector<SeqDist> > closest;  //closest[0] is a vector of sequence at are closest to queryseqs[0]...
		vector< vector<float> > averageRef;  //averageRef[0] is the average distance at each window for the references for querySeqs[0]
		vector< vector<float> > averageQuery;  //averageQuery[0] is the average distance at each winow for the query for querySeqs[0]
		vector< vector<float> >  sumRef;  //sumRef[0] is the sum of distances at each window for the references for querySeqs[0]
		vector< vector<float> >  sumSquaredRef;  //sumSquaredRef[0] is the sum of squared distances at each window for the references for querySeqs[0]
		vector< vector<float> > sumQuery;  //sumQuery[0] is the sum of distances at each window for the comparison of query to references for querySeqs[0]
		vector< vector<float> >  sumSquaredQuery;  //sumSquaredQuery[0] is the sum of squared distances at each window for the comparison of query to references for querySeqs[0]
		vector< vector<float> > varRef;  //varRef[0] is the variance among references seqs at each window for querySeqs[0]
		vector< vector<float> > varQuery;  //varQuery[0] is the variance among references and querySeqs[0] at each window
		vector< vector<float> > sdRef;  //sdRef[0] is the standard deviation of references seqs at each window for querySeqs[0]
		vector< vector<float> > sdQuery;  //sdQuery[0] is the standard deviation of references and querySeqs[0] at each window
		vector< vector<float> > anova;  //anova[0] is the vector of anova scores for each window for querySeqs[0]
		vector<int> refCombo;  //refCombo[0] is the number of reference sequences combinations for querySeqs[0]
		vector< vector<bool> > isChimericConfidence;  //isChimericConfidence[0] indicates whether querySeqs[0] is chimeric at a given window according to the confidence limits
		vector< vector<bool> > isChimericTStudent;  //isChimericConfidence[0] indicates whether querySeqs[0] is chimeric at a given window according to the confidence limits
		vector< vector<bool> > isChimericANOVA;  //isChimericConfidence[0] indicates whether querySeqs[0] is chimeric at a given window according to the confidence limits
		
		vector< vector<SeqDist> > findClosest(int, int, int); 
		void removeBadReferenceSeqs(vector<SeqDist>&, int);  //removes sequences from closest that are to different of too similar to eachother. 
		void trimSequences(int);
		vector<int> findWindows(int);
		void getAverageRef(vector<SeqDist>, int);		//fills sumRef[i], averageRef[i], sumSquaredRef[i] and refCombo[i].
		void getAverageQuery (vector<SeqDist>, int);	//fills sumQuery[i], averageQuery[i], sumSquaredQuery[i].
		void findVarianceRef (int);						//fills varRef[i] and sdRef[i] also sets minimum error rate to 0.001 to avoid divide by 0.
		void findVarianceQuery (int);					//fills varQuery[i] and sdQuery[i]
		void determineChimeras (int);					//fills anova, isChimericConfidence[i], isChimericTStudent[i] and isChimericANOVA[i].
		
		int getDiff(string, string);  //return number of mismatched bases, a gap to base is not counted as a mismatch
		float getT(int); 
		float getF(int); 
		
		void createProcessesClosest();
		void createProcessesRemoveBad();
		void createProcessesAverages();
		void createProcessesVariances();
		void createProcessesDetermine();
				
};

/***********************************************************/

#endif


