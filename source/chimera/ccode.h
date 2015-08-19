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
//This class was created using the algorithms described in the 
// "Evaluating putative chimeric sequences from PCR-amplified products" paper 
//by Juan M. Gonzalez, Johannes Zimmerman and Cesareo Saiz-Jimenez.

/***********************************************************/

class Ccode : public Chimera {
	
	public:
		Ccode(string, string, bool, string, int, int, string);	//fasta, template, filter, mask, window, numWanted, outputDir
		~Ccode();
		
		int getChimeras(Sequence* query);
		Sequence print(ostream&, ostream&);
    
	private:
	
		Dist* distCalc;
		DeCalculator* decalc;
		int iters, window, numWanted;
		string fastafile, mapInfo;
		
		Sequence* querySeq;
		
		map<int, int> spotMap;
		map<int, int>::iterator it;
		
		vector<int>  windows; //windows is the vector of window breaks for query
		int windowSizes;  //windowSizes is the size of the windows for query
		map<int, int> trim;  //trim is the map containing the starting and ending positions for query
		vector<SeqDist>  closest;  //closest is a vector of sequence at are closest to query
		vector<float>  averageRef;  //averageRef is the average distance at each window for the references for query
		vector<float>  averageQuery;  //averageQuery is the average distance at each winow for the query for query
		vector<float>   sumRef;  //sumRef is the sum of distances at each window for the references for query
		vector<float>   sumSquaredRef;  //sumSquaredRef is the sum of squared distances at each window for the references for query
		vector<float> sumQuery;  //sumQuery is the sum of distances at each window for the comparison of query to references for query
		vector<float>  sumSquaredQuery;  //sumSquaredQuery is the sum of squared distances at each window for the comparison of query to references for query
		vector<float> varRef;  //varRef is the variance among references seqs at each window for query
		vector<float> varQuery;  //varQuery is the variance among references and query at each window
		vector<float> sdRef;  //sdRef is the standard deviation of references seqs at each window for query
		vector<float> sdQuery;  //sdQuery is the standard deviation of references and query at each window
		vector<float> anova;  //anova is the vector of anova scores for each window for query
		int refCombo;  //refCombo is the number of reference sequences combinations for query
		vector<bool>  isChimericConfidence;  //isChimericConfidence indicates whether query is chimeric at a given window according to the confidence limits
		vector<bool>  isChimericTStudent;  //isChimericConfidence indicates whether query is chimeric at a given window according to the confidence limits
		vector<bool>  isChimericANOVA;  //isChimericConfidence indicates whether query is chimeric at a given window according to the confidence limits
		
		vector<SeqDist>  findClosest(Sequence*, int); 
		void removeBadReferenceSeqs(vector<SeqDist>&);  //removes sequences from closest that are to different of too similar to eachother. 
		void trimSequences(Sequence*);
		vector<int> findWindows();
		void getAverageRef(vector<SeqDist>);		//fills sumRef, averageRef, sumSquaredRef and refCombo.
		void getAverageQuery (vector<SeqDist>, Sequence*);	//fills sumQuery, averageQuery, sumSquaredQuery.
		void findVarianceRef ();						//fills varRef and sdRef also sets minimum error rate to 0.001 to avoid divide by 0.
		void findVarianceQuery ();					//fills varQuery and sdQuery
		void determineChimeras ();					//fills anova, isChimericConfidence, isChimericTStudent and isChimericANOVA.
		
		int getDiff(string, string);  //return number of mismatched bases, a gap to base is not counted as a mismatch
		float getT(int); 
		float getF(int); 
};

/***********************************************************/

#endif


