#ifndef READBLAST_H
#define READBLAST_H
/*
 *  readblast.h
 *  Mothur
 *
 *  Created by westcott on 12/10/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "sparsematrix.hpp"
#include "nameassignment.hpp"

/****************************************************************************************/
struct DistNode {
	int seq1;
	int seq2;
	float dist;
	DistNode(int s1, int s2, float d) : seq1(s1), seq2(s2), dist(d) {}
	DistNode() {}
	~DistNode() {}
};
/****************************************************************************************/

//Note: this class creates a sparsematrix and list if the read is executed, but does not delete them on deconstruction.
//the user of this object is responsible for deleting the matrix and list if they call the read or there will be a memory leak
//it is done this way so the read can be deleted and the information still used.

class ReadBlast {
	
public:
	ReadBlast(string, float, float, int, bool, bool); //blastfile, cutoff, penalty, length of overlap, min or max bsr, hclusterWanted
	~ReadBlast() {}
	
	void read(NameAssignment*);
	SparseMatrix* getDistMatrix()		{	return matrix;		}
	vector<DistNode> getOverlapMatrix()	{	return overlap;		}
	string getOverlapFile()				{	return overlapFile;	}
	string getDistFile()				{	return distFile;	}
	
private:
	string blastfile, overlapFile, distFile;
	int length;	//number of amino acids overlapped
	float penalty, cutoff;  //penalty is used to adjust error rate
	bool minWanted;  //if true choose min bsr, if false choose max bsr
	bool hclusterWanted;
	
	SparseMatrix* matrix;
	vector<DistNode> overlap;
	
	void readNames(NameAssignment*);
};

/*******************************************************************************************/

#endif

