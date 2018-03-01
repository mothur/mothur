#ifndef KNN_H
#define KNN_H

/*
 *  knn.h
 *  Mothur
 *
 *  Created by westcott on 11/4/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */
 
#include "mothur.h"
#include "classify.h"

/**************************************************************************************************/

class Knn : public Classify {
	
public:
	Knn(string, string, string, int, float, float, float, float, int, int, string);
	~Knn();
	
	void setDistName(string s);
	string getTaxonomy(Sequence*, string&, bool&);
	
private:
	int num;
	string findCommonTaxonomy(vector<string>);
	string search, outDistName;
	
};

/**************************************************************************************************/

#endif


