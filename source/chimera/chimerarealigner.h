#ifndef CHIMERAREALIGNER_H
#define CHIMERAREALIGNER_H

/*
 *  chimerarealigner.h
 *  Mothur
 *
 *  Created by westcott on 2/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimera.h"
#include "alignment.hpp"

/***********************************************************/

struct AlignCell {
	int score;
	char direction;
	AlignCell() : score(0), direction('x') {};
};

/***********************************************************/

struct  bases {
	int A, T, G, C, Gap, Chars;
	bases() : A(0), T(0), G(0), C(0), Gap(0), Chars(0){};
};

/***********************************************************/


class ChimeraReAligner  {
	
public:
	ChimeraReAligner();	 
	~ChimeraReAligner();
	
	void reAlign(Sequence*, vector<string>);
				
private:
	void buildTemplateProfile(vector<string>);
	void createAlignMatrix(int, int);
	void fillAlignMatrix(string);
	int calcMatchScore(bases, char);
	string getNewAlignment(string);

	int alignmentLength;
	vector<bases> profile;
	vector<vector<AlignCell> > alignMatrix;

	MothurOut* m;
};

/***********************************************************/

#endif

