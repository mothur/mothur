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

class ChimeraReAligner  {
	
	public:
		ChimeraReAligner();	 
		~ChimeraReAligner();
		
		void reAlign(Sequence*, vector<results>);
				
	private:
		Sequence* querySeq;
		MothurOut* m;
};
/***********************************************************/

#endif

