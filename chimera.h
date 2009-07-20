#ifndef CHIMERA_H
#define CHIMERA_H

/*
 *  chimera.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "mothur.h"
#include "sparsematrix.hpp"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;  //maps sequence to all distance for that seqeunce


/***********************************************************************/

class Chimera {

	public:
	
		Chimera(){};
		Chimera(string);
		Chimera(string, string);
		virtual ~Chimera(){};
		virtual void setFilter(bool f)			{	filter = f;			}
		virtual void setCorrection(bool c)		{	correction = c;		}
		virtual void setProcessors(int p)		{	processors = p;		}
		virtual void setWindow(int w)			{	window = w;			}
		virtual void setIncrement(int i)		{	increment = i;		}
		
		virtual void setCons(string) {};
		
		
		//pure functions
		virtual void getChimeras() = 0;	
		virtual void print(ostream&) = 0;	
		
	protected:
		
		bool filter, correction;
		int processors, window, increment;
			

};

/***********************************************************************/

#endif

