#ifndef TAXONOMYNODE
#define TAXONOMYNODE

/*
 *  taxonomynode.h
 *  
 *
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "mothurout.h"
/**************************************************************************************************/

class TaxonomyNode {
	
public:
	TaxonomyNode();
	TaxonomyNode(string, int);
    
	void setName(string n)  {    name = n;      }
    string getName()        { return name;      }
    void setParent(int p)   { parent = p;       }
    int getParent()         { return parent;    }
	
	void makeChild(string c, int i)     {    children[c] = i;               }
	map<string, int> getChildren()      {    return children;               }
	int	getNumKids()                    {    return (int)children.size();   }
	int getNumSeqs()                    {    return numSeqs;                }
	void setTotalSeqs(int n)            {    totalSeqs = n;                 }
	int getLevel()                      {    return level;                  }
    
    int getChildIndex(string);
	
private:
	int parent;
	map<string, int> children;
	int numChildren;
	int level;
	
protected:
    MothurOut* m;
	int numSeqs;
	int totalSeqs;
	string name;
};

/**************************************************************************************************/

#endif
