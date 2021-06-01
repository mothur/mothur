//
//  protein.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/24/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef protein_hpp
#define protein_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "writer.h"


/**************************************************************************************************/

class Protein {
    
#ifdef UNIT_TEST
    friend class TestProtein;
#endif
    
public:
    
    Protein();
    Protein(string, string);
    Protein(ifstream&);
    Protein(ifstream&, string&, bool);
    Protein(istringstream&);
    ~Protein() {}
    
    void setName(string);
    string getName();
    void setUnaligned(string);
    string getUnaligned();
    void setAligned(string);
    string getAligned();
    void setComment(string);
    string getComment();
    string getInlineProtein();
    
    int getNumBases();
    int getStartPos();
    int getEndPos();
    
    void reverseComplement();
    void trim(int);
    void padToPos(int);
    void padFromPos(int);
    int filterToPos(int); //any character before the pos is changed to . and aligned and unaligned strings changed
    int filterFromPos(int); //any character after the pos is changed to . and aligned and unaligned strings changed
    int getAlignLength();
    
    void printProtein(ostream&);
    void printProtein(OutputWriter*);
    void printUnAlignedProtein(ostream&);
    
protected:
    
    MothurOut* m;
    Utils util;
    
    void initialize();
    string getProteinString(ifstream&, int&);
    string getCommentString(ifstream&);
    string getProteinString(istringstream&, int&);
    string getCommentString(istringstream&);
    string getProteinName(ifstream&);
    string getProteinName(istringstream&);
    
    string name;
    string unaligned;
    string aligned;
    string comment;
    int numBases;
    int alignmentLength;
    int startPos, endPos;
    
    //string pairwise;
    //int longHomoPolymer;
    //int ambigBases;
   
    
};

/**************************************************************************************************/


#endif /* protein_hpp */
