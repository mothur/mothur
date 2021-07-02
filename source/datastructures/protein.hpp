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
#include "aminoacid.hpp"

/**************************************************************************************************/

class Protein {
    
#ifdef UNIT_TEST
    friend class TestProtein;
#endif
    
public:
    
    Protein();
    Protein(string, vector<AminoAcid>);
    Protein(ifstream&);
    Protein(ifstream&, string&, bool);
    Protein(istringstream&);
#ifdef USE_BOOST
    Protein(boost::iostreams::filtering_istream&);
#endif
    ~Protein() {}
    
    void setName(string);
    string getName();
    void setUnaligned(vector<AminoAcid>);
    vector<AminoAcid> getUnaligned();
    void setAligned(vector<AminoAcid>);
    vector<AminoAcid> getAligned();
    void setComment(string);
    string getComment();
    string getInlineProtein();
    void setPairwise(vector<AminoAcid>);
    vector<AminoAcid> getPairwise();
    
    int getNumBases();
    int getStartPos();
    int getEndPos();
    
    void trim(int);
    void padToPos(int);
    void padFromPos(int);
    void filterToPos(int); //any character before the pos is changed to . and aligned and unaligned strings changed
    void filterFromPos(int); //any character after the pos is changed to . and aligned and unaligned strings changed
    int getAlignLength();
    
    void printProtein(ostream&);
    void printProtein(OutputWriter*);
    void printUnAlignedProtein(ostream&);
    
protected:
    
    MothurOut* m;
    Utils util;
    
    void initialize();
    vector<AminoAcid> getProtein(ifstream&);
    vector<AminoAcid> getProtein(istringstream&);
    string getCommentString(ifstream&);
    string getCommentString(istringstream&);
    string getProteinName(ifstream&);
    string getProteinName(istringstream&);
    string getProteinString(vector<AminoAcid>);
    
#ifdef USE_BOOST
    string getCommentString(boost::iostreams::filtering_istream&);
    vector<AminoAcid> getProtein(boost::iostreams::filtering_istream&);
    string getSequenceName(boost::iostreams::filtering_istream&);
#endif
    
    string name;
    vector<AminoAcid> unaligned;
    vector<AminoAcid> aligned;
    string comment;
    int numBases;
    int alignmentLength;
    int startPos, endPos;
    
    vector<AminoAcid> pairwise;
    
   
    
};

/**************************************************************************************************/


#endif /* protein_hpp */
