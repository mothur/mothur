//
//  aminoacid.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/24/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef aminoacid_hpp
#define aminoacid_hpp

#include "mothurout.h"
#include "utils.hpp"

/*
 https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
 
 AminoAcids : A,R,N,D,B,C,Q,E,Z,G,H,I,L,K,M,F,P,S,T,W,Y,V
 
 AminoAcid      DNA codons                  Compressed
 
 Ala, A         GCT,GCC,GCA,GCG             GCN
 Arg, R         CGT,CGC,CGA,CGG; AGA,AGG    CGN,AGR; or CGY,MGR
 Asn, N         AAT,AAC                     AAY
 Asp, D         GAT,GAC                     GAY
 Asn or Asp, B  AAT,AAC; GAT,GAC            RAY
 Cys, C         TGT,TGC                     TGY
 Gln, Q         CAA,CAG                     CAR
 Glu, E         GAA,GAG                     GAR
 Gln or Glu, Z  CAA,CAG; GAA,GAG            SAR
 Gly, G         GGT,GGC,GGA,GGG             GGN
 His, H         CAT,CAC                     CAY
 Ile, I         ATT,ATC,ATA                 ATH
 Leu, L         CTT,CTC,CTA,CTG; TTA,TTG    CTN,TTR; or CTY,YTR
 Lys, K         AAA,AAG                     AAR
 Met, M         ATG                         ATG
 Phe, F         TTT,TTC                     TTY
 Pro, P         CCT,CCC,CCA,CCG             CCN
 Ser, S         TCT,TCC,TCA,TCG; AGT,AGC    TCN,AGY
 Thr, T         ACT,ACC,ACA,ACG             ACN
 Trp, W         TGG                         TGG
 Tyr, Y         TAT,TAC                     TAY
 Val, V         GTT,GTC,GTA,GTG             GTN
 
 START          ATG
 STOP           TAA,TGA,TAG                 TRA,TAR
 
 . and -        ./-           gap
*/

typedef enum {
  ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro,
  ser1, ser2, thr, trp, tyr, val, del, stop, asx, glx, ser, unk, quest
} aas;

/**************************************************************************************************/

class AminoAcid {
    
public:
    AminoAcid();
    AminoAcid(char);    //AminoAcid character
    AminoAcid(string);    //dna codon length of 3
    ~AminoAcid() = default;
    
    string getName();
    char getAmino()          { return aminoBase;    }
    int getNum()             { return aminoNum;     }
    void setAmino(char c);
        
protected:
    
    MothurOut* m;
    Utils util;
    
    char aminoBase;
    int aminoNum;
    
    map<char, int> indexes; //A -> 0, T -> 1, G -> 2, C -> 3, N -> 4, Gap -> 5
    map<char, int>::iterator it;

    char findAmino(string);
    char getAminoBase(string); //from name
    int getIndex(char x);
    
    
};

#endif /* aminoacid_hpp */
