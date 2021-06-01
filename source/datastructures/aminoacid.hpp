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
 
 
*/

/**************************************************************************************************/

class AminoAcid {
    
public:
    AminoAcid(char);    //AminoAcid character
    //AminoAcid(string);  //DNA codon
    ~AminoAcid() {}
    
    string getName();
    char get()          { return aminoBase; }
    
    //vector<string> getCodons;
    
protected:
    
    MothurOut* m;
    Utils util;
    
    char aminoBase;
    map<char, string> aminoNameMap;
    map<string, char> codonMap;
    
    void fillNameMap();
    void fillCodons();
    char getAminoBase(string);
    
};

#endif /* aminoacid_hpp */
