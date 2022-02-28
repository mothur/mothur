//
//  codon.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/24/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "aminoacid.hpp"

/******************************************************************************************************************/
AminoAcid::AminoAcid() {
    try {
        m = MothurOut::getInstance();
        
        indexes['A'] = 0; indexes['T'] = 1; indexes['G'] = 2; indexes['C'] = 3; indexes['N'] = 4; indexes['-'] = 5; indexes['.'] = 5;

        setAmino('?');
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "AminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************/
AminoAcid::AminoAcid(char c) {
    try {
        m = MothurOut::getInstance();
        
        indexes['A'] = 0; indexes['T'] = 1; indexes['G'] = 2; indexes['C'] = 3; indexes['N'] = 4; indexes['-'] = 5; indexes['.'] = 5;

        setAmino(c);
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "AminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************/
//requires the codon to be 3 characters long. Only valid characters are a,t,g,c,n.
AminoAcid::AminoAcid(string codon) {
    try {
        m = MothurOut::getInstance();
        
        indexes['A'] = 0; indexes['T'] = 1; indexes['G'] = 2; indexes['C'] = 3; indexes['N'] = 4; indexes['-'] = 5; indexes['.'] = 5;

        char amino = findAmino(codon);
        
        setAmino(amino);
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "AminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************/
void AminoAcid::setAmino(char c) {
    try {
        
        c = toupper(c);
        
        if (m->validAminoAcids.count(c) != 0) {
            aminoBase = c;
            getName(); //sets name, number and compressed dna
        }else { m->mothurOut("[ERROR]: " + toString(c) + " is an invalid amino acid, please correct.\n"); m->setControl_pressed(true); }
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "setAmino");
        exit(1);
    }
}
/******************************************************************************************************************/
char AminoAcid::findAmino(string codon) {
    try {
        
        char amino = '?';
        if (codon.length() != 3) { m->mothurOut("[ERROR]: " + codon + " is not the correct length. Codons must be 3 characters long, quitting.\n"); m->setControl_pressed(true); return amino; }
        
        
        int index1 = -1; int index2 = -1; int index3 = -1;
        it = indexes.find(codon[0]); if (it != indexes.end()) { index1 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[0]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return amino; }
        it = indexes.find(codon[1]); if (it != indexes.end()) { index2 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[1]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return amino; }
        it = indexes.find(codon[2]); if (it != indexes.end()) { index3 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[2]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return amino; }
        
        if ((index1 == 5) && (index2 == 5) && (index3 == 5)) { amino = '-';  return amino; }

        //if no N's then the set should contain one amino acid. if N's, then try all possible values for N in that position.
        //for example:ACN -> Threonine (T) because ACA,ACT,ACG,ACC all map to Threonine
        //   but      GAN -> could be Glutamate (E) (for N=A or G) or Aspartate (D) (for N=T or C)
        if ((index1 > 3) || (index2 > 3) || (index3 > 3)) { //any position of the codon is an N or gap
            set<char> possibleAminoAcids;
            
            if (((index1 > 3) && (index2 > 3)) || ((index1 > 3) && (index3 > 3)) || ((index3 > 3) && (index2 > 3))) { //2 N's or gaps in codon
            }else{ //only 1 N
                if (index1 > 3) {
                    possibleAminoAcids.insert(m->codons[0][index2][index3]); possibleAminoAcids.insert(m->codons[1][index2][index3]);
                    possibleAminoAcids.insert(m->codons[2][index2][index3]); possibleAminoAcids.insert(m->codons[3][index2][index3]);
                }else if (index2 > 3) {
                    possibleAminoAcids.insert(m->codons[index1][0][index3]); possibleAminoAcids.insert(m->codons[index1][1][index3]);
                    possibleAminoAcids.insert(m->codons[index1][2][index3]); possibleAminoAcids.insert(m->codons[index1][3][index3]);
                }else {
                    possibleAminoAcids.insert(m->codons[index1][index2][0]); possibleAminoAcids.insert(m->codons[index1][index2][1]);
                    possibleAminoAcids.insert(m->codons[index1][index2][2]); possibleAminoAcids.insert(m->codons[index1][index2][3]);
                }
                if (possibleAminoAcids.size() == 1) {
                    amino = (*possibleAminoAcids.begin());
                }
            }
        }else {
            amino = m->codons[index1][index2][index3];
        }
        
        return amino;
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "findAmino");
        exit(1);
    }
}
 /******************************************************************************************************************/
//ala(0), arg(1), asn(2), asp(3), cys(4), gln(5), glu(6), gly(7), his(8), ileu(9), leu(10), lys(11), met(12), phe(13), pro(14),
//ser1(15), ser2(16), thr(17), trp(18), tyr(19), val(20), del(21), stop(22), asx(23), glx(24), ser(25), unk(26), quest(27)
 string AminoAcid::getName() {
     try {
         string aminoName = "unknown"; aminoNum = unk;
         
         if (aminoBase == 'A')          { aminoName = "Alanine";        aminoNum = ala;  } //0
         else if (aminoBase == 'R')     { aminoName = "Arginine";       aminoNum = arg;  } //1
         else if (aminoBase == 'N')     { aminoName = "Asparagine";     aminoNum = asn; } //2
         else if (aminoBase == 'D')     { aminoName = "Aspartic";       aminoNum = asp;  } //3
         
         else if (aminoBase == 'B')     { aminoName = "Asparagine or Aspartic"; aminoNum = asx;  } //23
         else if (aminoBase == 'C')     { aminoName = "Cysteine";               aminoNum = cys;  } //4
         else if (aminoBase == 'Q')     { aminoName = "Glutamine";              aminoNum = gln;  } //5
         else if (aminoBase == 'E')     { aminoName = "Glutamic";               aminoNum = glu;  } //6
         
         else if (aminoBase == 'Z')     { aminoName = "Glutamine or Glutamic_Acid"; aminoNum = glx; } //24
         else if (aminoBase == 'G')     { aminoName = "Glycine";        aminoNum = gly;     } //7
         else if (aminoBase == 'H')     { aminoName = "Histidine";      aminoNum = his;    } //8
         else if (aminoBase == 'I')     { aminoName = "Isoleucine";     aminoNum = ileu;    } //9
         
         else if (aminoBase == 'L')     { aminoName = "Leucine";        aminoNum = leu;   } //10
         else if (aminoBase == 'K')     { aminoName = "Lysine";         aminoNum = lys;    } //11
         else if (aminoBase == 'M')     { aminoName = "Methionine";     aminoNum = met;    } //12
         else if (aminoBase == 'F')     { aminoName = "Phenylalanine";  aminoNum = phe;     } //13
         
         else if (aminoBase == 'P')     { aminoName = "Proline";        aminoNum = pro;      } //14
         else if (aminoBase == 'S')     { aminoName = "Serine";         aminoNum = ser1;    } //15
         else if (aminoBase == 'T')     { aminoName = "Threonine";      aminoNum = thr;       } //17
         else if (aminoBase == 'W')     { aminoName = "Tryptophan";     aminoNum = trp;     } //18
         
         else if (aminoBase == 'Y')     { aminoName = "Tyrosine";       aminoNum = tyr;       } //19
         else if (aminoBase == 'V')     { aminoName = "Valine";         aminoNum = val;     } //20
         else if ((aminoBase == '.') || (aminoBase == '-'))     { aminoName = "Gap"; aminoNum = del;    } //21
         else if ((aminoBase == '*') || (aminoBase == 'X'))     { aminoName = "STOP";           aminoNum = stop;     } //22
         else if (aminoBase == '?')     { aminoName = "QUESTION";       aminoNum = quest;   } //27
         
         return aminoName;
     }
     catch(exception& e) {
         m->errorOut(e, "AminoAcid", "getName");
         exit(1);
     }
 }
/******************************************************************************************************************/
