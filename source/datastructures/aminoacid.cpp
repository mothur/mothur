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
        fillValidAminoAcid();
        aminoBase = '0';
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
        fillValidAminoAcid();
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
        
        if (codon.length() != 3) { m->mothurOut("[ERROR]: " + codon + " is not the correct length. Codons must be 3 characters long, quitting.\n"); m->setControl_pressed(true); return; }
        
        fillValidAminoAcid();
        fillCodonsMap();
        
        int index1 = -1; int index2 = -1; int index3 = -1;
        it = indexes.find(codon[0]); if (it != indexes.end()) { index1 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[0]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return; }
        it = indexes.find(codon[1]); if (it != indexes.end()) { index2 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[1]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return; }
        it = indexes.find(codon[2]); if (it != indexes.end()) { index3 = it->second; } else { m->mothurOut("[ERROR]: " + toString(codon[2]) + " is not A, T, G, C, or N, quitting.\n"); m->setControl_pressed(true); return; }
        
        //if no N's then the set should contain one amino acid. if N's, then try all possible values for N in that position.
        //for example:ACN -> Threonine (T) because ACA,ACT,ACG,ACC all map to Threonine
        //   but      GAN -> could be Glutamate (E) (for N=A or G) or Aspartate (D) (for N=T or C)
        if ((index1 == 4) || (index2 == 4) || (index3 == 4)) { //any position of the codon is an N
            set<char> possibleAminoAcids;
            
            if (((index1 == 4) && (index2 == 4)) || ((index1 == 4) && (index3 == 4)) || ((index3 == 4) && (index2 == 4))) { //2 N's in codon
                setAmino('?');
            }else{ //only 1 N
                if (index1 == 4) {
                    possibleAminoAcids.insert(codons[0][index2][index3]); possibleAminoAcids.insert(codons[1][index2][index3]);
                    possibleAminoAcids.insert(codons[2][index2][index3]); possibleAminoAcids.insert(codons[3][index2][index3]);
                }else if (index2 == 4) {
                    possibleAminoAcids.insert(codons[index1][0][index3]); possibleAminoAcids.insert(codons[index1][1][index3]);
                    possibleAminoAcids.insert(codons[index1][2][index3]); possibleAminoAcids.insert(codons[index1][3][index3]);
                }else {
                    possibleAminoAcids.insert(codons[index1][index2][0]); possibleAminoAcids.insert(codons[index1][index2][1]);
                    possibleAminoAcids.insert(codons[index1][index2][2]); possibleAminoAcids.insert(codons[index1][index2][3]);
                }
                if (possibleAminoAcids.size() == 1) {
                    setAmino(*possibleAminoAcids.begin());
                }else { setAmino('?'); }
            }
        }else {
            char c = codons[index1][index2][index3];
            setAmino(c);
        }
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
        
        if (validAminoAcids.count(c) != 0) {
            aminoBase = c;
            getName(); //sets name, number and compressed dna
        }else { m->mothurOut("[ERROR]: " + toString(c) + " is an invalid amino acid, please correct.\n"); m->setControl_pressed(true); }
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "setAmino");
        exit(1);
    }
}
//******************************************************************************************************************
/*
 Codons are 3 characters long. The codons vector is [4][4][4].
 codons[0][0][0] -> AAA -> Lysine (K) -> K.
 codons[0][0][1] -> AAT -> Asparagine (N) -> N.
 */
void AminoAcid::fillCodonsMap() {
    try {
        codons.clear(); codons.resize(4);
        for (int i = 0; i < codons.size(); i++) {
            codons[i].resize(4);
            for (int j = 0; j < codons[i].size(); j++) {
                codons[i][j].resize(4);
            }
        }
        
        indexes['A'] = 0; indexes['T'] = 1; indexes['G'] = 2; indexes['C'] = 3; indexes['N'] = 4;
        
        //AAX
        codons[0][0][0] = 'K';   //AAA |  Lysine (K) -> 11. where 11 is the index into the aas enum.
        codons[0][0][1] = 'N';   //AAT |  Asparagine (N) -> 2.
        codons[0][0][2] = 'K';   //AAG |  Lysine (K) -> 11.
        codons[0][0][3] = 'N';   //AAC |  Asparagine (N) -> 2.
        
        //ATX
        codons[0][1][0] = 'I';   //ATA |  Isoleucine (I) -> 9.
        codons[0][1][1] = 'I';   //ATT |  Isoleucine (I) -> 9.
        codons[0][1][2] = 'M';   //ATG |  Methionine (M) -> 12.
        codons[0][1][3] = 'I';   //ATC |  Isoleucine (I) -> 9.
        
        //AGX
        codons[0][2][0] = 'R';   //AGA |  Arginine (R) -> 1.
        codons[0][2][1] = 'S';   //AGT |  Serine (S) -> 15.
        codons[0][2][2] = 'R';   //AGG |  Arginine (R) -> 1.
        codons[0][2][3] = 'S';   //AGC |  Serine (S) -> 15.
        
        //ACX
        codons[0][3][0] = 'T';    //ACA |  Threonine (T) -> 17.
        codons[0][3][1] = 'T';    //ACT |  Threonine (T) -> 17.
        codons[0][3][2] = 'T';    //ACG |  Threonine (T) -> 17.
        codons[0][3][3] = 'T';    //ACC |  Threonine (T) -> 17.
        
        
        //TAX
        codons[1][0][0] = '*';   //TAA | Termination (X) -> 22
        codons[1][0][1] = 'Y';   //TAT | Tyrosine (Y) -> 19
        codons[1][0][2] = '*';   //TAG | Termination (X) -> 22
        codons[1][0][3] = 'Y';   //TAC | Tyrosine (Y) -> 19
        
        //TTX
        codons[1][1][0] = 'L';    //TTA | Leucine (L) -> 10
        codons[1][1][1] = 'F';    //TTT | Phenylalanine (F) -> 13
        codons[1][1][2] = 'L';    //TTG | Leucine (L) -> 10
        codons[1][1][3] = 'F';    //TTC | Phenylalanine (F) -> 13
        
        //TGX
        codons[1][2][0] = '*';    //TGA | Termination (X) -> 22
        codons[1][2][1] = 'C';    //TGT | Cysteine (C) -> 4
        codons[1][2][2] = 'W';    //TGG | Tryptophan (W) -> 18
        codons[1][2][3] = 'C';    //TGC | Cysteine (C) -> 4
        
        //TCX
        codons[1][3][0] = 'S';    //TCA | Serine (S) -> 15
        codons[1][3][1] = 'S';    //TCT | Serine (S) -> 15
        codons[1][3][2] = 'S';    //TCG | Serine (S) -> 15
        codons[1][3][3] = 'S';    //TCC | Serine (S) -> 15
        
        //GAX
        codons[2][0][0] = 'E';   //GAA | Glutamate (E) -> 6
        codons[2][0][1] = 'D';   //GAT | Aspartate (D) -> 3
        codons[2][0][2] = 'E';   //GAG | Glutamate (E) -> 6
        codons[2][0][3] = 'D';   //GAC | Aspartate (D) -> 3
        
        //GTX
        codons[2][1][0] = 'V';    //GTA | Valine (V)
        codons[2][1][1] = 'V';    //GTT | Valine (V)
        codons[2][1][2] = 'V';    //GTG | Valine (V)
        codons[2][1][3] = 'V';    //GTC | Valine (V)
        
        //GGX
        codons[2][2][0] = 'G';    //GGA | Glycine (G)
        codons[2][2][1] = 'G';    //GGT | Glycine (G)
        codons[2][2][2] = 'G';    //GGG | Glycine (G)
        codons[2][2][3] = 'G';    //GGC | Glycine (G)
        
        //GCX
        codons[2][3][0] = 'A';    //GCA | Alanine (A)
        codons[2][3][1] = 'A';    //GCT | Alanine (A)
        codons[2][3][2] = 'A';    //GCG | Alanine (A)
        codons[2][3][3] = 'A';    //GCC | Alanine (A)
        
        //CAX
        codons[3][0][0] = 'Q';   //CAA | Glutamine (Q)
        codons[3][0][1] = 'H';   //CAT | Histidine (H)
        codons[3][0][2] = 'Q';   //CAG | Glutamine (Q)
        codons[3][0][3] = 'H';   //CAC | Histidine (H)
        
        //CTX
        codons[3][1][0] = 'L';    //CTA | Leucine (L)
        codons[3][1][1] = 'L';    //CTT | Leucine (L)
        codons[3][1][2] = 'L';    //CTG | Leucine (L)
        codons[3][1][3] = 'L';    //CTC | Leucine (L)
        
        //CGX
        codons[3][2][0] = 'R';    //CGA | Arginine (R)
        codons[3][2][1] = 'R';    //CGT | Arginine (R)
        codons[3][2][2] = 'R';    //CGG | Arginine (R)
        codons[3][2][3] = 'R';    //CGC | Arginine (R)
        
        //CCX
        codons[3][3][0] = 'P';    //CCA | Proline (P)
        codons[3][3][1] = 'P';    //CCT | Proline (P)
        codons[3][3][2] = 'P';    //CCG | Proline (P)
        codons[3][3][3] = 'P';    //CCC | Proline (P)
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "fillCodons");
        exit(1);
    }
}
 /******************************************************************************************************************/
//ala(0), arg(1), asn(2), asp(3), cys(4), gln(5), glu(6), gly(7), his(8), ileu(9), leu(10), lys(11), met(12), phe(13), pro(14),
//ser1(15), ser2(16), thr(17), trp(18), tyr(19), val(20), del(21), stop(22), asx(23), glx(24), ser(25), unk(26), quest(27)
 string AminoAcid::getName() {
     try {
         string aminoName = "unknown"; aminoNum = unk; dna.clear();
         
         if (aminoBase == 'A')          { aminoName = "Alanine";        aminoNum = ala; dna.push_back("GCN"); } //0
         else if (aminoBase == 'R')     { aminoName = "Arginine";       aminoNum = arg; dna.push_back("CGN"); dna.push_back("AGR"); dna.push_back("CGY"); dna.push_back("MGR"); } //1
         else if (aminoBase == 'N')     { aminoName = "Asparagine";     aminoNum = asn; dna.push_back("AAY"); } //2
         else if (aminoBase == 'D')     { aminoName = "Aspartic";       aminoNum = asp; dna.push_back("GAY"); } //3
         
         else if (aminoBase == 'B')     { aminoName = "Asparagine or Aspartic"; aminoNum = asx; dna.push_back("RAY"); } //23
         else if (aminoBase == 'C')     { aminoName = "Cysteine";               aminoNum = cys; dna.push_back("TGY"); } //4
         else if (aminoBase == 'Q')     { aminoName = "Glutamine";              aminoNum = gln; dna.push_back("CAR"); } //5
         else if (aminoBase == 'E')     { aminoName = "Glutamic";               aminoNum = glu; dna.push_back("GAR"); } //6
         
         else if (aminoBase == 'Z')     { aminoName = "Glutamine or Glutamic_Acid"; aminoNum = glx; dna.push_back("SAR"); } //24
         else if (aminoBase == 'G')     { aminoName = "Glycine";        aminoNum = gly;  dna.push_back("GGN");   } //7
         else if (aminoBase == 'H')     { aminoName = "Histidine";      aminoNum = his;  dna.push_back("CAY");   } //8
         else if (aminoBase == 'I')     { aminoName = "Isoleucine";     aminoNum = ileu; dna.push_back("ATH");   } //9
         
         else if (aminoBase == 'L')     { aminoName = "Leucine";        aminoNum = leu;  dna.push_back("CTN"); dna.push_back("TTR"); dna.push_back("CTY"); dna.push_back("YTR");  } //10
         else if (aminoBase == 'K')     { aminoName = "Lysine";         aminoNum = lys;  dna.push_back("AAR");   } //11
         else if (aminoBase == 'M')     { aminoName = "Methionine";     aminoNum = met;  dna.push_back("ATG");   } //12
         else if (aminoBase == 'F')     { aminoName = "Phenylalanine";  aminoNum = phe;   dna.push_back("TTY");  } //13
         
         else if (aminoBase == 'P')     { aminoName = "Proline";        aminoNum = pro;   dna.push_back("CCN");    } //14
         else if (aminoBase == 'S')     { aminoName = "Serine";         aminoNum = ser1;  dna.push_back("TCN"); dna.push_back("AGY");  } //15
         else if (aminoBase == 'T')     { aminoName = "Threonine";      aminoNum = thr;   dna.push_back("ACN");    } //17
         else if (aminoBase == 'W')     { aminoName = "Tryptophan";     aminoNum = trp;   dna.push_back("TGG");   } //18
         
         else if (aminoBase == 'Y')     { aminoName = "Tyrosine";       aminoNum = tyr;  dna.push_back("TAY");     } //19
         else if (aminoBase == 'V')     { aminoName = "Valine";         aminoNum = val;   dna.push_back("GTN");    } //20
         else if ((aminoBase == '.') || (aminoBase == '-'))     { aminoName = "Gap"; aminoNum = del; dna.push_back("---");    } //21
         else if ((aminoBase == '*') || (aminoBase == 'X'))     { aminoName = "STOP";           aminoNum = stop;   dna.push_back("TRA"); dna.push_back("TAR");   } //22
         else if (aminoBase == '?')     { aminoName = "QUESTION";       aminoNum = quest; dna.push_back("NNN");    } //27
         
         return aminoName;
     }
     catch(exception& e) {
         m->errorOut(e, "AminoAcid", "getName");
         exit(1);
     }
 }
/******************************************************************************************************************/
void AminoAcid::fillValidAminoAcid() {
    try {
        validAminoAcids.insert('A');
        validAminoAcids.insert('R');
        validAminoAcids.insert('N');
        validAminoAcids.insert('D');
        
        validAminoAcids.insert('B');
        validAminoAcids.insert('C');
        validAminoAcids.insert('Q');
        validAminoAcids.insert('E');
        
        validAminoAcids.insert('Z');
        validAminoAcids.insert('G');
        validAminoAcids.insert('H');
        validAminoAcids.insert('I');
        
        validAminoAcids.insert('L');
        validAminoAcids.insert('K');
        validAminoAcids.insert('M');
        validAminoAcids.insert('F');
        
        validAminoAcids.insert('P');
        validAminoAcids.insert('S');
        validAminoAcids.insert('T');
        validAminoAcids.insert('W');
        
        validAminoAcids.insert('Y');
        validAminoAcids.insert('V');
        validAminoAcids.insert('X');
        validAminoAcids.insert('-');
        validAminoAcids.insert('.');
        
        validAminoAcids.insert('*');
        validAminoAcids.insert('?');
        
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "fillValidAminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************/
