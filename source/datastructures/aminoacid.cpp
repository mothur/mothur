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
void AminoAcid::setAmino(char c) {
    try {
        
        c = toupper(c);
        
        if (validAminoAcids.count(c) != 0) { aminoBase = c; getName(); }
        else {
            m->mothurOut("[ERROR]: " + toString(c) + " is an invalid amino acid, please correct.\n"); m->setControl_pressed(true);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "setAmino");
        exit(1);
    }
}

/******************************************************************************************************************
void AminoAcid::fillCodons() {
    try {
        //Ala, A         GCT,GCC,GCA,GCG
        codonMap["GCT"] = 'A';
        codonMap["GCC"] = 'A';
        codonMap["GCA"] = 'A';
        codonMap["GCG"] = 'A';
        
        //Arg, R         CGT,CGC,CGA,CGG; AGA,AGG
        codonMap["CGT"] = 'R';
        codonMap["CGC"] = 'R';
        codonMap["CGA"] = 'R';
        codonMap["CGG"] = 'R';
        codonMap["AGA"] = 'R';
        codonMap["AGG"] = 'R';
        
        //Asn, N         AAT,AAC
        codonMap["AAT"] = 'N';
        codonMap["AAC"] = 'N';
        
        //Asp, D         GAT,GAC
        codonMap["GAT"] = 'D';
        codonMap["GAC"] = 'D';
        
        //Cys, C         TGT,TGC
        codonMap["TGT"] = 'C';
        codonMap["TGC"] = 'C';
        
        //Gln, Q         CAA,CAG
        codonMap["CAA"] = 'Q';
        codonMap["CAG"] = 'Q';
        
        //Glu, E         GAA,GAG
        codonMap["GAA"] = 'E';
        codonMap["GAG"] = 'E';
        
        //Gly, G         GGT,GGC,GGA,GGG
        codonMap["GGT"] = 'G';
        codonMap["GGC"] = 'G';
        codonMap["GGA"] = 'G';
        codonMap["GGG"] = 'G';
        
        //His, H         CAT,CAC
        codonMap["CAT"] = 'H';
        codonMap["CAC"] = 'H';
        
        //Ile, I         ATT,ATC,ATA
        codonMap["ATT"] = 'I';
        codonMap["ATC"] = 'I';
        codonMap["ATA"] = 'I';
        
        //Leu, L         CTT,CTC,CTA,CTG; TTA,TTG
        codonMap["CTT"] = 'L';
        codonMap["CTC"] = 'L';
        codonMap["CTA"] = 'L';
        codonMap["CTG"] = 'L';
        codonMap["TTA"] = 'L';
        codonMap["TTG"] = 'L';
        
        //Lys, K         AAA,AAG
        codonMap["AAA"] = 'K';
        codonMap["AAG"] = 'K';
        
        //Met, M         ATG
        codonMap["ATG"] = 'M';
        
        //Phe, F         TTT,TTC
        codonMap["TTT"] = 'F';
        codonMap["TTC"] = 'F';
        
        //Pro, P         CCT,CCC,CCA,CCG
        codonMap["CCT"] = 'P';
        codonMap["CCC"] = 'P';
        codonMap["CCA"] = 'P';
        codonMap["CCG"] = 'P';
        
        //Ser, S         TCT,TCC,TCA,TCG; AGT,AGC
        codonMap["TCT"] = 'S';
        codonMap["TCC"] = 'S';
        codonMap["TCA"] = 'S';
        codonMap["TCG"] = 'S';
        codonMap["AGT"] = 'S';
        codonMap["AGC"] = 'S';
        
        //Thr, T         ACT,ACC,ACA,ACG
        codonMap["ACT"] = 'T';
        codonMap["ACC"] = 'T';
        codonMap["ACA"] = 'T';
        codonMap["ACG"] = 'T';
        
        //Trp, W         TGG
        codonMap["TGG"] = 'W';
        
        //Tyr, Y         TAT,TAC
        codonMap["TAT"] = 'Y';
        codonMap["TAC"] = 'Y';
        
        //Val, V         GTT,GTC,GTA,GTG
        codonMap["GTT"] = 'V';
        codonMap["GTC"] = 'V';
        codonMap["GTA"] = 'V';
        codonMap["GTG"] = 'V';
        
        //TODO::resolve codons assigned to multiple aminoacids
        //TODO::start and stop ???
        
        //Gln or Glu, Z  CAA,CAG; GAA,GAG
        codonMap["CAA"] = 'Z';
        codonMap["CAG"] = 'Z';
        codonMap["GAA"] = 'Z';
        codonMap["GAG"] = 'Z';
        
        //Asn or Asp, B AAT,AAC; GAT,GAC
        codonMap["AAT"] = 'B';
        codonMap["AAC"] = 'B';
        codonMap["GAT"] = 'B';
        codonMap["GAC"] = 'B';
        condonMap["."] = '.';
        condonMap["-"] = '-';
        
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "fillCodons");
        exit(1);
    }
}
 /******************************************************************************************************************/
//TODO::resolve multi names issue
//TODO::start and stop ???
//ala(0), arg(1), asn(2), asp(3), cys(4), gln(5), glu(6), gly(7), his(8), ileu(9), leu(10), lys(11), met(12), phe(13), pro(14),
//ser1(15), ser2(16), thr(17), trp(18), tyr(19), val(20), del(21), stop(22), asx(23), glx(24), ser(25), unk(26), quest(27)
 string AminoAcid::getName() {
     try {
         string aminoName = "unknown"; aminoNum = unk;
         
         if (aminoBase == 'A')          { aminoName = "Alanine";        aminoNum = ala; } //0
         else if (aminoBase == 'R')     { aminoName = "Arginine";       aminoNum = arg; } //1
         else if (aminoBase == 'N')     { aminoName = "Asparagine";     aminoNum = asn; } //2
         else if (aminoBase == 'D')     { aminoName = "Aspartic";       aminoNum = asp; } //3
         
         else if (aminoBase == 'B')     { aminoName = "Asparagine or Aspartic"; aminoNum = asx; } //23
         else if (aminoBase == 'C')     { aminoName = "Cysteine";               aminoNum = cys; } //4
         else if (aminoBase == 'Q')     { aminoName = "Glutamine";              aminoNum = gln; }
         else if (aminoBase == 'E')     { aminoName = "Glutamic";               aminoNum = glu; }
         
         else if (aminoBase == 'Z')     { aminoName = "Glutamine or Glutamic_Acid"; aminoNum = glx; } //24
         else if (aminoBase == 'G')     { aminoName = "Glycine";        aminoNum = gly;     }
         else if (aminoBase == 'H')     { aminoName = "Histidine";      aminoNum = his;     }
         else if (aminoBase == 'I')     { aminoName = "Isoleucine";     aminoNum = ileu;    }
         
         else if (aminoBase == 'L')     { aminoName = "Leucine";        aminoNum = leu;     }
         else if (aminoBase == 'K')     { aminoName = "Lysine";         aminoNum = lys;     }
         else if (aminoBase == 'M')     { aminoName = "Methionine";     aminoNum = met;     }
         else if (aminoBase == 'F')     { aminoName = "Phenylalanine";  aminoNum = phe;     }
         
         else if (aminoBase == 'P')     { aminoName = "Proline";        aminoNum = pro;     }
         else if (aminoBase == 'S')     { aminoName = "Serine";         aminoNum = ser1;    }
         else if (aminoBase == 'T')     { aminoName = "Threonine";      aminoNum = thr;     }
         else if (aminoBase == 'W')     { aminoName = "Tryptophan";     aminoNum = trp;     }
         
         else if (aminoBase == 'Y')     { aminoName = "Tyrosine";       aminoNum = tyr;     }
         else if (aminoBase == 'V')     { aminoName = "Valine";         aminoNum = val;     }
         else if ((aminoBase == '.') || (aminoBase == '-'))     { aminoName = "Gap"; aminoNum = del;  } //21
         else if (aminoBase == '*')     { aminoName = "STOP";           aminoNum = stop;    } //22
         else if (aminoBase == '?')     { aminoName = "QUESTION";       aminoNum = quest;   } //27
         
         
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
