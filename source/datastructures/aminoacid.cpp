//
//  codon.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/24/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "aminoacid.hpp"

/******************************************************************************************************************/
AminoAcid::AminoAcid(char c) {
    try {
        m = MothurOut::getInstance();
        
        fillNameMap();
        fillCodons();
        
        aminoBase = toupper(c);
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "AminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************
AminoAcid::AminoAcid(string codon) {
    try {
        m = MothurOut::getInstance();
    
        fillNameMap();
        fillCodons();
        
        //trim in case its not a string of 3
        codon = codon.substr(3);
        
        aminoBase = getAminoBase(codon);
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "AminoAcid");
        exit(1);
    }
}
/******************************************************************************************************************/
char AminoAcid::getAminoBase(string codon) {
    try {
        char base = ' ';
        
        map<string, char>::iterator it = codonMap.find(codon);
        if (it != codonMap.end()) {
            base = it->second;
        }else {
            m->mothurOut("[ERROR]: " + codon + " is not a valid codon, please correct.\n"); m->setControl_pressed(true);
        }
        
        return base;
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "getAminoBase");
        exit(1);
    }
}
/******************************************************************************************************************/
string AminoAcid::getName() {
    try {
        string aminoName = "unknown";
        
        map<char, string>::iterator it = aminoNameMap.find(aminoBase);
        if (it != aminoNameMap.end()) {
            aminoName = it->second;
        }else {
            m->mothurOut("[ERROR]: " + toString(aminoBase) + " is not an amino acid base, please correct.\n"); m->setControl_pressed(true);
        }
        
        return aminoName;
        
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "getName");
        exit(1);
    }
}
/******************************************************************************************************************/
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
        
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "fillCodons");
        exit(1);
    }
}
/******************************************************************************************************************/
void AminoAcid::fillNameMap() {
    try {
        aminoNameMap['A'] = "Alanine";
        aminoNameMap['R'] = "Arginine";
        aminoNameMap['N'] = "Asparagine";
        aminoNameMap['D'] = "Aspartic_acid";
        
    //TODO::resolve multi names issue
    //TODO::start and stop ???
        
        aminoNameMap['B'] = "Asparagine or Aspartic";
        aminoNameMap['C'] = "Cysteine";
        aminoNameMap['Q'] = "Glutamine";
        aminoNameMap['E'] = "Glutamic_Acid";
        
        aminoNameMap['Z'] = "Glutamine or Glutamic_Acid";
        aminoNameMap['G'] = "Glycine";
        aminoNameMap['H'] = "Histidine";
        aminoNameMap['I'] = "Isoleucine";
        
        aminoNameMap['L'] = "Leucine";
        aminoNameMap['K'] = "Lysine";
        aminoNameMap['M'] = "Methionine";
        aminoNameMap['F'] = "Phenylalanine";
        
        aminoNameMap['P'] = "Proline";
        aminoNameMap['S'] = "Serine";
        aminoNameMap['T'] = "Threonine";
        aminoNameMap['W'] = "Tryptophan";
        
        aminoNameMap['Y'] = "Tyrosine";
        aminoNameMap['V'] = "Valine";
        
        
    }
    catch(exception& e) {
        m->errorOut(e, "AminoAcid", "fillNameMap");
        exit(1);
    }
}
/******************************************************************************************************************/
