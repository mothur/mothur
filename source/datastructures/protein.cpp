//
//  protein.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/24/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "protein.hpp"

/***********************************************************************/
Protein::Protein(){
    m = MothurOut::getInstance();
    initialize();
}
/***********************************************************************/
Protein::Protein(string newName, vector<AminoAcid> sequence) {
    try {
        m = MothurOut::getInstance();
        initialize();
        name = newName;
        
        util.checkName(name);
        
        setUnaligned(sequence); //setUnaligned removes any gap characters for us
        setAligned(sequence);
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}
/***********************************************************************/
Protein::Protein(string newName, string seq) {
    try {
        m = MothurOut::getInstance();
        initialize();
        name = newName;
        
        util.checkName(name);
        
        vector<AminoAcid> sequence;
        for (int i = 0; i < seq.size(); i++) {
            AminoAcid temp(seq[i]);
            sequence.push_back(temp);
        }
        
        setUnaligned(sequence); //setUnaligned removes any gap characters for us
        setAligned(sequence);
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Protein::Protein(istringstream& fastaString){
    try {
        m = MothurOut::getInstance();
    
        initialize();
        name = getProteinName(fastaString);
        
        if (!m->getControl_pressed()) {
            string proteinComment;
        
            //read comments
            while ((name[0] == '#') && fastaString) {
                while (!fastaString.eof())    {    char c = fastaString.get(); if (c == 10 || c == 13){    break;    }    } // get rest of line if there's any crap there
                proteinComment = getCommentString(fastaString);
                
                if (fastaString) {
                    fastaString >> name;
                    name = name.substr(1);
                }else {
                    name = "";
                    break;
                }
            }
            
            comment = getCommentString(fastaString);
            vector<AminoAcid> proteinSeq = getProtein(fastaString);
            
            setAligned(proteinSeq);
            setUnaligned(proteinSeq); //setUnaligned removes any gap characters for us
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}

//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Protein::Protein(ifstream& fastaFile){
    try {
        m = MothurOut::getInstance();
        initialize();
        name = getProteinName(fastaFile);
        
        if (!m->getControl_pressed()) {
            
            string proteinComment;
        
            //read comments
            while ((name[0] == '#') && fastaFile) {
                while (!fastaFile.eof())    {    char c = fastaFile.get(); if (c == 10 || c == 13){    break;    }    } // get rest of line if there's any crap there
                proteinComment = getCommentString(fastaFile);
                
                if (fastaFile) {
                    fastaFile >> name;
                    name = name.substr(1);
                }else {
                    name = "";
                    break;
                }
            }
            
            //while (!fastaFile.eof())    {    char c = fastaFile.get(); if (c == 10 || c == 13){  break;    }    } // get rest of line if there's any crap there
            comment = getCommentString(fastaFile);
            vector<AminoAcid> proteinSeq = getProtein(fastaFile);
            
            setAligned(proteinSeq);
            setUnaligned(proteinSeq); //setUnaligned removes any gap characters for us
        }

    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
#ifdef USE_BOOST
Protein::Protein(boost::iostreams::filtering_istream& fastaFile){
    try {
        m = MothurOut::getInstance();
        initialize();
        name = getSequenceName(fastaFile);
        
        if (!m->getControl_pressed()) {
            
            string sequence;
            
            //read comments
            while ((name[0] == '#') && fastaFile) {
                while (!fastaFile.eof())    {    char c = fastaFile.get(); if (c == 10 || c == 13){    break;    }    } // get rest of line if there's any crap there
                sequence = getCommentString(fastaFile);
                
                if (fastaFile) {
                    fastaFile >> name;
                    name = name.substr(1);
                }else {
                    name = "";
                    break;
                }
            }
            
            //while (!fastaFile.eof())    {    char c = fastaFile.get(); if (c == 10 || c == 13){  break;    }    } // get rest of line if there's any crap there
            comment = getCommentString(fastaFile);
            vector<AminoAcid> proteinSeq = getProtein(fastaFile);
            
            setAligned(proteinSeq);
            setUnaligned(proteinSeq); //setUnaligned removes any gap characters for us
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Protein::Protein(ifstream& fastaFile, string& extraInfo, bool getInfo){
    try {
        m = MothurOut::getInstance();
        initialize();
        extraInfo = "";
        
        name = getProteinName(fastaFile);
        
        if (!m->getControl_pressed()) {
            string sequence;
            
            //read comments
            while ((name[0] == '#') && fastaFile) {
                while (!fastaFile.eof())    {    char c = fastaFile.get(); if (c == 10 || c == 13){    break;    }    } // get rest of line if there's any crap there
                sequence = getCommentString(fastaFile);
                
                if (fastaFile) {
                    fastaFile >> name;
                    name = name.substr(1);
                }else {
                    name = "";
                    break;
                }
            }
            
            //read info after sequence name
            while (!fastaFile.eof())    {
                char c = fastaFile.get();
                if (c == 10 || c == 13 || c == -1){   break;    }
                extraInfo += c;
            }
            
            comment = extraInfo;
            
            vector<AminoAcid> proteinSeq = getProtein(fastaFile);

            setAligned(proteinSeq);
            setUnaligned(proteinSeq); //setUnaligned removes any gap characters for us
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "Protein");
        exit(1);
    }
}
//********************************************************************************************************************
string Protein::getProteinName(ifstream& fastaFile) {
    try {
        string name = "";
        
        fastaFile >> name;
        
        if (name.length() != 0) {
            
            name = name.substr(1);
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name.\n");  m->setControl_pressed(true);  } }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getProteinName");
        exit(1);
    }
}
//********************************************************************************************************************
#ifdef USE_BOOST
string Protein::getSequenceName(boost::iostreams::filtering_istream& fastaFile) {
    try {
        string name = "";
        
        fastaFile >> name;
        
        if (name.length() != 0) {
            
            name = name.substr(1);
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name.\n");  m->setControl_pressed(true);  }  }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getSequenceName");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
string Protein::getProteinName(istringstream& fastaFile) {
    try {
        string name = "";
        
        fastaFile >> name;
        
        if (name.length() != 0) {
            
            name = name.substr(1);
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name.\n");  m->setControl_pressed(true);  }  }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getProteinName");
        exit(1);
    }
}
//********************************************************************************************************************
vector<AminoAcid> Protein::getProtein(ifstream& fastaFile) {
    try {
        char letter;
        vector<AminoAcid> protein;
        
        while(!fastaFile.eof()){
            letter= fastaFile.get();
            if(letter == '>'){
                fastaFile.putback(letter);
                break;
            }else if (letter == ' ') {;}
            else if(isprint(letter)){
                letter = toupper(letter);
                if(letter == 'U'){letter = 'T';}
                AminoAcid amino(letter);
                protein.push_back(amino);
            }
        }
        
        return protein;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getProtein");
        exit(1);
    }
}
//********************************************************************************************************************
#ifdef USE_BOOST
vector<AminoAcid> Protein::getProtein(boost::iostreams::filtering_istream& fastaFile) {
    try {
        char letter;
        vector<AminoAcid> protein;
        
        while(fastaFile){
            letter= fastaFile.get();
            if(letter == '>'){
                fastaFile.putback(letter);
                break;
            }else if (letter == ' ') {;}
            else if(isprint(letter)){
                letter = toupper(letter);
                if(letter == 'U'){letter = 'T';}
                AminoAcid amino(letter);
                protein.push_back(amino);
            }
        }
        
        return protein;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getProtein");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Protein::getCommentString(ifstream& fastaFile) {
    try {
        char letter;
        string temp = "";
        
        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }
        
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getCommentString");
        exit(1);
    }
}
//********************************************************************************************************************
#ifdef USE_BOOST
//comment can contain '>' so we need to account for that
string Protein::getCommentString(boost::iostreams::filtering_istream& fastaFile) {
    try {
        char letter;
        string temp = "";
        
        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }
        
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getCommentString");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
vector<AminoAcid> Protein::getProtein(istringstream& fastaFile) {
    try {
        char letter;
        vector<AminoAcid> protein;
        
        while(!fastaFile.eof()){
            letter= fastaFile.get();
    
            if(letter == '>'){
                fastaFile.putback(letter);
                break;
            }else if (letter == ' ') {;}
            else if(isprint(letter)){
                letter = toupper(letter);
                if(letter == 'U'){letter = 'T';}
                AminoAcid amino(letter);
                protein.push_back(amino);
            }
        }
        
        return protein;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getProtein");
        exit(1);
    }
}
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Protein::getCommentString(istringstream& fastaFile) {
    try {
        char letter;
        string temp = "";
        
        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }
        
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "Protein", "getCommentString");
        exit(1);
    }
}
//********************************************************************************************************************

void Protein::initialize(){
    
    name = "";
    unaligned.clear();
    aligned.clear();
    pairwise.clear();
    comment = "";
    
    numBases = 0;
    alignmentLength = 0;
    startPos = -1;
    endPos = -1;
}

//********************************************************************************************************************

void Protein::setName(string seqName) {
    if(seqName[0] == '>')    {    name = seqName.substr(1);    }
    else                    {    name = seqName;                }
}

//********************************************************************************************************************

void Protein::setUnaligned(vector<AminoAcid> protein){
    unaligned.clear();
 
    for(int j=0;j<protein.size();j++) {
        if(isalpha(protein[j].getAmino()))    {    unaligned.push_back(protein[j]);    }
    }
       
    numBases = unaligned.size();
    startPos = -1;
    endPos = -1;
}
//********************************************************************************************************************

void Protein::setAligned(string seq){
    vector<AminoAcid> sequence;
    for (int i = 0; i < seq.size(); i++) {
        AminoAcid temp(seq[i]);
        sequence.push_back(temp);
    }
    
    setAligned(sequence);
}
//********************************************************************************************************************

void Protein::setAligned(vector<AminoAcid> sequence){
    
    //if the alignment starts or ends with a gap, replace it with a period to indicate missing data
    aligned.clear();
    aligned = sequence;
    alignmentLength = aligned.size();
    setUnaligned(sequence);

    if(aligned[0].getAmino() == '-'){ //convert ending gaps
        for(int i=0;i<alignmentLength;i++){
            if(aligned[i].getAmino() == '-'){ aligned[i].setAmino('.'); }
            else{ break; }
        }
        for(int i=alignmentLength-1;i>=0;i--){
            if(aligned[i].getAmino() == '-'){ aligned[i].setAmino('.'); }
            else{ break; }
        }
    }
}
//********************************************************************************************************************
bool Protein::isAligned(){
    
    for (int i = 0; i < aligned.size(); i++) {
        if ((aligned[i].getAmino() == '.') || (aligned[i].getAmino() == '-')) { return true; }
    }
    return false;
}
//********************************************************************************************************************

void Protein::setPairwise(vector<AminoAcid> sequence){ pairwise = sequence; }

//********************************************************************************************************************

string Protein::getName(){  return name; }

//********************************************************************************************************************

vector<AminoAcid> Protein::getAligned(){ return aligned; }

//********************************************************************************************************************

string Protein::getProteinString(vector<AminoAcid> prot){
    
    string inlinePro = "";
    
    for (int i = 0; i < prot.size(); i++) {
        inlinePro += prot[i].getAmino();
    }
    
    return inlinePro;
}

//********************************************************************************************************************

string Protein::getInlineProtein(){
    
    string inlinePro = name + '\t';
    
    inlinePro += getProteinString(aligned);
    
    return inlinePro;
}


//********************************************************************************************************************

vector<AminoAcid> Protein::getPairwise(){ return pairwise; }

//********************************************************************************************************************

vector<AminoAcid> Protein::getUnaligned(){ return unaligned; }

//********************************************************************************************************************

string Protein::getComment(){ return comment; }

//********************************************************************************************************************

int Protein::getNumBases(){ return numBases; }

//********************************************************************************************************************
void Protein::printProtein(OutputWriter* out){
    string seqOutput = ">";
    seqOutput += name + comment + '\n' + getProteinString(aligned) + '\n';
    out->write(seqOutput);
}
//********************************************************************************************************************

void Protein::printProtein(ostream& out){

    out << ">" << name << comment << endl;
    out << getProteinString(aligned) << endl;
}
//********************************************************************************************************************

void Protein::printUnAlignedProtein(ostream& out){
    
    out << ">" << name << comment << endl;
    out << getProteinString(unaligned) << endl;
}
//********************************************************************************************************************

int Protein::getAlignLength(){
    return alignmentLength;
}

//********************************************************************************************************************

int Protein::getStartPos(){
    bool isAligned = false;
    if(startPos == -1){
        for(int j = 0; j < alignmentLength; j++) {
            if((aligned[j].getAmino() != '.')&&(aligned[j].getAmino() != '-')){
                startPos = j + 1;
                break;
            }else { isAligned = true; }
        }
    }
    
    if(!isAligned){    startPos = 1;    }

    return startPos;
}

//********************************************************************************************************************

void Protein::filterToPos(int start){
    
    if (start > aligned.size()) { start = aligned.size(); m->mothurOut("[ERROR]: start to large.\n"); }
    
    for(int j = 0; j < start; j++) { aligned[j].setAmino('.'); }
    
    //things like ......----------AT become ................AT
    for(int j = start; j < aligned.size(); j++) {
        if (isalpha(aligned[j].getAmino())) { break; }
        else { aligned[j].setAmino('.'); }
    }
    setUnaligned(aligned);
    
}
//********************************************************************************************************************

void Protein::filterFromPos(int end){
    
    if (end > aligned.size()) { end = aligned.size(); m->mothurOut("[ERROR]: end to large.\n"); }
    
    for(int j = end; j < aligned.size(); j++) {
        aligned[j].setAmino('.');
    }
    
    for(int j = aligned.size()-1; j < 0; j--) {
        if (isalpha(aligned[j].getAmino())) { break; }
        else { aligned[j].setAmino('.'); }
    }
    
    setUnaligned(aligned);
}

//********************************************************************************************************************

int Protein::getEndPos(){
    bool isAligned = false;
    if (alignmentLength != numBases) { isAligned = true; }
    
    if(endPos == -1){
        for(int j=alignmentLength-1;j>=0;j--){
            if((aligned[j].getAmino() != '.')&&(aligned[j].getAmino() != '-')){
                endPos = j + 1;
                break;
            }else { isAligned = true; }
        }
    }
    if(!isAligned){    endPos = numBases;    }
    
    return endPos;
}
//********************************************************************************************************************

void Protein::padToPos(int start){
    
    for(int j = getStartPos()-1; j < start-1; j++) {
        aligned[j].setAmino('.');
    }
    startPos = start;
}

//********************************************************************************************************************

void Protein::padFromPos(int end){
    
    for(int j = end; j < getEndPos(); j++) {
        aligned[j].setAmino('.');
    }
    endPos = end;
    
}
//********************************************************************************************************************

void Protein::setComment(string c){
    comment = c;
}
//********************************************************************************************************************

void Protein::trim(int length){
    
    if(numBases > length){
        unaligned.resize(length);
        numBases = length;
        setAligned(unaligned);
    }
}

///**************************************************************************************************/
