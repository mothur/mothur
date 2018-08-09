/*
 *  sequence.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "sequence.hpp"

/***********************************************************************/
Sequence::Sequence(){
	m = MothurOut::getInstance();
	initialize();
}
/***********************************************************************/
Sequence::Sequence(string newName, string sequence) {
	try {
		m = MothurOut::getInstance();
		initialize();	
		name = newName;
        
        util.checkName(name);
		
		//setUnaligned removes any gap characters for us
		setUnaligned(sequence);
		setAligned(sequence);
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}			
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(istringstream& fastaString){
	try {
		m = MothurOut::getInstance();
	
		initialize();
        name = getSequenceName(fastaString);
		
		if (!m->getControl_pressed()) { 
			string sequence;
		
			//read comments
			while ((name[0] == '#') && fastaString) { 
				while (!fastaString.eof())	{	char c = fastaString.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				sequence = getCommentString(fastaString);
				
				if (fastaString) {  
					fastaString >> name;  
					name = name.substr(1);	
				}else { 
					name = "";
					break;
				}
			}
			
			//while (!fastaString.eof())	{	char c = fastaString.get();  if (c == 10 || c == 13){ break;	}	} // get rest of line if there's any crap there
            comment = getCommentString(fastaString);
			
			int numAmbig = 0;
			sequence = getSequenceString(fastaString, numAmbig);
			
			setAligned(sequence);	
			//setUnaligned removes any gap characters for us						
			setUnaligned(sequence);	
			
			if ((numAmbig / (float) numBases) > 0.25) { m->mothurOut("[WARNING]: We found more than 25% of the bases in sequence " + name + " to be ambiguous. Mothur is not setup to process protein sequences."); m->mothurOutEndLine(); }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}								
}

//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(ifstream& fastaFile){
	try {
		m = MothurOut::getInstance();
		initialize();
		name = getSequenceName(fastaFile);
		
		if (!m->getControl_pressed()) { 
			
			string sequence;
		
			//read comments
			while ((name[0] == '#') && fastaFile) { 
				while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				sequence = getCommentString(fastaFile);
				
				if (fastaFile) {  
					fastaFile >> name;  
					name = name.substr(1);	
				}else { 
					name = "";
					break;
				}
			}
			
			//while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){  break;	}	} // get rest of line if there's any crap there
            comment = getCommentString(fastaFile);
			
			int numAmbig = 0;
			sequence = getSequenceString(fastaFile, numAmbig);
			
			setAligned(sequence);	
			//setUnaligned removes any gap characters for us						
			setUnaligned(sequence);	
			
			if ((numAmbig / (float) numBases) > 0.25) { m->mothurOut("[WARNING]: We found more than 25% of the bases in sequence " + name + " to be ambiguous. Mothur is not setup to process protein sequences."); m->mothurOutEndLine(); }
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}							
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
#ifdef USE_BOOST
Sequence::Sequence(boost::iostreams::filtering_istream& fastaFile){
    try {
        m = MothurOut::getInstance();
        initialize();
        name = getSequenceName(fastaFile);
        
        if (!m->getControl_pressed()) {
            
            string sequence;
            
            //read comments
            while ((name[0] == '#') && fastaFile) {
                while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
                sequence = getCommentString(fastaFile);
                
                if (fastaFile) {
                    fastaFile >> name;
                    name = name.substr(1);
                }else {
                    name = "";
                    break;
                }
            }
            
            //while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){  break;	}	} // get rest of line if there's any crap there
            comment = getCommentString(fastaFile);
            
            int numAmbig = 0;
            sequence = getSequenceString(fastaFile, numAmbig);
            
            setAligned(sequence);
            //setUnaligned removes any gap characters for us
            setUnaligned(sequence);
            
            if ((numAmbig / (float) numBases) > 0.25) { m->mothurOut("[WARNING]: We found more than 25% of the bases in sequence " + name + " to be ambiguous. Mothur is not setup to process protein sequences."); m->mothurOutEndLine(); }
            
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Sequence", "Sequence");
        exit(1);
    }							
}
#endif
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(ifstream& fastaFile, string& extraInfo, bool getInfo){
	try {
		m = MothurOut::getInstance();
		initialize();
        extraInfo = "";
		
		name = getSequenceName(fastaFile);
		
		if (!m->getControl_pressed()) { 			
			string sequence;
            
			//read comments
			while ((name[0] == '#') && fastaFile) { 
				while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
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
			while (!fastaFile.eof())	{	
                char c = fastaFile.get(); 
                if (c == 10 || c == 13 || c == -1){   break;	}
                extraInfo += c;
            }
            
            comment = extraInfo;
			
			int numAmbig = 0;
			sequence = getSequenceString(fastaFile, numAmbig);
			
			setAligned(sequence);	
			//setUnaligned removes any gap characters for us						
			setUnaligned(sequence);	
			
			if ((numAmbig / (float) numBases) > 0.25) { m->mothurOut("[WARNING]: We found more than 25% of the bases in sequence " + name + " to be ambiguous. Mothur is not setup to process protein sequences."); m->mothurOutEndLine(); }
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}							
}
//********************************************************************************************************************
string Sequence::getSequenceName(ifstream& fastaFile) {
	try {
		string name = "";
		
        fastaFile >> name;
		
		if (name.length() != 0) { 
            
			name = name.substr(1); 
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->setControl_pressed(true);  } }
        
		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceName");
		exit(1);
	}
}
//********************************************************************************************************************
#ifdef USE_BOOST
string Sequence::getSequenceName(boost::iostreams::filtering_istream& fastaFile) {
    try {
        string name = "";
        
        fastaFile >> name;
        
        if (name.length() != 0) {
            
            name = name.substr(1);
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->setControl_pressed(true);  }  }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "Sequence", "getSequenceName");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
string Sequence::getSequenceName(istringstream& fastaFile) {
	try {
		string name = "";
		
        fastaFile >> name;
		
		if (name.length() != 0) { 
            
			name = name.substr(1); 
            
            util.checkName(name);
            
        }else{ if (!fastaFile.eof()) { m->mothurOut("Error in reading your fastafile, at position " + toString(fastaFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->setControl_pressed(true);  }  }
        
		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceName");
		exit(1);
	}
}
//********************************************************************************************************************
string Sequence::getSequenceString(ifstream& fastaFile, int& numAmbig) {
	try {
		char letter;
		string sequence = "";	
		numAmbig = 0;
		
		while(fastaFile){
			letter= fastaFile.get();
			if(letter == '>'){
				fastaFile.putback(letter);
				break;
			}else if (letter == ' ') {;}
			else if(isprint(letter)){
				letter = toupper(letter);
				if(letter == 'U'){letter = 'T';}
				if(letter != '.' && letter != '-' && letter != 'A' && letter != 'T' && letter != 'G'  && letter != 'C' && letter != 'N'){
					letter = 'N';
					numAmbig++;
				}
				sequence += letter;
			}
		}
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceString");
		exit(1);
	}
}
//********************************************************************************************************************
#ifdef USE_BOOST
string Sequence::getSequenceString(boost::iostreams::filtering_istream& fastaFile, int& numAmbig) {
    try {
        char letter;
        string sequence = "";
        numAmbig = 0;
        
        while(fastaFile){
            letter= fastaFile.get();
            if(letter == '>'){
                fastaFile.putback(letter);
                break;
            }else if (letter == ' ') {;}
            else if(isprint(letter)){
                letter = toupper(letter);
                if(letter == 'U'){letter = 'T';}
                if(letter != '.' && letter != '-' && letter != 'A' && letter != 'T' && letter != 'G'  && letter != 'C' && letter != 'N'){
                    letter = 'N';
                    numAmbig++;
                }
                sequence += letter;
            }
        }
        
        return sequence;
    }
    catch(exception& e) {
        m->errorOut(e, "Sequence", "getSequenceString");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Sequence::getCommentString(ifstream& fastaFile) {
	try {
		char letter;
		string temp = "";
		
		while(fastaFile){
			letter=fastaFile.get();
			if((letter == '\r') || (letter == '\n') || letter == -1){
				util.gobble(fastaFile);  //in case its a \r\n situation
				break;
			}else {
                temp += letter;
            }
		}
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getCommentString");
		exit(1);
	}
}
//********************************************************************************************************************
#ifdef USE_BOOST
//comment can contain '>' so we need to account for that
string Sequence::getCommentString(boost::iostreams::filtering_istream& fastaFile) {
    try {
        char letter;
        string temp = "";
        
        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                util.gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }
        
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "Sequence", "getCommentString");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
string Sequence::getSequenceString(istringstream& fastaFile, int& numAmbig) {
	try {
		char letter;
		string sequence = "";
		numAmbig = 0;
		
		while(!fastaFile.eof()){
			letter= fastaFile.get();
	
			if(letter == '>'){
				fastaFile.putback(letter);
				break;
			}else if (letter == ' ') {;}
			else if(isprint(letter)){
				letter = toupper(letter);
				if(letter == 'U'){letter = 'T';}
				if(letter != '.' && letter != '-' && letter != 'A' && letter != 'T' && letter != 'G'  && letter != 'C' && letter != 'N'){
					letter = 'N';
					numAmbig++;
				}
				sequence += letter;
			}
		}
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceString");
		exit(1);
	}
}
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Sequence::getCommentString(istringstream& fastaFile) {
	try {
		char letter;
		string temp = "";
		
		while(fastaFile){
			letter=fastaFile.get();
			if((letter == '\r') || (letter == '\n') || letter == -1){  
				util.gobble(fastaFile);  //in case its a \r\n situation
				break;
			}else {
                temp += letter;
            }
		}
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getCommentString");
		exit(1);
	}
}
//********************************************************************************************************************

void Sequence::initialize(){
	
	name = "";
	unaligned = "";
	aligned = "";
	pairwise = "";
    comment = "";
	
	numBases = 0;
	alignmentLength = 0;
	startPos = -1;
	endPos = -1;
	longHomoPolymer = -1;
	ambigBases = -1;
	
}	

//********************************************************************************************************************

void Sequence::setName(string seqName) {
	if(seqName[0] == '>')	{	name = seqName.substr(1);	}
	else					{	name = seqName;				}
}

//********************************************************************************************************************

void Sequence::setUnaligned(string sequence){
	
	if(sequence.find_first_of('.') != string::npos || sequence.find_first_of('-') != string::npos) {
		string temp = "";
		for(int j=0;j<sequence.length();j++) {
			if(isalpha(sequence[j]))	{	temp += sequence[j];	}
		}
		unaligned = temp;
	}
	else {
		unaligned = sequence;
	}
	numBases = unaligned.length();
    startPos = -1;
    endPos = -1;
    longHomoPolymer = -1;
    ambigBases = -1;
	
}

//********************************************************************************************************************

void Sequence::setAligned(string sequence){
	
	//if the alignment starts or ends with a gap, replace it with a period to indicate missing data
	aligned = sequence;
	alignmentLength = aligned.length();
	setUnaligned(sequence);	

	if(aligned[0] == '-'){
		for(int i=0;i<alignmentLength;i++){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
		for(int i=alignmentLength-1;i>=0;i--){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
	}
}

//********************************************************************************************************************

void Sequence::setPairwise(string sequence){
	pairwise = sequence;
}

//********************************************************************************************************************

string Sequence::convert2ints() {
	
	if(unaligned == "")	{	/* need to throw an error */	}
	
	string processed;
	
	for(int i=0;i<unaligned.length();i++) {
		if(toupper(unaligned[i]) == 'A')		{	processed += '0';	}
		else if(toupper(unaligned[i]) == 'C')	{	processed += '1';	}
		else if(toupper(unaligned[i]) == 'G')	{	processed += '2';	}
		else if(toupper(unaligned[i]) == 'T')	{	processed += '3';	}
		else if(toupper(unaligned[i]) == 'U')	{	processed += '3';	}
		else									{	processed += '4';	}
	}
	return processed;
}

//********************************************************************************************************************

string Sequence::getName(){
	return name;
}

//********************************************************************************************************************

string Sequence::getAligned(){
	 return aligned;
}

//********************************************************************************************************************

string Sequence::getInlineSeq(){
	return name + '\t' + aligned;	
}


//********************************************************************************************************************

string Sequence::getPairwise(){
	return pairwise;
}

//********************************************************************************************************************

string Sequence::getUnaligned(){
	return unaligned;
}
//********************************************************************************************************************

string Sequence::getComment(){
    return comment;
}
//********************************************************************************************************************

int Sequence::getNumBases(){
	return numBases;
}
//********************************************************************************************************************

int Sequence::getNumNs(){
    int numNs = 0;
	for (int i = 0; i < unaligned.length(); i++) {
        if(toupper(unaligned[i]) == 'N') { numNs++; }
    }
    return numNs;
}
//********************************************************************************************************************
void Sequence::printSequence(OutputWriter* out){
    string seqOutput = ">";
    seqOutput += name + '\t' + comment + '\n';
    seqOutput += aligned + '\n';
    out->write(seqOutput);
}
//********************************************************************************************************************

void Sequence::printSequence(ostream& out){

	out << ">" << name << comment << endl;
    out << aligned << endl;
}
//********************************************************************************************************************

void Sequence::printUnAlignedSequence(ostream& out){
    
    out << ">" << name << comment << endl;
    out << unaligned << endl;
}
//********************************************************************************************************************

int Sequence::getAlignLength(){
	return alignmentLength;
}

//********************************************************************************************************************

int Sequence::getAmbigBases(){
	if(ambigBases == -1){
		ambigBases = 0;
		for(int j=0;j<numBases;j++){
			if(unaligned[j] != 'A' && unaligned[j] != 'T' && unaligned[j] != 'G' && unaligned[j] != 'C'){
				ambigBases++;
			}
		}
	}	
	
	return ambigBases;
}

//********************************************************************************************************************

void Sequence::removeAmbigBases(){
	
	for(int j=getStartPos();j<getEndPos();j++){
		if(aligned[j] != 'A' && aligned[j] != 'T' && aligned[j] != 'G' && aligned[j] != 'C'){
			aligned[j] = '-';
		}
	}
	setUnaligned(aligned);
}
	
//********************************************************************************************************************

int Sequence::getLongHomoPolymer(){
	if(longHomoPolymer == -1){
		longHomoPolymer = 1;
		int homoPolymer = 1;
		for(int j=1;j<numBases;j++){
			if(unaligned[j] == unaligned[j-1]){
				homoPolymer++;
			}
			else{
				if(homoPolymer > longHomoPolymer){	longHomoPolymer = homoPolymer;	}
				homoPolymer = 1;
			}
		}
		if(homoPolymer > longHomoPolymer){	longHomoPolymer = homoPolymer;	}
	}
	return longHomoPolymer;
}

//********************************************************************************************************************

int Sequence::getStartPos(){
    bool isAligned = false;
	if(startPos == -1){
		for(int j = 0; j < alignmentLength; j++) {
			if((aligned[j] != '.')&&(aligned[j] != '-')){
				startPos = j + 1;
				break;
            }else { isAligned = true; }
		}
	}
    
	if(!isAligned){	startPos = 1;	}

	return startPos;
}

//********************************************************************************************************************

int Sequence::filterToPos(int start){
    
    if (start > aligned.length()) { start = aligned.length(); m->mothurOut("[ERROR]: start to large.\n"); }
    
	for(int j = 0; j < start; j++) { aligned[j] = '.'; }
	
    //things like ......----------AT become ................AT
    for(int j = start; j < aligned.length(); j++) {
        if (isalpha(aligned[j])) { break; }
        else { aligned[j] = '.'; }
    }
    setUnaligned(aligned);
    
    return 0;
    
}
//********************************************************************************************************************

int Sequence::filterFromPos(int end){
    
    if (end > aligned.length()) { end = aligned.length(); m->mothurOut("[ERROR]: end to large.\n"); }
    
	for(int j = end; j < aligned.length(); j++) {
		aligned[j] = '.';
	}
	
    for(int j = aligned.length()-1; j < 0; j--) {
        if (isalpha(aligned[j])) { break; }
        else { aligned[j] = '.'; }
    }
    
    setUnaligned(aligned);
    
    return 0;
}

//********************************************************************************************************************

int Sequence::getEndPos(){
    bool isAligned = false;
    if (alignmentLength != numBases) { isAligned = true; }
    
	if(endPos == -1){
		for(int j=alignmentLength-1;j>=0;j--){
			if((aligned[j] != '.')&&(aligned[j] != '-')){
				endPos = j + 1;
				break;
            }else { isAligned = true; }
		}
	}
	if(!isAligned){	endPos = numBases;	}
	
	return endPos;
}
//********************************************************************************************************************

void Sequence::padToPos(int start){
    
    for(int j = getStartPos()-1; j < start-1; j++) {
        aligned[j] = '.';
    }
    startPos = start;
    
}

//********************************************************************************************************************

void Sequence::padFromPos(int end){
	
	for(int j = end; j < getEndPos(); j++) {
		aligned[j] = '.';
	}
	endPos = end;
	
}
//********************************************************************************************************************

void Sequence::setComment(string c){
    comment = c;
}
//********************************************************************************************************************

void Sequence::reverseComplement(){

	string temp;
	for(int i=numBases-1;i>=0;i--){
		if(unaligned[i] == 'A')		{	temp += 'T';	}
		else if(unaligned[i] == 'T'){	temp += 'A';	}
		else if(unaligned[i] == 'G'){	temp += 'C';	}
		else if(unaligned[i] == 'C'){	temp += 'G';	}
		else						{	temp += 'N';	}
	}
	
	setAligned(temp);
	
}

//********************************************************************************************************************

void Sequence::trim(int length){
	
	if(numBases > length){
		unaligned = unaligned.substr(0,length);
		numBases = length;
        setAligned(unaligned);
	}
}

///**************************************************************************************************/
