/*
 *  chimera.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 8/11/09.
 *  Copyright 2009 Schloss Lab Umass Amherst. All rights reserved.
 *
 */

#include "chimera.h"
#include "referencedb.h"

//***************************************************************************************************************
//this is a vertical soft filter
string Chimera::createFilter(vector<Sequence*> seqs, float t) {
	try {
		filterString = "";
		int threshold = int (t * seqs.size());
//cout << "threshhold = " << threshold << endl;
		
		vector<int> gaps;	gaps.resize(seqs[0]->getAligned().length(), 0);
		vector<int> a;		a.resize(seqs[0]->getAligned().length(), 0);
		vector<int> t;		t.resize(seqs[0]->getAligned().length(), 0);
		vector<int> g;		g.resize(seqs[0]->getAligned().length(), 0);
		vector<int> c;		c.resize(seqs[0]->getAligned().length(), 0);
	
		filterString = (string(seqs[0]->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			if (m->control_pressed) { return filterString; }
		
			string seqAligned = seqs[i]->getAligned();
			
			if (seqAligned.length() != filterString.length()) {  m->mothurOut(seqs[i]->getName() + " is not the same length as the template sequences. Aborting!\n");  exit(1); }
		
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.'))	{	gaps[j]++;	}
				else if (toupper(seqAligned[j]) == 'A')					{	a[j]++;		}
				else if (toupper(seqAligned[j]) == 'T')					{	t[j]++;		}
				else if (toupper(seqAligned[j]) == 'G')					{	g[j]++;		}
				else if (toupper(seqAligned[j]) == 'C')					{	c[j]++;		}
			}
		}
		
		//zero out spot where all sequences have blanks
		int numColRemoved = 0;
		for(int i = 0;i < seqs[0]->getAligned().length(); i++){
		
			if (m->control_pressed) { return filterString; }
			
			if(gaps[i] == seqs.size())	{	filterString[i] = '0'; 	numColRemoved++;  }
			
			else if (((a[i] < threshold) && (t[i] < threshold) && (g[i] < threshold) && (c[i] < threshold))) {	filterString[i] = '0';	numColRemoved++;  }
			//cout << "a = " << a[i] <<  " t = " << t[i] <<  " g = " << g[i] <<  " c = " << c[i] << endl;
		}

		if (threshold != 0.0) {  m->mothurOut("Filter removed " + toString(numColRemoved) + " columns.");  m->mothurOutEndLine();  }
		
		return filterString;
	}
	catch(exception& e) {
		m->errorOut(e, "Chimera", "createFilter");
		exit(1);
	}
}
//***************************************************************************************************************
map<int, int> Chimera::runFilter(Sequence* seq) {
	try {
		map<int, int> maskMap;
		string seqAligned = seq->getAligned();
		string newAligned = "";
		int count = 0;
			
		for (int j = 0; j < seqAligned.length(); j++) {
			//if this spot is a gap
			if (filterString[j] == '1') { 
				newAligned += seqAligned[j]; 
				maskMap[count] = j;
				count++;
			}
		}
			
		seq->setAligned(newAligned);
		
		return maskMap;
	}
	catch(exception& e) {
		m->errorOut(e, "Chimera", "runFilter");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> Chimera::readSeqs(string file) {
	try {
		
		vector<Sequence*> container;
		int count = 0;
		length = 0;
		unaligned = false;
		ReferenceDB* rdb = ReferenceDB::getInstance();
		
		if (file == "saved") {
			
			
			m->mothurOutEndLine();  m->mothurOut("Using sequences from " + rdb->getSavedReference() + " that are saved in memory.");	m->mothurOutEndLine();
			
			for (int i = 0; i < rdb->referenceSeqs.size(); i++) {
				Sequence* temp = new Sequence(rdb->referenceSeqs[i].getName(), rdb->referenceSeqs[i].getAligned());
				
				if (count == 0) {  length = temp->getAligned().length();  count++;  } //gets first seqs length
				else if (length != temp->getAligned().length()) {	unaligned = true;	}
				
				if (temp->getName() != "") {  container.push_back(temp);  }
			}
			
			templateFileName = rdb->getSavedReference();
			
		}else {
			
			m->mothurOut("Reading sequences from " + file + "..."); cout.flush();
			
			
			ifstream in;
			m->openInputFile(file, in);
			
			//read in seqs and store in vector
			while(!in.eof()){
				
				if (m->control_pressed) { return container; }
				
				Sequence* current = new Sequence(in);  m->gobble(in);
				
				if (count == 0) {  length = current->getAligned().length();  count++;  } //gets first seqs length
				else if (length != current->getAligned().length()) {   unaligned = true;	}
							
				if (current->getName() != "") {  
					container.push_back(current);  
					if (rdb->save) { rdb->referenceSeqs.push_back(*current); }
				}
			}
			in.close();
		
			m->mothurOut("Done."); m->mothurOutEndLine();
			
			filterString = (string(container[0]->getAligned().length(), '1'));
		}
		
		return container;
	}
	catch(exception& e) {
		m->errorOut(e, "Chimera", "readSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
void Chimera::setMask(string filename) {
	try {
		
		if (filename == "default") {
			//default is from wigeon  236627 EU009184.1 Shigella dysenteriae str. FBD013
			seqMask = ".....................................................................................................AAATTGAAGAGTTT-GA--T-CA-T-G-GCTC-AG-AT-TGAA-C-GC--TGG-C--G-GC-A-GG--C----C-T--AACACA-T-GC-A-AGT-CGA-A-CG----------G-TAA-CA-G----------------------------GAAG-A-AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CT-TCTTT----------------G-CT--G--AC--G--AG-T-GG-C-GG-A--C-------------GGG-TGAGT-A--AT-GT-C-T-G-GG---A-A--A-CT-G--C-C-TGA--TG-G------------------------------------------------------------------A-GG----GGG-AT-AA-CTA-------------------------C-T-G-----------------------GAA-A---CGG-TAG-CTAA-TA---CC-G--C-AT-A----------A--------------------C-------------------------------------GT-C-----------------------------------------------------------------------------------------------------------------------G-CA-A--------------------------------------------------------------------------------------------------------------------------------------G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-A-G-GG-----G--GA-C-CT--------------------------------------------------------------------------------------------------------------------TCG-G----------------------------------------------------------------------------------------------------------------------G----CC-TC--T---T-G--------------C----C-A---T-CG-G---AT---G-T-----G-CCC-AGA--T-GGG--A------TT--A--G-CT-A----G---TAGG-T-G-GG-G-T----AAC-GG-C-T-C-ACCT--A-GG-C-G--A-CG-A------------TCC-C-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AC-C-AG-CCAC-A-CTGGA--A-C-TG-A-GA-C-AC-G-G-TCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GC-GC-A----A-G-CC-T-GA-TG-CA-GCCA-TGCC-G-CG-T---G-T-A--T--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-T-G-T--A---AA-G-TAC--------TT-TC-A-G--C-GGG----GA-G--G---AA-GGGA---GTAA-AG----T--T--AA-T---A----C-----CT-T-TGC-TCA-TT-GA-CG-TT-A-C-CC-G-CA-G---------AA-----------GAAGC-ACC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GGT-GCA-A-G-CG-TTAA-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GCGC-AC--G-CA-G-G-C-G------------G--T-TT-G-T-T-AA----G-T-C-A---G-ATG-TG-A-AA-TC--CC-CGG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-T-G-GG-AA-C----T-G-C-A-T-C--------T--GA-T-A-C-T-G-GCA--A-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-C-----T-CG--TA-G-A------------G-GG-G-GG-T----AG--AATT-CCA-G-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TC-T-GGA--GG-A-AT-A-CC-GG--T--G--GC-GAA-G--G-C---G----G--C-C-CCCTG------G-AC-GA--------------------------------------------------------------AG-A-C-T--GA--CG-----CT-CA-GG--T-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-TC--GA-CT---------T-GG--A--G-G-TT-G-TG-C--C--------------------------------------------------------------------------------------CTT-GA--------------------------------------------------------------------------------------------------------------------------------------------------G-G-C-GT--G-G-C-T-TC-C------GG--A----GC-TAA--CG-C-G-T--T--AA-GT--C----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGT-T--AAA-ACTC-AAA---------TGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-TGGTC-TT-G-AC-A-T-C--------------CAC-G-G-------------A-AG-T-T-T--TC--A-GA-G-A-T--G-A-G--A-A-T-G--T-G-----CC-------------------------------------T--TC-G------------------------------------------GG----A----A---CC-GTG---A--GA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TT-G--TGA-AA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TC--C-TTTG--T-T-G-C-C---AG-C-G-----G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCC------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG---C----C-G------------G----G---A-A--CT---------------C-A-A-A-G-GA-G--AC-T-G-CCA--G-T------------------------------------G-A---TAA----------------------------------A-C-T-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-G-C-C-CTT----AC-G--AC-C-A-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-CGCA-T-A--C-AAA-GA-GA--------------------------------------------------------------------------------------------------A-G-C-G-A--C-CTCG-C--G---------------------------------------A-GA-G-C-----------A--A-G-CG---G----------A--CCT-C------A-T-AAAGT-GC-G-T-C-G-TAG-TCC--------GGA-T-TGGAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTCC-A-T-G-AA-G-TC-GGAAT-CG-C-TA--G-TA-AT-C-G-T----GGA-TC-A-G--A------AT--GCC-AC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---TGG-G-TT-GC-AAA--A-GAA------G--T-AGG-TA-G-C-T-T-AA-C-C--------------------------------------------------------------TT----C-------------------------------------------------------------------------------------------------G--GG-A--GG-G--C---GC-TTA--CC--ACT-T----T-GTG-AT-TCA------------------------TG--ACT-GGGG-TG-AAG-TCGTAACAA-GGTAA-CCGT-AGGGGAA-CCTG-CGGT-TGGATCACCTCCTTA................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................";
		}else if (filename == "") {  //do nothing 
			seqMask = "";
		}else{
		
            ifstream infile;
            m->openInputFile(filename, infile);
            
            if (!infile.eof()) {
                Sequence temp(infile);
                seqMask = temp.getAligned();
            }else {
                m->mothurOut("Problem with mask."); m->mothurOutEndLine();
                seqMask = "";
            }
            infile.close();
	
	
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Chimera", "setMask");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* Chimera::getSequence(string name) {
	try{
		Sequence* temp;
		
		//look through templateSeqs til you find it
		int spot = -1;
		for (int i = 0; i < templateSeqs.size(); i++) {
			if (name == templateSeqs[i]->getName()) {  
				spot = i;
				break;
			}
		}
		
		if(spot == -1) { m->mothurOut("Error: Could not find sequence."); m->mothurOutEndLine(); return NULL; }
		
		temp = new Sequence(templateSeqs[spot]->getName(), templateSeqs[spot]->getAligned());
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "Chimera", "getSequence");
		exit(1);
	}
}
//***************************************************************************************************************




