/*
 *  taxonomyequalizer.cpp
 *  Mothur
 *
 *  Created by westcott on 11/20/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "taxonomyequalizer.h"

/**************************************************************************************************/
TaxEqualizer::TaxEqualizer(string tfile, int c) : cutoff(c) {
	try {
		m = MothurOut::getInstance();
		containsConfidence = false;
		
		ifstream inTax;
		openInputFile(tfile, inTax);
	
		highestLevel = getHighestLevel(inTax);
		
		if (!m->control_pressed) { 
			
			//if the user has specified a cutoff and it's smaller than the highest level
			if ((cutoff != -1) && (cutoff < highestLevel)) { 
				highestLevel = cutoff;
			}else if (cutoff > highestLevel) {
				m->mothurOut("The highest level taxonomy you have is " + toString(highestLevel) + " and your cutoff is " + toString(cutoff) + ". I will set the cutoff to " + toString(highestLevel));
				m->mothurOutEndLine();
			}
			
			inTax.close(); 
			ifstream in; 
			openInputFile(tfile, in);
			
			ofstream out;
			equalizedFile = getRootName(tfile) + "equalized.taxonomy";
			openOutputFile(equalizedFile, out);
			
	
			string name, tax;
			while (in) {
			
				if (m->control_pressed) {  break; }
				
				in >> name >> tax;   gobble(in);
				
				if (containsConfidence) {  removeConfidences(tax);	}
				
				//is this a taxonomy that needs to be extended?
				if (seqLevels[name] < highestLevel) {
					extendTaxonomy(name, tax, highestLevel);
				}else if (seqLevels[name] > highestLevel) { //this can happen if the user enters a cutoff
					truncateTaxonomy(name, tax, highestLevel);
				}
				
				out << name << '\t' << tax << endl;
			}
			
			in.close();
			out.close();
			
			if (m->control_pressed) { remove(equalizedFile.c_str());  }
		}else { inTax.close(); }
		
	}
	catch(exception& e) {
		m->errorOut(e, "TaxEqualizer", "TaxEqualizer");
		exit(1);
	}
}
/**************************************************************************************************/
int TaxEqualizer::getHighestLevel(ifstream& in) {
	try {
		
		int level = 0;
		string name, tax;
		
		while (in) {
			in >> name >> tax;   gobble(in);
		
			//count levels in this taxonomy
			int thisLevel = 0;
			for (int i = 0; i < tax.length(); i++) {
				if (tax[i] == ';') {  thisLevel++;	}
			}
		
			//save sequences level
			seqLevels[name] = thisLevel;
		
			//is this the longest taxonomy?
			if (thisLevel > level) {  
				level = thisLevel;  
				testTax = tax; //testTax is used to figure out if this file has confidences we need to strip out
			}  
		}
		
		int pos = testTax.find_first_of('(');

		//if there are '(' then there are confidences we need to take out
		if (pos != -1) {  containsConfidence = true;  }
		
		return level;
					
	}
	catch(exception& e) {
		m->errorOut(e, "TaxEqualizer", "getHighestLevel");
		exit(1);
	}
}
/**************************************************************************************************/
void TaxEqualizer::extendTaxonomy(string name, string& tax, int desiredLevel) {
	try {
			
		//get last taxon
		tax = tax.substr(0, tax.length()-1);  //take off final ";"
		int pos = tax.find_last_of(';');
		string lastTaxon = tax.substr(pos+1);  
		lastTaxon += ";"; //add back on delimiting char
		tax += ";";
		
		int currentLevel = seqLevels[name];
		
		//added last taxon until you get desired level
		for (int i = currentLevel; i < desiredLevel; i++) {
			tax += lastTaxon;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "TaxEqualizer", "extendTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
void TaxEqualizer::truncateTaxonomy(string name, string& tax, int desiredLevel) {
	try {
			
		int currentLevel = seqLevels[name];
		tax = tax.substr(0, tax.length()-1);  //take off final ";"
		
		//remove a taxon until you get to desired level
		for (int i = currentLevel; i > desiredLevel; i--) {
			tax = tax.substr(0,  tax.find_last_of(';'));
		}
		
		tax += ";";
	}
	catch(exception& e) {
		m->errorOut(e, "TaxEqualizer", "truncateTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
void TaxEqualizer::removeConfidences(string& tax) {
	try {
		
		string taxon;
		string newTax = "";
		
		while (tax.find_first_of(';') != -1) {
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
			taxon = taxon.substr(0, taxon.find_first_of('(')); //rip off confidence
			taxon += ";";
			
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
		}
		
		tax = newTax;
	}
	catch(exception& e) {
		m->errorOut(e, "TaxEqualizer", "removeConfidences");
		exit(1);
	}
}

/**************************************************************************************************/

