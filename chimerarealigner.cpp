/*
 *  chimerarealigner.cpp
 *  Mothur
 *
 *  Created by westcott on 2/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerarealigner.h"
#include "needlemanoverlap.hpp"
#include "nast.hpp"

//***************************************************************************************************************
ChimeraReAligner::ChimeraReAligner(vector<Sequence*> t, int ms, int mm) : match(ms), misMatch(mm) {  templateSeqs = t;   m = MothurOut::getInstance(); }
//***************************************************************************************************************
ChimeraReAligner::~ChimeraReAligner() {}	
//***************************************************************************************************************
void ChimeraReAligner::reAlign(Sequence* query, vector<results> parents) {
	try {
		if (parents.size() != 0) {
			vector<Sequence*> queryParts;
			vector<Sequence*> parentParts;  //queryParts[0] relates to parentParts[0]
			
			string qAligned = query->getAligned();
			string newQuery = "";
			
			//sort parents by region start
			sort(parents.begin(), parents.end(), compareRegionStart);

			//make sure you don't cutoff beginning of query 
			if (parents[0].nastRegionStart > 0) {  newQuery += qAligned.substr(0, parents[0].nastRegionStart);  }
			int longest = 0;

			//take query and break apart into pieces using breakpoints given by results of parents
			for (int i = 0; i < parents.size(); i++) {
			
				int length = parents[i].nastRegionEnd - parents[i].nastRegionStart+1;
				string q = qAligned.substr(parents[i].nastRegionStart, length);
				Sequence* queryFrag = new Sequence(query->getName(), q);

				queryParts.push_back(queryFrag);
			
				Sequence* parent = getSequence(parents[i].parent);
				string p = parent->getAligned();
		
				p = p.substr(parents[i].nastRegionStart, length);
				parent->setAligned(p);
				
				parentParts.push_back(parent);

				if (q.length() > longest)	{ longest = q.length(); }
				if (p.length() > longest)	{ longest = p.length();	}
			}

			//align each peice to correct parent from results
			for (int i = 0; i < queryParts.size(); i++) {
				alignment = new NeedlemanOverlap(-2.0, match, misMatch, longest+1); //default gapopen, match, mismatch, longestbase
				Nast nast(alignment, queryParts[i], parentParts[i]);
				delete alignment;
			}

			//recombine pieces to form new query sequence
			for (int i = 0; i < queryParts.size(); i++) {
				//sometimes the parent regions do not meet, for example region 1 may end at 1000 and region 2 starts at 1100.  
				//we don't want to loose length so in this case we will leave query alone
				if (i != 0) {
					int space = parents[i].nastRegionStart - parents[i-1].nastRegionEnd - 1;
					if (space > 0) { //they don't meet and we need to add query piece
						string q = qAligned.substr(parents[i-1].nastRegionEnd+1, space);
						newQuery += q;
					}
				}

				newQuery += queryParts[i]->getAligned();
			}
			
			//make sure you don't cutoff end of query 
			if (parents[parents.size()-1].nastRegionEnd < (qAligned.length()-1)) {  newQuery += qAligned.substr(parents[parents.size()-1].nastRegionEnd+1);  }
		
			//set query to new aligned string
			query->setAligned(newQuery);
			
			//free memory
			for (int i = 0; i < queryParts.size(); i++) { delete queryParts[i];  }
			for (int i = 0; i < parentParts.size(); i++) { delete parentParts[i];  }
			
		} //else leave query alone, you have bigger problems...
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "reAlign");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* ChimeraReAligner::getSequence(string name) {
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
		m->errorOut(e, "ChimeraReAligner", "getSequence");
		exit(1);
	}
}
//***************************************************************************************************************
