/*
 *  chimeraslayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Pschloss Lab. All rights reserved.
 *
 */

#include "chimeraslayer.h"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

ChimeraSlayer::~ChimeraSlayer() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)			{  delete querySeqs[i];			}
		for (int i = 0; i < templateSeqs.size(); i++)		{  delete templateSeqs[i];		}
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "~ChimeraSlayer");
		exit(1);
	}
}	
//***************************************************************************************************************
void ChimeraSlayer::print(ostream& out) {
	try {
		mothurOutEndLine();
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			if (chimeraFlags[i] == "yes") {	
				mothurOut(querySeqs[i]->getName() + "\tyes"); mothurOutEndLine();
			
			}else{
				out << querySeqs[i]->getName() << "\tno" << endl;
				mothurOut("no");
			}
		}
/*		
	
	my $div_ratio_QLA_QRB = $data_struct->{div_ratio_QLA_QRB};
	my $div_ratio_QRA_QLB = $data_struct->{div_ratio_QLB_QRA};
	
	my $per_id_QLA = $data_struct->{per_id_QLA};
	my $per_id_QRB = $data_struct->{per_id_QRB};
	my $per_id_AB = $data_struct->{per_id_AB};
	my $per_id_QA = $data_struct->{per_id_QA};
	my $per_id_QB = $data_struct->{per_id_QB}; 
	my $per_id_LAB = $data_struct->{per_id_LAB};
	my $per_id_RAB = $data_struct->{per_id_RAB};
	my $per_id_QRA = $data_struct->{per_id_QRA};
	my $per_id_QLB = $data_struct->{per_id_QLB};
	my $per_id_QLB_QRA = $data_struct->{per_id_QLB_QRA};
	my $per_id_QLA_QRB = $data_struct->{per_id_QLA_QRB};
	
	my $win_left_end5 = $data_struct->{win_left_end5};
	my $win_left_end3 = $data_struct->{win_left_end3};
	my $win_right_end5 = $data_struct->{win_right_end5};
	my $win_right_end3 = $data_struct->{win_right_end3};
	my $Q = $data_struct->{query_alignment};
	my $A = $data_struct->{parent_A_alignment};
	my $B = $data_struct->{parent_B_alignment}; 
	my $BS_A = $data_struct->{BS_A};
	my $BS_B = $data_struct->{BS_B};
	
	my @Q_chars = @{$Q->{align}};
	my @A_chars = @{$A->{align}};
	my @B_chars = @{$B->{align}};
	
	my $query_acc = $Q->{acc};
	my $A_acc = $A->{acc};
	my $B_acc = $B->{acc};
	
	my $break_left = $Q->{seqPos}->[$win_left_end3];
	my $break_right = $Q->{seqPos}->[$win_right_end5];
	
	
	cout << "//\n## CHIMERA\t" << querySeqs[i]->getName() << "\t" << $break_left-$break_right" << endl  
		<< "\tDIV_QLARB: ". sprintf("%.3f", $div_ratio_QLA_QRB)
		<< "\tBS_QLARB: " . sprintf("%.2f", $BS_A)
		<< "\tDIV_QRALB: " . sprintf("%.3f", $div_ratio_QRA_QLB)
		<< "\tBS_QRALB: " . sprintf("%.2f", $BS_B)
		<< "\t$A_acc\t$B_acc" 
		<< "\tbreakpoint: $break_left-$break_right\n\n";
	
	## draw illustration:

	print "            Per_id parents: " . sprintf("%.2f", $per_id_AB) . "\n\n";
	print "           Per_id(Q,A): " . sprintf("%.2f", $per_id_QA) . "\n";
	print "--------------------------------------------------- A: $A_acc\n"
		. " " . sprintf("%.2f", $per_id_QLA) . "                                " . sprintf("%.2f", $per_id_QRA) . "\n"
		. "~~~~~~~~~~~~~~~~~~~~~~~~\\ /~~~~~~~~~~~~~~~~~~~~~~~~ Q: $query_acc\n"
		. "DivR: " . sprintf("%.3f", $div_ratio_QLA_QRB) . " BS: " . sprintf("%.2f", $BS_A) . "     |\n"
		. "Per_id(QLA,QRB): " . sprintf("%.2f", $per_id_QLA_QRB) . "   |\n"
		. "                         |\n"
		. "   (L-AB: " . sprintf("%.2f", $per_id_LAB) . ")         |      (R-AB: " . sprintf("%.2f", $per_id_RAB) . ")\n"
		. "   WinL:$win_left_end5-$win_left_end3            |      WinR:$win_right_end5-$win_right_end3\n"
		. "                         |\n"
		. "Per_id(QLB,QRA): " . sprintf("%.2f", $per_id_QLB_QRA) . "   |\n"
		. "DivR: " . sprintf("%.3f", $div_ratio_QRA_QLB) . " BS: " . sprintf("%.2f", $BS_B) . "    |\n"
		. "~~~~~~~~~~~~~~~~~~~~~~~~/ \\~~~~~~~~~~~~~~~~~~~~~~~~~ Q: $query_acc\n"
		. " " . sprintf("%.2f", $per_id_QLB) . "                                " . sprintf("%.2f", $per_id_QRB) . "\n"
		. "---------------------------------------------------- B: $B_acc\n";
	print "            Per_id(Q,B): ". sprintf("%.2f", $per_id_QB) . "\n\n";
	
	my $deltaL = $per_id_QLA - $per_id_QLB;
	my $deltaR = $per_id_QRA - $per_id_QRB;

	print "DeltaL: " . sprintf("%.2f", $deltaL) . "                   DeltaR: " . sprintf("%.2f", $deltaR) . "\n\n";
	
	unless ($printAlignmentsFlag) { return; }
	
	
	## build the left windows:
	my @Q_left_win = @Q_chars[$win_left_end5..$win_left_end3];
	my @A_left_win = @A_chars[$win_left_end5..$win_left_end3];
	my @B_left_win = @B_chars[$win_left_end5..$win_left_end3];
	
	&print_alignment($A_acc, \@A_left_win, 
					 $query_acc, \@Q_left_win, 
					 $B_acc, \@B_left_win);
	
	print "\t\t** Breakpoint **\n\n";
	
	my @Q_right_win = @Q_chars[$win_right_end5..$win_right_end3];
	my @A_right_win = @A_chars[$win_right_end5..$win_right_end3];
	my @B_right_win = @B_chars[$win_right_end5..$win_right_end3];
	
	&print_alignment($A_acc, \@A_right_win, 
					 $query_acc, \@Q_right_win, 
					 $B_acc, \@B_right_win);
	
	return;
}


####
		
	*/	
				
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void ChimeraSlayer::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		chimeraResults.resize(numSeqs);
		chimeraFlags.resize(numSeqs, "no");
		
		//break up file if needed
		int linesPerProcess = numSeqs / processors ;
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			//find breakup of sequences for all times we will Parallelize
			if (processors == 1) {   lines.push_back(new linePair(0, numSeqs));  }
			else {
				//fill line pairs
				for (int i = 0; i < (processors-1); i++) {			
					lines.push_back(new linePair((i*linesPerProcess), ((i*linesPerProcess) + linesPerProcess)));
				}
				//this is necessary to get remainder of processors / numSeqs so you don't miss any lines at the end
				int i = processors - 1;
				lines.push_back(new linePair((i*linesPerProcess), numSeqs));
			}
		#else
			lines.push_back(new linePair(0, numSeqs));
		#endif
		
		if (seqMask != "") {	decalc = new DeCalculator();	} //to use below
		
		//referenceSeqs, numWanted, matchScore, misMatchPenalty, divR, minSimilarity
		maligner = new Maligner(templateSeqs, numWanted, match, misMatch, 1.01, minSim);
		slayer = new Slayer(window, increment, minSim, divR);
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			string chimeraFlag = maligner->getResults(querySeqs[i]);
			float percentIdentical = maligner->getPercentID();
			vector<results> Results = maligner->getOutput();
			
			//cout << querySeqs[i]->getName() << '\t' << chimeraFlag << '\t' << percentIdentical << endl;
			
			for (int j = 0; j < Results.size(); j++) {
				//cout << "regionStart = " << Results[j].regionStart << "\tRegionEnd = " << Results[j].regionEnd << "\tName = " << Results[j].parent << "\tPerQP = " << Results[j].queryToParent << "\tLocalPerQP = " << Results[j].queryToParentLocal << "\tdivR = " << Results[j].divR << endl;
			}
			
			if (chimeraFlag == "yes") {
			
				//get sequence that were given from maligner results
				vector<SeqDist> seqs;
				for (int j = 0; j < Results.size(); j++) {
					Sequence* seq = getSequence(Results[j].parent); //makes copy so you can filter and mask and not effect template
					
					//seq = NULL if error occurred in getSequence
					if (seq == NULL) {  break;	}
					else {	
						SeqDist member;
						member.seq = seq;
						member.dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
						seqs.push_back(member);	
					}
				}
			
				//limit number of parents to explore - default 5
				if (Results.size() > parents) {
					//sort by distance
					sort(seqs.begin(), seqs.end(), compareSeqDist);
					//prioritize larger more similiar sequence fragments
					reverse(seqs.begin(), seqs.end());
					
					for (int k = seqs.size()-1; k > (parents-1); k--)  {  
						delete seqs[k].seq;
						seqs.pop_back();	
					}
				}
		
				//put seqs into vector to send to slayer
				vector<Sequence*> seqsForSlayer;
				for (int k = 0; k < seqs.size(); k++) {  seqsForSlayer.push_back(seqs[k].seq);	}
			
				//mask then send to slayer...
				if (seqMask != "") {
					decalc->setMask(seqMask);

					//mask querys
					decalc->runMask(querySeqs[i]);
					
					//mask parents
					for (int k = 0; k < seqsForSlayer.size(); k++) {
						decalc->runMask(seqsForSlayer[k]);
					}
					
				}
				
				//send to slayer
				chimeraFlags[i] = slayer->getResults(querySeqs[i], seqsForSlayer);
				chimeraResults[i] = slayer->getOutput();
			
				//free memory
				for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
			}
			
		}	
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];	}
		
		if (seqMask != "") {
			delete decalc; 
		}

			
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* ChimeraSlayer::getSequence(string name) {
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
		
		if(spot == -1) { mothurOut("Error: Could not find sequence in chimeraSlayer."); mothurOutEndLine(); return NULL; }
		
		temp = new Sequence(templateSeqs[spot]->getName(), templateSeqs[spot]->getAligned());
		
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "getSequence");
		exit(1);
	}
}
//***************************************************************************************************************



