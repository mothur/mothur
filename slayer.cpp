/*
 *  slayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "slayer.h"

/***********************************************************************/
Slayer::Slayer(int win, int increment, int parentThreshold, float div, int i, int snp, int mi) :
		minBS(mi), windowSize(win), windowStep(increment), parentFragmentThreshold(parentThreshold), divRThreshold(div), iters(i), percentSNPSample(snp){ m = MothurOut::getInstance(); }
/***********************************************************************/
string Slayer::getResults(Sequence* query, vector<Sequence*> refSeqs) {
	try {
		vector<data_struct> all; all.clear();
		myQuery = *query;
		/*cout << "here" << endl;	
		query->setName("S000381746"); query->setAligned("...............................................................................................................................................A-C-GC--TGG-C--G-GC-A-GG--C----C-T--AACACA-T-GC-A-AGT-CGA-G-CG----------G-CAG-CG-G---------------------------GA-GG-A-AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTT-CCTC----------------G-CC--G--GC--G--AG-C-GG-C-GG-A--C-------------GGG-TGAGT-A--AT-GT-C-T-G-GG---G-A--T-CT-G--C-C-CGA--TG-G------------------------------------------------------------------A-GG----GGG-AT-AA-CCA-------------------------C-T-G-----------------------GAA-A---CGG-TGG-CTAA-TA---CC-G--C-AT-A----------A--------------------C-------------------------------------GT-C-----------------------------------------------------------------------------------------------------------------------G-CA-A--------------------------------------------------------------------------------------------------------------------------------------G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-T-G-GG-----G--GA-C--C--------------------------------------------------------------------------------------------------------------------TTC-G----------------------------------------------------------------------------------------------------------------------G-G--CC-TC--A---C-A--------------C----C-A---T-CG-G---AT---G-A-----A-CCC-AGA--T-GGG--A------TT--A--G-CT-A----G---TAGG-T-G-GG-G-T----AAT-GG-C-T-C-ACCT--A-GG-C-G--A-CG-A------------TCC-C-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AC-C-AG-CCAC-A-CTGGA--A-C-TG-A-GA-C-AC-G-G-TCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GC-GC-A----A-G-CC-T-GA-TG-CA-GCCA-TGCC-G-CG-T---G-T-G--T--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-T-G-T--A---AA-G-CAC--------TT-TC-A-G--C-GAG----GA-G--G---AA-GGTG---GTGA-GC----T--T--AA-T---A----------CG-CTCAT-CAA-TT-GA-CG-TT-A-C-TC-G-CA-G---------AA-----------GAAGC-ACC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GGT-GCA-A-G-CG-TTAA-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GCGC-AC--G-CA-G-G-C-G------------G--T-TT-G-T-T-AA----G-T-C-A---G-ATG-TG-A-AA-TC--CC-CGA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-T-G-GG-AA-C----T-G-C-A-T-T--------T--GA-A-A-C-T-G-GCA--A-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-T-C-----T-CG--TA-G-A------------G-GG-G-GG-T----AG--AATT-CCA-G-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TC-T-GGA--GG-A-AT-A-CC-GG--T--G--GC-GAA-G--G-C---G----G--C-C-CCCTG------G-AC-GA--------------------------------------------------------------AG-A-C-T--GA--CG-----CT-CA-GG--T-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-TC--GA-TT---------T-GG--A--G-G-TT-G-TG-C--C--------------------------------------------------------------------------------------CTT-GA--------------------------------------------------------------------------------------------------------------------------------------------------G-G-C-GT--G-G-C-T-TC-C------GG--A----GC-TAA--CG-C-G-T--T--AA-AT--C----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGT-T--AAA-ACTC-AAA---------TGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-TACTC-TT-G-AC-A-T-C--------------CAG-A-G-------------A-AC-T-T-T--CC--A-GA-G-A-T--G-G-A--T-T-G-G--T-G-----CC-------------------------------------T--TC-G------------------------------------------GG----A----A---CT-CTG---A--GA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TT-G--TGA-AA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TC--C-TTTG--T-T-G-C-C---AG-C-G--G--T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G----C----C-G------------G----G---A-A--CT---------------C-A-A-A-G-GA-G--AC-T-G-CCA--G-T------------------------------------G-A---TAA----------------------------------A-C-T-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-G-C-C-CTT----AC-G--AG-T-A-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-CGTA-T-A--C-AAA-GA-GA--------------------------------------------------------------------------------------------------A-G-C-G-A--ACCT-G-C--G---------------------------------------A-GG-G-C-----------A--A-G-CG---G----------A--CCT-C------A-T-AAAGT-AC-G-T-C-G-TAG-TCC--------GGA-T-TGGAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTCC-A-T-G-AA-G-TC-GGAAT-CG-C-TA--G-TA-AT-C-G-T----AGA-TC-A-G--A------AT--GCT-AC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---TGG-G-TT-GC-AAA--A-GAA------G--T-AGG-TA-G-C-T-T-AA-C-C-------------------------------------------------------------T-TC-G------------------------------------------------------------------------------------------------------GG-A--GG-G--C---GC-TTA--CC--ACT-T----T-GT..........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................");
		refSeqs.resize(2);
		refSeqs[1]->setName("S000381740"); refSeqs[1]->setAligned("...............................................................................................................................................a-c-gc--tgg-c--g-gc-a-gg--c----c-t--aacaca-t-gc-a-agt-cga-g-cg----------g-tag-ca-c----------------------------agga-g-ag----------------------------------------------------ctt-g----------------------------------------------------------------------------------ctc-tctg----------------g-gt--g--ac--g--ag-c-gg-c-gg-a--c-------------ggg-tgagt-a--at-gt-c-t-g-gg---a-a--a-ct-g--c-c-tga--tg-g------------------------------------------------------------------a-gg----ggg-at-aa-cta-------------------------c-t-g-----------------------gaa-a---cgg-tag-ctaa-ta---cc-g--c-at-a----------a--------------------c-------------------------------------gt-c-----------------------------------------------------------------------------------------------------------------------t-ac-g--------------------------------------------------------------------------------------------------------------------------------------g-a-c---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------caaa--g-t-g-gg-----g--ga-c--c--------------------------------------------------------------------------------------------------------------------ttc-g----------------------------------------------------------------------------------------------------------------------g-g--cc-tc--a---c-g--------------c----c-a---t-ca-g---at---g-t-----g-ccc-aga--t-ggg--a------tt--a--g-ct-a----g---tagg-t-g-gg-g-t----aat-gg-c-t-c-acct--a-gg-c-g--a-cg-a------------tcc-c-t------ag-ct-g-g-tct-g-ag----a--gg-at--g-ac-c-ag-ccac-a-ctgga--a-c-tg-a-ga-c-ac-g-g-tccaga-ctcc-tac-g--g-g-a-g-gc-a-gc-a-g-tg---gg-g-a-ata-ttgca-c-aa-t-gg--gc-gc-a----a-g-cc-t-ga-tg-ca-gcca-tgcc-g-cg-t---g-t-g--t--ga-a-g--a--a-g-g-cc-----tt-cg---------g-g-t-t-g-t--a---aa-g-cac--------tt-tc-a-g--c-gag----ga-g--g---aa-gggc---gatg-tc----t--t--aa-t---a----c-----gg-c-agc-gca-tt-ga-cg-tt-a-c-tc-g-ca-g---------aa-----------gaagc-acc-gg-c-taa---c--t-ccgt--gcca--g-c---a--gccg---c-gg--ta-at--ac---gg-ag-ggt-gca-a-g-cg-ttaa-t-cgg-aa-tt-a--c-t--gggc-gta----aa-gcgc-ac--g-ca-g-g-c-g------------g--t-tt-g-t-t-aa----g-t-c-a---g-atg-tg-a-aa-tc--cc-cgc-g--------------------------------------------------------------------ct-t-aa-------------------------------------------------------------------------cg-t-g-gg-aa-c----t-g-c-a-t-t--------t--ga-a-a-c-t-g-gca--a-g-c---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------t-a-g-a-g-t-c-----t-cg--ta-g-a------------g-gg-g-gg-t----ag--aatt-cca-g-gt--gt-a-gcg-gtgaaa-tg-cgt-agag-a-tc-t-gga--gg-a-at-a-cc-gg--t--g--gc-gaa-g--g-c---g----g--c-c-ccctg------g-ac-ga--------------------------------------------------------------ag-a-c-t--ga--cg-----ct-ca-gg--t-g-cga--aa-g-c--------------g-tggg-gag-c-a-aaca--gg-atta-g-ata-c-----cc-t-g-gta-g-t----c-ca--c-g-ctg-t-aaa--c-gatg-tc--ga-tt---------t-gg--a--g-g-tt-g-tg-c--c--------------------------------------------------------------------------------------ctt-ga--------------------------------------------------------------------------------------------------------------------------------------------------g-g-c-gt--g-g-c-t-tc-c------gg--a----gc-taa--cg-c-g-t--t--aa-at--c----g-acc-gcc-t-g-gg-gag-ta---cgg-----c-c--g-c-a-a-ggt-t--aaa-actc-aaa---------tgaa-ttg-acggg-g-g-cccg----c-a--c-a-a-gcg-gt-g--g--ag-ca-t--gt-ggt-tt-aatt-c-g-atg-caac-g-cg-a-ag-a-a-cc-tt-a-cc-tactc-tt-g-ac-a-t-c--------------cag-a-g-------------a-ac-t-t-t--cc--a-ga-g-a-t--g-g-a--t-t-g-g--t-g-----cc-------------------------------------t--tc-g------------------------------------------gg----a----a---ct-ctg---a--ga---------------------------------------------------c-a-g-g-t-gctg-ca-tgg-ct--gtc-gtc-a-gc-tc---g-tg-tt-g--tga-aa-tgt-t-gg-g-tt-aa-gt-cccgc-aa--------c-gag-cgc-a-acc-c-t-ta--tc--c-tttg--t-t-g-c-c---ag-c-g--a--t-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------tcg------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------g----t----c-g------------g----g---a-a--ct---------------c-a-a-a-g-ga-g--ac-t-g-ccg--g-t------------------------------------g-a---taa----------------------------------a-c-c-g--g-a-gg-a--agg-t--gggg-a-tgac-gtc--aagt-c---atc-a-t-g-g-c-c-ctt----ac-g--ag-t-a-gg-gc-ta-cac-acgtg-c--ta--caatg---g-cgta-t-a--c-aaa-ga-ga--------------------------------------------------------------------------------------------------a-g-c-g-a--a-ctcg-c--g---------------------------------------a-ga-g-c-----------a--a-g-cg---g----------a--cct-c------a-t-aaagt-ac-g-t-c-g-tag-tcc--------gga-t-tggag-tc--t-gcaa-ct-c-------------------------------------------------------------------------------------------------g-actcc-a-t-g-aa-g-tc-ggaat-cg-c-ta--g-ta-at-c-g-t----aga-tc-a-g--a------at--gct-ac-g-gt-g-aat-acgt-t-cccgggcct-tgta----cacaccg-ccc-gtc-----a---ca--cca-tg-gg-a--g---tgg-g-tt-gc-aaa--a-gaa------g--t-agg-ta-g-c-t-t-aa-c-c-------------------------------------------------------------t-tc-g------------------------------------------------------------------------------------------------------gg-a--gg-g--c---gc-tta--cc--act-t----t-gtg-at-tca------------------------t.........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................");
		refSeqs[0]->setName("7000004131500404"); refSeqs[0]->setAligned(".........................................................................................................................................AT-TGAA-C-GC--TGG-C--G-GC-A-GG--C----C-T--AACACA-T-GC-A-AGT-CGA-G-CG----------G-CAG-CG-G----------------------------AAAG-A-AG----------------------------------------------------CTT-G---------------------------------------------------------------------------------ACTT-CTTT----------------G-CC--G--GC--G--AG-C-GG-C-GG-A--C-------------GGG-TGAGT-A--AT-GT-C-T-G-GG---G-A--T-CT-G--C-C-CGA--TG-G------------------------------------------------------------------A-GG----GGG-AT-AA-CTA-------------------------C-T-G-----------------------GAA-A---CGG-TAG-CTAA-TA---CC-G--C-AT-A----------A--------------------C-------------------------------------GT-C-----------------------------------------------------------------------------------------------------------------------G-CA-A--------------------------------------------------------------------------------------------------------------------------------------G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-T-G-GG-----G--GA-C--C--------------------------------------------------------------------------------------------------------------------TTC-G----------------------------------------------------------------------------------------------------------------------G-G--CC-TC--A---C-A--------------C----C-A---T-CG-G---AT---G-A-----A-CCC-AGA--T-GGG--A------TT--A--G-CT-A----G---TAGG-T-G-GG-G-T----AAT-GG-C-T-C-ACCT--A-GG-C-G--A-CG-A------------TCC-C-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AC-C-AG-CCAC-A-CTGGA--A-C-TG-A-GA-C-AC-G-G-TCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GG-GA-A----A-C-CC-T-GA-TG-CA-GCCA-TGCC-G-CG-T---G-T-G--T--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-T-G-T--A---AA-G-CAC--------TT-TC-A-G--C-GGG----GA-A--G---AA-GGCG---TT-A-GC---GT--T--AA-C---A----G-----CG-C-TAT-CGA-TT-GA-CG-TT-A-C-CT-G-CA-G---------AA-----------GAAGC-ACC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GGT-GCA-A-G-CG-TTAA-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GCGT-AC--G-CA-G-G-C-G------------G--T-CT-G-T-T-AA----G-T-C-A---G-ATG-TG-A-AA-TC--CC-CGG-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CC-T-G-GG-AA-C----T-G-C-A-T-T--------T--GA-A-A-C-T-G-GCA--G-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-T-C-----T-CG--TA-G-A------------G-GG-G-GG-T----AG--AATT-CCA-G-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TC-T-GGA--GG-A-AT-A-CC-GG--T--G--GC-GAA-G--G-C---G----G--C-C-CCCTG------G-AC-GA--------------------------------------------------------------AG-A-C-T--GA--CG-----CT-CA-GG--T-A-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CTG-T-AAA--C-GATG-TC--GA-TT---------T-GA--A--G-G-TT-G-TG-G--C--------------------------------------------------------------------------------------CTT-GA--------------------------------------------------------------------------------------------------------------------------------------------------G-C-T-GT--G-G-C-T-TT-C------GG--A----GC-TAA--CG-C-G-T--T--AA-AT--C----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGT-T--AAA-ACTC-AAA---------TGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-TACTC-TT-G-AC-A-T-C--------------CAG-A-G-------------A-AC-T-T-G--GC--A-GA-G-A-T--G-C-C--T-T-G-G--T-G-----CC-------------------------------------T--TC-G------------------------------------------GG----A----G---CT-CTG---A--GA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TT-G--TGA-AA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TC--C-TTTG--T-T-G-C-C---AG-C-G--A--T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G----T----C-G------------G----G---A-A--CT---------------C-A-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---TAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-G-C-C-CTT----AC-G--AG-T-A-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-CGCA-T-A--C-AAA-GA-GA--------------------------------------------------------------------------------------------------A-G-C-G-A--T-CTCG-C--G---------------------------------------A-GA-G-T-----------C--A-G-CG---G----------A--CCT-C------A-C-AAAGT-GC-G-T-C-G-TAG-TCC--------GGA-T-TGGAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTCC-A-T-G-AA-G-TC-GGAAT-CG-C-TA--G-TA-AT-C-G-T----GGA-TC-A-G--A------AT--GCC-AC-G-GT-G-AAT-ACGT-T-CCTGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---TGG-G-TT-GC-AAA--A-GAA------G--T-AGG-TA-G-C-T-T-AA-C-C-------------------------------------------------------------T-TC-G------------------------------------------------------------------------------------------------------GG-A--GG-G--C---GC-TTA--CC--ACT-T----T-GTG-AT-TCA------------------------TG--ACT-GGGG-TG-AAG-TCGTAACAA-GGTAA-CCGT-AGGGGAA-CCT......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................");		
		for (int i = 0; i < refSeqs.size(); i++) {
			string newAligned = "";
			for (int j = 0; j < refSeqs[i]->getAligned().length(); j++) {
				newAligned += toupper(refSeqs[i]->getAligned()[j]);
			}
			refSeqs[i]->setAligned(newAligned);
		}*/
		
		for (int i = 0; i < refSeqs.size(); i++) {
		
			for (int j = i+1; j < refSeqs.size(); j++) {
			
				if (m->control_pressed) { return "no";  }
	
				//make copies of query and each parent because runBellerophon removes gaps and messes them up
				Sequence* q = new Sequence(query->getName(), query->getAligned());
				Sequence* leftParent = new Sequence(refSeqs[i]->getName(), refSeqs[i]->getAligned());
				Sequence* rightParent = new Sequence(refSeqs[j]->getName(), refSeqs[j]->getAligned());

				map<int, int> spots;  //map from spot in original sequence to spot in filtered sequence for query and both parents
				vector<data_struct> divs = runBellerophon(q, leftParent, rightParent, spots);
				//cout << divs.size() << endl;
				if (m->control_pressed) { 
					delete q;
					delete leftParent;
					delete rightParent;
					return "no"; 
				}
					
//				cout << divs.size() << endl;
				vector<data_struct> selectedDivs;
				for (int k = 0; k < divs.size(); k++) {
					
					vector<snps> snpsLeft = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winLStart, divs[k].winLEnd);
					vector<snps> snpsRight = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winRStart, divs[k].winREnd);
					
					if (m->control_pressed) { 
						delete q;
						delete leftParent;
						delete rightParent;
						return "no"; 
					}
					
					int numSNPSLeft = snpsLeft.size();
					int numSNPSRight = snpsRight.size();
					
					//require at least 3 SNPs on each side of the break
//					if ((numSNPSLeft >= 3) && (numSNPSRight >= 3)) {
					if ((numSNPSLeft >= 4) && (numSNPSRight >= 4)) {
					
						//removed in 12/09 version of chimeraSlayer
						//int winSizeLeft = divs[k].winLEnd - divs[k].winLStart + 1;
						//int winSizeRight = divs[k].winREnd - divs[k].winRStart + 1;
						
						//float snpRateLeft = numSNPSLeft / (float) winSizeLeft;
						//float snpRateRight = numSNPSRight / (float) winSizeRight;
						//float logR = log(snpRateLeft / snpRateRight) / log(2.0); 
						
						// do not accept excess snp ratio on either side of the break
						//if (abs(logR) < 1 ) {  
							
							float BS_A, BS_B;
							bootstrapSNPS(snpsLeft, snpsRight, BS_A, BS_B, iters);
							
							if (m->control_pressed) { 
								delete q;
								delete leftParent;
								delete rightParent;
								return "no"; 
							}

							divs[k].bsa = BS_A;
							divs[k].bsb = BS_B;
						
							divs[k].bsMax = max(BS_A, BS_B);
						
							divs[k].chimeraMax = max(divs[k].qla_qrb, divs[k].qlb_qra);
							
							//so results reflect orignal alignment
							divs[k].winLStart = spots[divs[k].winLStart];
							divs[k].winLEnd = spots[divs[k].winLEnd];  
							divs[k].winRStart = spots[divs[k].winRStart]; 
							divs[k].winREnd = spots[divs[k].winREnd]; 
						
							selectedDivs.push_back(divs[k]);
						//}
					}
				}
				
				//save selected
				for (int m = 0; m < selectedDivs.size(); m++) {  all.push_back(selectedDivs[m]);	}
				
				delete q;
				delete leftParent;
				delete rightParent;
			}
		}
		
		// compute bootstrap support
		if (all.size() > 0) {
			//sort them
			sort(all.begin(), all.end(), compareDataStruct);
			reverse(all.begin(), all.end());
			
			outputResults = all;
			return "yes"; 
		}
		else {
			outputResults = all;
			return "no";
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "getResults");
		exit(1);
	}
}
/***********************************************************************/
vector<data_struct> Slayer::runBellerophon(Sequence* q, Sequence* pA, Sequence* pB, map<int, int>& spots) {
	try{
		
		vector<data_struct> data;
		
		//vertical filter
		vector<Sequence*> temp;
		temp.push_back(q); temp.push_back(pA); temp.push_back(pB);
		
		//maps spot in new alignment to spot in alignment before filter
		spots = verticalFilter(temp);  //fills baseSpots
		
		//get these to avoid numerous function calls
		string query = q->getAligned();
		string parentA = pA->getAligned();
		string parentB = pB->getAligned();
		int length = query.length();
/*cout << q->getName() << endl << q->getAligned() << endl << endl;	
cout << pA->getName() << endl << pA->getAligned() << endl << endl;		
cout << pB->getName() << endl << pB->getAligned() << endl << endl;	
cout << " length = " << length << endl;
cout << q->getName() << endl;
cout << pA->getName() << '\t';
cout << pB->getName() << endl;*/
	
		//check window size
		if (length < (2*windowSize+windowStep)) { 
//			m->mothurOut("Your window size is too large for " + q->getName() + ". I will make the window size " + toString(length/4) + " which is 1/4 the filtered length."); m->mothurOutEndLine();	
			windowSize = length / 4;
		}
		
		for (int i = windowSize-1; i <= (length - windowSize); i += windowStep) {
		
			if (m->control_pressed) { return data; }
		
			int breakpoint = i;
			int leftLength = breakpoint + 1;
			int rightLength = length - leftLength;
			
			float QLA = computePercentID(query, parentA, 0, breakpoint);
			float QRB = computePercentID(query, parentB, breakpoint+1, length - 1);
		
			float QLB = computePercentID(query, parentB, 0, breakpoint);
			float QRA = computePercentID(query, parentA, breakpoint+1, length - 1);
		
			float LAB = computePercentID(parentA, parentB, 0, breakpoint);
			float RAB = computePercentID(parentA, parentB, breakpoint+1, length - 1);
		
			float AB = ((LAB*leftLength) + (RAB*rightLength)) / (float) length;
			float QA = ((QLA*leftLength) + (QRA*rightLength)) / (float) length;
			float QB = ((QLB*leftLength) + (QRB*rightLength)) / (float) length;
		
			float QLA_QRB = ((QLA*leftLength) + (QRB*rightLength)) / (float) length;
			float QLB_QRA = ((QLB*leftLength) + (QRA*rightLength)) / (float) length;
		
			//in original and not used
			//float avgQA_QB = ((QA*leftLength) + (QB*rightLength)) / (float) length;
		
			float divR_QLA_QRB = min((QLA_QRB/QA), (QLA_QRB/QB));
			float divR_QLB_QRA = min((QLB_QRA/QA), (QLB_QRA/QB));
//cout << leftLength << '\t' << rightLength << '\t' << QLA << '\t' << QRB << '\t' << QLB << '\t' << QRA  << '\t' << LAB << '\t' << RAB << '\t' << AB << '\t' << QA << '\t' << QB << '\t' << QLA_QRB << '\t' <<  QLB_QRA <<    endl;    		

//cout << divRThreshold << endl;
//cout << breakpoint << '\t' << divR_QLA_QRB << '\t' << divR_QLB_QRA << endl;
			//is one of them above the 
			if (divR_QLA_QRB >= divRThreshold || divR_QLB_QRA >= divRThreshold) {
				
				if (((QLA_QRB > QA) && (QLA_QRB > QB) && (QLA >= parentFragmentThreshold) && (QRB >= parentFragmentThreshold))  ||
					((QLB_QRA > QA) && (QLB_QRA > QB) && (QLB >=parentFragmentThreshold) && (QRA >= parentFragmentThreshold)))  {
					
					data_struct member;
					
					member.divr_qla_qrb = divR_QLA_QRB;
					member.divr_qlb_qra = divR_QLB_QRA;
					member.qla_qrb = QLA_QRB;
					member.qlb_qra = QLB_QRA;
					member.qla = QLA;
					member.qrb = QRB;
					member.ab = AB; 
					member.qa = QA;
					member.qb = QB; 
					member.lab = LAB; 
					member.rab = RAB; 
					member.qra = QRA; 
					member.qlb = QLB; 
					member.winLStart = 0;
					member.winLEnd = breakpoint;  
					member.winRStart = breakpoint+1; 
					member.winREnd = length-1; 
					member.querySeq = *(q); 
					member.parentA = *(pA);
					member.parentB = *(pB);
					member.bsa = 0;
					member.bsb = 0;
					member.bsMax = 0;
					member.chimeraMax = 0;
					
					data.push_back(member);
					
				}//if
			}//if
		}//for
		
		return data;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "runBellerophon");
		exit(1);
	}
}
/***********************************************************************/
vector<snps> Slayer::getSNPS(string parentA, string query, string parentB, int left, int right) {
	try {
	
		vector<snps> data;
//cout << left << '\t' << right << endl;
		for (int i = left; i <= right; i++) {
			
			char A = parentA[i];
			char Q = query[i];
			char B = parentB[i];
			
			if ((A != Q) || (B != Q)) {
//cout << "not equal " << Q << '\t' << A << '\t' << B << endl;
			
				//ensure not neighboring a gap. change to 12/09 release of chimeraSlayer - not sure what this adds, but it eliminates alot of SNPS
				if (
					//did query loose a base here during filter??
					( i == 0 || abs (baseSpots[0][i] - baseSpots[0][i-1]) == 1) &&
					( i == query.length() || abs (baseSpots[0][i] - baseSpots[0][i+1]) == 1)
					&&
					//did parentA loose a base here during filter??
					( i == 0 || abs (baseSpots[1][i] - baseSpots[1][i-1]) == 1) &&
					( i == parentA.length() || abs (baseSpots[1][i] - baseSpots[1][i+1]) == 1) 
					&&
					//did parentB loose a base here during filter??
					( i == 0 || abs (baseSpots[2][i] - baseSpots[2][i-1]) == 1) &&
					( i == parentB.length() || abs (baseSpots[2][i] - baseSpots[2][i+1]) == 1)
					) 
				{ 
				
					snps member;
					member.queryChar = Q;
					member.parentAChar = A;
					member.parentBChar = B;
//cout << "not neighboring a gap " << Q << '\t' << A << '\t' << B << '\t' << baseSpots[0][i] << '\t' << baseSpots[0][i+1] << '\t' << baseSpots[0][i-1] << '\t' << baseSpots[1][i] << '\t' << baseSpots[1][i+1] << '\t' << baseSpots[1][i-1] << '\t' << baseSpots[2][i] << '\t' << baseSpots[2][i+1] << '\t' << baseSpots[2][i-1] << endl;				
					data.push_back(member);
				}
			}
//			cout << i << '\t' << data.size() << endl;
		}
		
		return data;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "getSNPS");
		exit(1);
	}
}
/***********************************************************************/
int Slayer::bootstrapSNPS(vector<snps> left, vector<snps> right, float& BSA, float& BSB, int numIters) {
	try {

		srand((unsigned)time( NULL ));

		int count_A = 0; // sceneario QLA,QRB supported
		int count_B = 0; // sceneario QLB,QRA supported
	
		int numLeft = max(1, int(left.size() * percentSNPSample/(float)100 + 0.5));
		int numRight = max(1, int(right.size() * percentSNPSample/(float)100 + 0.5));

		for (int i = 0; i < numIters; i++) {
			//random sampling with replacement.
		
			if (m->control_pressed) { return 0;  }
			
			vector<snps> selectedLeft;

			for (int j = 0; j < numLeft; j++) {
				int index = int(rand() % left.size());
				selectedLeft.push_back(left[index]);
			}

			vector<snps> selectedRight;
			for (int j = 0; j < numRight; j++) {
				int index = int(rand() % right.size());
				selectedRight.push_back(right[index]);
			}
		
			/* A  ------------------------------------------
			#       QLA                     QRA
			# Q  ------------------------------------------
			#                      |
			#                      |
			# Q  ------------------------------------------
			#       QLB                     QRB
			# B  ------------------------------------------ */
		
		
			float QLA = snpQA(selectedLeft);
			float QRA = snpQA(selectedRight);
		
			float QLB = snpQB(selectedLeft);
			float QRB = snpQB(selectedRight);
	
			//in original - not used - not sure why?
			//float ALB = snpAB(selectedLeft);
			//float ARB = snpAB(selectedRight);
		
			if ((QLA > QLB) && (QRB > QRA)) {
				count_A++;
			}
		
			if ((QLB > QLA) && (QRA > QRB)) {
				count_B++;
			}
			
//cout << "selected left snp: \n";
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].parentAChar;  } 
//cout << endl;
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].queryChar;  }
//cout << endl;
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].parentBChar;  }
//cout << endl;
//cout << "selected right snp: \n";
//for (int j = 0; j < selectedRight.size(); j++) {  cout << selectedRight[j].parentAChar;  } 
//cout << endl;
//for (int i = 0; i < selectedRight.size(); i++) {  cout << selectedRight[i].queryChar;  }
//cout << endl;
//for (int i = 0; i < selectedRight.size(); i++) {  cout << selectedRight[i].parentBChar;  }
//cout << endl;		
		}




		BSA = (float) count_A / (float) iters * 100;
		BSB = (float) count_B / (float) iters * 100;
//cout << "bsa = " << BSA << " bsb = " << BSB << endl;
		
		//run borderline bootstrap values longer
		//if (numIters < 1000) {
			//are you within 10 points of min bootstrap value cutoff
		//	if (((abs((double)(BSA - minBS))) <= 5) || ((abs((double)(BSB - minBS))) <= 5)) {
		//		m->mothurOut("extending bootstrap for " + myQuery.getName()); m->mothurOutEndLine();
		//		bootstrapSNPS(left, right, BSA, BSB, 1000);
		//	}
		//}

		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "bootstrapSNPS");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpQA(vector<snps> data) {
	try {
	
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentAChar == data[i].queryChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpQA");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpQB(vector<snps> data) {
	try {
	
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentBChar == data[i].queryChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpQB");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpAB(vector<snps> data) {
	try {
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentAChar == data[i].parentBChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpAB");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::computePercentID(string queryFrag, string parent, int left, int right) {
	try {
		int total = 0;
		int matches = 0;

		for (int i = left; i <= right; i++) {
			total++;
			if (queryFrag[i] == parent[i]) {
				matches++;
			}
		}

		float percentID =( matches/(float)total) * 100;
		
		return percentID;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "computePercentID");
		exit(1);
	}
}
/***********************************************************************/
//remove columns that contain any gaps
map<int, int> Slayer::verticalFilter(vector<Sequence*> seqs) {
	try {
		//find baseSpots
		baseSpots.clear(); 
		baseSpots.resize(3);  //query, parentA, parentB
	
		vector<int> gaps;	gaps.resize(seqs[0]->getAligned().length(), 0);
		
		string filterString = (string(seqs[0]->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{   gaps[j]++;	}
			}
		}
		
		//zero out spot where any sequences have blanks
		int numColRemoved = 0;
		int count = 0;
		map<int, int> maskMap; maskMap.clear();

		for(int i = 0; i < seqs[0]->getAligned().length(); i++){
			if(gaps[i] != 0)	{	filterString[i] = '0'; 	numColRemoved++;  }
			else {
				maskMap[count] = i;
				count++;
			}
		}

		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			string newAligned = "";
			
			int baseCount = 0;
			int count = 0;
			for (int j = 0; j < seqAligned.length(); j++) {
				//are you a base
				if ((seqAligned[j] != '-') && (seqAligned[j] != '.') && (toupper(seqAligned[j]) != 'N'))	{ baseCount++; }
			
				//if this spot is not a gap
				if (filterString[j] == '1') { 
					newAligned += seqAligned[j]; 
					baseSpots[i][count] = baseCount;
					count++;
				}
			}
			
			seqs[i]->setAligned(newAligned);
		}
		
		return maskMap;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "verticalFilter");
		exit(1);
	}
}
/***********************************************************************/
