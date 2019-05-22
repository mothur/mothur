/*
 *  validcalculator.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validcalculator.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "simpsoneven.h"
#include "invsimpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "smithwilson.h"
#include "heip.h"
#include "shannoneven.h"
#include "jackknife.h"
#include "geom.h"
#include "qstat.h"
#include "logsd.h"
#include "bergerparker.h"
#include "bstick.h"
#include "goodscoverage.h"
#include "efron.h"
#include "boneh.h"
#include "solow.h"
#include "shen.h"
#include "coverage.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharednseqs.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"
#include "sharedsobs.h"
#include "sharednseqs.h"
#include "sharedjsd.h"
#include "sharedrjsd.h"
#include "shannonrange.h"



/********************************************************************/
ValidCalculators::ValidCalculators() {
	try {
		m = MothurOut::getInstance();
		
		initialSingle();  
		initialShared();
		initialRarefaction();
		initialSharedRarefact();
		initialSummary();
		initialSharedSummary();
		initialVennSingle();
		initialVennShared();
		initialTreeGroups();
		initialDistance();
		initialMatrix();
		initialHeat();
        initialEstimators();
		
		for(it = single.begin(); it != single.end(); it++) { allCalcs.insert(*it); }
		for(it = shared.begin(); it != shared.end(); it++) { allCalcs.insert(*it); }
		for(it = rarefaction.begin(); it != rarefaction.end(); it++) { allCalcs.insert(*it); }
		for(it = summary.begin(); it != summary.end(); it++) { allCalcs.insert(*it); }
		for(it = sharedrarefaction.begin(); it != sharedrarefaction.end(); it++) { allCalcs.insert(*it); }
		for(it = sharedsummary.begin(); it != sharedsummary.end(); it++) { allCalcs.insert(*it); }
		for(it = vennsingle.begin(); it != vennsingle.end(); it++) { allCalcs.insert(*it); }
		for(it = vennshared.begin(); it != vennshared.end(); it++) { allCalcs.insert(*it); }
		for(it = treegroup.begin(); it != treegroup.end(); it++) { allCalcs.insert(*it); }
		for(it = matrix.begin(); it != matrix.end(); it++) { allCalcs.insert(*it); }
		for(it = heat.begin(); it != heat.end(); it++) { allCalcs.insert(*it); }
		for(it = distance.begin(); it != distance.end(); it++) { allCalcs.insert(*it); }
        for(it = estimators.begin(); it != estimators.end(); it++) { allCalcs.insert(*it); }
		
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "ValidCalculator");
		exit(1);
	}
}
/********************************************************************/
ValidCalculators::~ValidCalculators() {}
/********************************************************************/
void ValidCalculators::printCitations(vector<string> Estimators) {
	try {
		
		for (int i = 0; i < Estimators.size(); i++) {
			//is this citation, do nothing
			if ((Estimators[i] == "citation") || (Estimators[i] == "default") || (Estimators[i] == "eachgap") || (Estimators[i] == "nogaps") || (Estimators[i] == "onegap")) {}
			//is this a valid calculator
			else if (allCalcs.count(Estimators[i]) != 0) {
				if (Estimators[i] == "sobs") { Calculator* temp = new Sobs(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "chao") { Calculator* temp = new Chao1(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "nseqs") { Calculator* temp = new NSeqs(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "coverage") { Calculator* temp = new Coverage(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "ace") { Calculator* temp = new Ace(10); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				
				}else if (Estimators[i] == "jack") { Calculator* temp = new Jackknife(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "shannon") { Calculator* temp = new Shannon(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "shannoneven") { Calculator* temp = new ShannonEven(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
                }else if (Estimators[i] == "shannonrange") { Calculator* temp = new RangeShannon(0); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "npshannon") { Calculator* temp = new NPShannon(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "heip") { Calculator* temp = new Heip(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "smithwilson") { Calculator* temp = new SmithWilson(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "simpson") { Calculator* temp = new Simpson(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "simpsoneven") { Calculator* temp = new SimpsonEven(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "invsimpson") { Calculator* temp = new InvSimpson(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "bootstrap") { Calculator* temp = new Bootstrap(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "geometric") { Calculator* temp = new Geom(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "qstat") { Calculator* temp = new QStat(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "logseries") { Calculator* temp = new LogSD(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "bergerparker") { Calculator* temp = new BergerParker(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "bstick") { Calculator* temp = new BStick(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "goodscoverage") { Calculator* temp = new GoodsCoverage(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "efron") { Calculator* temp = new Efron(10); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "boneh") { Calculator* temp = new Boneh(10); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;  
				}else if (Estimators[i] == "solow") { Calculator* temp = new Solow(10); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "shen") { Calculator* temp = new Shen(10, 10); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "sharedchao") { Calculator* temp = new SharedChao1(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sharedsobs") { Calculator* temp = new SharedSobsCS(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sharedace") { Calculator* temp = new SharedAce(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "jabund") { 	Calculator* temp = new JAbund(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sorabund") { Calculator* temp = new SorAbund(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "jclass") { Calculator* temp = new Jclass(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sorclass") { Calculator* temp = new SorClass(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;  
				}else if (Estimators[i] == "jest") { Calculator* temp = new Jest(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sorest") { Calculator* temp = new SorEst(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "thetayc") { Calculator* temp = new ThetaYC(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "thetan") { Calculator* temp = new ThetaN(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "kstest") { Calculator* temp = new KSTest(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "whittaker") { Calculator* temp = new Whittaker(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "sharednseqs") { Calculator* temp = new SharedNSeqs(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;  
				}else if (Estimators[i] == "ochiai") { Calculator* temp = new Ochiai(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				
				}else if (Estimators[i] == "anderberg") { Calculator* temp = new Anderberg(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "kulczynski") { Calculator* temp = new Kulczynski(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "kulczynskicody") { Calculator* temp = new KulczynskiCody(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp; 
				}else if (Estimators[i] == "lennon") { Calculator* temp = new Lennon(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "morisitahorn") { Calculator* temp = new MorHorn(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				
				}else if (Estimators[i] == "braycurtis") { Calculator* temp = new BrayCurtis(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "odum") { Calculator* temp = new Odum(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "canberra") { Calculator* temp = new Canberra(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "structeuclidean") { Calculator* temp = new StructEuclidean(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "structchord") { Calculator* temp = new StructChord(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				
				}else if (Estimators[i] == "hellinger") { Calculator* temp = new Hellinger(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "manhattan") { Calculator* temp = new Manhattan(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "structpearson") { Calculator* temp = new StructPearson(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "soergel") { Calculator* temp = new Soergel(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "spearman") { Calculator* temp = new Spearman(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				
				}else if (Estimators[i] == "structkulczynski") { Calculator* temp = new StructKulczynski(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "speciesprofile") { Calculator* temp = new SpeciesProfile(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "hamming") { Calculator* temp = new Hamming(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "structchi2") { Calculator* temp = new StructChi2(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "gower") { Calculator* temp = new Gower(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				
				}else if (Estimators[i] == "memchi2") { Calculator* temp = new MemChi2(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "memchord") { Calculator* temp = new MemChord(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "memeuclidean") { Calculator* temp = new MemEuclidean(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "mempearson") { Calculator* temp = new MemPearson(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "sharedobserved") { Calculator* temp = new SharedSobs(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else if (Estimators[i] == "kulczynski") { Calculator* temp = new Kulczynski(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
                }else if (Estimators[i] == "jsd") { Calculator* temp = new JSD(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
                }else if (Estimators[i] == "rjsd") { Calculator* temp = new RJSD(); m->mothurOut(temp->getName() + ": "); temp->citation(); delete temp;
				}else { m->mothurOut("[ERROR]: Missing else if for " + Estimators[i] + " in printCitations."); m->mothurOutEndLine(); }
			}else { m->mothurOut(Estimators[i] + " is not a valid calculator, no citation will be given."); m->mothurOutEndLine(); }
		}
			
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "printCitations");
		exit(1);
	}
}
/********************************************************************/

bool ValidCalculators::isValidCalculator(string parameter, string calculator) {
	try {
        Utils util;
		//are you looking for a calculator for a single parameter
		if (parameter == "single") {
			//is it valid
			if ((single.find(calculator)) != (single.end())) { return true; }
            else {
				m->mothurOut(calculator + " is not a valid estimator for the collect.single command and will be disregarded. Valid estimators are " + util.getStringFromSet(single, ",") + ".\n");
                return false;
            }
		//are you looking for a calculator for a shared parameter
        }else if (parameter == "estimator") {
            //is it valid
            if ((estimators.find(calculator)) != (estimators.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the estimator.single command and will be disregarded. Valid estimators are " + util.getStringFromSet(estimators, ",") + ".\n");
                return false;
            }
		}else if (parameter == "shared") {
			//is it valid
            if ((shared.find(calculator)) != (shared.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the collect.shared command and will be disregarded. Valid estimators are " + util.getStringFromSet(shared, ",") + ".\n");
                return false;
            }
		//are you looking for a calculator for a rarefaction parameter
		}else if (parameter == "rarefaction") {
			//is it valid
            if ((rarefaction.find(calculator)) != (rarefaction.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the rarefaction.single command and will be disregarded. Valid estimators are " + util.getStringFromSet(rarefaction, ",") + ".\n");
                return false;
            }
		//are you looking for a calculator for a summary parameter
		}else if (parameter == "summary") {
			//is it valid
            if ((summary.find(calculator)) != (summary.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the summary.single command and will be disregarded. Valid estimators are " + util.getStringFromSet(summary, ",") + ".\n");
                return false;
            }
		//are you looking for a calculator for a sharedsummary parameter
		}else if (parameter == "sharedsummary") {
            //is it valid
            if ((sharedsummary.find(calculator)) != (sharedsummary.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the summary.shared command and will be disregarded. Valid estimators are " + util.getStringFromSet(sharedsummary, ",") + ".\n");
                return false;
            }
		}else if (parameter == "sharedrarefaction") {
            //is it valid
            if ((sharedrarefaction.find(calculator)) != (sharedrarefaction.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the rarefaction.shared command and will be disregarded. Valid estimators are " + util.getStringFromSet(sharedrarefaction, ",") + ".\n");
                return false;
            }
		}else if (parameter == "vennsingle") {
            //is it valid
            if ((vennsingle.find(calculator)) != (vennsingle.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the venn command and will be disregarded. Valid estimators are " + util.getStringFromSet(vennsingle, ",") + ".\n");
                return false;
            }
		}else if (parameter == "vennshared") {
            //is it valid
            if ((vennshared.find(calculator)) != (vennshared.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the venn command in shared mode and will be disregarded. Valid estimators are " + util.getStringFromSet(vennshared, ",") + ".\n");
                return false;
            }
		}else if (parameter == "treegroup") {
            //is it valid
            if ((treegroup.find(calculator)) != (treegroup.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the tree.shared command and will be disregarded. Valid estimators are " + util.getStringFromSet(treegroup, ",") + ".\n");
                return false;
            }
		}else if (parameter == "matrix") {
            //is it valid
            if ((matrix.find(calculator)) != (matrix.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the dist.shared command and will be disregarded. Valid estimators are " + util.getStringFromSet(matrix, ",") + ".\n");
                return false;
            }
		}else if (parameter == "heat") {
            //is it valid
            if ((heat.find(calculator)) != (heat.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the heatmap.sim command and will be disregarded. Valid estimators are " + util.getStringFromSet(heat, ",") + ".\n");
                return false;
            }
		}else if (parameter == "distance") {
            //is it valid
            if ((distance.find(calculator)) != (distance.end())) { return true; }
            else {
                m->mothurOut(calculator + " is not a valid estimator for the distance command and will be disregarded. Valid estimators are " + util.getStringFromSet(distance, ",") + ".\n");
                return false;
            }
		//not a valid parameter
		}else { return false; }
		
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "isValidCalculator");
		exit(1);
	}
}
/********************************************************************/

void ValidCalculators::initialEstimators() {
    try {
        
        estimators.insert("erarefact");
        estimators.insert("metroig");
        estimators.insert("metroln");
        estimators.insert("metrols");
        estimators.insert("metrosichel");
        estimators.insert("igabund");
        estimators.insert("lnabund");
        estimators.insert("igrarefact");
        estimators.insert("lnrarefact");
        estimators.insert("lnshift");
        estimators.insert("lsabund");
        estimators.insert("lsrarefact");
        estimators.insert("default");
    }
    catch(exception& e) {
        m->errorOut(e, "ValidCalculator", "initialEstimators");
        exit(1);
    }
}
/********************************************************************/
void ValidCalculators::initialSingle() {
	try {	
		single.insert("sobs");
		single.insert("chao");
		single.insert("ace");
		single.insert("jack");
		single.insert("shannon");
		single.insert("npshannon");
		single.insert("shannoneven");
        single.insert("shannonrange");
		single.insert("smithwilson");
		single.insert("heip");
		single.insert("simpson");
		single.insert("simpsoneven");
		single.insert("invsimpson");
		single.insert("bergerparker");
		single.insert("bootstrap");
		single.insert("geometric");
		single.insert("logseries");
		single.insert("qstat");
		single.insert("bstick");
		single.insert("goodscoverage");
		single.insert("nseqs");
		single.insert("coverage");
		single.insert("efron");
		single.insert("boneh");
		single.insert("solow");
		single.insert("shen");
		single.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSingle");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialShared() {
	try {	
		shared.insert("sharedsobs");
		shared.insert("sharedchao");
		shared.insert("sharedace");
		shared.insert("jabund");
		shared.insert("sorabund");
		shared.insert("jclass");
		shared.insert("sorclass");
		shared.insert("jest");
		shared.insert("sorest");
		shared.insert("thetayc");
		shared.insert("thetan");
		shared.insert("kstest");
		shared.insert("whittaker");
		shared.insert("sharednseqs");
		shared.insert("ochiai");
		shared.insert("anderberg");
		shared.insert("kulczynski");
		shared.insert("kulczynskicody");
		shared.insert("lennon");
		shared.insert("morisitahorn");
		shared.insert("braycurtis");
		shared.insert("odum");
		shared.insert("canberra");
		shared.insert("structeuclidean");
		shared.insert("structchord");
		shared.insert("hellinger");
		shared.insert("manhattan");
		shared.insert("structpearson");
		shared.insert("soergel");
		shared.insert("spearman");
		shared.insert("structkulczynski");
		shared.insert("structchi2");
		shared.insert("speciesprofile");
		shared.insert("hamming");
		shared.insert("gower");
		shared.insert("memchi2");
		shared.insert("memchord");
		shared.insert("memeuclidean");
		shared.insert("mempearson");
        shared.insert("jsd");
        shared.insert("rjsd");
		shared.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialShared");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialRarefaction() {
	try {	
		rarefaction.insert("sobs");
		rarefaction.insert("chao");
		rarefaction.insert("ace");
		rarefaction.insert("jack");
		rarefaction.insert("shannon");
		rarefaction.insert("smithwilson");
		rarefaction.insert("heip");
		rarefaction.insert("npshannon");
		rarefaction.insert("shannoneven");
        rarefaction.insert("shannonrange");
		rarefaction.insert("simpson");
		rarefaction.insert("invsimpson");
		rarefaction.insert("simpsoneven");
		rarefaction.insert("bootstrap");
		rarefaction.insert("nseqs");
		rarefaction.insert("coverage");
		rarefaction.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialRarefaction");
		exit(1);
	}
}

/********************************************************************/

void ValidCalculators::initialSummary() {
	try {	
		summary.insert("sobs");
		summary.insert("chao");
		summary.insert("ace");
		summary.insert("jack");
		summary.insert("shannon");
		summary.insert("heip");
		summary.insert("shannoneven");
		summary.insert("smithwilson");
		summary.insert("invsimpson");
		summary.insert("npshannon");
        summary.insert("shannonrange");
		summary.insert("simpson");
		summary.insert("simpsoneven");
		summary.insert("bergerparker");
		summary.insert("geometric");
		summary.insert("bootstrap");
		summary.insert("logseries");
		summary.insert("qstat");
		summary.insert("bstick");
		summary.insert("nseqs");
		summary.insert("goodscoverage");
		summary.insert("coverage");
		summary.insert("efron");
		summary.insert("boneh");
		summary.insert("solow");
		summary.insert("shen");
		summary.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSummary");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialSharedSummary() {
	try {	
		sharedsummary.insert("sharedsobs");
		sharedsummary.insert("sharedchao");
		sharedsummary.insert("sharedace");
		sharedsummary.insert("jabund");
		sharedsummary.insert("sorabund");
		sharedsummary.insert("jclass");
		sharedsummary.insert("sorclass");
		sharedsummary.insert("jest");
		sharedsummary.insert("sorest");
		sharedsummary.insert("thetayc");
		sharedsummary.insert("thetan");
		sharedsummary.insert("kstest");
		sharedsummary.insert("whittaker");
		sharedsummary.insert("sharednseqs");
		sharedsummary.insert("ochiai");
		sharedsummary.insert("anderberg");
		sharedsummary.insert("kulczynski");
		sharedsummary.insert("kulczynskicody");
		sharedsummary.insert("lennon");
		sharedsummary.insert("morisitahorn");
		sharedsummary.insert("braycurtis");
		sharedsummary.insert("odum");
		sharedsummary.insert("canberra");
		sharedsummary.insert("structeuclidean");
		sharedsummary.insert("structchord");
		sharedsummary.insert("hellinger");
		sharedsummary.insert("manhattan");
		sharedsummary.insert("structpearson");
		sharedsummary.insert("structkulczynski");
		sharedsummary.insert("structchi2");
		sharedsummary.insert("soergel");
		sharedsummary.insert("spearman");
		sharedsummary.insert("speciesprofile");
		sharedsummary.insert("hamming");
		sharedsummary.insert("gower");
		sharedsummary.insert("memchi2");
		sharedsummary.insert("memchord");
		sharedsummary.insert("memeuclidean");
		sharedsummary.insert("mempearson");
        sharedsummary.insert("jsd");
        sharedsummary.insert("rjsd");
		sharedsummary.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSharedSummary");
		exit(1);
	}
}


/********************************************************************/

void ValidCalculators::initialSharedRarefact() {
	try {	
		sharedrarefaction.insert("sharedobserved");
		sharedrarefaction.insert("sharednseqs");
		sharedrarefaction.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSharedRarefact");
		exit(1);
	}
}


/********************************************************************/
void ValidCalculators::initialVennSingle() {
	try {
		vennsingle.insert("sobs");
		vennsingle.insert("chao");
		vennsingle.insert("ace");
		vennsingle.insert("jack");
		vennsingle.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialVennSingle");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialVennShared() {
	try {
		vennshared.insert("sharedsobs");
		vennshared.insert("sharedchao");
		vennshared.insert("sharedace");
		vennshared.insert("default");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialVennShared");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialTreeGroups() {
	try {	
		treegroup.insert("sharedsobs");
		treegroup.insert("sharedchao");
		treegroup.insert("sharedace");
		treegroup.insert("jabund");
		treegroup.insert("sorabund");
		treegroup.insert("jclass");
		treegroup.insert("sorclass");
		treegroup.insert("jest");
		treegroup.insert("sorest");
		treegroup.insert("thetayc");
		treegroup.insert("thetan");
		treegroup.insert("kstest");
		treegroup.insert("whittaker");
		treegroup.insert("sharednseqs");
		treegroup.insert("ochiai");
		treegroup.insert("anderberg");
		treegroup.insert("kulczynski");
		treegroup.insert("kulczynskicody");
		treegroup.insert("lennon");
		treegroup.insert("morisitahorn");
		treegroup.insert("braycurtis");
		treegroup.insert("odum");
		treegroup.insert("canberra");
		treegroup.insert("structeuclidean");
		treegroup.insert("structchord");
		treegroup.insert("hellinger");
		treegroup.insert("manhattan");
		treegroup.insert("structpearson");
		treegroup.insert("structkulczynski");
		treegroup.insert("structchi2");
		treegroup.insert("soergel");
		treegroup.insert("spearman");
		treegroup.insert("speciesprofile");
		treegroup.insert("hamming");
		treegroup.insert("gower");
		treegroup.insert("memchi2");
		treegroup.insert("memchord");
        treegroup.insert("jsd");
        treegroup.insert("rjsd");
		treegroup.insert("memeuclidean");
		treegroup.insert("mempearson");
		
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialTreeGroups");
		exit(1);
	}
}
/********************************************************************/
void ValidCalculators::initialHeat() {
	try {	
		heat.insert("jabund");
		heat.insert("sorabund");
		heat.insert("jclass");
		heat.insert("sorclass");
		heat.insert("jest");
		heat.insert("sorest");
		heat.insert("thetayc");
		heat.insert("thetan");
		heat.insert("morisitahorn");
		heat.insert("braycurtis");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialHeat");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialMatrix() {
	try {	
		matrix.insert("sharedsobs");
		matrix.insert("sharedchao");
		matrix.insert("sharedace");
		matrix.insert("jabund");
		matrix.insert("sorabund");
		matrix.insert("jclass");
		matrix.insert("sorclass");
		matrix.insert("jest");
		matrix.insert("sorest");
		matrix.insert("thetayc");
		matrix.insert("thetan");
		matrix.insert("kstest");
		matrix.insert("whittaker");
		matrix.insert("sharednseqs");
		matrix.insert("ochiai");
		matrix.insert("anderberg");
		matrix.insert("kulczynski");
		matrix.insert("kulczynskicody");
		matrix.insert("lennon");
		matrix.insert("morisitahorn");
		matrix.insert("braycurtis");
		matrix.insert("odum");
		matrix.insert("canberra");
		matrix.insert("structeuclidean");
		matrix.insert("structchord");
		matrix.insert("hellinger");
		matrix.insert("manhattan");
		matrix.insert("structpearson");
		matrix.insert("structkulczynski");
		matrix.insert("structchi2");
		matrix.insert("soergel");
		matrix.insert("spearman");
		matrix.insert("speciesprofile");
		matrix.insert("hamming");
		matrix.insert("gower");
		matrix.insert("memchi2");
		matrix.insert("memchord");
		matrix.insert("memeuclidean");
		matrix.insert("mempearson");
        matrix.insert("rjsd");
        matrix.insert("jsd");
		
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialMatrix");
		exit(1);
	}
}
/********************************************************************/
void ValidCalculators::initialDistance() {
	try {	
		distance.insert("nogaps");
		distance.insert("eachgap");
		distance.insert("onegap");
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialDistance");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::printCalc(string parameter, ostream& out) {
	try{
        out << printCalc(parameter);
    }
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "printCalc");
		exit(1);
	}
}
/********************************************************************/
string ValidCalculators::printCalc(string parameter) {
	try{
        Utils util;
        
		string output = "The available estimators for calc are ";
		
        if (parameter == "single")                  {  output += util.getStringFromSet(single, ", ");               }
        else if (parameter == "shared")             {  output += util.getStringFromSet(shared, ", ");               }
        else if (parameter == "rarefaction")        {  output += util.getStringFromSet(rarefaction, ", ");          }
        else if (parameter == "summary")            {  output += util.getStringFromSet(summary, ", ");              }
        else if (parameter == "sharedsummary")      {  output += util.getStringFromSet(sharedsummary, ", ");        }
        else if (parameter == "sharedrarefaction")  {  output += util.getStringFromSet(sharedrarefaction, ", ");    }
        else if (parameter == "vennsingle")         {  output += util.getStringFromSet(vennsingle, ", ");           }
        else if (parameter == "vennshared")         {  output += util.getStringFromSet(vennshared, ", ");           }
        else if (parameter == "treegroup")          {  output += util.getStringFromSet(treegroup, ", ");            }
        else if (parameter == "matrix")             {  output += util.getStringFromSet(matrix, ", ");               }
        else if (parameter == "heat")               {  output += util.getStringFromSet(heat, ", ");                 }
        else if (parameter == "distance")           {  output += util.getStringFromSet(distance, ", ");             }
		output += "\n";
		
		return output;
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "printCalc");
		exit(1);
	}
}
/********************************************************************/


