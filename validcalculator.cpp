/*
 *  validcalculator.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validcalculator.h"

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
		 initialBoot();
		 initialDistance();
		 initialMatrix();
		 initialHeat();
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "ValidCalculator");
		exit(1);
	}
}

/********************************************************************/

ValidCalculators::~ValidCalculators() {}

/********************************************************************/

bool ValidCalculators::isValidCalculator(string parameter, string calculator) {
	try {	
		//are you looking for a calculator for a single parameter
		if (parameter == "single") {
			//is it valid
			if ((single.find(calculator)) != (single.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the collect.single command and will be disregarded. Valid estimators are ");
				for (it = single.begin(); it != single.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		//are you looking for a calculator for a shared parameter
		}else if (parameter == "shared") {
			//is it valid
			if ((shared.find(calculator)) != (shared.end())) {
				return true;
			}else { 
				m->mothurOut(calculator +  " is not a valid estimator for the collect.shared command and will be disregarded.  Valid estimators are ");
				for (it = shared.begin(); it != shared.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		//are you looking for a calculator for a rarefaction parameter
		}else if (parameter == "rarefaction") {
			//is it valid
			if ((rarefaction.find(calculator)) != (rarefaction.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the rarefaction.single command and will be disregarded. Valid estimators are ");
				for (it = rarefaction.begin(); it != rarefaction.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		//are you looking for a calculator for a summary parameter
		}else if (parameter == "summary") {
			//is it valid
			if ((summary.find(calculator)) != (summary.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the summary.shared command and will be disregarded. Valid estimators are ");
				for (it = summary.begin(); it != summary.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		//are you looking for a calculator for a sharedsummary parameter
		}else if (parameter == "sharedsummary") {
			//is it valid
			if ((sharedsummary.find(calculator)) != (sharedsummary.end())) {
				return true;
			}else { 
				m->mothurOut(calculator +  " is not a valid estimator for the summary.shared command and will be disregarded. Valid estimators are ");
				for (it = sharedsummary.begin(); it != sharedsummary.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "sharedrarefaction") {
			//is it valid
			if ((sharedrarefaction.find(calculator)) != (sharedrarefaction.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the rarefaction.shared command and will be disregarded. Valid estimator is ");
				for (it = sharedrarefaction.begin(); it != sharedrarefaction.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "vennsingle") {
			//is it valid
			if ((vennsingle.find(calculator)) != (vennsingle.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the venn command in single mode and will be disregarded. Valid estimators are ");
				for (it = vennsingle.begin(); it != vennsingle.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "vennshared") {
			//is it valid
			if ((vennshared.find(calculator)) != (vennshared.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the venn command in shared mode and will be disregarded. Valid estimators are ");
				for (it = vennshared.begin(); it != vennshared.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "treegroup") {
			//is it valid
			if ((treegroup.find(calculator)) != (treegroup.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the tree.shared command and will be disregarded. Valid estimators are ");
				for (it = treegroup.begin(); it != treegroup.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "matrix") {
			//is it valid
			if ((matrix.find(calculator)) != (matrix.end())) {
				return true;
			}else { 
				m->mothurOut(calculator +  " is not a valid estimator for the matrix.output command and will be disregarded. Valid estimators are ");
				for (it = matrix.begin(); it != matrix.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "heat") {
			//is it valid
			if ((heat.find(calculator)) != (heat.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the heatmap.sim command and will be disregarded. Valid estimators are ");
				for (it = heat.begin(); it != heat.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "boot") {
			//is it valid
			if ((boot.find(calculator)) != (boot.end())) {
				return true;
			}else { 
				m->mothurOut(calculator + " is not a valid estimator for the bootstrap.shared command and will be disregarded. Valid estimators are ");
				for (it = boot.begin(); it != boot.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		}else if (parameter == "distance") {
			//is it valid
			if ((distance.find(calculator)) != (distance.end())) {
				return true;
			}else { 
				m->mothurOut(calculator +  " is not a valid estimator for the distance command and will be disregarded. Valid calculators are ");
				for (it = distance.begin(); it != distance.end(); it++) {
					m->mothurOut(it->first + ", ");
				}
				m->mothurOutEndLine();
				return false; }
		//not a valid parameter
		}else { return false; }
		
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "isValidCalculator");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialSingle() {
	try {	
		single["sobs"]	        = "sobs";
		single["chao"]		    = "chao";
		single["ace"]		    = "ace";
		single["jack"]		    = "jack";
		single["shannon"]	    = "shannon";
		single["npshannon"]    	= "npshannon";
		single["shannoneven"]	= "shannoneven";
		single["smithwilson"]	= "smithwilson";
		single["heip"]			= "heip";
		single["simpson"]	    = "simpson";
		single["simpsoneven"]   = "simpsoneven";
		single["invsimpson"]	= "invsimpson";
		single["bergerparker"]  = "bergerparker";
		single["bootstrap"]    	= "bootstrap";
		single["geometric"]     = "geometric";
		single["logseries"]		= "logseries";
		single["qstat"]         = "qstat";
		single["bstick"]        = "bstick";
		single["goodscoverage"] = "goodscoverage";
		single["nseqs"]			= "nseqs";
		single["coverage"]		= "coverage";
		single["efron"]         = "efron";
		single["boneh"]         = "boneh";
		single["solow"]         = "solow";
		single["shen"]          = "shen";
		single["default"]	    = "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSingle");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialShared() {
	try {	
		shared["sharedsobs"]			= "sharedsobs";
		shared["sharedchao"]			= "sharedchao";
		shared["sharedace"]				= "sharedace";
		shared["jabund"]				= "jabund";
		shared["sorabund"]				= "sorabund";
		shared["jclass"]				= "jclass";
		shared["sorclass"]				= "sorclass";
		shared["jest"]					= "jest";
		shared["sorest"]				= "sorest";
		shared["thetayc"]				= "thetayc";
		shared["thetan"]				= "thetan";
		shared["kstest"]				= "kstest";
		shared["whittaker"]			    = "whittaker";
		shared["sharednseqs"]			= "sharednseqs";
		shared["ochiai"]				= "ochiai";
		shared["anderberg"]				= "anderberg";
		shared["kulczynski"]			= "kulczynski";
		shared["kulczynskicody"]		= "kulczynskicody";
		shared["lennon"]				= "lennon";
		shared["morisitahorn"]			= "morisitahorn";
		shared["braycurtis"]			= "braycurtis";
		shared["default"]	            = "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialShared");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialRarefaction() {
	try {	
		rarefaction["sobs"]			= "sobs";
		rarefaction["chao"]			= "chao";
		rarefaction["ace"]			= "ace";
		rarefaction["jack"]			= "jack";
		rarefaction["shannon"]		= "shannon";
		rarefaction["smithwilson"]	= "smithwilson";
		rarefaction["heip"]			= "heip";
		rarefaction["npshannon"]	= "npshannon";
		rarefaction["shannoneven"]	= "shannoneven";
		rarefaction["simpson"]		= "simpson";
		rarefaction["invsimpson"]	= "invsimpson";
		rarefaction["simpsoneven"]	= "simpsoneven";
		rarefaction["bootstrap"]	= "bootstrap";
		rarefaction["nseqs"]		= "nseqs";
		rarefaction["coverage"]		= "coverage";
		rarefaction["default"]	    = "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialRarefaction");
		exit(1);
	}
}

/********************************************************************/

void ValidCalculators::initialSummary() {
	try {	
		summary["sobs"]			= "sobs";
		summary["chao"]			= "chao";
		summary["ace"]			= "ace";
		summary["jack"]			= "jack";
		summary["shannon"]		= "shannon";
		summary["heip"]			= "heip";
		summary["shannoneven"]	= "shannoneven";
		summary["smithwilson"]	= "smithwilson";
		summary["invsimpson"]	= "invsimpson";
		summary["npshannon"]	= "npshannon";
		summary["simpson"]		= "simpson";
		summary["simpsoneven"]	= "simpsoneven";
		summary["bergerparker"] = "bergerparker";
		summary["geometric"]    = "geometric";
		summary["bootstrap"]	= "bootstrap";
		summary["logseries"]    = "logseries";
		summary["qstat"]        = "qstat";
		summary["bstick"]       = "bstick";
		summary["nseqs"]		= "nseqs";
		summary["goodscoverage"]= "goodscoverage";
		summary["coverage"]		= "coverage";
		summary["efron"]        = "efron";
		summary["boneh"]        = "boneh";
		summary["solow"]        = "solow";
		summary["shen"]         = "shen";
		summary["default"]	    = "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSummary");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialSharedSummary() {
	try {	
		sharedsummary["sharedsobs"]				= "sharedsobs";
		sharedsummary["sharedchao"]				= "sharedchao";
		sharedsummary["sharedace"]				= "sharedace";
		sharedsummary["jabund"]					= "jabund";
		sharedsummary["sorabund"]			= "sorabund";
		sharedsummary["jclass"]					= "jclass";
		sharedsummary["sorclass"]				= "sorclass";
		sharedsummary["jest"]					= "jest";
		sharedsummary["sorest"]					= "sorest";
		sharedsummary["thetayc"]				= "thetayc";
		sharedsummary["thetan"]					= "thetan";
		sharedsummary["kstest"]					= "kstest";
		sharedsummary["whittaker"]				= "whittaker";
		sharedsummary["sharednseqs"]			= "sharednseqs";
		sharedsummary["ochiai"]					= "ochiai";
		sharedsummary["anderberg"]				= "anderberg";
		sharedsummary["kulczynski"]				= "kulczynski";
		sharedsummary["kulczynskicody"]			= "kulczynskicody";
		sharedsummary["lennon"]					= "lennon";
		sharedsummary["morisitahorn"]			= "morisitahorn";
		sharedsummary["braycurtis"]				= "braycurtis";
		sharedsummary["default"]				= "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSharedSummary");
		exit(1);
	}
}


/********************************************************************/

void ValidCalculators::initialSharedRarefact() {
	try {	
		sharedrarefaction["sharedobserved"]	= "sharedobserved";
		sharedrarefaction["sharednseqs"]	= "sharednseqs";
		sharedrarefaction["default"]	    = "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialSharedRarefact");
		exit(1);
	}
}


/********************************************************************/
void ValidCalculators::initialVennSingle() {
	try {
		vennsingle["sobs"]	        = "sobs";
		vennsingle["chao"]		    = "chao";
		vennsingle["ace"]			= "ace";
		vennsingle["jack"]		    = "jack";
		vennsingle["default"]		= "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialVennSingle");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialVennShared() {
	try {
		vennshared["sharedsobs"]	= "sharedsobs";
		vennshared["sharedchao"]	= "sharedchao";
		vennshared["sharedace"]		= "sharedace";
		vennshared["default"]		= "default";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialVennShared");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialTreeGroups() {
	try {	
		treegroup["jabund"]					= "jabund";
		treegroup["sorabund"]				= "sorabund";
		treegroup["jclass"]					= "jclass";
		treegroup["sorclass"]				= "sorclass";
		treegroup["jest"]					= "jest";
		treegroup["sorest"]					= "sorest";
		treegroup["thetayc"]				= "thetayc";
		treegroup["thetan"]					= "thetan";
		treegroup["morisitahorn"]			= "morisitahorn";
		treegroup["braycurtis"]				= "braycurtis";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialTreeGroups");
		exit(1);
	}
}
/********************************************************************/
void ValidCalculators::initialHeat() {
	try {	
		heat["jabund"]				= "jabund";
		heat["sorabund"]			= "sorabund";
		heat["jclass"]				= "jclass";
		heat["sorclass"]			= "sorclass";
		heat["jest"]				= "jest";
		heat["sorest"]				= "sorest";
		heat["thetayc"]				= "thetayc";
		heat["thetan"]				= "thetan";
		heat["morisitahorn"]		= "morisitahorn";
		heat["braycurtis"]			= "braycurtis";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialHeat");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialMatrix() {
	try {	
		matrix["jabund"]				= "jabund";
		matrix["sorabund"]				= "sorabund";
		matrix["jclass"]				= "jclass";
		matrix["sorclass"]				= "sorclass";
		matrix["jest"]					= "jest";
		matrix["sorest"]				= "sorest";
		matrix["thetayc"]				= "thetayc";
		matrix["thetan"]				= "thetan";
		matrix["morisitahorn"]			= "morisitahorn";
		matrix["braycurtis"]			= "braycurtis";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialMatrix");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::initialBoot() {
	try {	
		boot["jabund"]				= "jabund";
		boot["sorabund"]			= "sorabund";
		boot["jclass"]				= "jclass";
		boot["sorclass"]			= "orclass";
		boot["jest"]				= "jest";
		boot["sorest"]				= "sorest";
		boot["thetayc"]				= "thetayc";
		boot["thetan"]				= "thetan";
		boot["morisitahorn"]		= "morisitahorn";
		boot["braycurtis"]			= "braycurtis";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialBoot");
		exit(1);
	}
}
/********************************************************************/
void ValidCalculators::initialDistance() {
	try {	
		distance["nogaps"]		= "nogaps";
		distance["eachgap"]		= "eachgap";
		distance["onegap"]		= "onegap";
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "initialDistance");
		exit(1);
	}
}

/********************************************************************/
void ValidCalculators::printCalc(string parameter, ostream& out) {
	try{
		out << "The available estimators for calc are ";
		//are you looking for a calculator for a single parameter
		if (parameter == "single") {
			for (it = single.begin(); it != single.end(); it++) {
				out << it->first << ", ";
			}
		//are you looking for a calculator for a shared parameter
		}else if (parameter == "shared") {
			for (it = shared.begin(); it != shared.end(); it++) {
				out << it->first << ", ";
			}
		//are you looking for a calculator for a rarefaction parameter
		}else if (parameter == "rarefaction") {
			for (it = rarefaction.begin(); it != rarefaction.end(); it++) {
				out << it->first << ", ";
			}
		//are you looking for a calculator for a summary parameter
		}else if (parameter == "summary") {
			for (it = summary.begin(); it != summary.end(); it++) {
				out << it->first << ", ";
			}
		//are you looking for a calculator for a sharedsummary parameter
		}else if (parameter == "sharedsummary") {
			for (it = sharedsummary.begin(); it != sharedsummary.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "sharedrarefaction") {
			for (it = sharedrarefaction.begin(); it != sharedrarefaction.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "vennsingle") {
			for (it = vennsingle.begin(); it != vennsingle.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "vennshared") {
			for (it = vennshared.begin(); it != vennshared.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "treegroup") {
			for (it = treegroup.begin(); it != treegroup.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "matrix") {
			for (it = matrix.begin(); it != matrix.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "heat") {
			for (it = heat.begin(); it != heat.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "boot") {
			for (it = boot.begin(); it != boot.end(); it++) {
				out << it->first << ", ";
			}
		}else if (parameter == "distance") {
			for (it = distance.begin(); it != distance.end(); it++) {
				out << it->first << ", ";
			}
		}

		out << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "ValidCalculator", "printCalc");
		exit(1);
	}
}
/********************************************************************/

