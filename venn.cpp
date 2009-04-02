/*
 *  venn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venn.h"

//**********************************************************************************************************************
Venn::Venn(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function Venn. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function Venn. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void Venn::getPic(OrderVector* order) {
	try {
		
		rabund = order->getRAbundVector();
		
		string filenamesvg = globaldata->inputFileName + ".venn." + order->getLabel() + ".svg";
		
		openOutputFile(filenamesvg, outsvg);
		
			
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
		outsvg << "<g>\n";
				
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + order->getLabel() + "</text>\n";
		outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"350\" cy=\"200\" r=\"100\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(rabund.getNumBins()).length() / 2)) + "\" y=\"195\">" + toString(rabund.getNumBins()) + "</text>\n";  
		outsvg << "</g>\n</svg>\n";
		
		outsvg.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void Venn::getPic(SharedOrderVector* sharedorder) {
	try {
		
		//fills vector of sharedsabunds - lookup
		getSharedVectors(sharedorder);
				
		string filenamesvg = globaldata->inputFileName + ".venn." + sharedorder->getLabel() + "." + groupComb + ".svg";
		openOutputFile(filenamesvg, outsvg);
		
		//image window
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
		outsvg << "<g>\n";
				
		if (lookup.size() == 1) {
			int numA = 0;
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				//are they only in one
				if (lookup[0]->getAbundance(i) != 0)  { numA++; }
			}

			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"350\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(numA).length() / 2)) + "\" y=\"195\">" + toString(numA) + "</text>\n";  
			
		}else if (lookup.size() == 2) {
			//calc the shared otu
			int shared = 0;
			int numA = 0;
			int numB = 0;
			
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				//are they only in one
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0)) { numA++; }
				if ((lookup[1]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0)) { numB++; }
				//are they shared
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0)) { shared++; }
			}
			
			//draw circles
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"250\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"435\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA).length() / 2)) + "\" y=\"195\">" + toString(numA) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)toString(numB).length() / 2)) + "\" y=\"195\">" + toString(numB) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(shared).length() / 2)) + "\" y=\"195\">" + toString(shared) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"580\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString((shared / (float)(numA + numB + shared))) + "</text>\n";
		
		}else if (lookup.size() == 3) {
			//calc the shared otu
			int sharedABC = 0;
			int numA = 0; int numB = 0; int numC = 0;
			int sharedAB = 0; int sharedAC = 0; int sharedBC = 0;
			
			//float scalerB;
			
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				//are they only in one
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[2]->getAbundance(i) == 0)) { numA++; }
				if ((lookup[1]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0) && (lookup[2]->getAbundance(i) == 0)) { numB++; }
				if ((lookup[2]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[0]->getAbundance(i) == 0)) { numC++; }
				//are they shared by 2
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) == 0)) { sharedAB++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0)) { sharedAC++; }
				if ((lookup[2]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0)) { sharedBC++; }
				
				//are they shared by all
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0)) { sharedABC++; }
			}
						
			//draw circles
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"230\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"455\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"343\" cy=\"400\" r=\"150\"/>"; 
			
			//place labels within overlaps
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA).length() / 2)) + "\" y=\"170\">" + toString(numA) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"170\">" + toString(sharedAB) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)toString(numB).length() / 2)) + "\"  y=\"170\">" + toString(numB) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(268 - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"305\">" + toString(sharedAC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)toString(numC).length() / 2)) + "\"   y=\"430\">" + toString(numC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(408 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"305\">" + toString(sharedBC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"280\">" + toString(sharedABC) + "</text>\n"; 
			
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"580\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString(((sharedAB + sharedABC) / (float)(numA + numB + sharedAB + sharedABC))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"610\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[2] + " is " + toString(((sharedAC + sharedABC) / (float)(numA + numC + sharedAC + sharedABC))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"640\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(((sharedBC + sharedABC) / (float)(numB + numC + sharedBC + sharedABC))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"670\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString((sharedABC / (float)(numA + numB + numC + sharedAB + sharedAC + sharedBC + sharedABC))) + "</text>\n";
		
		}else if (lookup.size() == 4) {
			//calc the shared otu
			int sharedABCD = 0;
			int numA = 0; int numB = 0; int numC = 0; int numD = 0;
			int sharedAB = 0; int sharedAC = 0; int sharedBC = 0; int sharedAD = 0; int sharedBD = 0; int sharedCD = 0;
			int sharedABC = 0; int sharedACD = 0; int sharedBCD = 0; int sharedABD = 0;
			
			//A = red, B = green, C = blue, D = yellow
			
			//float scalerB;
			
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				//are they only in one
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[2]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { numA++; }
				if ((lookup[1]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0) && (lookup[2]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { numB++; }
				if ((lookup[2]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { numC++; }
				if ((lookup[3]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[2]->getAbundance(i) == 0)) { numD++; }
				//are they shared by 2
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { sharedAB++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { sharedAC++; }
				if ((lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0) && (lookup[3]->getAbundance(i) == 0)) { sharedBC++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) == 0) && (lookup[1]->getAbundance(i) == 0)) { sharedAD++; }
				if ((lookup[3]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) == 0) && (lookup[0]->getAbundance(i) == 0)) { sharedBD++; }
				if ((lookup[2]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0) && (lookup[0]->getAbundance(i) == 0)) { sharedCD++; }
				//are they shared by 3
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) == 0)) { sharedABC++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) == 0)) { sharedACD++; }
				if ((lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0) && (lookup[0]->getAbundance(i) == 0)) { sharedBCD++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) == 0)) { sharedABD++; }
				//are they shared by all
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0) && (lookup[3]->getAbundance(i) != 0)) { sharedABCD++; }
			}
						
			//draw circles
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n";
			outsvg << "<ellipse fill=\"red\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n "; 
			outsvg << "<ellipse fill=\"green\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n ";
			outsvg << "<ellipse fill=\"blue\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-40 440 315) \" cx=\"440\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
			outsvg << "<ellipse fill=\"yellow\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+40 270 315) \" cx=\"270\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
			
			//A = red, B = green, C = blue, D = yellow
			
			//place labels within overlaps
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(460 - ((int)toString(numA).length() / 2)) + "\" y=\"110\">" + toString(numA) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"160\">" + toString(sharedAB) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(250 - ((int)toString(numB).length() / 2)) + "\"  y=\"110\">" + toString(numB) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"190\">" + toString(sharedAC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(550 - ((int)toString(numC).length() / 2)) + "\"   y=\"230\">" + toString(numC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(215 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"190\">" + toString(sharedBC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(150 - ((int)toString(numC).length() / 2)) + "\"   y=\"230\">" + toString(numD) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(240 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"325\">" + toString(sharedAD) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(470 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"325\">" + toString(sharedBD) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedCD).length() / 2)) + "\" y=\"430\">" + toString(sharedCD) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(275 - ((int)toString(sharedABD).length() / 2)) + "\" y=\"240\">" + toString(sharedABD) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(400 - ((int)toString(sharedBCD).length() / 2)) + "\" y=\"360\">" + toString(sharedBCD) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(305 - ((int)toString(sharedACD).length() / 2)) + "\" y=\"360\">" + toString(sharedACD) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(440 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"240\">" + toString(sharedABC) + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedABCD).length() / 2)) + "\"  y=\"320\">" + toString(sharedABCD) + "</text>\n"; 
			
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"490\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString(((sharedAB + sharedABD + sharedABC + sharedABCD) / (float)(numA + numB + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"510\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[2] + " is " + toString(((sharedAC + sharedACD + sharedABC + sharedABCD) / (float)(numA + numC + sharedAB + sharedAC + sharedAD + sharedBC + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"530\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[3] + " is " + toString(((sharedAD + sharedACD + sharedABD + sharedABCD) / (float)(numA + numD + sharedAB + sharedAC + sharedAD + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"550\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(((sharedBC + sharedABC + sharedBCD + sharedABCD) / (float)(numB + numC + sharedAB + sharedAC + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"570\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[3] + " is " + toString(((sharedBD + sharedABD + sharedBCD + sharedABCD) / (float)(numB + numD + sharedAB + sharedAD + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"590\">Percentage of species that are shared in groups " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedCD + sharedBCD + sharedACD + sharedABCD) / (float)(numC + numD + sharedAC + sharedAD + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"610\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(((sharedABC + sharedABCD) / (float)(numA + numB + numC + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"630\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[3] + " is " + toString(((sharedABD + sharedABCD) / (float)(numA + numB + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"650\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedACD + sharedABCD) / (float)(numA + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"670\">Percentage of species that are shared in groups " + globaldata->Groups[1] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedBCD + sharedABCD) / (float)(numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"690\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString((sharedABCD / (float)(numA + numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
		}

		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();

		
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void Venn::getSharedVectors(SharedOrderVector* order){
	try {
	
		//delete lookup
		for (int j = 0; j < lookup.size(); j++) {
			delete lookup[j];
		}

		lookup.clear();
		
		groupComb = "";
		
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < globaldata->Groups.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(globaldata->Groups[i]);
			groupComb += globaldata->Groups[i];
			lookup.push_back(temp);
		}
		
		int numSeqs = order->size();
		//sample all the members
		for(int i=0;i<numSeqs;i++){
			//get first sample
			individual chosen = order->get(i);
			int abundance; 
					
			//set info for sharedvector in chosens group
			for (int j = 0; j < lookup.size(); j++) { 
				if (chosen.group == lookup[j]->getGroup()) {
					 abundance = lookup[j]->getAbundance(chosen.bin);
					 lookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
					 break;
				}
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************



