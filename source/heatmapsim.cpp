/*
 *  heatmapsim.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapsim.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"

//**********************************************************************************************************************
HeatMapSim::HeatMapSim(string dir, string i, int f) : outputDir(dir), inputfile(i), fontSize(f) {
		m = MothurOut::getInstance();
}
//**********************************************************************************************************************
vector<string> HeatMapSim::getPic(vector<SharedRAbundVector*> lookup, vector<Calculator*> calcs, vector<string> groups) {
	try {
		EstOutput data;
		vector<double> sims;
		vector<string> outputNames;
				
		//make file for each calculator selected
		for (int k = 0; k < calcs.size(); k++) {
		
			if (m->control_pressed) { return outputNames; }
		
			string filenamesvg = outputDir + m->getRootName(m->getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + calcs[k]->getName() + ".heatmap.sim.svg";
			m->openOutputFile(filenamesvg, outsvg);
			outputNames.push_back(filenamesvg);
			
			//svg image
			outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString((lookup.size() * 150) + 160) + " " + toString((lookup.size() * 150) + 160)  + "\">\n";
			outsvg << "<g>\n";
		
			//white backround
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString((lookup.size() * 150) + 160) + "\" height=\"" + toString((lookup.size() * 150) + 160)  + "\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString((lookup.size() * 75) - 40) + "\" y=\"25\">Heatmap at distance " + lookup[0]->getLabel() + "</text>\n";
		
			//column labels
			for (int h = 0; h < lookup.size(); h++) {
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(((150 * (h+1)) ) - ((int)groups[h].length() / 2)) + "\" y=\"50\">" + groups[h] + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" y=\"" + toString(((150 * (h+1)) ) - ((int)groups[h].length() / 2)) + "\" x=\"50\">" + groups[h] + "</text>\n";
			}
			
			sims.clear();
			double biggest = 0;
			double smallest = 10000000;
			//float scaler;

			//get sim for each comparison and save them so you can find the relative similairity
			for(int i = 0; i < (lookup.size()-1); i++){
				for(int j = (i+1); j < lookup.size(); j++){
						
						if (m->control_pressed) { outsvg.close(); return outputNames; }
						
						vector<SharedRAbundVector*> subset;
						subset.push_back(lookup[i]);  subset.push_back(lookup[j]); 
					
						//get similairity between groups
						data = calcs[k]->getValues(subset);
						sims.push_back(1.0 - data[0]);
                
						//save biggest similairity to set relative sim
						if ((1.0 - data[0]) > biggest) { biggest = (1.0 - data[0]); }
                    
                        //save smalllest similairity to set relative sim
                        if ((1.0 - data[0]) < smallest) { smallest = (1.0 - data[0]); }

				}
			}
            
			//map biggest similairity found to red
			float scalerBig = 255.0 / biggest;
    
			int count = 0;
			//output similairites to file
			for(int i = 0; i < (lookup.size()-1); i++){
				for(int j = (i+1); j < lookup.size(); j++){
				
						//find relative color
						int color = scalerBig * sims[count];
						//draw box
						outsvg << "<rect fill=\"rgb(" + toString(color) + ",0,0)\" stroke=\"rgb(" + toString(color) + ",0,0)\" x=\"" + toString((i*150)+80) + "\" y=\"" + toString((j*150)+75) + "\" width=\"150\" height=\"150\"/>\n";
						count++;
				}
			}
			
			int y = ((lookup.size() * 150) + 120);
			printLegend(y, biggest, smallest);
		
			outsvg << "</g>\n</svg>\n";
			outsvg.close();

		}
		
		return outputNames;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSim", "getPic");
		exit(1);
	}
}
//**********************************************************************************************************************
string HeatMapSim::getPic(vector< vector<double> > dists, vector<string> groups) {
	try {
		
		vector<double> sims;
		
		string filenamesvg = outputDir + m->getRootName(m->getSimpleName(inputfile)) + "heatmap.sim.svg";
		m->openOutputFile(filenamesvg, outsvg);
			
		//svg image
		outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString((dists.size() * 150) + 160) + " " + toString((dists.size() * 150) + 160)  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString((dists.size() * 150) + 160) + "\" height=\"" + toString((dists.size() * 150) + 160)  + "\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString((dists.size() * 75) - 40) + "\" y=\"25\">Heatmap for " + inputfile + "</text>\n";
		
		//column labels
		for (int h = 0; h < groups.size(); h++) {
			outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(((150 * (h+1)) ) - ((int)groups[h].length() / 2)) + "\" y=\"50\">" + groups[h] + "</text>\n"; 
			outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" y=\"" + toString(((150 * (h+1)) ) - ((int)groups[h].length() / 2)) + "\" x=\"50\">" + groups[h] + "</text>\n";
		}
			
		double biggest = -1;
        double smallest = 10000000;
		float scaler;

		//get sim for each comparison and save them so you can find the relative similairity
		for(int i = 0; i < (dists.size()-1); i++){
			for(int j = (i+1); j < dists.size(); j++){
			
				if (m->control_pressed) { outsvg.close(); return filenamesvg; }
				
				float sim = 1.0 - dists[i][j];
				sims.push_back(sim);
					
				//save biggest similairity to set relative sim
				if (sim > biggest) { biggest = sim; }
                
                //save smalllest similairity to set relative sim
                if (sim < smallest) { smallest = sim; }
			}
		}
			
		//map biggest similairity found to red
		scaler = 255.0 / biggest;
			
		int count = 0;
		//output similairites to file
		for(int i = 0; i < (dists.size()-1); i++){
			for(int j = (i+1); j < dists.size(); j++){
				
					//find relative color
					int color = scaler * sims[count];					
					//draw box
					outsvg << "<rect fill=\"rgb(" + toString(color) + ",0,0)\" stroke=\"rgb(" + toString(color) + ",0,0)\" x=\"" + toString((i*150)+80) + "\" y=\"" + toString((j*150)+75) + "\" width=\"150\" height=\"150\"/>\n";
					count++;
			}
		}
			
		int y = ((dists.size() * 150) + 120);
		printLegend(y, biggest, smallest);
		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();
		
		return filenamesvg;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSim", "getPic");
		exit(1);
	}
}

//**********************************************************************************************************************

void HeatMapSim::printLegend(int y, float maxSim, float minSim) {
	try {
		//output legend and color labels
		//go through map and give each score a color value
		string color;
		int x = 10;
		
		//prints legend
		for (int i = 1; i < 255; i++) {
			color = toHex(int((float)(i)));
			outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"3\" height=\"30\"/>\n";
			x += 3;
		}
		
        float scaler = (maxSim-minSim) / 5.0; 
		
		//prints legend labels
		x = 0;
		for (int i = 0; i<=5; i++) {
			float label = scaler*i;
			label = int(label * 1000 + 0.5);
			label /= 1000.0;
			string text = toString(label, 3);
			
			outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(x) + "\" y=\"" + toString(y-3) + "\">" + text + "</text>\n";
			x += 153;
		}
	}
	
	catch(exception& e) {
		m->errorOut(e, "HeatMapSim", "printLegend");
		exit(1);
	}
}

//**********************************************************************************************************************




