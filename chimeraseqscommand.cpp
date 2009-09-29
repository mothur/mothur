/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"
#include "bellerophon.h"
#include "pintail.h"
#include "ccode.h"
#include "chimeracheckrdp.h"
#include "chimeraslayer.h"


//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "correction", "processors", "method", "window", "increment", "template", "conservation", "quantile", "mask", "numwanted", "ksize", "svg", "name", "match","mismatch", "divergence", "minsim", "parents", "printall" };
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; mothurOut("fasta is a required parameter for the chimera.seqs command."); mothurOutEndLine(); abort = true;  }	
			
			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  }	
			
			consfile = validParameter.validFile(parameters, "conservation", true);
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { consfile = "";  }	
			
			quanfile = validParameter.validFile(parameters, "quantile", true);
			if (quanfile == "not open") { abort = true; }
			else if (quanfile == "not found") { quanfile = "";  }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = "";  }

			maskfile = validParameter.validFile(parameters, "mask", false);
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				ifstream in;
				int	ableToOpen = openInputFile(maskfile, in);
				if (ableToOpen == 1) { abort = true; }
				in.close();
			}
			
			method = validParameter.validFile(parameters, "method", false);			if (method == "not found") { method = "pintail"; }
			
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "printall", false);			if (temp == "not found") { temp = "F"; }
			printAll = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
			
			temp = validParameter.validFile(parameters, "svg", false);				if (temp == "not found") { temp = "F"; }
			svg = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "window", false);	
			if ((temp == "not found") && (method == "chimeraslayer")) { temp = "100"; }			
			else if (temp == "not found") { temp = "0"; }
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "match", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, match);
			
			temp = validParameter.validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-4"; }
			convert(temp, mismatch);
			
			temp = validParameter.validFile(parameters, "divergence", false);		if (temp == "not found") { temp = "1.0"; }
			convert(temp, divR);
			
			temp = validParameter.validFile(parameters, "minsim", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minSimilarity);
			
			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, parents); 
			 
			temp = validParameter.validFile(parameters, "increment", false);		
			if ((temp == "not found") && ((method == "chimeracheck") || (method == "chimeraslayer"))) { temp = "10"; }
			else if (temp == "not found") { temp = "25"; }
			convert(temp, increment);
			
			temp = validParameter.validFile(parameters, "numwanted", false);
			if ((temp == "not found") && (method == "chimeraslayer")) { temp = "10"; }		
			else if (temp == "not found") { temp = "20"; }
			convert(temp, numwanted);

			
			
			if (((method != "bellerophon")) && (templatefile == "")) { mothurOut("You must provide a template file with the pintail, ccode or chimeracheck methods."); mothurOutEndLine(); abort = true;  }
			

		}
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "ChimeraSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraSeqsCommand::help(){
	try {
	
		//"fasta", "filter", "correction", "processors", "method", "window", "increment", "template", "conservation", "quantile", "mask", "numwanted", "ksize", "svg", "name"
		//mothurOut("chimera.seqs ASSUMES that your sequences are ALIGNED and if using a template that the template file sequences are the same length as the fasta file sequences.\n\n");
		mothurOut("The chimera.seqs command reads a fastafile and creates list of potentially chimeric sequences.\n");
		mothurOut("The chimera.seqs command parameters are fasta, filter, correction, processors, mask, method, window, increment, template, conservation, quantile, numwanted, ksize, svg, name.\n");
		mothurOut("The fasta parameter is always required and template is required if using pintail, ccode or chimeracheck.\n");
		mothurOut("The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n");
		mothurOut("The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs.\n");
		mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		mothurOut("The method parameter allows you to specify the method for finding chimeric sequences.  The default is pintail. Options include bellerophon, ccode and chimeracheck \n");
		mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences. \n");
		mothurOut("The window parameter allows you to specify the window size for searching for chimeras. \n");
		mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences.\n");
		mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences. \n");
		mothurOut("The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment.\n");
		mothurOut("The quantile parameter allows you to enter a file containing quantiles for a template files sequences.\n");
		mothurOut("The numwanted parameter allows you to specify how many sequences you would each query sequence compared with.\n");
		mothurOut("The ksize parameter allows you to input kmersize. \n");
		mothurOut("The svg parameter allows you to specify whether or not you would like a svg file outputted for each query sequence.\n");
		mothurOut("The name parameter allows you to enter a file containing names of sequences you would like .svg files for.\n");
		mothurOut("NOT ALL PARAMETERS ARE USED BY ALL METHODS. Please look below for method specifics.\n\n");
		mothurOut("Details for each method: \n"); 
		mothurOut("\tpintail: \n"); 
		mothurOut("\t\tparameters: fasta=required, template=required, filter=F, mask=no mask, processors=1, window=300, increment=25, conservation=not required, but will improve speed, quantile=not required, but will greatly improve speed. \n"); 
		mothurOut("\t\tIf you have run chimera.seqs using pintail a .quan and .freq file will be created for your template, if you have not provided them for use in future command executions.\n");
		mothurOut("\tbellerophon: \n"); 
		mothurOut("\t\tparameters: fasta=required, filter=F, processors=1, window=1/4 length of seq, increment=25, correction=T. \n"); 
		mothurOut("\tccode: \n"); 
		mothurOut("\t\tparameters: fasta=required, template=required, filter=F, mask=no mask, processors=1, window=10% of length, numwanted=20\n"); 
		mothurOut("\tchimeracheck: \n"); 
		mothurOut("\t\tparameters: fasta=required, template=required, processors=1, increment=10, ksize=7, svg=F, name=none\n\n"); 
		mothurOut("The chimera.seqs command should be in the following format: \n");
		mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, method=bellerophon, window=200) \n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraSeqsCommand::~ChimeraSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		if (method == "bellerophon")			{		chimera = new Bellerophon(fastafile);						}
		else if (method == "pintail")			{		chimera = new Pintail(fastafile, templatefile);				}
		else if (method == "ccode")				{		chimera = new Ccode(fastafile, templatefile);				}
		else if (method == "chimeracheck")		{		chimera = new ChimeraCheckRDP(fastafile, templatefile);		}
		else if (method == "chimeraslayer")		{		chimera = new ChimeraSlayer(fastafile, templatefile);		}
		else { mothurOut("Not a valid method."); mothurOutEndLine(); return 0;		}
		
		//set user options
		if (maskfile == "default") { mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); mothurOutEndLine();  }
		
		//saves time to avoid generating it
		chimera->setCons(consfile);	
		
		//saves time to avoid generating it
		chimera->setQuantiles(quanfile);				
		
		chimera->setMask(maskfile);
		chimera->setFilter(filter);
		chimera->setCorrection(correction);
		chimera->setProcessors(processors);
		chimera->setWindow(window);
		chimera->setIncrement(increment);
		chimera->setNumWanted(numwanted);
		chimera->setKmerSize(ksize);
		chimera->setSVG(svg);
		chimera->setName(namefile);
		chimera->setMatch(match);
		chimera->setMisMatch(mismatch);
		chimera->setDivR(divR);
		chimera->setParents(parents);
		chimera->setMinSim(minSimilarity);
		chimera->setPrint(printAll);
		
				
		//find chimeras
		chimera->getChimeras();
		
		string outputFileName = getRootName(fastafile) + method + maskfile + ".chimeras";
		ofstream out;
		openOutputFile(outputFileName, out);
		
		//print results
		chimera->print(out);
		
		out.close();
		
		delete chimera;
		
		return 0;
		
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/

