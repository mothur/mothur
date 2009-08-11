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

//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "correction", "processors", "method", "window", "increment", "template", "conservation", "quantile", "mask" };
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
				
			maskfile = validParameter.validFile(parameters, "mask", false);
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				ifstream in;
				int	ableToOpen = openInputFile(maskfile, in);
				if (ableToOpen == 1) { abort = true; }
				in.close();
			}
		

			

			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, window);
					
			temp = validParameter.validFile(parameters, "increment", false);			if (temp == "not found") { temp = "25"; }
			convert(temp, increment);
				
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found") { method = "pintail"; }
			
			if ((method == "pintail") && (templatefile == "")) { mothurOut("You must provide a template file with the pintail method."); mothurOutEndLine(); abort = true;  }
			

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
		mothurOut("chimera.seqs ASSUMES that your sequences are ALIGNED.\n");
		mothurOut("The chimera.seqs command reads a fastafile and creates a sorted priority score list of potentially chimeric sequences (ideally, the sequences should already be aligned).\n");
		mothurOut("The chimera.seqs command parameters are fasta, filter, correction, processors, mask, method, window, increment, template, conservation and quantile.\n");
		mothurOut("The filter parameter allows you to specify if you would like to apply a 50% soft filter this only applies when the method is bellerphon.  The default is false. \n");
		mothurOut("The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs.   The default is true. This only applies when the method is bellerphon.\n");
		mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		mothurOut("The method parameter allows you to specify the method for finding chimeric sequences.  The default is pintail. Options include..... \n");
		mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the pintail and mallard method.  The default is no mask.  If you enter mask=default, then the mask is 236627 EU009184.1 Shigella dysenteriae str. FBD013. \n");
		mothurOut("The window parameter allows you to specify the window size for searching for chimeras.  The default is 300 is method is pintail unless the sequence length is less than 300, and 1/4 sequence length for bellerphon.\n");
		mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences.  The default is 25.\n");
		mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences for use by the pintail algorythm. It is a required parameter if using pintail.\n");
		mothurOut("The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment for use by the pintail algorythm. It is not required, but will speed up the pintail method.\n");
		mothurOut("The quantile parameter allows you to enter a file containing quantiles for a template files sequences for use by the pintail algorythm. It is not required, but will speed up the pintail method.\n");
		mothurOut("If you have run chimera.seqs using pintail a .quan and .freq file will be created if you have not provided them for use in future command executions.");
		mothurOut("The chimera.seqs command should be in the following format: \n");
		mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, processors=2, method=yourMethod) \n");
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
		
		if (method == "bellerophon")	{		chimera = new Bellerophon(fastafile);			}
		else if (method == "pintail")	{		chimera = new Pintail(fastafile, templatefile);	
			//saves time to avoid generating it
			if (consfile != "")			{		chimera->setCons(consfile);						}
			else						{		chimera->setCons("");							}
			
			//saves time to avoid generating it
			if (quanfile != "")			{		chimera->setQuantiles(quanfile);				}
			else						{		chimera->setQuantiles("");						}
			
			if (maskfile == "default") { mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); mothurOutEndLine();  }
			chimera->setMask(maskfile);
						
		}else { mothurOut("Not a valid method."); mothurOutEndLine(); return 0;		}
		
		//set user options
		chimera->setFilter(filter);
		chimera->setCorrection(correction);
		chimera->setProcessors(processors);
		chimera->setWindow(window);
		chimera->setIncrement(increment);
				
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

