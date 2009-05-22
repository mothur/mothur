#include "globaldata.hpp"
#include "tree.h"
#include "sparsematrix.hpp"

/*******************************************************/

/******************************************************/
GlobalData* GlobalData::getInstance() {
	if( _uniqueInstance == 0 ) {
		_uniqueInstance = new GlobalData();
	}
	return _uniqueInstance;
}
/*******************************************************/

/******************************************************/
//This function parses through the option string of the command to remove its parameters
void GlobalData::parseGlobalData(string commandString, string optionText){
	try {
		commandName = commandString; //save command name to be used by other classes
		
		//set all non filename paramters to default
		reset();
		
		//clears out data from previous read
		if ((commandName == "read.dist") || (commandName == "read.otu") || (commandName == "read.tree")) { 
			clear();
			gGroupmap = NULL;
			gTree.clear();
			Treenames.clear();
			labels.clear(); lines.clear(); Groups.clear();
			allLines = 1;
		}
		
		//saves help request
		if (commandName =="help") {
			helpRequest = optionText;
		}
		
		if (commandName == "libshuff") {
			iters = "10000";
			cutoff = "1.0";
		}
		
		//set default value for cutoff
		if (commandName == "dist.seqs") {	cutoff = "1.0";		}

		string key, value;		
		//reads in parameters and values
		if((optionText != "") && (commandName != "help")){
			while((optionText.find_first_of(',') != -1)) {  //while there are parameters
				splitAtComma(value, optionText);
				splitAtEquals(key, value);
				
				if (key == "phylip" )	{ phylipfile = value; inputFileName = value; fileroot = value; format = "phylip";	}
				if (key == "column" )	{ columnfile = value; inputFileName = value; fileroot = value; format = "column";	}
				if (key == "list" )		{ listfile = value; inputFileName = value; fileroot = value; format = "list";		}
				if (key == "rabund" )	{ rabundfile = value; inputFileName = value; fileroot = value; format = "rabund";	}
				if (key == "sabund" )	{ sabundfile = value; inputFileName = value; fileroot = value; format = "sabund";	} 
				if (key == "fasta" )	{ fastafile = value; inputFileName = value; fileroot = value; format = "fasta";		}
				if (key == "nexus" )	{ nexusfile = value; inputFileName = value; fileroot = value; format = "nexus";		} 
				if (key == "clustal" )	{ clustalfile = value; inputFileName = value; fileroot = value; format = "clustal"; }
				if (key == "tree" )		{ treefile = value; inputFileName = value; fileroot = value; format = "tree";		}
				if (key == "shared" )	{ sharedfile = value; inputFileName = value; fileroot = value; format = "sharedfile";	}
				if (key == "name" )		{ namefile = value;		}
				if (key == "order" )	{ orderfile = value;	}
				if (key == "group" )	{ groupfile = value;	}
				if (key == "cutoff" )		{ cutoff = value;		}
				if (key == "precision" )	{ precision = value;	}
				if (key == "iters" )		{ iters = value;		}
				if (key == "jumble" )		{ jumble = value;		}
				if (key == "freq" )			{ freq = value;			}
				if (key == "method" )		{ method = value;		}
				if (key == "fileroot" )		{ fileroot = value;		}
				if (key == "abund" )        { abund = value;        }
				if (key == "random" )		{ randomtree = value;	}
				if (key == "calc")			{ calc = value;			}
				if (key == "step")			{ step = value;			}
				if (key == "form")			{ form = value;			}
				if (key == "sorted")		{ sorted = value;		}
				if (key == "vertical")		{ vertical = value;		}
				if (key == "trump")		    { trump = value;		}
				if (key == "filter")		{ filter = value;		}
				if (key == "soft")		    { soft = value;		    }
				if (key == "scale")			{ scale = value;		}
				if (key == "ends" )			{ ends = value;			}
				if (key == "processors" )	{ processors = value;	}
				if (key == "size" )         { size = value;         }
				if (key == "candidate")		{ candidatefile = value;	}
				if (key == "search")		{ search = value;		}
				if (key == "ksize")			{ ksize = value;		}
				if (key == "align")		    { align = value;		}
				if (key == "match")			{ match = value;		}
				if (key == "mismatch")		{ mismatch = value;	    }
				if (key == "gapopen")		{ gapopen = value;		}
				if (key == "gapextend" )	{ gapextend = value;	}
				
				if (key == "line") {//stores lines to be used in a set
					lines.clear();
					labels.clear();
					line = value;
					label = "";
					splitAtDash(value, lines);
					allLines = 0;
				}
				if (key == "label") {//stores labels to be used in a set
					labels.clear();
					lines.clear();
					label = value;
					line = "";
					splitAtDash(value, labels);
					allLines = 0;
				}

				if (key == "groups") {//stores groups to be used in a vector
					Groups.clear();
					groups = value;
					splitAtDash(value, Groups);
				}

			}
			
			//saves the last parameter
			value = optionText;
			splitAtEquals(key, value);
			if (key == "phylip" )	{ phylipfile = value; inputFileName = value; fileroot = value; format = "phylip";	}
			if (key == "column" )	{ columnfile = value; inputFileName = value; fileroot = value; format = "column";	}
			if (key == "list" )		{ listfile = value; inputFileName = value; fileroot = value; format = "list";		}
			if (key == "rabund" )	{ rabundfile = value; inputFileName = value; fileroot = value; format = "rabund";	}
			if (key == "sabund" )	{ sabundfile = value; inputFileName = value; fileroot = value; format = "sabund";	}
			if (key == "fasta" )	{ fastafile = value; inputFileName = value; fileroot = value; format = "fasta";		}
			if (key == "nexus" )	{ nexusfile = value; inputFileName = value; fileroot = value; format = "nexus";		}
			if (key == "clustal" )	{ clustalfile = value; inputFileName = value; fileroot = value; format = "clustal"; } 
			if (key == "tree" )		{ treefile = value; inputFileName = value; fileroot = value; format = "tree";		} 
			if (key == "shared" )	{ sharedfile = value; inputFileName = value; fileroot = value; format = "sharedfile";	} 
			if (key == "name" )		{ namefile = value;		}
			if (key == "order" )	{ orderfile = value;	}
			if (key == "group" )	{ groupfile = value;	}
			if (key == "cutoff" )		{ cutoff = value;		}
			if (key == "precision" )	{ precision = value;	}
			if (key == "iters" )		{ iters = value;		}
			if (key == "jumble" )		{ jumble = value;		}
			if (key == "freq" )			{ freq = value;			}
			if (key == "method" )		{ method = value;		}
			if (key == "fileroot" )		{ fileroot = value;		}
			if (key == "abund" )        { abund = value;        }
			if (key == "random" )		{ randomtree = value;	}
			if (key == "calc")			{ calc = value;			}
			if (key == "step")			{ step = value;			}
			if (key == "form")			{ form = value;			}
			if (key == "sorted")		{ sorted = value;		}
			if (key == "vertical")		{ vertical = value;		}
			if (key == "trump")		    { trump = value;		}
			if (key == "filter")		{ filter = value;		}
			if (key == "soft")		    { soft = value;		    }
			if (key == "scale")			{ scale = value;		}
			if (key == "ends" )			{ ends = value;			}
			if (key == "processors" )	{ processors = value;	}
			if (key == "size" )         { size = value;         }
			if (key == "candidate")		{ candidatefile = value;	}
			if (key == "search")		{ search = value;		}
			if (key == "ksize")			{ ksize = value;		}
			if (key == "align")		    { align = value;		}
			if (key == "match")			{ match = value;		}
			if (key == "mismatch")		{ mismatch = value;	    }
			if (key == "gapopen")		{ gapopen = value;		}
			if (key == "gapextend" )	{ gapextend = value;	}

			if (key == "line") {//stores lines to be used in a vector
				lines.clear();
				labels.clear();
				line = value;
				label = "";
				if (line != "all") {  splitAtDash(value, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if (key == "label") {//stores lines to be used in a vector
				labels.clear();
				lines.clear();
				label = value;
				line = "";
				if (label != "all") {  splitAtDash(value, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if (key == "groups") {//stores groups to be used in a vector
					Groups.clear();
					groups = value;
					splitAtDash(value, Groups);
			}
		}
		
		//set format for shared
		if ((listfile != "") && (groupfile != "")) { format = "shared"; }
		if ((phylipfile != "") && (groupfile != "")) { format = "matrix"; }
				
		//input defaults for calculators
		if (commandName == "collect.single") {

			if ((calc == "default") || (calc == "")) { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "rarefaction.single") {
			if ((calc == "default") || (calc == "")) { calc = "sobs"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "collect.shared") {

			if ((calc == "default") || (calc == "")) { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "summary.single") {
			if ((calc == "default") || (calc == "")) { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "summary.shared") {
			if ((calc == "default") || (calc == "")) { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "rarefaction.shared") {
			if ((calc == "default") || (calc == "")) { calc = "sharedobserved"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "dist.seqs") {
			if ((calc == "default") || (calc == "")) {  calc = "onegap";  }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "venn") {
			if ((calc == "default") || (calc == "")) { 
				if (format == "list") { calc = "sobs"; }
				else { calc = "sharedsobs"; }
			}
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if ((commandName == "tree.shared") || (commandName == "bootstrap.shared") || (commandName == "dist.shared")) {
			if (calc != "") { 
				Estimators.clear();
				splitAtDash(calc, Estimators);			
			}else { cout << "You have not specified any calculators." << endl; }
		}


		//if you have done a read.otu with a groupfile but don't want to use it anymore because you want to do single commands
		if ((commandName == "collect.single") || (commandName == "rarefaction.single") || (commandName == "summary.single")) {
			if (listfile != "") { format = "list"; }
			else if (sabundfile != "") { format = "sabund"; }
			else if (rabundfile != "") { format = "rabund"; }
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function parseGlobalData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function parseGlobalData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************/

/******************************************************/
// These functions give you the option parameters of the commands
string GlobalData::getPhylipFile()		{	return phylipfile;	}
string GlobalData::getColumnFile()		{	return columnfile;	}
string GlobalData::getListFile()		{	return listfile;	}
string GlobalData::getRabundFile()		{	return rabundfile;	}
string GlobalData::getSabundFile()		{	return sabundfile;	}
string GlobalData::getNameFile()		{	return namefile;	}
string GlobalData::getGroupFile()		{	return groupfile;	}
string GlobalData::getOrderFile()		{	return orderfile;	}
string GlobalData::getTreeFile()		{	return treefile;	}
string GlobalData::getSharedFile()		{	return sharedfile;	}
string GlobalData::getFastaFile()		{	return fastafile;	}
string GlobalData::getNexusFile()		{	return nexusfile;	}
string GlobalData::getClustalFile()     {   return clustalfile; }
string GlobalData::getCutOff()			{	return cutoff;		}
string GlobalData::getFormat()			{	return format;		}
string GlobalData::getPrecision()		{	return precision;	}
string GlobalData::getMethod()			{	return method;		}
string GlobalData::getFileRoot()		{	return fileroot;	}
string GlobalData::getIters()			{	return iters;		}
string GlobalData::getJumble()			{	return jumble;		}
string GlobalData::getFreq()			{	return freq;		}
string GlobalData::getAbund()           {   return abund;       }
string GlobalData::getRandomTree()		{	return randomtree;	}
string GlobalData::getGroups()			{	return groups;		}
string GlobalData::getStep()			{	return step;		}
string GlobalData::getForm()			{	return form;		}
string GlobalData::getSorted()			{	return sorted;		}
string GlobalData::getTrump()			{   return trump;       }
string GlobalData::getSoft()			{   return soft;		}
string GlobalData::getFilter()			{   return filter;		}
string GlobalData::getScale()			{	return scale;		}
string GlobalData::getEnds()			{   return ends;		}
string GlobalData::getProcessors()		{	return processors;	}
string GlobalData::getSize()            {   return size;        }
string GlobalData::getCandidateFile()	{	return candidatefile;}
string GlobalData::getSearch()			{	return search;		}
string GlobalData::getKSize()			{	return ksize;		}
string GlobalData::getAlign()			{	return align;		}
string GlobalData::getMatch()			{	return match;		}
string GlobalData::getMismatch()		{	return mismatch;	}
string GlobalData::getGapopen()			{	return gapopen;		}
string GlobalData::getGapextend()		{	return gapextend;	}


void GlobalData::setListFile(string file)	{	listfile = file;	inputFileName = file;}
void GlobalData::setRabundFile(string file)	{	rabundfile = file;	inputFileName = file;}
void GlobalData::setSabundFile(string file)	{	sabundfile = file;	inputFileName = file;}
void GlobalData::setPhylipFile(string file)	{	phylipfile = file;    inputFileName = file;}
void GlobalData::setColumnFile(string file)	{	columnfile = file;    inputFileName = file;}
void GlobalData::setGroupFile(string file)		{	groupfile = file;	}
void GlobalData::setSharedFile(string file)		{	sharedfile = file;	inputFileName = file; fileroot = file;}
void GlobalData::setNameFile(string file)		{	namefile = file;		}
void GlobalData::setFormat(string Format)		{	format = Format;		}
void GlobalData::setRandomTree(string Random)	{	randomtree = Random;	}
void GlobalData::setGroups(string g)			{	groups = g;				}
void GlobalData::setCalc(string Calc)			{	calc = Calc;			}
void GlobalData::setEnds(string e)				{   ends = e;				}
void GlobalData::setProcessors(string p)		{	processors = p;			}


/*******************************************************/

/******************************************************/
GlobalData::GlobalData() {
	//option definitions should go here...
	helpRequest = "";
	clear();
	gListVector == NULL;		
	gSparseMatrix == NULL;	
}
/*******************************************************/

/******************************************************/
void GlobalData::clear() {
	//option definitions should go here...
	phylipfile		=	"";
	columnfile		=	"";
	listfile		=	"";
	rabundfile		=	"";
	sabundfile		=	"";
	namefile		=	"";
	groupfile		=	""; 
	orderfile		=	"";
	fastafile		=   "";
	nexusfile		=   "";
	clustalfile		=   "";
	treefile		=	"";
	sharedfile		=	"";
	candidatefile	=	"";
	cutoff			=	"10.00";
	format			=	"";
	precision		=	"100";
	iters			=	"1000"; 
	line			=   "";
	label			=	"";
	groups			=	"";
	jumble			=	"1";	//0 means don't jumble, 1 means jumble.
	randomtree		=	"";  //"" means user will enter some user trees, "outputfile" means they just want the random tree distribution to be outputted to outputfile.
	freq			=	"100";
	method			=	"furthest";
	fileroot		=	"";
	abund           =   "10";
	step			=	"0.01";
	form			=	"integral";
	sorted			=	"T";  //F means don't sort, T means sort.
	vertical        =   "";		
	trump           =   "";		
	filter          =   "";		
	soft            =   "";	
	scale			=	"log10";
	ends			=   "T";  //yes
	processors		=	"1";
	size            =   "1000";
	search			=	"kmer";
	ksize			=	"7";
	align			=	"needleman";
	match			=	"1.0";
	mismatch		=	"-1.0";
	gapopen			=	"-1.0";
	gapextend		=	"-2.0";
}

//*******************************************************/

/******************************************************/
void GlobalData::reset() {
	cutoff			=	"10.00";
	precision		=	"100";
	iters			=	"1000"; 
	groups			=	"";
	jumble			=	"1";	//0 means don't jumble, 1 means jumble.
	sorted			=	"T";  //F means don't sort, T means sort.
	randomtree		=	"";  //"" means user will enter some user trees, "outputfile" means they just want the random tree distribution to be outputted to outputfile.
	freq			=	"100";
	method			=	"furthest";
	calc			=	"";
	abund			=   "10";
	step			=	"0.01";
	form			=	"integral";
	ends			=   "T";
	processors		=	"1";
	size            =   "1000";
	search			=	"kmer";
	ksize			=	"7";
	align			=	"needleman";
	match			=	"1.0";
	mismatch		=	"-1.0";
	gapopen			=	"-1.0";
	gapextend		=	"-2.0";
}
/*******************************************************/

/******************************************************/
GlobalData::~GlobalData() {
	_uniqueInstance = 0;
	if(gListVector != NULL)		{	delete gListVector;		}
	if(gSparseMatrix != NULL)	{	delete gSparseMatrix;	}
	if(gorder != NULL)			{	delete gorder;		}
}
/*******************************************************/

/*******************************************************/
void GlobalData::parseTreeFile() {
	//only takes names from the first tree and assumes that all trees use the same names.
	try {
		string filename = treefile;
		ifstream filehandle;
		openInputFile(filename, filehandle);
		int c, comment;
		comment = 0;
		
		//if you are not a nexus file 
		if ((c = filehandle.peek()) != '#') {  
			while((c = filehandle.peek()) != ';') { 
				while ((c = filehandle.peek()) != ';') {
					// get past comments
					if(c == '[') {
						comment = 1;
					}
					if(c == ']'){
						comment = 0;
					}
					if((c == '(') && (comment != 1)){ break; }
					filehandle.get();
				}

				readTreeString(filehandle); 
			}
		//if you are a nexus file
		}else if ((c = filehandle.peek()) == '#') {
			string holder = "";
					
			// get past comments
			while(holder != "translate" && holder != "Translate"){	
				if(holder == "[" || holder == "[!"){
					comment = 1;
				}
				if(holder == "]"){
					comment = 0;
				}
				filehandle >> holder; 
	
				//if there is no translate then you must read tree string otherwise use translate to get names
				if(holder == "tree" && comment != 1){	
					//pass over the "tree rep.6878900 = "
					while (((c = filehandle.get()) != '(') && ((c = filehandle.peek()) != EOF) ) {;}

					if (c == EOF ) { break; }
					filehandle.putback(c);  //put back first ( of tree.
					readTreeString(filehandle);	
					break;
				}
			}
			
			//use nexus translation rather than parsing tree to save time
			if ((holder == "translate") || (holder == "Translate")) {

				string number, name, h;
				h = ""; // so it enters the loop the first time
				while((h != ";") && (number != ";")) { 
					filehandle >> number;
					filehandle >> name;
	
					//c = , until done with translation then c = ;
					h = name.substr(name.length()-1, name.length()); 
					name.erase(name.end()-1);  //erase the comma
					Treenames.push_back(number);
				}
				if (number == ";") { Treenames.pop_back(); }  //in case ';' from translation is on next line instead of next to last name
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function parseTreeFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function parseTreeFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/*******************************************************/

/*******************************************************/
void GlobalData::readTreeString(ifstream& filehandle)	{
	try {
		int c;
		string name; //k
		
		while((c = filehandle.peek()) != ';') { 
				//if you are a name
			if ((c != '(') && (c != ')') && (c != ',') && (c != ':') && (c != '\n') && (c != '\t') && (c != 32)) { //32 is space
				name = "";
				c = filehandle.get();
	//		k = c;
//cout << k << endl;
				while ((c != '(') && (c != ')') && (c != ',') && (c != ':') && (c != '\n') && (c != 32) && (c != '\t')) {			
					name += c;
					c = filehandle.get();
		//	k = c;
//cout << " in name while " << k << endl;
				}
				
//cout << "name = " << name << endl;
				Treenames.push_back(name);
				filehandle.putback(c);
//k = c;
//cout << " after putback" <<  k << endl;
			} 
			
			if (c  == ':') { //read until you reach the end of the branch length
				while ((c != '(') && (c != ')') && (c != ',') && (c != ';') && (c != '\n') && (c != '\t') && (c != 32)) {
					c = filehandle.get();
				//	k = c;
	//cout << " in branch while " << k << endl;
				}
				filehandle.putback(c);
			}
			c = filehandle.get();
			if (c == ';') { break; }
		//	k = c;
//cout << k << endl;

		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function parseTreeFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function parseTreeFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}	

/*******************************************************/

/*******************************************************/


