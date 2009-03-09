/*
 *  readtree.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readtree.h"

/***********************************************************************/
ReadTree::ReadTree() {
	try {
		globaldata = GlobalData::getInstance();
		globaldata->gTree.clear();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************************/
int ReadTree::readSpecialChar(istream& f, char c, string name) {
    try {
	
		gobble(f);
		char d = f.get();
	
		if(d == EOF){
			cerr << "Error: Input file ends prematurely, expecting a " << name << "\n";
			exit(1);
		}
		if(d != c){
			cerr << "Error: Expected " << name << " in input file.  Found " << d << ".\n";
			exit(1);
		}
		if(d == ')' && f.peek() == '\n'){
			gobble(f);
		}	
		return d;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readSpecialChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readSpecialChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/

int ReadTree::readNodeChar(istream& f) {
	try {
//		while(isspace(d=f.get()))		{;}
		gobble(f);
		char d = f.get();

		if(d == EOF){
			cerr << "Error: Input file ends prematurely, expecting a left parenthesis\n";
			exit(1);
		}
		return d;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readNodeChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readNodeChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/**************************************************************************************************/

float ReadTree::readBranchLength(istream& f) {
    try {
		float b;
	
		if(!(f >> b)){
			cerr << "Error: Missing branch length in input tree.\n";
			exit(1);
		}
		gobble(f);
		return b;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readBranchLength. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readBranchLength. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}


/***********************************************************************/
/***********************************************************************/


//Child Classes Below

/***********************************************************************/
/***********************************************************************/
//This class reads a file in Newick form and stores it in a tree.

int ReadNewickTree::read() {
	try {
		int c, error;
		int comment = 0;
		
		//if you are not a nexus file 
		if ((c = filehandle.peek()) != '#') {  
			while((c = filehandle.peek()) != EOF) { 
				//make new tree
				T = new Tree(); 
				numNodes = T->getNumNodes();
				numLeaves = T->getNumLeaves();
				
				error = readTreeString(); 
				
				//save trees for later commands
				globaldata->gTree.push_back(T); 
				gobble(filehandle);
			}
		//if you are a nexus file
		}else if ((c = filehandle.peek()) == '#') {
			nexusTranslation();  //reads file through the translation and updates treemap
			while((c = filehandle.peek()) != EOF) { 
				// get past comments
				while ((c = filehandle.peek()) != EOF) {	
					if(holder == "[" || holder == "[!"){
						comment = 1;
					}
					if(holder == "]"){
						comment = 0;
					}
					if((holder == "tree" || holder == "end;") && comment != 1){ holder = ""; comment = 0; break;}
					filehandle >> holder;
				}
			
				//pass over the "tree rep.6878900 = "
				while (((c = filehandle.get()) != '(') && ((c = filehandle.peek()) != EOF) ) {;}
					
				if (c == EOF ) { break; }
				filehandle.putback(c);  //put back first ( of tree.
				
				//make new tree
				T = new Tree(); 
				numNodes = T->getNumNodes();
				numLeaves = T->getNumLeaves();
				
				//read tree info
				error = readTreeString(); 
				 
				//save trees for later commands
				globaldata->gTree.push_back(T); 
			}
		}
		return readOk;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/
//This function read the file through the translation of the sequences names and updates treemap.
void ReadNewickTree::nexusTranslation() {
	try {
		
		holder = "";
		int numSeqs = globaldata->gTreemap->getNumSeqs(); //must save this some when we clear old names we can still know how many sequences there were
		int comment = 0;
		
		// get past comments
		while(holder != "translate" && holder != "Translate"){	
			if(holder == "[" || holder == "[!"){
				comment = 1;
			}
			if(holder == "]"){
				comment = 0;
			}
			filehandle >> holder; 
			if(holder == "tree" && comment != 1){return;}
		}
		
		//update treemap
		globaldata->gTreemap->namesOfSeqs.clear();
		for(int i=0;i<numSeqs;i++){
			string number, name;
			filehandle >> number;
			filehandle >> name;
			name.erase(name.end()-1);  //erase the comma
			//insert new one with new name
			globaldata->gTreemap->treemap[toString(number)].groupname = globaldata->gTreemap->treemap[name].groupname;
			globaldata->gTreemap->treemap[toString(number)].vectorIndex = globaldata->gTreemap->treemap[name].vectorIndex;
			//erase old one.  so treemap[sarah].groupnumber is now treemap[1].groupnumber. if number is 1 and name is sarah.
			globaldata->gTreemap->treemap.erase(name);
			globaldata->gTreemap->namesOfSeqs.push_back(number);
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function nexus. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function nexus. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/**************************************************************************************************/
int ReadNewickTree::readTreeString() {
	try {
		
		int n = 0;
		int lc, rc; 
		
		int rooted = 0;
	
		int ch = filehandle.peek();	
		
		if(ch == '('){
			n = numLeaves;  //number of leaves / sequences, we want node 1 to start where the leaves left off

			lc = readNewickInt(filehandle, n, T);
			if (lc == -1) { return -1; } //reports an error in reading
		
			if(filehandle.peek()==','){							
				readSpecialChar(filehandle,',',"comma");
			}
			// ';' means end of tree.												
			else if((ch=filehandle.peek())==';' || ch=='['){		
				rooted = 1;									
			}												
			if(rooted != 1){								
				rc = readNewickInt(filehandle, n, T);
				if (rc == -1) { return -1; } //reports an error in reading
				if(filehandle.peek() == ')'){					
					readSpecialChar(filehandle,')',"right parenthesis");
				}											
			}												
		}
		//note: treeclimber had the code below added - not sure why?
		else{
			filehandle.putback(ch);
			char name[MAX_LINE];
			filehandle.get(name, MAX_LINE,'\n');
			SKIPLINE(filehandle, ch);
		
			n = T->getIndex(name);

			if(n!=0){
				cerr << "Internal error: The only taxon is not taxon 0.\n";
				//exit(1);
				readOk = -1; return -1;
			}
			lc = rc = -1;
		} 
		
		while((ch=filehandle.get())!=';'){;}						
		if(rooted != 1){									
			T->tree[n].setChildren(lc,rc);
			T->tree[n].setBranchLength(0);
			T->tree[n].setParent(-1);
			if(lc!=-1){		T->tree[lc].setParent(n);		}
			if(rc!=-1){		T->tree[rc].setParent(n);		}
		}
		return 0;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function readTreeString. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function readTreeString. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/**************************************************************************************************/

int ReadNewickTree::readNewickInt(istream& f, int& n, Tree* T) {
	try {
		int c = readNodeChar(f);
    
		if(c == '('){
			int lc = readNewickInt(f, n, T);
			if (lc == -1) { return -1; } //reports an error in reading
			readSpecialChar(f,',',"comma");

			int rc = readNewickInt(f, n, T);
			if (rc == -1) { return -1; }  //reports an error in reading	
			if(f.peek()==')'){	
				readSpecialChar(f,')',"right parenthesis");	
			}			
		
			if(f.peek() == ':'){									      
				readSpecialChar(f,':',"colon");	
										
				if(n >= numNodes){	cerr << "Error: Too many nodes in input tree\n";  readOk = -1; return -1; }
				
				T->tree[n].setBranchLength(readBranchLength(f));
			}else{T->tree[n].setBranchLength(0.0); }						
		
			T->tree[n].setChildren(lc,rc);
			T->tree[lc].setParent(n);
			T->tree[rc].setParent(n);
		
			return n++;
		}else{
			f.putback(c);
			string name = "";
			char d=f.get();
			while(d != ':' && d != ',' && d!=')' && d!='\n'){					
				name += d;
				d=f.get();
			}
		
			int blen = 0;
			if(d == ':')	{		blen = 1;			}		
		
			f.putback(d);
		
			//set group info
			string group = globaldata->gTreemap->getGroup(name);
			
			//find index in tree of name
			int n1 = T->getIndex(name);
			
			//adds sequence names that are not in group file to the "xxx" group
			if(n1 == -1) {
				cerr << "Name: " << name << " not found in your groupfile. \n"; readOk = -1; return n1;
				
				//globaldata->gTreemap->namesOfSeqs.push_back(name);
				//globaldata->gTreemap->treemap[name].groupname = "xxx";
				//globaldata->gTreemap->treemap[name].vectorIndex = (globaldata->gTreemap->namesOfSeqs.size() - 1);
				
				//map<string, int>::iterator it;
				//it = globaldata->gTreemap->seqsPerGroup.find("xxx");
				//if (it == globaldata->gTreemap->seqsPerGroup.end()) { //its a new group
				//	globaldata->gTreemap->namesOfGroups.push_back("xxx");
				//	globaldata->gTreemap->seqsPerGroup["xxx"] = 1;
				//}else {
				//	globaldata->gTreemap->seqsPerGroup["xxx"]++;
				//}
				
				//find index in tree of name
				//n1 = T->getIndex(name);
				//group = "xxx";
				//numLeaves++;
				//numNodes = 2*numLeaves - 1;
			}
			
			T->tree[n1].setGroup(group);
			T->tree[n1].setChildren(-1,-1);
		
			if(blen == 1){	
				f.get();
				T->tree[n1].setBranchLength(readBranchLength(f));
			}else{
				T->tree[n1].setBranchLength(0.0);
			}
		
			while((c=f.get())!=0 && (c != ':' && c != ',' && c!=')') )		{;}		
			f.putback(c);
		
			return n1;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function readNewickInt. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function readNewickInt. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/
/**************************************************************************************************/

