/*
 *  readtree.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readtree.h"


/* Special characters to trees:
 
 , ) ( ; [ ] :
 
 
 */
/***********************************************************************/
ReadTree::ReadTree() {
	try {
		m = MothurOut::getInstance();
	}
	catch(exception& e) {
		m->errorOut(e, "ReadTree", "ReadTree");
		exit(1);
	}
}
/***********************************************************************/
int ReadTree::AssembleTrees() {
	 try {
		 //assemble users trees
		 for (int i = 0; i < Trees.size(); i++) {
			 if (m->getControl_pressed()) { return 0;  }
			 Trees[i]->assembleTree();
		 }
		 return 0;
	 }
	catch(exception& e) {
		m->errorOut(e, "ReadTree", "AssembleTrees");
		exit(1);
	}
}
/***********************************************************************/
int ReadTree::readSpecialChar(istream& f, char c, string name) {
    try {
        Utils util;
		util.gobble(f);
		char d = f.get();
	
		if(d == EOF){ m->mothurOut("Error: Input file ends prematurely, expecting a " + name + "\n"); exit(1); }
        
		if(d != c){ m->mothurOut("Error: Expected " + name + " in input file.  Found " + toString(d) + ".\n"); exit(1); }
        
		if(d == ')' && f.peek() == '\n'){ util.gobble(f); }
		return d;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadTree", "readSpecialChar");
		exit(1);
	}
}
/**************************************************************************************************/

int ReadTree::readNodeChar(istream& f) {
	try {
        Utils util;
		util.gobble(f);
		char d = f.get();

		if(d == EOF){ m->mothurOut("Error: Input file ends prematurely, expecting a left parenthesis\n"); exit(1); }
        
		return d;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadTree", "readNodeChar");
		exit(1);
	}
}

/**************************************************************************************************/

float ReadTree::readBranchLength(istream& f) {
    try {
		float b;
	
		if(!(f >> b)){
			m->mothurOut("Error: Missing branch length in input tree.\n");
			exit(1);
		}
		util.gobble(f);
		return b;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadTree", "readBranchLength");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/

//Child Classes Below

/***********************************************************************/
/***********************************************************************/
//This class reads a file in Newick form and stores it in a tree.

int ReadNewickTree::read(CountTable* ct) {
	try {
		holder = "";
		int c, error;
		int comment = 0;
		
		//if you are not a nexus file 
		if ((c = filehandle.peek()) != '#') {  
			while((c = filehandle.peek()) != EOF) {
                if (m->getControl_pressed()) {  filehandle.close(); return 0; }
				while ((c = filehandle.peek()) != EOF) {
                    if (m->getControl_pressed()) {  filehandle.close(); return 0; }
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

				//make new tree
				T = new Tree(ct, Treenames);

				numNodes = T->getNumNodes();
				numLeaves = T->getNumLeaves();
				
				error = readTreeString(ct); 
				
				//save trees for later commands
				Trees.push_back(T); 
				util.gobble(filehandle);
			}
		//if you are a nexus file
		}else if ((c = filehandle.peek()) == '#') {
			//get right number of seqs from nexus file.
			Tree* temp = new Tree(ct, Treenames);  delete temp;
			
			nexusTranslation(ct);  //reads file through the translation and updates treemap
			while((c = filehandle.peek()) != EOF) {
                if (m->getControl_pressed()) {  filehandle.close(); return 0; }
				// get past comments
				while ((c = filehandle.peek()) != EOF) {
                    if (m->getControl_pressed()) {  filehandle.close(); return 0; }
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
				T = new Tree(ct, Treenames);
				numNodes = T->getNumNodes();
				numLeaves = T->getNumLeaves();
				
				//read tree info
				error = readTreeString(ct); 
				 
				//save trees for later commands
				Trees.push_back(T); 
			}
		}
		
		if (error != 0) { readOk = error; } 
		
		filehandle.close();

		return readOk;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadNewickTree", "read");
		exit(1);
	}
}
/**************************************************************************************************/
//This function read the file through the translation of the sequences names and updates treemap.
string ReadNewickTree::nexusTranslation(CountTable* ct) {
	try {
		
		holder = "";
		int numSeqs = Treenames.size(); //must save this some when we clear old names we can still know how many sequences there were
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
			if(holder == "tree" && comment != 1){return holder;}
		}
    
		string number, name;
		for(int i=0;i<numSeqs;i++){
			
			filehandle >> number;
			filehandle >> name;
			name.erase(name.end()-1);  //erase the comma
			ct->renameSeq(name, toString(number));
		}
		
		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadNewickTree", "nexusTranslation");
		exit(1);
	}
}

/**************************************************************************************************/
int ReadNewickTree::readTreeString(CountTable* ct) {
	try {
		
		int n = 0;
		int lc, rc; 
		
		int rooted = 0;
	
		int ch = filehandle.peek();	
		
		if(ch == '('){
			n = numLeaves;  //number of leaves / sequences, we want node 1 to start where the leaves left off

			lc = readNewickInt(filehandle, n, T, ct);
			if (lc == -1) { m->mothurOut("error with lc\n");  m->setControl_pressed(true); return -1; } //reports an error in reading
	
			if(filehandle.peek()==','){							
				readSpecialChar(filehandle,',',"comma");
			}
			// ';' means end of tree.												
			else if((ch=filehandle.peek())==';' || ch=='['){		
				rooted = 1;									
			}	
		
			if(rooted != 1){								
				rc = readNewickInt(filehandle, n, T, ct);
				if (rc == -1) { m->mothurOut("error with rc\n");  m->setControl_pressed(true); return -1; } //reports an error in reading
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
				m->mothurOut("Internal error: The only taxon is not taxon 0.\n");
				//exit(1);
				readOk = -1; return -1;
			}
			lc = rc = -1;
		} 
		
		while(((ch=filehandle.get())!=';') && (filehandle.eof() != true)){;}  	
							
		if(rooted != 1){									
			T->tree[n].setChildren(lc,rc);
			T->tree[n].setBranchLength(0);
			T->tree[n].setParent(-1);
			if(lc!=-1){		T->tree[lc].setParent(n);		}
			if(rc!=-1){		T->tree[rc].setParent(n);		}
		}
		
		//T->printTree(); cout << endl;
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ReadNewickTree", "readTreeString");
		exit(1);
	}
}
/**************************************************************************************************/

int ReadNewickTree::readNewickInt(istream& f, int& n, Tree* T, CountTable* ct) {
	try {
		
		if (m->getControl_pressed()) { return -1; } 
		
		int c = readNodeChar(f);

		if(c == '('){
		
			//to account for multifurcating trees generated by fasttree, we are forcing them to be bifurcating
			//read all children
			vector<int> childrenNodes;
			while(f.peek() != ')'){
				int child = readNewickInt(f, n, T, ct);
				if (child == -1) { return -1; } //reports an error in reading
		//cout << "child = " << child << endl;		
				childrenNodes.push_back(child);
				
				//after a child you either have , or ), check for both
				if(f.peek()==')'){  break;  }
				else if (f.peek()==',') {   readSpecialChar(f,',',"comma");  }
				else {;}
			}
	//cout << childrenNodes.size() << endl;		
			if (childrenNodes.size() < 2) {  m->mothurOut("Error in tree, please correct."); m->mothurOutEndLine(); return -1; }
			
			//then force into 2 node structure
			for (int i = 1; i < childrenNodes.size(); i++) {
			
				int lc, rc;
				if (i == 1) { lc = childrenNodes[i-1]; rc = childrenNodes[i]; }
				else { lc = n-1; rc = childrenNodes[i]; }
			//cout << i << '\t' << lc << '\t' << rc << endl;	
				T->tree[n].setChildren(lc,rc);
				T->tree[lc].setParent(n);
				T->tree[rc].setParent(n);
				
				//T->printTree(); cout << endl;
				n++;
			}
			
			//to account for extra ++ in looping
			n--;
			
			if(f.peek()==')'){	
				readSpecialChar(f,')',"right parenthesis");	
				//to pass over labels in trees
				c=filehandle.get();
				while((c!=',') && (c != -1) && (c!= ':') && (c!=';')&& (c!=')')){ c=filehandle.get(); }
				filehandle.putback(c);
			}			
		
			if(f.peek() == ':'){									      
				readSpecialChar(f,':',"colon");	
										
				if(n >= numNodes){ m->mothurOut("Error: Too many nodes in input tree\n");  readOk = -1; return -1; }
				
				T->tree[n].setBranchLength(readBranchLength(f));
			}else{
				T->tree[n].setBranchLength(0.0); 
			}						
						
			return n++;
		
		}else{
			f.putback(c);
			string name = "";
			char d=f.get();
			while(d != ':' && d != ',' && d!=')' && d!='\n'){					
				name += d;
				d=f.get();
			}
//cout << name << endl;
			int blen = 0;
			if(d == ':')	{		blen = 1;	}		
		
			f.putback(d);
		
			//set group info
			vector<string> group = ct->getGroups(name);
            //cout << name << endl;	
			//find index in tree of name
			int n1 = T->getIndex(name);
			
			//adds sequence names that are not in group file to the "xxx" group
			if(group.size() == 0) {
				m->mothurOut("Name: " + name + " is not in your groupfile, and will be disregarded. \n");  //readOk = -1; return n1;
				
                vector<string> currentGroups = ct->getNamesOfGroups();
                Utils util; if (!util.inUsersGroups("xxx", currentGroups)) {  ct->addGroup("xxx");  }
                currentGroups = ct->getNamesOfGroups();
                vector<int> thisCounts; thisCounts.resize(currentGroups.size(), 0);
                for (int h = 0; h < currentGroups.size(); h++) {  
                    if (currentGroups[h] == "xxx") {  thisCounts[h] = 1;  break; }
                }
                ct->push_back(name, thisCounts);
                
				group.push_back("xxx");
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
		m->errorOut(e, "ReadNewickTree", "readNewickInt");
		exit(1);
	}
}
/**************************************************************************************************/
/**************************************************************************************************/

