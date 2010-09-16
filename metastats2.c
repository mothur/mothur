#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "fisher2.h"

void testp(double *permuted_ttests,int *B,double *permuted,double 
	   *Imatrix,int *nc,int *nr,int *g,double *Tinitial,double *ps);
void permute_matrix(double *Imatrix,int *nc,int *nr,double 
		    *permuted,int *g,double *trial_ts,double *Tinitial,double 
		    *counter);
void permute_array(int *array, int n);
void calc_twosample_ts(double *Pmatrix,int *g,int *nc,int *nr,double 
		       *Ts,double *Tinitial,double *counter1);
void meanvar(double *pmatrix,int *g,int *nr,int *nc,double *storage);
void start(double *Imatrix,int *g,int *nr,int *nc,double *testing,
			double storage[][9]);

int main (int argc, char *argv[]){
  
  int col=-1, row=-1, size,g=0,c=0,i=0,j=0,k,counter=0,lines=0, B=10000;
  double placeholder=0,min=0, thresh=0.05;

  char arr[10000], str[51], a;
  char location[41]="jobj.txt", output[41]="out.txt";
  
  for(i=0;i<10000;i++){
  	arr[i]='q';
  }
  
int u,rflag=0,cflag=0, bflag=0;
char *filename;
int numbers;
double numb;
extern char *optarg;
extern int optind, optopt, opterr;

while ((u = getopt(argc, argv, ":r:c:g:b:t:f:o:")) != -1) {
    switch(u) {
    case 'r':
        numbers = atoi(optarg);
        printf("The number of features/rows is %d.\n", numbers);
        row = numbers;
        rflag = 1;
        break;
    case 'c':
        numbers = atoi(optarg);
        printf("The number of samples/columns is %d.\n", numbers);
        col = numbers;
        cflag = 1;
        break;
    case 'g':
        numbers = atoi(optarg);
        printf("Your g-value is %d.\n", numbers);
        g = numbers;
        break;
    case 'b':
        numbers = atoi(optarg);
        printf("The number of permutations is %d\n", numbers);
        B = numbers;
        break;
    case 't':
        numb = atof(optarg);
        printf("Threshold is is %lf\n", numb);
        thresh = numb;
        break;        
    case 'f':
        filename = optarg;
        printf("filename input is %s\n", filename);
        strcpy(location,filename);
        break;
    case 'o':
        filename = optarg;
        printf("filename output %s\n", filename);
        strcpy(output,filename);
        break;
    case ':':
        printf("-%c without filename\n", optopt);
        break;
    case '?':
        printf("unknown arg %c\n", optopt);
        break;
    }
}
  
  FILE *jobj, *out;
  jobj=fopen(location,"r");
  
  if(jobj == NULL){
    printf("Please don't forget to save your matrix in the active");
    printf(" directory as \"%s\".\n",location);
    return 0;
  }
 
  // Gets the first line of samples names and checks for user error.
  fgets(arr,10000,jobj);
  
  for(i=0;i<10000;i++){
    if(isspace(arr[i])){
      counter++; }
  }

  if (cflag == 0) {
      printf("You didn't tell us how many subjects there are!\n");
      printf("But we'll still do the analysis as if there are %d subjects.\n\n",col=counter-1);
  }
  if (cflag == 1) {
  	if (col != counter-1){
  	  printf("We would expect %d subjects, but you said %d.\n",counter-1,col);
    } 
  }
  
  while((a = fgetc(jobj)) != EOF){
    if(a == '\n'){
      lines++; }
  }
  
  if (rflag == 0) {
      printf("You didn't specify the number of features!\n");
      printf("But we'll still do the analysis assuming %d features.\n\n", row=lines-1);   
  }
  if (rflag == 1) {
    if ( lines-1 != row ){
     printf("We would expect %d features, but you said %d.\n",lines-1,row);
    }
  }
 
  if (g>=col || g<=0){
  	printf("Check your g value\n");	
  }
 
  // Initialize the matrices
  size = row*col;
  double matrix[row][col];
  double pmatrix[size],pmatrix2[size],permuted[size];  
  double storage[row][9];
  
  for (i=0;i<row;i++){
  	for (j =0;j<9;j++){
      storage[i][j]=0; 		
	}	
  }
  // Reset file below and create a separate matrix.
  rewind(jobj);
  fgets(arr,10000,jobj);
  
  for(i=0; i<row; i++){
    fscanf(jobj,"%s",str); 	
    for(j=0; j<col;j++){
      fscanf(jobj,"%lf",&placeholder);
      matrix[i][j]=placeholder;
      if(isalnum(placeholder)!=0){ // check for ""
      	printf("Your matrix isn't set up properly!\n");
      	return 0;
		}
      pmatrix[c]=0; // initializing to zero
      permuted[c]=0;
      c++;
    }
  }
  
  fclose(jobj);


  // Produces the sum of each column
  double total[col],total1=0,total2=0;
  double ratio[col];
  
  for(i=0;i<col;i++){
    total[i]=0;
    ratio[i]=0; }
  
  for(i=0; i<col; i++){
    for(j=0;j<row;j++){
      total[i]=total[i]+matrix[j][i];  		
    }
  }
  
  for(i=0;i<g-1;i++){
  	total1=total1+total[i];}
  	
  for(i=g-1;i<col;i++){
    total2=total2+total[i];}
  

  // Creates the ratios by first finding the minimum of totals
  min = total[0];
  if (col==2){
  if (total[0]<total[1]){
  	min = total[1];}	
  }
  if (col >2){
  for(i=1;i<col;i++){
    if (min > total[i]){
      min = total[i];}
  }
  }
  if (min<=0){
    printf("Error, the sum of one of the columns <= 0.");
    return 0;	
  }

  
  // Ratio time...
  for(i=0;i<col;i++){
    ratio[i]=total[i]/min;
  }
  
  // Change matrix into an array as received by R for compatibility.
  
  c=0;
  for(i=0;i<col;i++){
    for(j=0;j<row;j++){
      pmatrix[c]=matrix[j][i];
      c++;
    }
  }
  
  if(row == 1){
  	for (i =0; i<col;i++){
  		pmatrix[i]=pmatrix[i]/ratio[i];
	}
  }
  else {
    counter = 0;
    j=-1;
    for (i=0; i<size; i++) {
      if (counter % row == 0) {
        j++;
      }
      pmatrix[i]=pmatrix[i]/ratio[j];
      counter++; 
    }   
  }
  // pass everything to the rest of the code using pointers. then 
  // write to output file. below pointers for most of the values are 
  // created to send everything by reference.
  
  int ptt_size, *permutes,*nc,*nr,*gvalue;
  
  nc = &col;
  nr = &row;
  gvalue = &g;
  
  permutes = &B;
  ptt_size = B*row;
  
  //changing ptt_size to row
  double permuted_ttests[row], pvalues[row], tinitial[row];
  
  for(i=0;i<row;i++){
    permuted_ttests[i]=0;}
  
  for(i=0;i<row;i++){
    pvalues[i]=0;
    tinitial[i]=0; }
  
  // Find the initial values for the matrix.
  start(pmatrix,gvalue,nr,nc,tinitial,storage);
  
  // Start the calculations.
 
  if ( (col==2) || ((g-1)<8) || ((col-g+1) < 8) ){  

  double fish[row], fish2[row];
  for(i=0;i<row;i++){
  	fish[i]=0;
  	fish2[i]=0;}
 
  for(i=0;i<row;i++){
  	
  	for(j=0;j<g-1;j++){
  		fish[i]=fish[i]+matrix[i][j];
	}

	for(j=g-1;j<col;j++){ 
  		fish2[i]=fish2[i]+matrix[i][j];
	}
	
	double  f11,f12,f21,f22;

	f11=fish[i];
	f12=fish2[i];

	f21=total1-f11;
	f22=total2-f12;
	
	double data[] = {f11, f12, f21, f22};
	
	// CONTINGENGCY TABLE:
	//   f11   f12
	//   f21   f22
	
	int *nr, *nc, *ldtabl, *work;
	int nrow=2, ncol=2, ldtable=2, workspace=100000;
	double *expect, *prc, *emin,*prt,*pre;
	double e=0, prc1=0, emin1=0, prt1=0, pre1=0;

	nr = &nrow;
	nc = &ncol;
	ldtabl=&ldtable;
	work = &workspace;
	
	expect = &e;
	prc = &prc1;
	emin=&emin1;
	prt=&prt1;
	pre=&pre1;
	
	fexact(nr,nc,data, ldtabl,expect,prc,emin,prt,pre,work);
	
	if (*pre>.999999999){
		*pre=1;
	}
	storage[i][8] = *pre;
	pvalues[i]=*pre;
    }
  }
  else{
  	 
  testp(permuted_ttests, permutes, permuted,pmatrix, nc, nr, gvalue,tinitial,pvalues);
  	 
  	       // Checks to make sure the matrix isn't sparse.
  double sparse[row], sparse2[row];
  for(i=0;i<row;i++){
  	sparse[i]=0;
  	sparse2[i]=0;}
  
  c=0;	
  for(i=0;i<row;i++){
  	
  	for(j=0;j<g-1;j++){
  		sparse[i]=sparse[i]+matrix[i][j];
	}
	
	if(sparse[i] < (double)(g-1)){
		c++;
	}
	for(j=g-1;j<col;j++){ // ?<= for col
  		sparse2[i]=sparse2[i]+matrix[i][j];
	}
	
	if( (sparse2[i] <(double)(col-g+1))) {
		c++;
	}
		
	if (c==2){
		c=0;
	
	double  f11,f12,f21,f22;

	f11=sparse[i];
	sparse[i]=0;
	
	f12=sparse2[i];
	sparse2[i]=0;
	
	f21=total1-f11;
	f22=total2-f12;
	
	double data[] = {f11, f12, f21, f22};

	int *nr, *nc, *ldtabl, *work;
	int nrow=2, ncol=2, ldtable=2, workspace=10000000; // I added two zeros for larger data sets
	double *expect, *prc, *emin,*prt,*pre;
	double e=0, prc1=0, emin1=0, prt1=0, pre1=0;

	nr = &nrow;
	nc = &ncol;
	ldtabl=&ldtable;
	work = &workspace;
	
	expect = &e;
	prc = &prc1;
	emin=&emin1;
	prt=&prt1;
	pre=&pre1;
	
	fexact(nr,nc,data, ldtabl,expect,prc,emin,prt,pre,work);
	
	if (*pre>.999999999){
		*pre=1;
	}
	storage[i][8] = *pre;
	pvalues[i]=*pre;
    }
  }  
  // End of else statement
  bflag = 1;
  }
  
  // Calculates the mean of counts (not normalized)
  double temp[row][2];
  
  for(j=0;j<row;j++){
  	for(i=0;i<2;i++){
  		temp[j][i]=0;
	}
  }
  
  for (j=0;j<row;j++){
  	for (i=1; i<=(g-1); i++){
  		temp[j][0]=temp[j][0]+matrix[j][i-1];
  	}
  	temp[j][0]= (double) temp[j][0]/(g-1);
  	for(i=g;i<=col;i++){
  		temp[j][1]=temp[j][1]+matrix[j][i-1];
	}
	temp[j][1]= (double) temp[j][1]/(col-g+1);
  }

  for(i=0;i<row;i++){
  	storage[i][3]=temp[i][0];
  	storage[i][7]=temp[i][1];
  	storage[i][8]=pvalues[i];
  }
  
// BACKUP checks
  
  for (i=0;i<row;i++){
    if(pvalues[i]<thresh){
    	printf("Feature %d is significant, p = %.10lf \n",i+1,pvalues[i]);
	}	
  }
  	  
  // And now we write the files to a text file.
  struct tm *local;
  time_t t;
  t = time(NULL);
  local = localtime(&t);
  
  jobj= fopen(location,"r");
  fgets(arr,10000,jobj);
  
  out = fopen(output,"a+");
  
  fprintf(out,"Local time and date of test: %s\n", asctime(local));
  fprintf(out,"# rows = %d, # col = %d, g = %d\n\n",row,col,g);
  if (bflag == 1){
    fprintf(out,"%d permutations\n\n",B);	
  }
  for(i=0; i<row; i++){
    fscanf(jobj,"%s",str); 	
	fprintf(out,"%s",str);
	for(k=0;k<col;k++){
		fscanf(jobj,"%*lf",&placeholder);
	}
    for(j=0; j<9;j++){
      fprintf(out,"\t%.12lf",storage[i][j]);
    }
    fprintf(out,"\n");
  }  
  
  fprintf(out,"\n \n");
  
  fclose(jobj);
  fclose(out);
  
  return 0;
}

void testp(double *permuted_ttests,int *B,double *permuted,
	   double *Imatrix,int *nc,int *nr,int *g,double *Tinitial,double 
	   *ps) {
  
  double Tvalues[*nr];
  int a, b, n, i, j,k=0;
  
  a = *B;
  b = *nr;
  n = a*b;
  
  double counter[b];
  
  for(j=0;j<b;j++){
    counter[j]=0;
  }    

  for (j=1; j<=*B; j++){
    permute_matrix(Imatrix,nc,nr,permuted,g,Tvalues,Tinitial,counter);
   // for(i=0;i<*nr;i++){
   //   permuted_ttests[k]=fabs(Tvalues[i]);
	//    k++;
    }
  
  
  for(j=0;j<*nr;j++){
    ps[j]=((counter[j]+1)/(double)(a+1));
  }
}	

void permute_matrix(double *Imatrix,int *nc,int *nr,double *permuted,
		    int *g,double *trial_ts,double *Tinitial,double *counter1){
		    	
  int i=0,j=0,n=0,a=0,b=0,f=0,c=0,k=0;
  
  a = *nr; // number of rows
  b = *nc;
  n = a*b;
  
  int y[b];
  
  for (i=1; i<=*nc; i++){
    y[i-1] = i;
  }
  
  permute_array(y, b); 
  
  for (i=0; i<*nc; i++){
    f = y[i]; //column number
    c=1;
    c*=(f-1);
    c*=a;
    if (f == 1){
      c = 0;
    } // starting value position in the Imatrix
    for(j=1; j<=*nr; j++){
      permuted[k] = Imatrix[c];
      c++;
      k++;
    }
  }
  
  calc_twosample_ts(permuted,g,nc,nr,trial_ts,Tinitial,counter1);
}

void permute_array(int *array, int n) {
  static int seeded = 0;
  int i;
  
  if (! seeded) {
    seeded = 1;
    srandom(time(NULL));
  }
  
  for (i = 0; i < n; i++) {
    int selection = random() % (n - i);
    int tmp = array[i + selection];
    array[i + selection] = array[i];
    array[i] = tmp;
  }
}

void calc_twosample_ts(double *Pmatrix,int *g,int *nc,int *nr,
		       double *Ts,double *Tinitial,double *counter) {
  int i,a;
  a = *nr;
  a*=4;
  
  double C1[*nr][3], C2[*nr][3], storage[a],tool[a];
  double nrows,ncols,gvalue, xbardiff=0, denom=0;
  
  nrows = (double) *nr;
  ncols = (double) *nc;
  gvalue= (double) *g;
  
    meanvar(Pmatrix,g,nr,nc,storage);
    for(i=0;i<=a-1;i++){
      tool[i]=storage[i];
    }
    for (i=0; i<*nr;i++){
      C1[i][0]=tool[i];
      C1[i][1]=tool[i+*nr+*nr];
      C1[i][2]=C1[i][1]/(gvalue-1);
      
      C2[i][0]=tool[i+*nr];
      C2[i][1]=tool[i+*nr+*nr+*nr]; // var group 2 
      C2[i][2]=C2[i][1]/(ncols-gvalue+1);
    }
    
    for (i=0; i<*nr; i++){
      xbardiff = C1[i][0]-C2[i][0];
      denom = sqrt(C1[i][2]+C2[i][2]);
      Ts[i]=fabs(xbardiff/denom);
      if (fabs(Ts[i])>(fabs(Tinitial[i])+.0000000000001)){ //13th place
	    counter[i]++;
      }
    }	
}

void meanvar(double *pmatrix,int *g,int *nr,int *nc,double *store){
  double temp[*nr], temp2[*nr],var[*nr],var2[*nr],a,b;
  
  int i,m,k,l,n;
  
  a = (double) *g-1;          
  b = (double) (*nc-a);
  
  for (i = 0; i<*nr; i++){
    temp[i]=0;
    temp2[i]=0;
    var[i]=0;
    var2[i]=0;
  }
  
  k = *nr; // number of rows 
  l = *nc;
  n = k*l;	
  
    m=0;
    m=*g-1;
    k=*nr;
    m*=k; // m = g * nr now
    for (i=0;i<m;i++){
      temp[i%k]=temp[i%k]+pmatrix[i];
    }
    for (i=0;i<n;i++){
      temp2[i%k]=temp2[i%k]+pmatrix[i];
    }
    for (i=0;i<*nr;i++){
      temp2[i]=temp2[i]-temp[i];
    }
    for (i=0;i<=*nr-1;i++){
      store[i]=temp[i]/a;
      store[i+*nr]=temp2[i]/b;
    }
    
    // That completes the mean calculations.
    
    for (i=0;i<m;i++){
      var[i%k]=var[i%k]+pow((pmatrix[i]-store[i%k]),2);
    }
    for (i=m;i<n;i++){
      var2[i%k]=var2[i%k]+pow((pmatrix[i]-store[(i%k)+*nr]),2);
    }
    
    for (i=0;i<=*nr-1;i++){
      store[i+2*k]=var[i]/(a-1);
      store[i+3*k]=var2[i]/(b-1);
    }
    // That completes var calculations.
}

void start(double *Imatrix,int *g,int *nr,int *nc,double *initial,
		double storage[][9]){
  int i, a = *nr;
  a*=4;
  
  double store[a], tool[a], C1[*nr][3], C2[*nr][3];
  double nrows,ncols,gvalue, xbardiff=0, denom=0;
  
  nrows = (double) *nr;
  ncols = (double) *nc;
  gvalue= (double) *g;
  
  meanvar(Imatrix,g,nr,nc,store);
  
  for(i=0;i<=a-1;i++){
    tool[i]=store[i];
  }
  for (i=0; i<*nr;i++){
    C1[i][0]=tool[i]; //mean group 1
    storage[i][0]=C1[i][0];
    C1[i][1]=tool[i+*nr+*nr]; // var group 1
    storage[i][1]=C1[i][1];
    C1[i][2]=C1[i][1]/(gvalue-1);
    storage[i][2]=sqrt(C1[i][2]);
    
    C2[i][0]=tool[i+*nr]; // mean group 2
    storage[i][4]=C2[i][0];    
    C2[i][1]=tool[i+*nr+*nr+*nr]; // var group 2 
    storage[i][5]=C2[i][1];        
    C2[i][2]=C2[i][1]/(ncols-gvalue+1);
    storage[i][6]=sqrt(C2[i][2]);   
  }
  for (i=0; i<*nr; i++){
    xbardiff = C1[i][0]-C2[i][0];
    denom = sqrt(C1[i][2]+C2[i][2]);
    initial[i]=fabs(xbardiff/denom);
  }											
}
