#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_MATRIX_LINE_LENGTH 102400
#define DELIM ","
#define DELIM2 ",\n"

int allocateMatrix2di(int ***paanInt, int m, int n)
{
  int ** aanTempMatrix;
  
  aanTempMatrix = (int **) malloc(m*sizeof(int*));
  if(aanTempMatrix){
    int i = 0;
    int j = 0;
    
    for(i = 0; i < m; i++){
      aanTempMatrix[i] = (int *) malloc(n*sizeof(int));
      
      if(!aanTempMatrix[i]){
	for(j = 0 ; j < i; j++){
	  free(aanTempMatrix[j]);
	  aanTempMatrix[i] = NULL;
	  free(aanTempMatrix);
	  aanTempMatrix = NULL;
	  fprintf(stderr,"Failed to allocate memory in allocateMatrix2di\n");
	  fflush(stderr);
	  return 0;
	}
      }
    }
  }
  else{
    fprintf(stderr,"Failed to allocate memory in allocateMatrix2di\n");
    fflush(stderr);
    return 0;
  }  
  
  (*paanInt) = aanTempMatrix;

  return 1;
}

int reallocateMatrix2di(int ***paanInt, int m1, int n1,int m2,int n2)
{
  int ** aanTempMatrix = *paanInt;
  
  aanTempMatrix = (int **) realloc(aanTempMatrix,m2*sizeof(int*));
  if(aanTempMatrix){
    int i = 0;
    int j = 0;
    
    for(i = 0; i < m1; i++){
      aanTempMatrix[i] = (int *) realloc(aanTempMatrix[i],n2*sizeof(int));
      
      if(!aanTempMatrix[i]){
	for(j = 0 ; j < i; j++){
	  free(aanTempMatrix[j]);
	  aanTempMatrix[i] = NULL;
	  free(aanTempMatrix);
	  aanTempMatrix = NULL;
	  fprintf(stderr,"Failed to reallocate memory in allocateMatrix2di\n");
	  fflush(stderr);
	  return 0;
	}
      }
    }
    
    for(i = m1; i < m2; i++){
      aanTempMatrix[i] = (int *) malloc(n2*sizeof(int));
      
      if(!aanTempMatrix[i]){
	for(j = 0 ; j < i; j++){
	  free(aanTempMatrix[j]);
	  aanTempMatrix[i] = NULL;
	  free(aanTempMatrix);
	  aanTempMatrix = NULL;
	  fprintf(stderr,"Failed to reallocate memory in allocateMatrix2di\n");
	  fflush(stderr);
	  return 0;
	}
      }
    }
  }
  else{
    fprintf(stderr,"Failed to reallocate memory in allocateMatrix2di\n");
    fflush(stderr);
    return 0;
  }  
  
  (*paanInt) = aanTempMatrix;
  return 1;
}

int callocateMatrix2di(int ***paanInt, int m, int n)
{
  int ** aanTempMatrix = NULL;
  int    i             = 0;
  int    j             = 0;

  allocateMatrix2di(&aanTempMatrix,m,n);
  for(i = 0; i < m;i++){
    for(j = 0; j < n;j++){
      aanTempMatrix[i][j] = 0;
    }
  }
  
  (*paanInt) = aanTempMatrix;

  return 1;
}

int recallocateMatrix2di(int ***paanInt,int m1,int n1,int m2,int n2)
{
  int **aanTempMatrix = *paanInt;
  
  aanTempMatrix = (int **) realloc(aanTempMatrix,m2*sizeof(int*));
  if(aanTempMatrix){
    int i = 0;
    int j = 0;
    
    for(i = 0; i < m1; i++){
      aanTempMatrix[i] = (int *) realloc(aanTempMatrix[i],n2*sizeof(int));
      
      if(!aanTempMatrix[i]){
	for(j = 0 ; j < i; j++){
	  free(aanTempMatrix[j]);
	  aanTempMatrix[i] = NULL;
	  free(aanTempMatrix);
	  aanTempMatrix = NULL;
	  fprintf(stderr,"Failed to reallocate memory in recallocateMatrix2di\n");
	  fflush(stderr);
	  return 0;
	}
      }
      else{
	for(j = n1 ; j < n2;j++){
	  aanTempMatrix[i][j] = 0;
	}
      }
    }
    
    for(i = m1 ; i < m2; i++){
      
      aanTempMatrix[i] = (int *) malloc(n2*sizeof(int));
      
      if(!aanTempMatrix[i]){
	for(j = 0 ; j < i; j++){
	  free(aanTempMatrix[j]);
	  aanTempMatrix[i] = NULL;
	  free(aanTempMatrix);
	  aanTempMatrix = NULL;
	  fprintf(stderr,"Failed to reallocate memory in recallocateMatrix2di\n");
	  fflush(stderr);
	  return 0;
	}
      }
      else{
	for(j = 0 ; j < n2;j++){
	  aanTempMatrix[i][j] = 0;
	}
      }
    }
  }
  else{
    fprintf(stderr,"Failed to reallocate memory in recallocateMatrix2di\n");
    fflush(stderr);
    return 0;
  }  
  
  (*paanInt) = aanTempMatrix;  

  return 1;
}

/*reads 2di matrix*/
/* assumes comma column and new line row*/
int readMatrix2di(FILE *ifp, int **aanMatrix, int m, int n)
{
  int i,j;
  char* szBuffer = (char *)malloc(MAX_MATRIX_LINE_LENGTH*sizeof(char));
  char* tok = NULL;
  char* cError;
  if(m <=0 || n <= 0)
    return 0;

  for(i = 0; i < m; i++){
    if(fgets(szBuffer, MAX_MATRIX_LINE_LENGTH, ifp)){
      tok = strtok(szBuffer,DELIM);
      aanMatrix[i][0] = strtol(tok,&cError,10);
      if(*cError!='\0')
	return 0;
      
      for(j = 1; j < (n - 1); j++){
	tok = strtok(NULL,DELIM);
	aanMatrix[i][j] = strtol(tok,&cError,10);
	if(*cError!='\0')
	  return 0;
      }

      tok = strtok(NULL,DELIM);
      aanMatrix[i][j] = strtol(tok,&cError,10);
      if(*cError!='\n')
	  return 0;      

    }
    else{
      return 0;
    }
  }
  
  free(szBuffer);
  return 1;
}

void writeMatrix2di(FILE *ofp, int **aanMatrix, int m, int n)
{
  int i,j;

  for(i = 0; i < m; i++){
    for(j = 0; j < (n -1) ; j++){
      fprintf(ofp, "%d,",aanMatrix[i][j]);
    }    
    fprintf(ofp,"%d\n",aanMatrix[i][j]);
  }
}

/*reads 2df matrix*/
/* assumes comma column and new line row*/
int readMatrix2df(FILE *ifp, float **aafMatrix, int m, int n)
{
  int i,j;
  char* szBuffer = (char *)malloc(MAX_MATRIX_LINE_LENGTH*sizeof(char));
  char* tok = NULL;
  char* cError;
  if(m <=0 || n <= 0)
    return 0;

  for(i = 0; i < m; i++){
    if(fgets(szBuffer, MAX_MATRIX_LINE_LENGTH, ifp)){
      tok = strtok(szBuffer,DELIM);
      aafMatrix[i][0] = strtod(tok,&cError);
      if(*cError!='\0')
	return 0;
      
      for(j = 1; j < (n - 1); j++){
	tok = strtok(NULL,DELIM);
	aafMatrix[i][j] = strtod(tok,&cError);
	if(*cError!='\0')
	  return 0;
      }

      tok = strtok(NULL,DELIM);
      aafMatrix[i][j] = strtod(tok,&cError);
      if(*cError!='\n')
	  return 0;      

    }
    else{
      return 0;
    }
  }
  
  free(szBuffer);
  return 1;
}

void writeMatrix2df(FILE *ofp, float **aafMatrix, int m, int n)
{
  int i,j;

  for(i = 0; i < m; i++){
    for(j = 0; j < (n -1) ; j++){
      fprintf(ofp, "%e,",aafMatrix[i][j]);
    }    
    fprintf(ofp,"%f\n",aafMatrix[i][j]);
  }
}

/* destroys square int matrix*/
int destroyMatrix2di(int ***paanMatrix,int m)
{
  int           i = 0;
  int **aanMatrix = (*paanMatrix);

  for(i = 0; i < m; i++)
  {
    free(aanMatrix[i]);
    aanMatrix[i] = NULL;
  }
  
  free(aanMatrix);
  aanMatrix = NULL;

  return 1;
}

int allocateMatrix2df(float ***paafMatrix,int n, int m)
{
  float **aafTempMatrix = NULL;

  aafTempMatrix = (float **) malloc(n*sizeof(float*));

  if(aafTempMatrix){
    int i = 0;
    int j = 0;

    for(i = 0; i < n; i++){
      aafTempMatrix[i] = (float *) malloc(sizeof(float)*m);
      if(!aafTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aafTempMatrix[j]);
	  aafTempMatrix[j] = NULL;
	}
	
	free(aafTempMatrix);
	aafTempMatrix = NULL;
	fprintf(stderr,"Failed in allocateMatrix2df\n");
	return 0;
      
      }
    }
  }
  else{
    fprintf(stderr,"Failed to allocate memory in allocateMatrix2df\n");
    return 0;
  }
  
  (*paafMatrix) = aafTempMatrix;
  return 1;
}

/* reallocates 2d float matrix*/
int reallocateMatrix2df(float ***paafMatrix,int m1, int n1,int m2,int n2)
{
  float **aafTempMatrix = *paafMatrix;

  aafTempMatrix = (float **) realloc(aafTempMatrix,m2*sizeof(float*));

  if(aafTempMatrix){
    int i = 0;
    int j = 0;

    for(i = 0; i < m1; i++){
      aafTempMatrix[i] = (float *) realloc(aafTempMatrix[i],sizeof(float)*n2);
      if(!aafTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aafTempMatrix[j]);
	  aafTempMatrix[j] = NULL;
	}
	
	free(aafTempMatrix);
	aafTempMatrix = NULL;
	fprintf(stderr,"Failed in reallocateMatrix2df\n");
	return 0;
      }
    }
    
    for(i = m1; i < m2; i++){
      aafTempMatrix[i] = (float *) malloc(sizeof(float)*n2);
      if(!aafTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aafTempMatrix[j]);
	  aafTempMatrix[j] = NULL;
	}
	
	free(aafTempMatrix);
	aafTempMatrix = NULL;
	fprintf(stderr,"Failed in reallocateMatrix2df\n");
	return 0;
      }
    }
  }
  else{
    fprintf(stderr,"Failed to allocate memory in allocateMatrix2df\n");
    return 0;
  }
  
  (*paafMatrix) = aafTempMatrix;
  return 1;
}

/* destroys a square float matrix*/
int destroyMatrix2df(float ***paafMatrix, int m)
{
  int   i           = 0;
  float **aafMatrix = (*paafMatrix);

  for(i = 0; i < m;i++)
  {
    free(aafMatrix[i]);
    aafMatrix[i] = NULL;
  }

  free(aafMatrix);
  aafMatrix = NULL;
  
  return 1;
}

int allocateMatrix2dd(double  ***paadMatrix,int n, int m)
{
  double **aadTempMatrix = NULL;

  aadTempMatrix = (double **) malloc(n*sizeof(double*));

  if(aadTempMatrix){
    int i = 0;
    int j = 0;

    for(i = 0; i < n; i++){
      aadTempMatrix[i] = (double *) malloc(sizeof(double)*m);
      if(!aadTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aadTempMatrix[j]);
	  aadTempMatrix[j] = NULL;
	  free(aadTempMatrix);
	  aadTempMatrix = NULL;
	  fprintf(stderr,"Failed in allocateMatrix2df\n");
	  return 0;
	}
      }
    }
  }
  else{
    fprintf(stderr,"Failed to allocate memory in allocateMatrix2df\n");
    return 0;
  }
  
  (*paadMatrix) = aadTempMatrix;
  return 1;
}

int callocateMatrix2dd(double  ***paadMatrix,int n, int m)
{
  int i = 0, j = 0;
  int ret = allocateMatrix2dd(paadMatrix,n,m);

  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      (*paadMatrix)[i][j] = 0.0;
    }
  }

  return ret;
}

/* reallocates 2d double matrix*/
int reallocateMatrix2dd(double ***paadMatrix,int m1, int n1,int m2,int n2)
{
  double **aadTempMatrix = *paadMatrix;

  aadTempMatrix = (double **) realloc(aadTempMatrix,m2*sizeof(double*));

  if(aadTempMatrix){
    int i = 0;
    int j = 0;

    for(i = 0; i < m1; i++){
      aadTempMatrix[i] = (double *) realloc(aadTempMatrix[i],sizeof(double)*n2);
      if(!aadTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aadTempMatrix[j]);
	  aadTempMatrix[j] = NULL;
	}
	
	free(aadTempMatrix);
	aadTempMatrix = NULL;
	fprintf(stderr,"Failed in reallocateMatrix2dd\n");
	return 0;
      }
    }
    
    for(i = m1; i < m2; i++){
      aadTempMatrix[i] = (double *) malloc(sizeof(double)*n2);
      if(!aadTempMatrix[i]){
	for(j = 0; j < i ; j++){
	  free(aadTempMatrix[j]);
	  aadTempMatrix[j] = NULL;
	}
	
	free(aadTempMatrix);
	aadTempMatrix = NULL;
	fprintf(stderr,"Failed in reallocateMatrix2df\n");
	return 0;
      }
    }
  }
  else{
    fprintf(stderr,"Failed to allocate memory in allocateMatrix2df\n");
    return 0;
  }
  
  (*paadMatrix) = aadTempMatrix;
  return 1;
}

int readMatrix2dd(FILE *ifp, double **aadMatrix, int m, int n)
{
  int i,j;
  char* szBuffer = (char *)malloc(MAX_MATRIX_LINE_LENGTH*sizeof(char));
  char* tok = NULL;
  char* cError;
  if(m <=0 || n <= 0)
    return 0;

  for(i = 0; i < m; i++){
    if(fgets(szBuffer, MAX_MATRIX_LINE_LENGTH, ifp)){
      tok = strtok(szBuffer,DELIM2);
      aadMatrix[i][0] = strtod(tok,&cError);
      if(*cError!='\0')
	return 0;
      
      for(j = 1; j < n; j++){
	tok = strtok(NULL,DELIM2);
	aadMatrix[i][j] = strtod(tok,&cError);
	if(*cError!='\0')
	  return 0;
      }

    }
    else{
      return 0;
    }
  }
  
  free(szBuffer);
  return 1;
}

void writeMatrix2dd(FILE *ofp,double **aadMatrix, int m, int n)
{
  int i,j;

  for(i = 0; i < m; i++){
    for(j = 0; j < (n -1) ; j++){
      fprintf(ofp, "%e,",aadMatrix[i][j]);
    }    
    fprintf(ofp,"%e\n",aadMatrix[i][j]);
  }
}

/* destroys square double matrix*/
int destroyMatrix2dd(double ***paadMatrix, int m)
{
  int   i           = 0;
  double **aadMatrix = (*paadMatrix);

  for(i = 0; i < m;i++)
  {
    free(aadMatrix[i]);
    aadMatrix[i] = NULL;
  }

  free(aadMatrix);
  aadMatrix = NULL;
  
  return 1;
}

/*  allocates 1d float matrix */
int allocateMatrix1df(float **pafMatrix,int m)
{
  float *afMatrix;

  afMatrix = (float *) malloc(m*sizeof(float));
  
  if(!afMatrix){
    fprintf(stderr,"Failed to allocate memory in allocateMatrix1df");
    fflush(stderr);
    return 0;
  }
  
  (*pafMatrix) = afMatrix;
  return 1;
}

/*  allocates 1d float matrix */
int allocateMatrix1dd(double **padMatrix,int m)
{
  double *adMatrix;

  adMatrix = (double *) malloc(m*sizeof(double));
  
  if(!adMatrix){
    fprintf(stderr,"Failed to allocate memory in allocateMatrix1df");
    fflush(stderr);
    return 0;
  }
  
  (*padMatrix) = adMatrix;
  return 1;
}

int callocateMatrix1dd(double **padMatrix,int m)
{
  int i = 0;
  int ret = allocateMatrix1dd(padMatrix, m);

  for(i = 0; i < m; i++){
    (*padMatrix)[i] = 0.0;
  }

  return ret;
}

/*  reallocates 1d float matrix */
int reallocateMatrix1df(float **pafMatrix,int m)
{
  float *afMatrix = *pafMatrix;

  afMatrix = (float *) realloc(afMatrix,m*sizeof(float));
  
  if(!afMatrix)
  {
    fprintf(stderr,"Failed to allocate memory in allocateMatrix1df");
    fflush(stderr);
    return 0;
  }
  
  (*pafMatrix) = afMatrix;
  return 1;
}

/*  reallocates 1d double matrix */
int reallocateMatrix1dd(double **padMatrix,int m)
{
  double *adMatrix = *padMatrix;

  adMatrix = (double *) realloc(adMatrix,m*sizeof(double));
  
  if(!adMatrix)
  {
    fprintf(stderr,"Failed to allocate memory in allocateMatrix1dd");
    fflush(stderr);
    return 0;
  }
  
  (*padMatrix) = adMatrix;
  return 1;
}

int destroyMatrix1df(float **pafMatrix)
{
  free((*pafMatrix));
  (*pafMatrix) = NULL;
  return 1;
}

/* allocate 1d int matrix*/
int allocateMatrix1di(int **panMatrix,int m)
{
  int  *anMatrix;

  anMatrix = (int *) malloc(m*sizeof(int));
  
  if(!anMatrix)
  {
    fprintf(stderr,"Failed to allocate memory in allocateMatrix1di");
    fflush(stderr);
    return 0;
  }
  
  (*panMatrix) = anMatrix;
  return 1;
}

/* clear and allocate 1d int matrix*/
int callocateMatrix1di(int **panMatrix,int m)
{
  int  *anMatrix;
  int  i = 0;

  anMatrix = (int *) malloc(m*sizeof(int));
  
  if(!anMatrix)
  {
    fprintf(stderr,"Failed to allocate memory in callocateMatrix1di");
    fflush(stderr);
    return 0;
  }
  
  for(i = 0; i < m; i++){
    anMatrix[i] = 0;
  }

  (*panMatrix) = anMatrix;
  return 1;
}

/* clear and reallocate 1d int matrix*/
int recallocateMatrix1di(int **panMatrix,int m1, int m2)
{
  int  *anMatrix = *panMatrix;
  int  i = 0;

  anMatrix = (int *) realloc(anMatrix, m2*sizeof(int));
  
  if(!anMatrix)
  {
    fprintf(stderr,"Failed to allocate memory in callocateMatrix1di");
    fflush(stderr);
    return 0;
  }
  
  for(i = m1; i < m2; i++){
    anMatrix[i] = 0;
  }

  (*panMatrix) = anMatrix;
  return 1;
}

