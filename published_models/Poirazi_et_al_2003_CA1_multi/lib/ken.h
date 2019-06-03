#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
/*#include "/home/niebur/tools/graphics/xview/header.h"*/

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif
#define SQU(a) ((a)*(a))
#define SQUARE(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define RINT(x) ( ((x) > 0) ?  \
( ((x) - ((int)(x))) < 0.5 ? ((int)(x)) : ((int)(x)+1) ) :  \
( ((x) - ((int)(x))) > (-0.5) ? ((int)(x)) : ((int)(x)-1) ) )
#ifndef PI
#define PI 3.14159265359
#endif
#define ABS(x) ( (x) > 0. ? (x) : (-(x)) )
#define EXP(x) (exp((double)(x)))
#define SQRT(x) ( ((x)<0.) ? FUNCERROR("sqrt",((double)x)) : sqrt((double)(x)))
#define LOG(x) ( ((x)<0.) ? FUNCERROR("log",((double)x)) : log((double)(x)))
#define LOG10(x) ( ((x)<0.) ? FUNCERROR("log10",((double)x)) : log10((double)(x)))
#define POW(x,y) ( ((x)<0.) ? ( ((y) == (int)(y)) ? \
 pow((double)(x),(double)(y)) :  FUNC2ERROR("pow",(double)(x),(double)(y)) ) \
 : pow((double)(x),(double)(y)) )
#define COS(x) (cos((double)(x)))
#define SIN(x) (sin((double)(x)))
#define TAN(x) (tan((double)(x)))
#define APRECISION 0.000000000001
#define ALIMIT(x,y) ( ABS(x) > (y) ? ((x) > 0 ? (y) : (-(y))) : (x) )
#define ACOS(x) ( (ABS(x) > (1. + APRECISION)) ? \
 FUNCERROR("acos",((double)(x))) : acos(ALIMIT((x),(1.-APRECISION))))
#define ASIN(x) ( (ABS(x) > (1. + APRECISION)) ? \
 FUNCERROR("asin",((double)(x))) : asin(ALIMIT((x),(1.-APRECISION))))
#define ATAN(x) (atan((double)(x)))
#define ATAN2(y,x) (atan2((double)y,(double)x))

#define FUNCERROR(func,dub) (fprintf(stderr,"Warning! Function %s failed: argument = %g; source file %s, line %d;\n", \
 func,dub,__FILE__,__LINE__))
#define FUNC2ERROR(func,dub1,dub2) (fprintf(stderr,"Warning: Function %s failed: arguments = %g, %g; source file %s, line %d\n", \
 func,dub1,dub2,__FILE__,__LINE__))

#ifdef DEBUG
#define DPRINTF(x) printf x ; fflush(stdout)  /* use as DPRINTF((" ", )), i.e. two parens */
#else
#define DPRINTF(x)  /* do nothing */
#endif

#define OPENFILE(fpr,fname,fmode,action){\
if((fpr = fopen(fname,fmode)) == NULL) {\
	char errorstring[120];\
	sprintf(errorstring, "Couldn't open file %s, mode %s; source file %s, line %d", fname, fmode, __FILE__, __LINE__);\
	perror(errorstring);\
	action;\
}\
}

#define READFILE(pointer, ptype, numread, fname, action){\
if (1) {\
	int readcount;\
	if ( (readcount = fread(pointer, sizeof(ptype), (unsigned)(numread), fname)) < numread ) { \
		char errorstring[120];\
		sprintf(errorstring, "fread fails; read %d, wanted %d; source file %s, line %d", \
		readcount, numread, __FILE__, __LINE__);\
		perror(errorstring);\
		action;\
	}\
}\
}

#define WRITEFILE(pointer, ptype, numwrite, fname, action){\
if (1) {\
	int writecount;\
	if ( (writecount = fwrite(pointer, sizeof(ptype), (unsigned)(numwrite), fname)) < numwrite ) { \
		char errorstring[120];\
		sprintf(errorstring, "fwrite fails; wrote %d, wanted %d; source file %s, line %d", \
		writecount, numwrite, __FILE__, __LINE__);\
		perror(errorstring);\
		action;\
	}\
}\
}

#define ALLOCATE(pointer, ptype, numelements, action){\
if ( (pointer = \
(ptype *)calloc((unsigned)numelements,(unsigned)sizeof(ptype))) == NULL) { \
	char errorstring[120];\
	sprintf(errorstring, "calloc fails; source file %s, line %d", \
	__FILE__, __LINE__);\
	perror(errorstring);\
	action;\
}\
}

/*  CONMATRIX: contiguous matrix pmat[nrow][ncol] 
    mtype is int, double, etc; declare: mtype **pmat, *parray 
    pmat[i][j] will then refer to correct matrix element,
    which is parray[i*ncol + j] 
*/

#define CONMATRIX(pmat, parray, mtype, nrow, ncol){\
if (1) {\
	int rowcount;\
	ALLOCATE(parray,mtype,((nrow)*(ncol)), exit(-1));\
	ALLOCATE(pmat,mtype *,(nrow),exit(-1));\
	for(rowcount=0; rowcount< (nrow); rowcount++) {\
		pmat[rowcount]= &parray[rowcount*(ncol)];\
	}\
}\
}

/* CONTHREE: contiguous three-index object pmat[n1][n2][n3] 
   mtype ***pmat, **ptemp, *parray 
   pmat[i][j][k] will then refer to correct matrix element,
   which is ptemp[i*n2 + j][k] = parray[i*n2*n3 + j*n3 + k]
*/
   
#define CONTHREE(pmat, ptemp, parray, mtype, n1, n2, n3){\
if (1) {\
	int count1,count2;\
	ALLOCATE(parray,mtype,((n1)*(n2)*(n3)), exit(-1));\
	ALLOCATE(ptemp,mtype *,((n1)*(n2)),exit(-1));\
	ALLOCATE(pmat,mtype **,(n1),exit(-1));\
	count2 = (n1)*(n2);\
	for(count1=0; count1< count2; count1++) {\
		ptemp[count1]= &parray[count1*(n3)];\
	}\
	for(count1=0; count1< (n1); count1++) {\
		pmat[count1]= &ptemp[count1*(n2)];\
	}\
}\
}

/* CONFOUR: contiguous four-index object pmat[n1][n2][n3][n4] 
   mtype ****pmat, ***pthree, **ptwo, *parray 
Then
   pmat[i][j][k][l] = pthree[i*n2 + j][k][l] =
   ptwo[i*n2*n3 + j*n3 + k][l] = parray[i*n2*n3*n4 + j*n3*n4 + k*n4 + l]
*/
   
#define CONFOUR(pmat, pthree, ptwo, parray, mtype, n1, n2, n3, n4){\
if (1) {\
	int count1,count2;\
	ALLOCATE(parray,mtype,((n1)*(n2)*(n3)*(n4)), exit(-1));\
	ALLOCATE(ptwo,mtype *,((n1)*(n2)*(n3)),exit(-1));\
	ALLOCATE(pthree,mtype **,((n1)*(n2)),exit(-1));\
	ALLOCATE(pmat,mtype ***,(n1),exit(-1));\
	count2 = (n1)*(n2)*(n3);\
	for(count1=0; count1< count2; count1++) {\
		ptwo[count1]= &parray[count1*(n4)];\
	}\
	count2 = (n1)*(n2);\
	for(count1=0; count1< count2; count1++) {\
		pthree[count1]= &ptwo[count1*(n3)];\
	}\
	for(count1=0; count1< (n1); count1++) {\
		pmat[count1]= &pthree[count1*(n2)];\
	}\
}\
}




