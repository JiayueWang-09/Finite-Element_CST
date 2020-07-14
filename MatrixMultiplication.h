#ifndef _MATRIXMULTIPLICATION_ 
#define _MATRIXMULTIPLICATION_ 

void MatrixMultiply(double **A,int m1,int n1,double **B,int m2,int n2,double **Result);
void Matrix_vector(double **A,int m,int n,double *b,double *result);
void vectorT_vector(double *a,double *b,int n,double result);

#endif