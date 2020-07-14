#include <iostream>
using namespace std;

void MatrixMultiply(double **A,int m1,int n1,double **B,int m2,int n2,double **Result)
{
	double sum=0.0;
	if(n1!=m2)
	{
		cout<<"The given matrix can't mutliplied together"<<endl;
	}
	else
	{
		for(int i=0;i<m1;i++)
		{
			for(int j=0;j<n2;j++)
			{
				for(int k=0;k<n1;k++)
				{
					sum+=A[i][k]*B[k][j];
				}
				Result[i][j]=sum;
				sum=0.0;
			}
		}
	}
}

void Matrix_vector(double **A,int m,int n,double *b,double *result)
{
	for(int i=0;i<m;i++)
	{
		double sum=0.0;
		for(int j=0;j<n;j++)
		{
			sum+=A[i][j]*b[j];
		}
		result[i]=sum;
	}
}

void vectorT_vector(double *a,double *b,int n,double result)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum+=a[i]*b[i];
	}
	result=sum;
}
