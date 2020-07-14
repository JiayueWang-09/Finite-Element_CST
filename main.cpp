#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Mesh.h"
#include <iostream>
using namespace std;

void main()
{
	//Geometry//
	double t=1,L=100,W=100;
	double NumL=10,NumW=10,detJ=L*W/(NumL*NumW);
	
	//Material//
	double E=200000,v=0.33;
	double **D;
	D=(double**)malloc(3*sizeof(double*));
	for(int i=0;i<3;i++)
	{
		*(D+i)=(double*)malloc(3*sizeof(double));
		for(int j=0;j<3;j++)
		{
			D[i][j]=0.0;
		}
	}
	D[0][0]=D[1][1]=E/(1-v*v);
	D[0][1]=D[1][0]=v*E/(1-v*v);
	D[2][2]=E/(2*(1+v));
	
	cout<<"Elastic Matrix (MPa)"<<endl;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout<<D[i][j]<<" ";
		}
		cout<<endl;
	}

	//Mesh and Stiffness Matrix Generation//
	double ***J,***B,***Ke;
	double **Kg;
	int **Nodes;
	Nodes=(int**)malloc(2*NumL*NumW*sizeof(int*));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		*(Nodes+i)=(int*)malloc(3*sizeof(int));
	}
	
	J = (double***)malloc(2*NumL*NumW*sizeof(double**));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		*(J+i) = (double**)malloc(2*sizeof(double*));
		for(int j=0;j<2;j++)
		{
			*(*(J+i)+j) = (double*)malloc(2*sizeof(double));
			for(int k=0;k<2;k++)
			{
				J[i][j][k] = 0.0;
			}
		}
	}

	B = (double***)malloc(2*sizeof(double**));
	for(int i=0;i<2;i++)
	{
		*(B+i) = (double**)malloc(3*sizeof(double*));
		for(int j=0;j<3;j++)
		{
			*(*(B+i)+j) = (double*)malloc(6*sizeof(double));
			for(int k=0;k<6;k++)
			{
				B[i][j][k] = 0.0;
			}
		}
	}

	Ke = (double***)malloc(2*sizeof(double**));
	for(int i=0;i<2;i++)
	{
		*(Ke+i) = (double**)malloc(6*sizeof(double*));
		for(int j=0;j<6;j++)
		{
			*(*(Ke+i)+j) = (double*)malloc(6*sizeof(double));
			for(int k=0;k<6;k++)
			{
				Ke[i][j][k] = 0.0;
			}
		}
	}

	Kg=(double**)malloc(2*(NumL+1)*(NumW+1)*sizeof(double*));
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		*(Kg+i)=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			Kg[i][j]=0.0;
		}
	}

	CSTgeneration(L,W,NumL,NumW,J,detJ,Nodes);
	B_Matrix(B,J,detJ);
	ElementalStiffnessMatrix(Ke,D,t,B,detJ);
	AssemblyGlobalMatrix(Ke,NumW,NumL,Kg,Nodes);

	//Displacement and force//
	double Ta=1,Tb=2;
	double *Force,*Disp,*Force_Reduced,*Disp_Reduced;
	double **Kg_Reduced;
	
	Force=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
	Disp=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		Force[i]=0.0;
		Disp[i]=0.0;
	}
	Force_Reduced=(double*)malloc(2*NumL*(NumW+1)*sizeof(double));
	Disp_Reduced=(double*)malloc(2*NumL*(NumW+1)*sizeof(double));
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		Force_Reduced[i]=0.0;
		Disp_Reduced[i]=0.0;
	}
	Kg_Reduced=(double**)malloc(2*NumL*(NumW+1)*sizeof(double*));
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		*(Kg_Reduced+i)=(double*)malloc(2*NumL*(NumW+1)*sizeof(double));
		for(int j=0;j<2*NumL*(NumW+1);j++)
		{
			Kg_Reduced[i][j]=0.0;
		}
	}

    EquivalentNodalForce(Ta,Tb,L,W,NumW,t,Force);
	ReduceLinearSystem(Kg,Kg_Reduced,Disp,Force,Disp_Reduced,Force_Reduced,NumW,NumL);
	Gauss(Kg_Reduced,Disp_Reduced,Force_Reduced,2*NumL*(NumW+1));
	Displacement(Disp,Disp_Reduced,NumW,NumL);
	NodeForce(Force,Kg,Disp,NumL,NumW);


	//Stress and Strain//
	double *Von_Mises,*StrainEnergy;
	Von_Mises=(double*)malloc(2*NumL*NumW*sizeof(double));
	StrainEnergy=(double*)malloc(2*NumL*NumW*sizeof(double));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		Von_Mises[i]=0.0;
		StrainEnergy[i]=0.0;
	}
	double **Strain,**Stress;
	Strain=(double**)malloc(2*NumL*NumW*sizeof(double*));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		*(Strain+i)=(double*)malloc(3*sizeof(double));
		for(int j=0;j<3;j++)
		{
			Strain[i][j]=0.0;
		}
	}
	Stress=(double**)malloc(2*NumL*NumW*sizeof(double*));
	for(int i=0;i<2*NumL*NumW;i++)
	{
		*(Stress+i)=(double*)malloc(3*sizeof(double));
		for(int j=0;j<3;j++)
		{
			Stress[i][j]=0.0;
		}
	}
	
	StressStrain(Strain,Stress,StrainEnergy,Ke,B,Disp,Nodes,NumW,NumL,D);
	VonMisesStress(Von_Mises,Stress,NumW,NumL);


	system("Pause");
}