#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include<math.h>
#include "MatrixMultiplication.h"

using namespace std;

void CSTgeneration(double L,double W,int NumL,int NumW,double ***J,double detJ,int **Node)
{
	for(int i=0;i<NumL*NumW;i++)
	{
		
			J[2*i][0][0]=0;
			J[2*i][0][1]=-W/NumW;
			J[2*i][1][0]=L/NumL;
			J[2*i][1][1]=0;
			
			J[2*i+1][0][0]=-L/NumL;
			J[2*i+1][0][1]=-W/NumW;
			J[2*i+1][1][0]=0;
			J[2*i+1][1][1]=-W/NumW;

			detJ=L*W/(NumL*NumW);
	}
	cout<<"detJ = "<<detJ<<endl;

	cout<<"i j k"<<endl;
	for(int i=0;i<NumW;i++)
	{
		for(int j=0;j<NumL;j++)
		{
			Node[2*(i*NumW+j)][0]=i*(NumW+1)+j;
			Node[2*(i*NumW+j)][1]=(i+1)*(NumW+1)+j+1;
			Node[2*(i*NumW+j)][2]=(i+1)*(NumW+1)+j;

			Node[2*(i*NumL+j)+1][0]=i*(NumW+1)+j;
			Node[2*(i*NumL+j)+1][1]=i*(NumW+1)+j+1;
			Node[2*(i*NumL+j)+1][2]=(i+1)*(NumW+1)+j+1;

			cout<<Node[2*(i*NumL+j)][0]<<" "<<Node[2*(i*NumL+j)][1]<<" "<<Node[2*(i*NumL+j)][2]<<endl;
			cout<<Node[2*(i*NumL+j)+1][0]<<" "<<Node[2*(i*NumL+j)+1][1]<<" "<<Node[2*(i*NumL+j)+1][2]<<endl;
		}
	}
	cout<<2*NumW*NumL<<" elements are seperated."<<endl;
}

void B_Matrix(double ***B,double ***J,double detJ)
{
	for(int i=0;i<2;i++)
	{
		B[i][0][0]=J[i][1][1]/detJ;
		B[i][0][2]=-J[i][0][1]/detJ;
		B[i][0][4]=(J[i][0][1]-J[i][1][1])/detJ;
		B[i][1][1]=-J[i][1][0]/detJ;
		B[i][1][3]=J[i][0][0]/detJ;
		B[i][1][5]=(J[i][1][0]-J[i][0][0])/detJ;
		B[i][2][0]=-J[i][1][0]/detJ;
		B[i][2][1]=J[i][1][1]/detJ;
		B[i][2][2]=J[i][0][0]/detJ;
		B[i][2][3]=-J[i][0][1]/detJ;
		B[i][2][4]=(J[i][1][0]-J[i][0][0])/detJ;
		B[i][2][5]=(J[i][0][1]-J[i][1][1])/detJ;
	}

	for(int i=0;i<2;i++)
	{
		cout<<"B matrix for element"<<i+1<<endl;
		for(int m=0;m<3;m++)
		{
			for(int n=0;n<6;n++)
			{
				cout<<B[i][m][n]<<" ";
			}
			cout<<endl;
		}
	}

}

void ElementalStiffnessMatrix(double ***Ke,double **D,double t,double ***B,double detJ)
{
	double **Temp;
	Temp=(double**)malloc(3*sizeof(double*));
	for(int i=0;i<3;i++)
	{
		*(Temp+i)=(double*)malloc(6*sizeof(double));
		for(int j=0;j<6;j++)
		{
			Temp[i][j]=0.0;
		}
	}
	

	for(int i=0;i<2;i++)
	{
		MatrixMultiply(D,3,3,B[i],3,6,Temp);
			for(int j=0;j<6;j++)
			{		
				for(int k=0;k<6;k++)
				{
					double sum=0.0;
					for(int m=0;m<3;m++)
					{
						sum+=0.5*detJ*t*B[i][m][j]*Temp[m][k];
					}
					Ke[i][j][k]=sum;		
				}
			}
	}

	for(int i=0;i<2;i++)
	{
		cout<<"The elemental stiffness matrix for element"<<i+1<<endl;
		for(int m=0;m<6;m++)
		{
			for(int n=0;n<6;n++)
			{
				cout<<Ke[i][m][n]<<" ";
			}
			cout<<endl;
		}
	}

}

void EquivalentNodalForce(double Ta,double Tb,double L,double W,int NumW,int t,double *Fe)
{
	double s=0,c=1;
	double Tx1,Tx2,Ty1,Ty2;
	for(int i=0;i<NumW;i++)
	{
	    Tx1=c*(Tb+(Ta-Tb)*i/NumW);
		Tx2=c*(Tb+(Ta-Tb)*(i+1)/NumW);
		Ty1=s*(Tb+(Ta-Tb)*i/NumW);
		Ty2=s*(Tb+(Ta-Tb)*(i+1)/NumW);

		Fe[2*(NumW+1)*i+2*NumW]+=t*(W/NumW)*(Tx1+Tx2)/4;
		Fe[2*(NumW+1)*i+2*NumW+1]+=t*(W/NumW)*(Ty1+Ty2)/4;
		Fe[2*(NumW+1)*i+4*NumW+2]+=t*(W/NumW)*(Tx1+Tx2)/4;
		Fe[2*(NumW+1)*i+4*NumW+3]+=t*(W/NumW)*(Ty1+Ty2)/4;
	}
}
	   
void AssemblyGlobalMatrix(double ***Ke,int NumW,int NumL,double **Kg,int **Nodes)

{
	for(int k=0;k<NumL*NumW;k++)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int m=0;m<2;m++)
				{
					for(int n=0;n<2;n++)
					{
						Kg[2*Nodes[2*k][i]+m][2*Nodes[2*k][j]+n]+=Ke[0][2*i+m][2*j+n];
						Kg[2*Nodes[2*k+1][i]+m][2*Nodes[2*k+1][j]+n]+=Ke[1][2*i+m][2*j+n];
					}
				}
				
			}
		}
	}
	cout<<"The Global stiffness matrix is"<<endl;
	for(int i=0;i<2*(NumL+1)*(NumW+1);i++)
	{
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			cout<<Kg[i][j]<<" ";
		}
		cout<<endl;
	}
}



void ReduceLinearSystem(double **Kg,double **Kg_Reduced,double *Disp,double *F,double *Disp_Reduced,double *Force_Reduced,int NumW,int NumL)
{
	double **Temp;
	Temp=(double**)malloc(2*NumL*(NumW+1)*sizeof(double*));
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		*(Temp+i)=(double*)malloc(2*(NumL+1)*(NumW+1)*sizeof(double));
		for(int j=0;j<2*(NumL+1)*(NumW+1);j++)
		{
			Temp[i][j]=0.0;
		}
	}
	for(int k=0;k<NumW+1;k++)
	{
		for(int i=0;i<2*NumL;i++)
		{
			for(int j=0;j<2*(NumL+1)*(NumW+1);j++)	
			{
				Temp[k*2*NumL+i][j]=Kg[k*2*(NumL+1)+2+i][j];//Remove Rows with unknown forces
			}
			Disp_Reduced[k*2*NumL+i]=Disp[k*2*(NumL+1)+2+i];
			Force_Reduced[k*2*NumL+i]=F[k*2*(NumL+1)+2+i];
		}
	}
	
	for(int k=0;k<NumW+1;k++)
	{
		for(int j=0;j<2*NumL;j++)
		{
			for(int i=0;i<2*NumL*(NumW+1);i++)
			{
				Kg_Reduced[i][k*2*NumL+j]=Temp[i][k*2*(NumL+1)+2+j];//Remove columns with 0 displacements
			}
		}
	}

	cout<<"The Reduced Matrix is"<<endl;
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		for(int j=0;j<2*NumL*(NumW+1);j++)
		{
			cout<<Kg_Reduced[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<"Displacement  Force "<<endl;
	for(int i=0;i<2*NumL*(NumW+1);i++)
	{
		cout<<Disp_Reduced[i]<<" "<<Force_Reduced[i]<<endl;
	}
}

void Gauss(double **Kg,double *Disp,double *Force,int n)
{
	double c;
	for(int i=0;i<n-1;i++)
	{
		for(int k=i+1;k<n;k++)
		{
			c=-Kg[k][i]/Kg[i][i];
			for(int j=0;j<n;j++)
			{
				Kg[k][j]+=c*Kg[i][j];
			}
			Force[k]+=c*Force[i];
		}	
	}

	Disp[n-1]=Force[n-1]/Kg[n-1][n-1];

	for(int i=n-2;i>=0;i--)
	{
		double sum=0.0;
		for(int j=i+1;j<n;j++)
		{
			sum+=Kg[i][j]*Disp[j];
		}
		Disp[i]=(Force[i]-sum)/Kg[i][i];
	}

	cout<<"The answer for the linear system"<<endl;
	for(int i=0;i<n;i++)
	{
		cout<<Disp[i]<<endl;
	}
}

void Displacement(double *Disp,double *Disp_Reduced,int NumW,int NumL)
{
	for(int i=0;i<NumW+1;i++)
	{
		for(int j=0;j<NumL+1;j++)
		{
			if(j==0)
			{
				Disp[2*((NumW+1)*i+j)]=0.0;
				Disp[2*((NumW+1)*i+j)+1]=0.0;
			}
			else
			{
				Disp[2*((NumW+1)*i+j)]=Disp_Reduced[2*(NumW*i+j-1)];
				Disp[2*((NumW+1)*i+j)+1]=Disp_Reduced[2*(NumW*i+j)-1];
			}

			
		}
	}
	
	cout<<" Node "<<"      Ux      "<<"      Uy     "<<endl;
	for(int i=0;i<(NumL+1)*(NumW+1);i++)
	{
		cout<<"Node"<<i<<" "<<Disp[2*i]<<" "<<Disp[2*i+1]<<endl;
	}
}

void NodeForce(double *Force,double **Kg,double *Disp,int NumL,int NumW)
{
	for(int i=0;i<2*(NumW+1)*(NumL+1);i++)
		{
			double sum=0.0;
			for(int j=0;j<2*(NumW+1)*(NumL+1);j++)
			{
				sum+=Kg[i][j]*Disp[j];
			}
			Force[i]=sum;
		}

	cout<<"Node          Fx            Fy"<<endl;
	for(int i=0;i<(NumW+1)*(NumL+1);i++)
	{
		cout<<"Node"<<i<<" "<<Force[2*i]<<" "<<Force[2*i+1]<<endl;
	}
}

void StressStrain(double **strain,double **stress,double *StrainEnergy,double ***Ke,double ***B,double *Disp,int **Nodes,int NumW,int NumL,double **D)
{
	double *temp1,*temp2,*temp3,*temp4;
	temp1=(double*)malloc(6*sizeof(double));
	temp2=(double*)malloc(6*sizeof(double));
	temp3=(double*)malloc(6*sizeof(double));
	temp4=(double*)malloc(6*sizeof(double));
	for(int i=0;i<6;i++)
	{
		temp1[i]=0.0;
		temp2[i]=0.0;
		temp3[i]=0.0;
		temp4[i]=0.0;
	}
	for(int i=0;i<NumW*NumL;i++)
	{
		for(int k=0;k<3;k++)
		{
			temp1[2*k]=Disp[2*Nodes[2*i][k]];
			temp1[2*k+1]=Disp[2*Nodes[2*i][k]+1];
			temp2[2*k]=Disp[2*Nodes[2*i+1][k]];
			temp2[2*k+1]=Disp[2*Nodes[2*i+1][k]+1];
		}
		
		for(int m=0;m<6;m++)
		{
			double sum1=0.0,sum2=0.0;
			for(int n=0;n<6;n++)
			{
				sum1+=Ke[0][m][n]*temp1[n];
				sum2+=Ke[1][m][n]*temp2[n];
			}
			temp3[m]=sum1;
			temp4[m]=sum2;
		}
		
		for(int k=0;k<6;k++)
		{
			StrainEnergy[2*i]+=0.5*temp3[k]*temp1[k];
			StrainEnergy[2*i+1]+=0.5*temp4[k]*temp2[k];
		}

		for(int m=0;m<3;m++)
		{
			double sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0;
			for(int n=0;n<6;n++)
			{
				sum1+=B[0][m][n]*temp1[n];
				sum2+=B[1][m][n]*temp2[n];
			}
			strain[2*i][m]=sum1;
			strain[2*i+1][m]=sum2;
		}
	}

	for(int i=0;i<NumW*NumL;i++)
	{
		for(int m=0;m<3;m++)
		{
			double sum1=0.0,sum2=0.0;
			for(int n=0;n<3;n++)
			{
				sum1+=D[m][n]*strain[2*i][n];
				sum2+=D[m][n]*strain[2*i+1][n];
			}
			stress[2*i][m]=sum1;
			stress[2*i+1][m]=sum2;
		}
	}

	cout<<"Element   ex    ey   gxy"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		cout<<"Elem"<<2*i<<" "<<strain[2*i][0]<<" "<<strain[2*i][1]<<" "<<strain[2*i][2]<<endl;
		cout<<"Elem"<<2*i+1<<" "<<strain[2*i+1][0]<<" "<<strain[2*i+1][1]<<" "<<strain[2*i+1][2]<<endl;
	}
	cout<<endl;
	cout<<"Element   sx    sy   sxy"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		cout<<"Elem"<<2*i<<" "<<stress[2*i][0]<<"MPa "<<strain[2*i][1]<<"MPa "<<strain[2*i][2]<<"MPa"<<endl;
		cout<<"Elem"<<2*i+1<<" "<<stress[2*i+1][0]<<"MPa "<<strain[2*i+1][1]<<"MPa "<<strain[2*i+1][2]<<"MPa"<<endl;
	}
	cout<<endl;

	cout<<"Element   strain energy"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		cout<<"Elem"<<2*i<<" "<<StrainEnergy[2*i]<<"mJ"<<endl;
		cout<<"Elem"<<2*i+1<<" "<<StrainEnergy[2*i+1]<<"mJ"<<endl;
	}
	cout<<endl;
}

void VonMisesStress(double *Von_Mises,double **stress,int NumW,int NumL)
{
	for(int i=0;i<NumW*NumL;i++)
	{
		Von_Mises[2*i]=sqrt(((stress[2*i][0]-stress[2*i][1])*(stress[2*i][0]-stress[2*i][1])+stress[2*i][0]*stress[2*i][0]+stress[2*i][1]*stress[2*i][1]+6*stress[2*i][2]*stress[2*i][2])/2);
		Von_Mises[2*i+1]=sqrt(((stress[2*i+1][0]-stress[2*i+1][1])*(stress[2*i+1][0]-stress[2*i+1][1])+stress[2*i+1][0]*stress[2*i+1][0]+stress[2*i+1][1]*stress[2*i+1][1]+6*stress[2*i+1][2]*stress[2*i+1][2])/2);
	}


	cout<<"Element Von_Mises Stress"<<endl;
	for(int i=0;i<NumW*NumL;i++)
	{
		cout<<"Elem"<<2*i<<" "<<Von_Mises[2*i]<<endl;
		cout<<"Elem"<<2*i+1<<" "<<Von_Mises[2*i+1]<<endl;
	}
	cout<<endl;
}