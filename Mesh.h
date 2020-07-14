#ifndef _MESH_ 
#define _MESH_ 

void CSTgeneration(double L,double W,int NumL,int NumW,double ***J,double detJ,int **Node);

void B_Matrix(double ***B,double ***J,double detJ);

void ElementalStiffnessMatrix(double ***Ke,double **D,double t,double ***B,double detJ);

void EquivalentNodalForce(double Ta,double Tb,double L,double W,int NumW,int t,double *Fe);

void AssemblyGlobalMatrix(double ***Ke,int NumW,int NumL,double **Kg,int **Nodes);

void ReduceLinearSystem(double **Kg,double **Kg_Reduced,double *Disp,double *F,double *Disp_Reduced,double *Force_Reduced,int NumW,int NumL);

void Gauss(double **Kg,double *Disp,double *Force,int n);

void Displacement(double *Disp,double *Disp_Reduced,int NumW,int NumL);

void NodeForce(double *Force,double **Kg,double *Disp,int NumL,int NumW);

void StressStrain(double **strain,double **stress,double *StrainEnergy,double ***Ke,double ***B,double *Disp,int **Nodes,int NumW,int NumL,double **D);

void VonMisesStress(double *Von_Mises,double **stress,int NumW,int NumL);

#endif