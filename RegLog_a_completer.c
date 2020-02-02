#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define EPS 0.001

#define p 40 /* Nombre de points dans chacun des fichiers */
#define n 3  /* Nombre de dimension + 1 dimension pour la constante de theta */


/* Resolution du systeme       
       C.x + C0 = 0 
   pour une matrice C de taille nxn et un vecteur C0 de taille n
   Retourne 1 si le systeme a une unique solution et 0 sinon 
   Dans le cas vraie, le vecteur x contient l'unique solution */
int gauss(double **C, double *C0, double *x) {
     int i, j, k;
     int imin;
     double pivot;
     double sum, valmin, temp1, temp2;
     
     for(k=0;k<n-1;k++) {
        /* Recherche de l'element de valeur absolue minimum non nulle
           dans la colonne k et d'indice i>=k.*/
        
        valmin = C[k][k] ; imin = k ;
        for(i = k+1 ; i < n ; i++){
           if (valmin != 0){
              if (abs(C[i][k]) < abs(valmin) && C[i][k] != 0){
                 valmin = C[i][k] ;
                 imin = i ;
              }
           }
           else {
                 valmin = C[i][k] ;
                 imin = i ;
           }     
        }
        
        /* Si l'element minimum est nul, on peut en déduire
           que le systeme n'a pas de solution ou n'a pas de solution
           unique.*/
        
        if (valmin == 0.) return 0;
        
        /* Sinon inversion des elements de la ligne imax avec les elements
           de la ligne k. On fait de meme avec le vecteur C0. */
        
        for(j=0;j<n;j++){
           temp1 = C[imin][j] ;
           C[imin][j] = C[k][j] ;
           C[k][j] = temp1 ;
        }
        
        temp2 = C0[imin] ;
        C0[imin] = C0[k] ;
        C0[k] = temp2 ;
                
        /* Reduction de la matrice par la methode d'élimination de Gauss */
        
        for(i=k+1;i<n;i++) {
           pivot = C[i][k]/C[k][k];           
           for(j=0;j<n;j++)
	     C[i][j] = C[i][j] - pivot*C[k][j]; 
           C0[i] = C0[i] - pivot*C0[k]; 
        }
     }   
          
     if (C[n-1][n-1] == 0) return 0;
     
     /* Deduction de l'unique solution */
     
     x[n-1] = -C0[n-1]/C[n-1][n-1] ;
     
     for(i=n-2;i>-1;i--){
           sum = 0 ;
           for(j=n-1;j>i;j--)
              sum = sum + C[i][j]*x[j] ;
           x[i] = (-C0[i] - sum)/C[i][i] ;
     }
     return 1;
}


int main(){

  int i,k,l;
 
  double **A,**B;
  A=(double**) malloc(sizeof(double*)*p);
  for (i=0;i<p;i++) A[i]=(double*) malloc(sizeof(double)*2);
  B=(double**) malloc(sizeof(double*)*p);
  for (i=0;i<p;i++) B[i]=(double*) malloc(sizeof(double)*2);

  double **C;
  double *C0;
  C=(double**) malloc(sizeof(double*)*3);
  for (i=0;i<3;i++) C[i]=(double*) malloc(sizeof(double)*3);
  C0=(double*) malloc(sizeof(double)*3);

  double *pT,*qT;
  pT=(double*) malloc(sizeof(double)*3);
  qT=(double*) malloc(sizeof(double)*3);
	  
  FILE *f1=fopen("Base_apprentissage_Notes_Admis.dat","r");
  if (!f1){printf("Pbm fichier\n"); return 0;}
  for (i=0;i<p;i++)
    for (k=0;k<n-1;k++){
      fscanf(f1,"%lf",&(A[i][k]));
      A[i][n-1]=1;  /* Permet de creer une constante dans le produit scalaire theta^T.A */
    }
  fclose(f1);
  FILE *f2=fopen("Base_apprentissage_Notes_NonAdmis.dat","r");
  for (i=0;i<p;i++)
    for (k=0;k<n-1;k++){
      fscanf(f2,"%lf",&(B[i][k]));
      B[i][n-1]=1;  /* Permet de creer une constante dans le produit scalaire theta^T.B */
    }
  fclose(f2);



  /* .... */

  
  
  return 0;
}
