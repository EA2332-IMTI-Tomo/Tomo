#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>



#include <vMatrix.h>  // template
#include <vParse.h>


using namespace std;


// =============================================================================



static void
usage(int argc, char **argv) 
{
  if (argc - 1 != 4)
    {
      printf("Calcul pour diane\n");
      printf("Usage: %s <fichier W> <fichier Theta> <fichier Megamat> \n", argv[0]);
      exit(1);
    }  
}
// typiquement: 800 et 400


// =============================================================================

int
main(int argc, char** argv)
{
  usage(argc, argv);

  char *W_file, *Theta_file, *MegaMat_file; 
  double mu;
  W_file = argv[1];
  Theta_file = argv[2];
  MegaMat_file = argv[3];
  mu = atof(argv[4]);
  
  size_t Wnlig, Wncol;
  size_t Mnlig, Mncol;
  
  
  if (!retrieveFileDimensions(W_file, 1000, Wnlig, Wncol, 1))
    {cerr << endl << "fichier de données non carrées"; return 1;}

  cout << endl << Wnlig << " " << Wncol << " ";
    
  if (!retrieveFileDimensions(MegaMat_file, 1000, Mnlig, Mncol, 1))
    {cerr << endl << "fichier de données non carrées"; return 1;}
  
  
  cout << endl << "Dimensions de megamat: " << Mnlig << " " << Mncol << " ";


  double** w;
  w = readMatrixFromFile(W_file, 1000, 7, 7, 1, 1);
  vMatrix<double> W(7, 7, w);
  
  
  // pour voir 16 décimales
  cout.precision(16);
  cout << endl << "W: " << W;

  
  double** megamat;
  Mncol--; Mnlig--;// we don't want 1st line and 1st tokens
  megamat = readMatrixFromFile(MegaMat_file, 1000, Mnlig, Mncol, 1, 1); 


  // 1: skip first token of each line
  vMatrix<double> M(Mnlig, Mncol, megamat);
    

  double** theta;
  theta = readMatrixFromFile(Theta_file, 1000, 7, 1, 1, 1);
  vMatrix<double> Theta(7, 1, theta);

  //vMatrix<double> tTheta(7, 1);
  //tTheta.transpose(Theta);

  cout << endl << "Theta: " << Theta;
  cout << endl << "mu: " << mu;
  
  // --------------------------------------------------
  // TA = mA * mB * M
  // TA: {1,7}
  // mA: {1, 30k} que des 1
  // mB: {30k, 30k}
  // mC: {30k, 7}
  
  // --------------------------------------------------
  // M déjà prête, de taille {Mnlig Mncol}
  vMatrix<double> mA(1, Mnlig, 1);
  // mB: trop grosse et creuse: on ne l'alloue pas

  // --------------------------------------------------
  //mB = diag(exp(M * Theta)
  //M * tTheta: {30k, 1}
  vMatrix<double> B(Mnlig, 1);
  vMatrix<double> tB(1, Mnlig);
  B.multiply(M, Theta);
  // visu provisoire
  tB.transpose(B);
  tB.display_limited = true;
  cout << endl << "megamat * theta (t) " << tB << endl;

  B.map(&exp);
  tB.transpose(B);
  

  // --------------------------------------------------
  //mB * M calculée directement => mBM
  vMatrix<double> mBM(Mnlig, 7);
  mBM.multiply_from_diag(tB, M);
  
  // --------------------------------------------------
  // TA = mA * mBM
  vMatrix<double> ma(1, Mnlig, 1);
  vMatrix<double> TA(1, 7);
  TA.multiply(ma, mBM);
  
  // --------------------------------------------------
  // varmu = TA * W * tTA  // {1, 7} * {7, 7} * {7, 1} // scalaire
  vMatrix<double> TAW(1, 7);
  TAW.multiply(TA, W);
  vMatrix<double> tTA(7, 1);
  tTA.transpose(TA);
  vMatrix<double> vM(1, 1);
  vM.multiply(TAW, tTA);
  
  double varmu = vM[0][0];
  //double mu = 1849437;
  double ICl = mu - 2*sqrt(varmu);
  double ICh = mu + 2*sqrt(varmu);
  

  // --------------------------------------------------
  // --------------------------------------------------

  cout << endl <<  "--------------------------------------------------" << endl <<  "--------------------------------------------------";

  tB.display_limited = true;
  cout.precision(22);
  //cout << endl << tB << endl;

  M.display_limited = true;
  //cout << M;

  cout << endl << "t1: " << TA;
  cout << endl << "varmu: " << varmu;
  cout << endl << "IC: [" << ICl << " ; " << ICh << " ]" << endl;

  return 0;
}


